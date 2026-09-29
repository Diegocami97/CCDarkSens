// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ResponseFold.hh -- Header for the generic "dR/dE density spectrum -> per-
//  bin expected counts" contract shared by every analysis space.
// ===========================================================================

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

class TH1D;

namespace ccdarksens {

/**
 * Abstract "fold a dR/dE density spectrum into per-bin expected counts"
 * contract. One implementation per analysis_space (pattern, cluster_energy,
 * ...). Built ONCE per scan run, before the grid loop -- the expensive
 * components each implementation wraps (EfficiencyMC's pattern table,
 * ClusterFitMC's kernel) do not depend on (mass, coupling), so this object
 * is constructed once and its Fold() is called once per grid point.
 *
 * This is deliberately NOT DetectorResponsePipeline::Apply() -- reading the
 * reference app (ccdarksens_scan_dmelectron_pattern.cc) shows the actual
 * signal-path fold bypasses that class entirely (explicit code comment:
 * "Use S_true, not S_obs -- folding with S_obs would double-count") and
 * calls ChargeIonization::FoldToNe + FoldNeToPatternRates directly. That
 * two-step call, generalized, is exactly this interface's contract.
 *
 * Background is NOT a method here -- it's the same Fold() call applied to a
 * flat dummy spectrum (see MakeFlatDrdeSpectrum), matching how the
 * reference app already treats background: same fold path as signal, just
 * a different input histogram.
 */
class ResponseFold {
 public:
  virtual ~ResponseFold() = default;

  /**
   * dRdE_density: events/(kg*year*eV) histogram (e.g. from
   * ModelFactory::MakeSignalSpectrumE or RateTable::MakeTH1D).
   * Returns expected counts per output bin. Same length and same bin order
   * for every call on a given instance.
   */
  virtual std::vector<double> Fold(const TH1D& dRdE_density,
                                    double exposure_kg_year) const = 0;

  /**
   * Number of output bins (pattern_roi.size(), or Ereco_nbins). Invariant
   * for the lifetime of the object -- callers allocate B(bin) and call
   * ProfileLikelihood::SetBTemplate/SetData against this size exactly once,
   * before the grid loop.
   */
  virtual std::size_t NumBins() const = 0;

  /// Human-readable output-bin labels for diagnostics/dump histograms
  /// (pattern codes as strings, or "Ereco_bin_k"). Purely cosmetic.
  virtual std::vector<std::string> BinLabels() const = 0;

  /**
   * Fold a dark-current background already expressed as an n_e-space Asimov
   * histogram (BackgroundBuilder::BuildBkgAsimov's output) into this fold's
   * output bins. This does NOT fit the Fold() contract above -- dark current
   * is synthesized directly in n_e space from detector geometry/timing/lambda,
   * with no dR/dE spectrum involved, so it needs a separate entry point.
   *
   * Default throws: not every analysis space has a dark-current model
   * (cluster-energy has no pixel/readout concept), and a config asking for
   * background_source="dc_flat_migration" there is a real user error, not
   * something to silently no-op.
   */
  virtual std::vector<double> FoldNeBackground(const TH1D& B_ne_asimov) const {
    (void)B_ne_asimov;
    throw std::runtime_error(
        "FoldNeBackground: dark-current background not supported for this analysis space");
  }

  /// Whether FoldNeBackground is implemented -- callers (BackgroundFactory)
  /// check this before attempting the dark-current fold, so
  /// background_source="dc_flat_migration" degrades gracefully to
  /// flat-only for analysis spaces with no dark-current/pixel model
  /// (cluster_energy) instead of throwing.
  virtual bool SupportsDarkCurrentBackground() const { return false; }

  /**
   * Output-bin edges in eV, for analysis spaces where bins are a genuine
   * energy axis (cluster_energy's Ereco bins) rather than pattern/n_e
   * labels. Needed by BackgroundFactory's background-efficiency-curve path
   * (MakeClusterEnergyFlatBackgroundWithEfficiency), which evaluates an
   * energy-dependent background efficiency per bin directly -- see
   * docs/ClusterFitMC_Design.md Sec. 6.9. Default throws, same pattern as
   * FoldNeBackground: not every analysis space has an energy-axis concept.
   */
  virtual std::vector<double> ErecoEdgesEV() const {
    throw std::runtime_error("ErecoEdgesEV: no energy-bin-edge concept for this analysis space");
  }
};

}  // namespace ccdarksens
