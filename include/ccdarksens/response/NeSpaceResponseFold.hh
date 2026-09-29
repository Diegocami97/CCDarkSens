// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  NeSpaceResponseFold.hh -- Header for the single-pixel n_e-space
//  ResponseFold implementation, wrapping ChargeIonization::FoldToNe + a
//  precomputed epsilon(n_e) efficiency curve.
// ===========================================================================

#pragma once

#include <map>
#include <memory>
#include <vector>

#include "ccdarksens/response/ResponseFold.hh"

namespace ccdarksens {

class ChargeIonization;

// ----------------------------------------------------------------------------
// NeSpaceResponseFoldConfig
//   Inputs of the n_e-space fold, built once by the caller.
// ----------------------------------------------------------------------------
struct NeSpaceResponseFoldConfig {
  int ne_min = 0;  // lowest n_e of the ionization fold
  int ne_max = 0;  // highest n_e of the ionization fold
  std::vector<int> roi_bins;      ///< n_e values the likelihood is defined over (experiment.roi_bins)
  std::map<int, double> eps_ne;   ///< n_e -> efficiency, flattened from EfficiencyMC::PrecomputeEpsilon's TH1D
};

/**
 * Wraps the reference app's n_e-space ("single-pixel") fold sequence
 * (ccdarksens_scan_dmelectron_pattern.cc, observable_bins == "n_e" branch):
 *   S_true = ion->FoldToNe(dRdE, exposure, ne_min, ne_max)
 *   S_obs[ne] = S_true[ne] * clamp(eps_ne[ne], 0, 1)     for ne in roi_bins
 * A genuinely different observable space from pattern-space: no multi-pixel
 * classification, just per-n_e detection efficiency. Distinct from
 * PatternResponseFold, selected by experiment.observable_bins == "n_e"
 * (default) rather than response.analysis_space.
 */
class NeSpaceResponseFold : public ResponseFold {
 public:
  // Constructor: ionization model and the fold configuration.
  NeSpaceResponseFold(std::shared_ptr<ChargeIonization> ion, NeSpaceResponseFoldConfig cfg);

  std::vector<double> Fold(const TH1D& dRdE_density, double exposure_kg_year) const override;
  std::size_t NumBins() const override { return cfg_.roi_bins.size(); }
  std::vector<std::string> BinLabels() const override;
  std::vector<double> FoldNeBackground(const TH1D& B_ne_asimov) const override;
  bool SupportsDarkCurrentBackground() const override { return true; }

 private:
  std::shared_ptr<ChargeIonization> ion_;  // P(n_e | E) folding
  NeSpaceResponseFoldConfig cfg_;  // configuration
};

}  // namespace ccdarksens
