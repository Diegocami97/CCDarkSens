// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  BackgroundFactory.hh -- Header for building the background vector(s) a
//  generic scan app needs, mirroring the reference app's
//  background_source/background_model wiring.
// ===========================================================================

#pragma once

#include <vector>

#include "ccdarksens/experiment/ExperimentSetup.hh"
#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/response/BackgroundEfficiencyTable.hh"
#include "ccdarksens/response/ResponseFactory.hh"

namespace ccdarksens {

class ResponseFold;

/**
 * B_pat — nominal background per output bin (Bp + Br at theta=1). Always
 *         populated regardless of background_source; feed this to
 *         ProfileLikelihood::SetBTemplate when run.background_model=="scale".
 * Bp, Br — always populated too, so background_model=="Bp_theta_Br" works
 *          regardless of background_source: when background_source ==
 *          "bp_br_template", these are the config arrays directly; when
 *          "dc_flat_migration", Bp is the zero vector and Br == B_pat (the
 *          same degenerate trick the reference app's own n_e-space branch
 *          uses to reuse the theta-constrained-prior/2D-minimizer machinery
 *          for a plain scale-equivalent model).
 */
struct BackgroundResult {
  std::vector<double> B_pat;
  std::vector<double> Bp;
  std::vector<double> Br;
};

/**
 * Builds background per run.background_source:
 *  - "dc_flat_migration" (default): dark current via BackgroundBuilder ::
 *    BuildBkgAsimov() folded through fold.FoldNeBackground(...), plus flat
 *    d.r.u. via MakeFlatDrdeSpectrum(...) folded through fold.Fold(...),
 *    summed elementwise. Only valid where FoldNeBackground is implemented
 *    (pattern/n_e space) -- ResponseFold's default throws a clear error for
 *    cluster_energy, which has no dark-current model.
 *  - "bp_br_template": run.background_Bp/run.background_Br loaded directly
 *    from config, size-checked against fold.NumBins() (generically, unlike
 *    the reference app which hardcodes the check against pattern_roi.size()).
 *
 * response must be the ResponseFactoryResult that built `fold` (its `ion`,
 * `ne_min_bkg`, `ne_max` are needed for the dark-current/flat fold path).
 */
BackgroundResult MakeBackground(const ConfigManager& cfg,
                                 const ExperimentSummary& summary,
                                 const ResponseFold& fold,
                                 const ResponseFactoryResult& response);

/**
 * Flat-Compton-only background for one cluster_energy channel of a
 * joint-likelihood config (response.channels[]): folds a flat dR/dE(E)
 * spectrum at flat_rate_dru through fold.Fold(...), same construction as
 * MakeBackground's dc_flat_migration path but with an explicit per-channel
 * rate instead of cfg.backgrounds().flat_background, and no dark-current
 * term (cluster_energy never supports one -- see ResponseFold::
 * SupportsDarkCurrentBackground()).
 */
BackgroundResult MakeClusterEnergyFlatBackground(const ResponseFold& fold,
                                                  double exposure_kg_year,
                                                  double Emin_eV, double Emax_eV, int nbins,
                                                  double flat_rate_dru);

/**
 * Flat-Compton background for one cluster_energy channel, using the
 * background's OWN energy-dependent detection efficiency curve (e.g.
 * PhysRevD.94.082006 Fig. 9's dashed lines) instead of folding through the
 * signal's kernel (which is what MakeClusterEnergyFlatBackground /
 * MakeBackground's dc_flat_migration path do). Matches the paper's own
 * stated construction (Sec. IV.3: "shape is given by a flat Compton
 * scattering energy spectrum multiplied by the background efficiency") --
 * evaluated directly per E_reco bin, no signal-kernel smearing step (a flat
 * spectrum is invariant under resolution convolution, so that step is a
 * no-op for background specifically; see docs/ClusterFitMC_Design.md Sec.
 * 6.9). fold must support ErecoEdgesEV() (cluster_energy only).
 */
BackgroundResult MakeClusterEnergyFlatBackgroundWithEfficiency(
    const ResponseFold& fold,
    double exposure_kg_year,
    double flat_rate_dru,
    const BackgroundEfficiencyTable& eff_table);

/**
 * Low-energy excess (LEE) background for one cluster_energy channel --
 * DAMIC SNOLAB's unexplained bulk ionization excess, dR/dE =
 * rate_per_kg_day/day * (1/epsilon_eV)*exp(-E/epsilon_eV) (Tier A: fixed
 * shape, no free nuisance -- see docs/LowEnergyExcess_Design.md). Evaluated
 * as the exact closed-form integral of the normalized exponential over each
 * fold.ErecoEdgesEV() bin -- unlike MakeClusterEnergyFlatBackgroundWithEfficiency's
 * digitized efficiency curve, the exponential has a closed form, so no
 * Simpson's-rule sub-sampling is needed. cluster_energy only.
 */
BackgroundResult MakeClusterEnergyLEEBackground(const ResponseFold& fold,
                                                 double exposure_kg_year,
                                                 double rate_per_kg_day,
                                                 double epsilon_eV);

}  // namespace ccdarksens
