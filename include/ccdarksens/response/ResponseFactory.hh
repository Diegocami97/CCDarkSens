// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ResponseFactory.hh -- Header for building a ResponseFold from config,
//  mirroring the reference app's own
//  EfficiencyMC/ChargeIonization/ClusterFitMC construction.
// ===========================================================================

#pragma once

#include <memory>
#include <string>

#include "ccdarksens/experiment/ExperimentSetup.hh"
#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/response/ResponseFold.hh"

namespace ccdarksens {

class ChargeIonization;

/**
 * fold           — the constructed ResponseFold (pattern | n_e | cluster_energy).
 * ion            — the ChargeIonization used to build `fold`, exposed because
 *                   Phase-A background construction (background_source ==
 *                   "dc_flat_migration") folds the flat-Compton background
 *                   through the same ion->FoldToNe(...) path signal uses, for
 *                   pattern/n_e space. Null for cluster_energy (no such path
 *                   there — that channel's flat background goes through
 *                   fold->Fold() directly, see MakeFlatDrdeSpectrum).
 * ne_min_bkg     — n_e lower bound to use for background construction.
 *                   Pattern space includes n_e=0 (dark current has n_e=0);
 *                   n_e space does not. 0 for cluster_energy (unused there).
 * ne_max         — n_e upper bound for background construction (mirrors the
 *                   reference app's ne_max, the same for signal and background).
 */
struct ResponseFactoryResult {
  std::unique_ptr<ResponseFold> fold;
  std::shared_ptr<ChargeIonization> ion;
  int ne_min_bkg = 0;
  int ne_max = 0;
};

/**
 * Builds the ResponseFold selected by response.analysis_space
 * ("cluster_energy") first, else experiment.observable_bins ("pattern" |
 * "n_e"), replicating the reference app's own construction sequence
 * (ccdarksens_scan_dmelectron_pattern.cc, detector-response-component setup
 * through pattern_eff_map construction) so a generic scan app can be
 * config-driven the same way. Built ONCE per run, before the grid loop --
 * this is the expensive step (EfficiencyMC's pattern-table MC, or
 * ClusterFitMC's kernel).
 *
 * config_path is needed to resolve response.efficiency_mc.efficiency_csv /
 * efficiency_csv_reference paths relative to the config file's directory,
 * exactly as the reference app does.
 *
 * Diagnostic-only paths from the reference app (use_2d_image_efficiency,
 * include_dc_pileup) are intentionally not replicated here -- the project documentation
 * documents both as non-production flags; the 1D-row EfficiencyMC path is
 * the production default and the only one this factory builds.
 *
 * analysis_space == "pcd" throws std::runtime_error: the reference app's PCD
 * wiring exists but its own signal path bypasses it too, so there is no
 * working PCD fold to mirror.
 */
ResponseFactoryResult MakeResponseFold(const ConfigManager& cfg,
                                        const ExperimentSummary& summary,
                                        const std::string& config_path);

/**
 * Builds a cluster_energy ResponseFold (NoiseTailCalibrator -> ChargeTransport
 * -> ClusterFitMC::BuildKernel) directly from a ClusterFitMCJSON, without a
 * full ConfigManager -- the per-channel building block used by joint-channel
 * configs (response.channels[], e.g. DAMIC's 1x1 + 1x100 readout modes) to
 * build one independent ResponseFold per channel. MakeResponseFold's
 * cluster_energy branch calls this with cfg.response().cluster_fit_mc, so the
 * single-channel path is unchanged.
 */
ResponseFactoryResult MakeClusterEnergyFoldFromJSON(const ClusterFitMCJSON& jcfg);

}  // namespace ccdarksens
