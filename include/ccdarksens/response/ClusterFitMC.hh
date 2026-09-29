// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ClusterFitMC.hh -- Phase 4 of the WIMP-nucleus SI channel: forward-
//  simulates real nuclear-recoil-like pixel events across a grid of true
//  (electron-equivalent) energies, scores each with the Slice 1-2 fit
//  engine, applies the Slice 3 noise-tail detection cut and a fiducial
//  (surface-rejection) cut on the fitted width, and bins the survivors into
//  a kernel K[E_true, E_reco] -- the continuous-energy analogue of
//  EfficiencyMC's P(pattern|n_e) table, used where n_e reaches into the
//  hundreds and a discrete pattern classifier no longer applies.
// ===========================================================================

#pragma once

#include <memory>
#include <random>
#include <vector>

#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/ClusterFitEngine.hh"
#include "ccdarksens/response/PixelSimulator.hh"

namespace ccdarksens {

// ----------------------------------------------------------------------------
// KernelMatrix
//   Detector response for the cluster analysis: fraction of events at each true energy that are accepted and reconstructed in each E_reco bin.
// ----------------------------------------------------------------------------
struct KernelMatrix {
  std::vector<double> Etrue_grid_eV;   ///< the E_true grid points BuildKernel was called with
  std::vector<double> Ereco_edges_eV;  ///< size NEreco+1
  /// K[iEtrue][iEreco] = fraction of trials at Etrue_grid_eV[iEtrue] that
  /// were accepted (beat the ΔLL cut and the fiducial cut) and reconstructed
  /// into Ereco bin iEreco. Rows do NOT sum to 1 -- the shortfall from 1 is
  /// exactly the fraction rejected (missed entirely).
  std::vector<std::vector<double>> K;
  /// Convenience: total accepted fraction per Etrue point (row sums of K).
  /// This is the detection-efficiency curve, the primary validation target.
  std::vector<double> efficiency;
};

// ----------------------------------------------------------------------------
// ClusterFitMCConfig
//   Settings of the forward simulation that builds the kernel.
// ----------------------------------------------------------------------------
struct ClusterFitMCConfig {
  int ne_trials_per_point = 5000;  ///< toy trials per E_true grid point; mirrors EfficiencyMCConfig::ne_trials

  /// Pixel window and readout noise. MUST match the config the delta_ll_cut
  /// below was calibrated with (see NoiseTailCalibrator) -- a mismatched
  /// window size or sigma_readout_e means delta_ll_cut no longer describes
  /// this window's actual noise behavior.
  PixelSimulatorConfig pix_cfg;

  /// Fit engine configuration. sigma_pix_e here must match
  /// pix_cfg.sigma_readout_e, same requirement as NoiseTailCalibratorConfig.
  ClusterFitConfig fit_cfg;

  double delta_ll_cut = 0.0;  ///< from NoiseTailCalibrationResult (Slice 3)

  /// Fiducial (surface-rejection) cut on the fitted width, in pixel units --
  /// PhysRevD.94.082006's stated range.
  double sigma_xy_fid_min_px = 0.35;
  double sigma_xy_fid_max_px = 1.22;

  /// Charge-generation statistics for sampling a true electron count from
  /// E_true (see docs/ClusterFitMC_Design.md for the derivation and the
  /// paper citation): n_e ~ round(Gaussian(mean = E_true/eh_pair_eV,
  /// variance = fano_factor * mean)). Values are PhysRevD.94.082006's own:
  /// eh_pair_eV = 3.77 eV, fano_factor = 0.133 +- 0.005 (their measured
  /// value for ionizing particles -- the paper states this is UNKNOWN for
  /// nuclear recoils specifically and varies it from 0.13 to 1.0 as a
  /// systematic check, which is why this is a real config knob and not a
  /// hardcoded constant).
  double eh_pair_eV = 3.77;
  double fano_factor = 0.133;

  std::uint64_t rng_seed = 987654321ULL;  ///< seeds the n_e sampling step specifically
};

// ----------------------------------------------------------------------------
// ClusterFitMC
//   I simulate nuclear-recoil-like events at a grid of true energies (Fano-
//   smeared electron number, depth and diffusion from ChargeTransport, pixel
//   deposition, noise), fit each one, apply the Delta LL cut and the fiducial
//   cut on the fitted width, and histogram the reconstructed energies into
//   KernelMatrix.
// ----------------------------------------------------------------------------
class ClusterFitMC {
 public:
  /// ct is shared (not owned) so the same ChargeTransport instance used
  /// elsewhere in a pipeline can be reused here, mirroring EfficiencyMC's
  /// constructor convention.
  ClusterFitMC(const ClusterFitMCConfig& cfg, std::shared_ptr<ChargeTransport> ct);

  /// Builds the kernel across Etrue_grid_eV (point values, not bin edges --
  /// each is simulated independently) and Ereco_edges_eV (bin edges).
  /// Deterministic given cfg.rng_seed / pix_cfg.rng_seed / ct's own seed.
  KernelMatrix BuildKernel(const std::vector<double>& Etrue_grid_eV,
                            const std::vector<double>& Ereco_edges_eV) const;

 private:
  ClusterFitMCConfig cfg_;  // settings
  std::shared_ptr<ChargeTransport> ct_;  // shared depth/diffusion sampler
  mutable std::mt19937_64 rng_;  // RNG for the electron-number sampling
  mutable std::normal_distribution<double> gaus_{0.0, 1.0};  // unit Gaussian for the Fano smearing
};

}  // namespace ccdarksens
