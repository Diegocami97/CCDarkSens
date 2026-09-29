// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  NoiseTailCalibrator.hh -- Phase 3 of the WIMP-nucleus SI channel:
//  generates a pure-noise toy pixel-window ensemble (no signal deposited at
//  all), scores each with ClusterFitEngine, and calibrates a detection-cut
//  threshold on ΔLL from the resulting distribution's tail -- the analog of
//  PhysRevD.94.082006's Fig. 6. Empirical, not an asymptotic chi^2 formula:
//  see docs/ClusterFitMC_Design.md for why.
// ===========================================================================

#pragma once

#include <vector>
#include "ccdarksens/response/ClusterFitEngine.hh"
#include "ccdarksens/response/PixelSimulator.hh"

namespace ccdarksens {

// ----------------------------------------------------------------------------
// NoiseTailCalibratorConfig
//   Toy-ensemble settings for calibrating the pure-noise Delta LL cut.
// ----------------------------------------------------------------------------
struct NoiseTailCalibratorConfig {
  /// Pixel window geometry and readout noise for the toy ensemble.
  /// lambda_dc should normally be left at its PixelSimulatorConfig default
  /// (0) -- dark current is handled as a separate background term
  /// elsewhere in the pipeline, not folded into this noise model (same
  /// convention EfficiencyMC uses; see the "Efficiency MC" notes in the project documentation).
  /// A nonzero value adds Poisson DC to every pixel of every toy: a
  /// diagnostic knob for testing whether DC perturbs the ΔLL tail, not a
  /// production setting.
  /// pix_cfg.rng_seed controls the toy sequence -- fix it for a
  /// reproducible calibration.
  PixelSimulatorConfig pix_cfg;

  /// Fit engine configuration used to score every toy. sigma_pix_e here
  /// must match pix_cfg.sigma_readout_e.
  ClusterFitConfig fit_cfg;

  int n_toys = 100000;  // number of pure-noise toy windows
  double target_tail_prob = 1e-3;  ///< fraction of pure-noise toys allowed to beat the cut
};

// ----------------------------------------------------------------------------
// NoiseTailCalibrationResult
//   Output of the calibration: every toy's Delta LL and the resulting cut.
// ----------------------------------------------------------------------------
struct NoiseTailCalibrationResult {
  std::vector<double> delta_ll_toys;  ///< raw ΔLL per toy, sorted ascending (most signal-like first)
  double delta_ll_cut = 0.0;          ///< calibrated threshold at target_tail_prob
  int n_toys_used = 0;  // number of toys actually run
};

// ----------------------------------------------------------------------------
// NoiseTailCalibrator
//   I generate pure-noise pixel windows (no signal), fit each one, and take the
//   target_tail_prob quantile of the resulting Delta LL distribution as the
//   detection cut. This is my analogue of the noise-tail fit in the DAMIC
//   publication and is fully empirical.
// ----------------------------------------------------------------------------
class NoiseTailCalibrator {
 public:
  // Constructor: keep the calibration settings.
  explicit NoiseTailCalibrator(const NoiseTailCalibratorConfig& cfg);

  /// Runs the full toy ensemble and calibrates the cut. Deterministic given
  /// pix_cfg.rng_seed.
  NoiseTailCalibrationResult Run() const;

 private:
  NoiseTailCalibratorConfig cfg_;  // calibration settings
};

}  // namespace ccdarksens
