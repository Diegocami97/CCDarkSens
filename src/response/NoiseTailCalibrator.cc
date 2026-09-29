// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  NoiseTailCalibrator.cc -- Pure-noise toy ensemble + empirical ΔLL tail
//  calibration.
// ===========================================================================

#include "ccdarksens/response/NoiseTailCalibrator.hh"

#include <algorithm>

namespace ccdarksens {

// Constructor: keep the calibration settings.
NoiseTailCalibrator::NoiseTailCalibrator(const NoiseTailCalibratorConfig& cfg) : cfg_(cfg) {}

// ----------------------------------------------------------------------------
// NoiseTailCalibrator::Run
//   I run n_toys pure-noise windows (dark current if lambda_dc > 0, then
//   Gaussian readout noise; no deposited charge), fit each, sort the Delta LL
//   values ascending (most signal-like first) and pick the entry at index
//   target_tail_prob * n_toys as the cut. Deterministic for a given
//   pix_cfg.rng_seed.
// ----------------------------------------------------------------------------
NoiseTailCalibrationResult NoiseTailCalibrator::Run() const {
  PixelSimulator pixSim(cfg_.pix_cfg);
  ClusterFitEngine engine(cfg_.fit_cfg);

  NoiseTailCalibrationResult result;
  result.delta_ll_toys.reserve(static_cast<std::size_t>(cfg_.n_toys));

  for (int it = 0; it < cfg_.n_toys; ++it) {
    pixSim.Reset();
    pixSim.AddDarkCurrent();   // no-op (no RNG draws) unless pix_cfg.lambda_dc > 0
    pixSim.AddReadoutNoise();  // no DepositElectron() calls -- pure noise, no signal
    const auto r = engine.Fit(pixSim.PixelCharges(), pixSim.NX(), pixSim.NY());
    result.delta_ll_toys.push_back(r.delta_ll);
  }
  result.n_toys_used = cfg_.n_toys;

  // ΔLL <= 0, more negative = more signal-like. Sort ascending (most
  // signal-like first) so the target_tail_prob quantile from the front of
  // the list is the cut: only that fraction of pure-noise toys scored at
  // least this well by chance.
  std::sort(result.delta_ll_toys.begin(), result.delta_ll_toys.end());
  int idx = static_cast<int>(cfg_.target_tail_prob * cfg_.n_toys);
  idx = std::max(0, std::min(cfg_.n_toys - 1, idx));
  result.delta_ll_cut = result.delta_ll_toys[static_cast<std::size_t>(idx)];

  return result;
}

}  // namespace ccdarksens
