// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: ccdarksens_calibrate_noise_tail.cc
//  Slice 3 validation: runs the pure-noise toy ensemble at increasing
//  n_toys to check the calibrated ΔLL cut converges, re-runs at a fixed
//  n_toys with several independent seeds to quantify the calibration's own
//  statistical uncertainty, and prints a shape summary of the ΔLL
//  distribution (the Fig. 6 analog) -- no external reference dataset
//  exists for this, so this checks self-consistency/stability rather than
//  parity to a formula. See docs/ClusterFitMC_Design.md.
//
//  Usage:
//    ccdarksens_calibrate_noise_tail [--n_toys_full N]
// ===========================================================================

#include "ccdarksens/response/NoiseTailCalibrator.hh"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace ccdarksens;

namespace {

// ----------------------------------------------------------------------------
// MakeConfig
//   Calibration settings for the validation: 15x15 window, 0.16 e- readout noise, no
//   dark current, Nelder-Mead at full precision, target tail probability 1e-3.
// ----------------------------------------------------------------------------
NoiseTailCalibratorConfig MakeConfig(int n_toys, std::uint64_t seed) {
  NoiseTailCalibratorConfig cfg;
  cfg.pix_cfg.nx = 15;
  cfg.pix_cfg.ny = 15;
  cfg.pix_cfg.pixel_size_um = 15.0;
  cfg.pix_cfg.sigma_readout_e = 0.16;
  cfg.pix_cfg.lambda_dc = 0.0;  // DC handled as a separate background term, not here
  cfg.pix_cfg.rng_seed = seed;
  cfg.fit_cfg.sigma_pix_e = 0.16;
  cfg.fit_cfg.method = ClusterFitConfig::Method::kNelderMead;  // Minuit2 already cross-checked in Slice 2
  // Deliberately NOT lowering max_iterations below ClusterFitEngine's
  // accuracy-first default (300) here: measured directly (see
  // docs/ClusterFitMC_Design.md §3) that doing so biases the CALIBRATED
  // CUT itself -- not just individual toy imprecision -- by ~11% toward
  // being too lenient (an under-converged toy ensemble underestimates how
  // good noise can look, so the resulting threshold lets more real
  // false-positives through than the target_tail_prob promises). This
  // makes pure-noise toy fitting the expensive step (no real minimum to
  // converge toward, so every toy burns the full iteration budget) but
  // that cost is accepted deliberately -- see the runtime note printed
  // below.
  cfg.n_toys = n_toys;
  cfg.target_tail_prob = 1e-3;
  return cfg;
}

}  // namespace

// ----------------------------------------------------------------------------
// main
//   Validation of the noise-tail calibration:
//     1) convergence: the calibrated cut must stabilize as the number of toys grows;
//     2) stability: the cut must not depend too strongly on the seed;
//     3) shape: the Delta LL distribution must be sorted and one-sided (<= 0).
//   Option: --n_toys_full N (default 8000; the production default is 100000).
//   Exit code 0 if all checks pass.
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
  // Measured cost at full (accuracy-preserving) precision: ~11 ms/toy --
  // pure noise has no real minimum to converge toward, so every toy burns
  // the full iteration budget (see the comment above and
  // docs/ClusterFitMC_Design.md §3). 8000 keeps this validation app's
  // total runtime to a few minutes; the NoiseTailCalibratorConfig
  // PRODUCTION default is n_toys=100000 (~100000*11ms =~ 18 minutes,
  // expected and accepted as a one-time per-detector-config cost, similar
  // in spirit to EfficiencyMC's own ne_trials MC loops) -- pass a larger
  // --n_toys_full to actually exercise that scale here.
  int n_toys_full = 8000;  // toys of the largest run
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--n_toys_full" && i + 1 < argc) { n_toys_full = std::stoi(argv[++i]); continue; }
  }
  std::printf("(NoiseTailCalibrator fit precision is deliberately NOT reduced for speed -- "
              "see source comments. n_toys_full=%d; pass --n_toys_full 100000 to validate at "
              "the production default, ~18 min.)\n\n", n_toys_full);

  bool ok = true;  // stays true while every check passes

  // --- 1. Convergence: does the calibrated cut stabilize as n_toys grows? ---
  std::printf("=== Convergence across n_toys (fixed seed) ===\n");
  std::vector<int> n_toys_sweep = {500, 2000, 8000};
  // Only append n_toys_full if it's a genuinely new point -- otherwise the
  // "last two points" convergence check below would trivially compare a
  // value against itself.
  if (n_toys_full != n_toys_sweep.back()) n_toys_sweep.push_back(n_toys_full);
  std::vector<double> cuts;  // calibrated cut at each n_toys of the sweep
  for (int n : n_toys_sweep) {
    auto cfg = MakeConfig(n, 20260825ULL);
    NoiseTailCalibrator calib(cfg);
    auto result = calib.Run();
    cuts.push_back(result.delta_ll_cut);
    std::printf("  n_toys=%7d  ΔLL_cut=%.4f\n", n, result.delta_ll_cut);
  }
  {
    const double last = cuts.back();
    const double prev = cuts[cuts.size() - 2];
    const double rel_change = (last != 0.0) ? std::abs(last - prev) / std::abs(last) : 0.0;
    std::printf("  relative change, last two points: %.1f%%  (threshold: 25%%)\n", rel_change * 100.0);
    if (rel_change > 0.25) { std::printf("  FAIL -- cut has not converged\n"); ok = false; }
    else std::printf("  ok\n");
  }

  // --- 2. Stability across independent seeds, at a fixed, moderate n_toys ---
  std::printf("\n=== Stability across seeds (n_toys=2000) ===\n");
  std::vector<double> seed_cuts;  // calibrated cut for each seed
  for (std::uint64_t seed : {111ULL, 222ULL, 333ULL, 444ULL, 555ULL}) {
    auto cfg = MakeConfig(2000, seed);
    NoiseTailCalibrator calib(cfg);
    auto result = calib.Run();
    seed_cuts.push_back(result.delta_ll_cut);
    std::printf("  seed=%3llu  ΔLL_cut=%.4f\n", static_cast<unsigned long long>(seed), result.delta_ll_cut);
  }
  {
    double mean = 0.0;
    for (double c : seed_cuts) mean += c;
    mean /= static_cast<double>(seed_cuts.size());
    double ss = 0.0;
    for (double c : seed_cuts) ss += (c - mean) * (c - mean);
    const double sd = std::sqrt(ss / static_cast<double>(seed_cuts.size() - 1));
    const double rel_spread = (mean != 0.0) ? std::abs(sd / mean) : 0.0;
    std::printf("  mean=%.4f  sd=%.4f  relative spread=%.1f%%  (threshold: 30%%)\n",
                mean, sd, rel_spread * 100.0);
    if (rel_spread > 0.30) { std::printf("  FAIL -- calibration too seed-sensitive at this n_toys\n"); ok = false; }
    else std::printf("  ok\n");
  }

  // --- 3. Shape summary of the full ΔLL distribution (Fig. 6 analog) ---
  std::printf("\n=== ΔLL distribution shape (n_toys=%d) ===\n", n_toys_full);
  {
    auto cfg = MakeConfig(n_toys_full, 20260825ULL);
    NoiseTailCalibrator calib(cfg);
    auto result = calib.Run();
    const auto& d = result.delta_ll_toys;  // already sorted ascending by NoiseTailCalibrator
    auto pct = [&](double p) { return d[static_cast<std::size_t>(p * (d.size() - 1))]; };
    std::printf("  min (best-looking noise) = %.4f\n", d.front());
    std::printf("  0.1%% percentile          = %.4f  (this is ΔLL_cut)\n", pct(0.001));
    std::printf("  1%%   percentile          = %.4f\n", pct(0.01));
    std::printf("  5%%   percentile          = %.4f\n", pct(0.05));
    std::printf("  50%%  percentile (median) = %.4f\n", pct(0.50));
    std::printf("  max (worst-looking noise) = %.4f\n", d.back());
    // Sanity: the distribution should be monotonically non-decreasing (it's
    // sorted) and one-sided (all values <= 0, since ΔLL is <= 0 by
    // construction) -- a basic structural check, not a shape-matching one.
    bool monotone = true, one_sided = true;
    for (std::size_t i = 1; i < d.size(); ++i) if (d[i] < d[i - 1]) monotone = false;
    for (double v : d) if (v > 1e-9) one_sided = false;
    std::printf("  monotone sorted: %s   one-sided (<=0): %s\n",
                monotone ? "ok" : "FAIL", one_sided ? "ok" : "FAIL");
    if (!monotone || !one_sided) ok = false;
  }

  std::printf("\n%s\n", ok ? "ALL CHECKS PASS" : "SOME CHECKS FAILED");
  return ok ? 0 : 1;
}
