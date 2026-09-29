// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: ccdarksens_build_cluster_fit_kernel.cc
//  Slice 4 validation: calibrates a noise-tail ΔLL cut (Slice 3), forward-
//  simulates real nuclear-recoil-like events across a grid of true energies
//  (Slice 4), and reports the resulting kernel K[E_true,E_reco] and
//  detection-efficiency curve -- the Fig. 9/11 analog. No exact external
//  reference exists for this specific configuration (detector params here
//  are not tuned to reproduce the literal PhysRevD.94.082006 dataset) --
//  this checks physical sanity/shape, not numerical parity. See
//  docs/ClusterFitMC_Design.md.
//
//  Usage:
//    ccdarksens_build_cluster_fit_kernel [--n_toys_calib N] [--ne_trials N]
// ===========================================================================

#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/ClusterFitMC.hh"
#include "ccdarksens/response/NoiseTailCalibrator.hh"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace ccdarksens;

namespace {

// ----------------------------------------------------------------------------
// LogSpace
//   n log-spaced values from lo to hi inclusive (needs n >= 2 and lo, hi > 0).
// ----------------------------------------------------------------------------
std::vector<double> LogSpace(double lo, double hi, int n) {
  std::vector<double> v(static_cast<std::size_t>(n));
  const double log_lo = std::log(lo), log_hi = std::log(hi);
  for (int i = 0; i < n; ++i) {
    const double t = static_cast<double>(i) / (n - 1);
    v[static_cast<std::size_t>(i)] = std::exp(log_lo + t * (log_hi - log_lo));
  }
  return v;
}

}  // namespace

// ----------------------------------------------------------------------------
// main
//   Slice-4 validation of the cluster-fit kernel:
//     1) calibrate the Delta LL noise-tail cut with pure-noise toys;
//     2) forward-simulate events on a log-spaced grid of true energies and build
//        K[E_true][E_reco] and the efficiency curve;
//     3) sanity checks: efficiencies lie in [0,1], the efficiency rises from the
//        lowest to the highest energy, and every K row sums to its efficiency.
//   Options: --n_toys_calib N and --ne_trials N. Exit code 0 if all checks pass.
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
  int n_toys_calib = 5000;  // pure-noise toys for the calibration
  int ne_trials = 1500;  // simulated events per true-energy point
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--n_toys_calib" && i + 1 < argc) { n_toys_calib = std::stoi(argv[++i]); continue; }
    if (arg == "--ne_trials" && i + 1 < argc) { ne_trials = std::stoi(argv[++i]); continue; }
  }

  const int nx = 15, ny = 15;  // fit window [pixels]
  const double sigma_pix_e = 0.16;  // modern skipper-CCD projection -- see docs/ClusterFitMC_Design.md

  // --- Step 1: calibrate the detection cut (Slice 3), same window/fit
  // config that will be used for the real events below -- required for the
  // cut to describe this window's actual noise behavior. ---
  NoiseTailCalibratorConfig calib_cfg;
  calib_cfg.pix_cfg.nx = nx; calib_cfg.pix_cfg.ny = ny;
  calib_cfg.pix_cfg.sigma_readout_e = sigma_pix_e;
  calib_cfg.pix_cfg.rng_seed = 20260825ULL;
  calib_cfg.fit_cfg.sigma_pix_e = sigma_pix_e;
  calib_cfg.fit_cfg.method = ClusterFitConfig::Method::kNelderMead;  // full precision, see Slice 3 finding
  calib_cfg.n_toys = n_toys_calib;
  calib_cfg.target_tail_prob = 1e-3;

  std::printf("=== Step 1: noise-tail calibration (n_toys=%d) ===\n", n_toys_calib);
  NoiseTailCalibrator calibrator(calib_cfg);
  const auto calib_result = calibrator.Run();
  std::printf("  DeltaLL_cut = %.4f\n\n", calib_result.delta_ll_cut);

  // --- Step 2: build the kernel over a grid of true energies. ---
  ChargeTransportConfig ct_cfg;  // defaults match the rest of the pipeline (A_um2=803.25, b_umInv=6.5e-4, ...)
  ct_cfg.rng_seed = 13579ULL;
  auto ct = std::make_shared<ChargeTransport>(ct_cfg);

  ClusterFitMCConfig cfg;
  cfg.pix_cfg = calib_cfg.pix_cfg;   // MUST match the calibration
  cfg.fit_cfg = calib_cfg.fit_cfg;   // MUST match the calibration
  cfg.delta_ll_cut = calib_result.delta_ll_cut;
  cfg.sigma_xy_fid_min_px = 0.35;
  cfg.sigma_xy_fid_max_px = 1.22;
  cfg.eh_pair_eV = 3.77;
  cfg.fano_factor = 0.133;
  cfg.ne_trials_per_point = ne_trials;
  cfg.rng_seed = 24680ULL;

  ClusterFitMC cluster_mc(cfg, ct);

  const auto Etrue_grid = LogSpace(10.0, 300.0, 12);
  std::vector<double> Ereco_edges;  // 0..400 eV in 10 eV bins -- comfortably above Etrue's max (300 eV)
  for (int i = 0; i <= 40; ++i) Ereco_edges.push_back(i * 10.0);

  std::printf("=== Step 2: building kernel (%d E_true points x %d trials) ===\n",
              static_cast<int>(Etrue_grid.size()), ne_trials);
  const auto kernel = cluster_mc.BuildKernel(Etrue_grid, Ereco_edges);

  std::printf("  %10s  %12s\n", "E_true[eV]", "efficiency");
  for (std::size_t i = 0; i < kernel.Etrue_grid_eV.size(); ++i) {
    std::printf("  %10.2f  %12.4f\n", kernel.Etrue_grid_eV[i], kernel.efficiency[i]);
  }

  // --- Sanity checks ---
  bool ok = true;  // stays true while every check passes
  std::printf("\n=== Sanity checks ===\n");

  bool all_in_range = true;
  for (double e : kernel.efficiency) if (e < 0.0 || e > 1.0) all_in_range = false;
  std::printf("  efficiency values in [0,1]: %s\n", all_in_range ? "ok" : "FAIL");
  if (!all_in_range) ok = false;

  const double eff_low = kernel.efficiency.front();
  const double eff_high = kernel.efficiency.back();
  std::printf("  efficiency at lowest E_true (%.1f eV) = %.4f  (expect near 0)\n",
              kernel.Etrue_grid_eV.front(), eff_low);
  std::printf("  efficiency at highest E_true (%.1f eV) = %.4f  (expect meaningfully higher)\n",
              kernel.Etrue_grid_eV.back(), eff_high);
  const bool turns_on = (eff_high > eff_low + 0.1);
  std::printf("  efficiency rises from low to high E_true: %s\n", turns_on ? "ok" : "FAIL");
  if (!turns_on) ok = false;

  // K rows should sum to exactly the reported efficiency (by construction --
  // see docs/ClusterFitMC_Design.md for why these can differ if Ereco_edges
  // is too narrow; here it's chosen wide enough that they should match).
  bool k_consistent = true;
  for (std::size_t i = 0; i < kernel.K.size(); ++i) {
    double sum = 0.0;
    for (double v : kernel.K[i]) sum += v;
    if (std::abs(sum - kernel.efficiency[i]) > 1e-9) k_consistent = false;
  }
  std::printf("  K rows sum to reported efficiency: %s\n", k_consistent ? "ok" : "FAIL");
  if (!k_consistent) ok = false;

  std::printf("\n%s\n", ok ? "ALL CHECKS PASS" : "SOME CHECKS FAILED");
  return ok ? 0 : 1;
}
