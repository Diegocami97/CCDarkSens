// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: ccdarksens_example_one_point_cluster.cc
//  Slice 5 integration example: the WIMP-nucleus SI channel, cluster-fit
//  (continuous E_reco) reconstruction, wired end-to-end through the same
//  JSON config machinery every other channel uses. Builds the noise-tail
//  calibration and the ClusterFitMC kernel ONCE from
//  response.cluster_fit_mc, then folds two different rate spectra (two
//  different WIMP masses) through that SAME kernel -- demonstrating that
//  the kernel is a pure detector-response object, independent of the DM
//  physics parameters, and should not be rebuilt per grid point.
//
//  Usage:
//    ccdarksens_example_one_point_cluster [config.json]
// ===========================================================================

#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/io/RateTable.hh"
#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/ClusterEnergyRates.hh"
#include "ccdarksens/response/ClusterFitMC.hh"
#include "ccdarksens/response/NoiseTailCalibrator.hh"
#include "ccdarksens/stats/ProfileLikelihood.hh"
#include <TH1D.h>

#include <cmath>
#include <cstdio>
#include <memory>
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
//   End-to-end single-point run of the WIMP-nucleon (cluster_energy) channel from a
//   config: calibrate the noise-tail cut, build the kernel once, fold two different
//   rate spectra through that same kernel, and feed each into ProfileLikelihood
//   (Asimov). Usage: <program> [config.json]. Exit code 0 if all checks pass.
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
  const std::string config_path = (argc > 1) ? argv[1] : "configs/wimp_nucleon_cluster_example.json";

  ConfigManager cfg(config_path);
  cfg.parse();
  const auto& jcfg = cfg.response().cluster_fit_mc;  // cluster-fit settings from the config
  const double exposure_kg_year = cfg.run().exposure_kg_year;  // exposure [kg*year] from the run block

  std::printf("=== Config loaded: %s ===\n", config_path.c_str());
  std::printf("  analysis_space = %s\n", cfg.response().analysis_space.c_str());
  std::printf("  exposure_kg_year = %.4f\n\n", exposure_kg_year);

  // --- Step 1: calibrate the detection cut, once, from config. ---
  NoiseTailCalibratorConfig calib_cfg;
  calib_cfg.pix_cfg.nx = jcfg.window_nx;
  calib_cfg.pix_cfg.ny = jcfg.window_ny;
  calib_cfg.pix_cfg.pixel_size_um = jcfg.pixel_size_um;
  calib_cfg.pix_cfg.sigma_readout_e = jcfg.sigma_readout_e;
  calib_cfg.pix_cfg.rng_seed = jcfg.calib_rng_seed;
  calib_cfg.fit_cfg.sigma_pix_e = jcfg.sigma_readout_e;
  calib_cfg.fit_cfg.method = (jcfg.fit_method == "minuit2")
                                 ? ClusterFitConfig::Method::kMinuit2
                                 : ClusterFitConfig::Method::kNelderMead;
  calib_cfg.fit_cfg.sigma_xy_lo_px = jcfg.sigma_xy_lo_px;
  calib_cfg.fit_cfg.sigma_xy_hi_px = jcfg.sigma_xy_hi_px;
  calib_cfg.n_toys = jcfg.n_toys;
  calib_cfg.target_tail_prob = jcfg.target_tail_prob;

  std::printf("=== Step 1: noise-tail calibration (n_toys=%d) ===\n", jcfg.n_toys);
  NoiseTailCalibrator calibrator(calib_cfg);
  const auto calib_result = calibrator.Run();
  std::printf("  DeltaLL_cut = %.4f\n\n", calib_result.delta_ll_cut);

  // --- Step 2: build the kernel ONCE. It depends only on the detector
  // config above, not on any WIMP mass/cross-section -- reused for every
  // rate spectrum folded below. ---
  ChargeTransportConfig ct_cfg;
  ct_cfg.thickness_um = jcfg.thickness_um;
  ct_cfg.A_um2 = jcfg.A_um2;
  ct_cfg.b_umInv = jcfg.b_umInv;
  ct_cfg.alpha = jcfg.alpha;
  ct_cfg.beta_per_keV = jcfg.beta_per_keV;
  ct_cfg.rng_seed = jcfg.rng_seed;
  auto ct = std::make_shared<ChargeTransport>(ct_cfg);

  ClusterFitMCConfig mc_cfg;
  mc_cfg.pix_cfg = calib_cfg.pix_cfg;
  mc_cfg.fit_cfg = calib_cfg.fit_cfg;
  mc_cfg.delta_ll_cut = calib_result.delta_ll_cut;
  mc_cfg.sigma_xy_fid_min_px = jcfg.sigma_xy_fid_min_px;
  mc_cfg.sigma_xy_fid_max_px = jcfg.sigma_xy_fid_max_px;
  mc_cfg.eh_pair_eV = jcfg.eh_pair_eV;
  mc_cfg.fano_factor = jcfg.fano_factor;
  mc_cfg.ne_trials_per_point = jcfg.ne_trials_per_point;
  mc_cfg.rng_seed = jcfg.rng_seed;

  ClusterFitMC cluster_mc(mc_cfg, ct);

  const auto Etrue_grid = LogSpace(jcfg.Etrue_min_eV, jcfg.Etrue_max_eV, jcfg.Etrue_npoints);
  std::vector<double> Ereco_edges;
  for (int i = 0; i <= jcfg.Ereco_nbins; ++i) {
    Ereco_edges.push_back(jcfg.Ereco_min_eV +
        i * (jcfg.Ereco_max_eV - jcfg.Ereco_min_eV) / jcfg.Ereco_nbins);
  }

  std::printf("=== Step 2: building kernel ONCE (%d E_true points x %d trials) ===\n",
              jcfg.Etrue_npoints, jcfg.ne_trials_per_point);
  const auto kernel = cluster_mc.BuildKernel(Etrue_grid, Ereco_edges);
  std::printf("  kernel built (%d E_true points, %d E_reco bins)\n\n",
              static_cast<int>(kernel.Etrue_grid_eV.size()),
              static_cast<int>(kernel.Ereco_edges_eV.size()) - 1);

  // --- Step 3: fold two DIFFERENT rate spectra through the SAME kernel. ---
  struct Point { const char* label; std::string csv; };  // one rate table to fold
  const std::vector<Point> points = {
      {"m_chi=3000 MeV", "data/wimp_nucleon_rates_ee/Si/heavy/dRdE_ee_Si28_heavy_m3000.000000_s1.0e-38.csv"},
      {"m_chi=5000 MeV", "data/wimp_nucleon_rates_ee/Si/heavy/dRdE_ee_Si28_heavy_m5000.000000_s1.0e-38.csv"},
  };

  bool ok = true;  // stays true while every check passes
  for (const auto& pt : points) {
    std::printf("=== Step 3: folding %s (%s) ===\n", pt.label, pt.csv.c_str());
    RateTable table;
    if (!table.LoadCSV(pt.csv)) {
      std::fprintf(stderr, "  FAIL: could not load %s\n", pt.csv.c_str());
      ok = false;
      continue;
    }
    const auto h_Etrue = table.MakeTH1D("h_Etrue", jcfg.Etrue_min_eV, jcfg.Etrue_max_eV, 400);
    const auto S_reco = FoldEtrueToErecoRates(*h_Etrue, exposure_kg_year, kernel);

    double total = 0.0;
    bool all_finite_nonneg = true;
    for (double v : S_reco) {
      total += v;
      if (!std::isfinite(v) || v < 0.0) all_finite_nonneg = false;
    }
    std::printf("  total expected counts (S_reco, exposure=%.3f kg-yr) = %.6g\n", exposure_kg_year, total);
    std::printf("  all bins finite and non-negative: %s\n", all_finite_nonneg ? "ok" : "FAIL");
    if (!all_finite_nonneg) ok = false;

    // Feed into the existing, UNCHANGED ProfileLikelihood -- same call
    // pattern as ccdarksens_example_one_point_pattern.cc's Asimov path.
    ccdarksens::stats::ProfileLikelihood pl;
    pl.SetData(S_reco);
    pl.SetBTemplate(S_reco);
    std::vector<double> S_null(S_reco.size(), 0.0);
    const double q = pl.EvaluateRatio(S_null, S_reco, 0.01, 10.0);
    std::printf("  ProfileLikelihood q (S vs null, Asimov) = %.6g  %s\n\n",
                q, std::isfinite(q) ? "(finite, ok)" : "(FAIL: not finite)");
    if (!std::isfinite(q)) ok = false;
  }

  std::printf("%s\n", ok ? "ALL CHECKS PASS" : "SOME CHECKS FAILED");
  return ok ? 0 : 1;
}
