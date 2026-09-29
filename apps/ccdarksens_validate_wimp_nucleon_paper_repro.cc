// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: ccdarksens_validate_wimp_nucleon_paper_repro.cc
//  Phase 7 (WIMP-nucleon channel plan): config-driven reproduction check
//  against PhysRevD.94.082006 (DAMIC 0.6 kg-day, SNOLAB), Figs. 6/9/11, using
//  the paper's own detector parameters (sigma_pix=1.8 e-, 11x11 window,
//  675 um thickness) instead of the modern-projection defaults used
//  elsewhere. Prints the noise-tail calibration and the K[E_true,E_reco]
//  efficiency curve for direct comparison against the paper's Figs. 6/9.
//  See docs/ClusterFitMC_Design.md.
//
//  Usage: ccdarksens_validate_wimp_nucleon_paper_repro <config.json>
//             [--dump-toys=<path.csv>] [--dump-kernel=<path.csv>]
//
//  --dump-toys writes the full sorted per-toy DeltaLL array from Step 1
//  (one value per line) -- for plotting the noise-tail distribution
//  itself, not just its calibrated cut.
//  --dump-kernel writes the full K[E_true,E_reco] matrix from Step 2 (one
//  row per E_true grid point, one column per Ereco bin, plus the marginal
//  efficiency) -- for plotting the kernel as a heatmap.
// ===========================================================================

#include <cstdio>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/ClusterFitEngine.hh"
#include "ccdarksens/response/ClusterFitMC.hh"
#include "ccdarksens/response/NoiseTailCalibrator.hh"

using namespace ccdarksens;

namespace {

// ----------------------------------------------------------------------------
// LogSpace
//   n log-spaced values from lo to hi inclusive (a single point returns lo).
// ----------------------------------------------------------------------------
std::vector<double> LogSpace(double lo, double hi, int n) {
  std::vector<double> v(static_cast<std::size_t>(n));
  const double log_lo = std::log(lo), log_hi = std::log(hi);
  for (int i = 0; i < n; ++i) {
    const double t = (n > 1) ? static_cast<double>(i) / (n - 1) : 0.0;
    v[static_cast<std::size_t>(i)] = std::exp(log_lo + t * (log_hi - log_lo));
  }
  return v;
}

// ----------------------------------------------------------------------------
// FlagValue
//   Value of the first command-line argument that starts with prefix (e.g. "--dump-toys="), or an empty string.
// ----------------------------------------------------------------------------
std::string FlagValue(int argc, char** argv, const std::string& prefix) {
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg.rfind(prefix, 0) == 0) return arg.substr(prefix.size());
  }
  return "";
}

}  // namespace

// ----------------------------------------------------------------------------
// main
//   Config-driven reproduction check for the WIMP-nucleon channel:
//     Step 1: noise-tail calibration with the config's window, noise and toys;
//             prints the Delta LL cut next to the paper's Fig. 6 value;
//     Step 2: kernel build, printing the efficiency curve to compare with Fig. 9.
//   Options: --dump-toys=<csv> (sorted per-toy Delta LL) and --dump-kernel=<csv>
//   (K[E_true][E_reco] plus efficiency). Returns 0 on success, 1 on error.
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: %s <config.json> [--dump-toys=<path.csv>] [--dump-kernel=<path.csv>]\n", argv[0]);
    return 1;
  }
  const std::string dump_toys_path = FlagValue(argc, argv, "--dump-toys=");
  const std::string dump_kernel_path = FlagValue(argc, argv, "--dump-kernel=");
  try {
    ConfigManager cfg(argv[1]);
    cfg.parse();
    const auto& jcfg = cfg.response().cluster_fit_mc;  // cluster-fit settings from the config

    std::printf("=== Step 1: noise-tail calibration (window=%dx%d, sigma_pix=%.2f e-, n_toys=%d, target_tail_prob=%.1e) ===\n",
                jcfg.window_nx, jcfg.window_ny, jcfg.sigma_readout_e, jcfg.n_toys, jcfg.target_tail_prob);
    NoiseTailCalibratorConfig calib_cfg;
    calib_cfg.pix_cfg.nx = jcfg.window_nx;
    calib_cfg.pix_cfg.ny = jcfg.window_ny;
    calib_cfg.pix_cfg.pixel_size_um = jcfg.pixel_size_um;
    calib_cfg.pix_cfg.sigma_readout_e = jcfg.sigma_readout_e;
    calib_cfg.pix_cfg.rng_seed = jcfg.calib_rng_seed;
    calib_cfg.pix_cfg.collapse_y = jcfg.one_dimensional;
    calib_cfg.fit_cfg.sigma_pix_e = jcfg.sigma_readout_e;
    calib_cfg.fit_cfg.one_dimensional = jcfg.one_dimensional;
    calib_cfg.fit_cfg.method = (jcfg.fit_method == "minuit2") ? ClusterFitConfig::Method::kMinuit2
                                                               : ClusterFitConfig::Method::kNelderMead;
    calib_cfg.fit_cfg.sigma_xy_lo_px = jcfg.sigma_xy_lo_px;
    calib_cfg.fit_cfg.sigma_xy_hi_px = jcfg.sigma_xy_hi_px;
    calib_cfg.n_toys = jcfg.n_toys;
    calib_cfg.target_tail_prob = jcfg.target_tail_prob;

    NoiseTailCalibrator calibrator(calib_cfg);
    const auto calib_result = calibrator.Run();
    std::printf("  DeltaLL_cut = %.4f   (paper Fig. 6: -28 for 1x1, -25 for 1x100)\n\n", calib_result.delta_ll_cut);

    if (!dump_toys_path.empty()) {
      std::ofstream out(dump_toys_path);
      out << "# per-toy DeltaLL from the pure-noise calibration ensemble, sorted ascending (most signal-like first)\n";
      out << "# n_toys_used=" << calib_result.n_toys_used << " delta_ll_cut=" << calib_result.delta_ll_cut << "\n";
      out << "delta_ll\n";
      for (double v : calib_result.delta_ll_toys) out << v << "\n";
      std::printf("  wrote %s (%zu toys)\n\n", dump_toys_path.c_str(), calib_result.delta_ll_toys.size());
    }

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
    Ereco_edges.reserve(static_cast<std::size_t>(jcfg.Ereco_nbins) + 1);
    for (int i = 0; i <= jcfg.Ereco_nbins; ++i) {
      Ereco_edges.push_back(jcfg.Ereco_min_eV + i * (jcfg.Ereco_max_eV - jcfg.Ereco_min_eV) / jcfg.Ereco_nbins);
    }

    std::printf("=== Step 2: building kernel (%d E_true points x %d trials, %.0f-%.0f eV) ===\n",
                jcfg.Etrue_npoints, jcfg.ne_trials_per_point, jcfg.Etrue_min_eV, jcfg.Etrue_max_eV);
    const auto kernel = cluster_mc.BuildKernel(Etrue_grid, Ereco_edges);

    std::printf("\n  %12s  %12s   (paper Fig. 9: ~9-25%% at 60-75 eVee -> ~100%% by ~150-400 eVee -> plateaus ~75%% above ~1 keV)\n",
                "E_true[eV]", "efficiency");
    for (std::size_t i = 0; i < kernel.Etrue_grid_eV.size(); ++i) {
      std::printf("  %12.2f  %12.4f\n", kernel.Etrue_grid_eV[i], kernel.efficiency[i]);
    }

    if (!dump_kernel_path.empty()) {
      std::ofstream out(dump_kernel_path);
      out << "# K[E_true,E_reco] matrix. First row: Ereco bin edges (eV), size NEreco+1.\n";
      out << "# Each following row: Etrue_eV, then K[iEtrue][0..NEreco-1], then efficiency (row sum).\n";
      out << "Ereco_edges_eV";
      for (double e : kernel.Ereco_edges_eV) out << "," << e;
      out << "\n";
      for (std::size_t i = 0; i < kernel.Etrue_grid_eV.size(); ++i) {
        out << kernel.Etrue_grid_eV[i];
        for (double k : kernel.K[i]) out << "," << k;
        out << "," << kernel.efficiency[i] << "\n";
      }
      std::printf("\n  wrote %s (%zu Etrue rows x %zu Ereco bins)\n", dump_kernel_path.c_str(),
                  kernel.Etrue_grid_eV.size(), kernel.Ereco_edges_eV.size() - 1);
    }
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
}
