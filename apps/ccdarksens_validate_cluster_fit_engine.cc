// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: ccdarksens_validate_cluster_fit_engine.cc
//  Slice 2 validation: inject known (mux,muy,sigma_xy,I) onto synthetic
//  noisy pixel windows, fit with ClusterFitEngine, and check (a) recovery
//  is unbiased across repeated noise realizations at several SNR/width
//  points, and (b) the homemade Nelder-Mead and Minuit2 paths agree with
//  each other on the same windows.
//
//  No external reference dataset exists for this (unlike the Phase 1/2
//  parity checks against WIMPyCCD/Chavarria) -- this validates the fit
//  against known ground truth it was given itself, per the strategy in
//  docs/ClusterFitMC_Design.md.
//
//  Usage:
//    ccdarksens_validate_cluster_fit_engine [--fit_method nelder_mead|minuit2|both] [--n_repeat N]
//
//  Default: both methods, cross-checked against each other, N=200 repeats
//  per test point.
// ===========================================================================

#include "ccdarksens/response/ClusterFitEngine.hh"
#include "ccdarksens/response/ClusterFitModel.hh"

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace ccdarksens;

namespace {

// One 2D injection point: true centre, width [pixels] and total charge [e-].
struct TestPoint {
  const char* label;
  double mux_px, muy_px, sigma_xy_px, I_e;
};

// Fit minus truth for every repeat of one 2D test point.
struct Residuals {
  std::vector<double> d_mux, d_muy, d_sigma, d_I;
};

// mean and standard error of the mean (std/sqrt(N)) for a residual sample.
// ----------------------------------------------------------------------------
// MeanAndSem
//   Mean of the residuals and its standard error, std/sqrt(N).
// ----------------------------------------------------------------------------
void MeanAndSem(const std::vector<double>& v, double& mean, double& sem) {
  const double n = static_cast<double>(v.size());
  double sum = 0.0;
  for (double x : v) sum += x;
  mean = sum / n;
  double ss = 0.0;
  for (double x : v) ss += (x - mean) * (x - mean);
  const double sd = std::sqrt(ss / std::max(1.0, n - 1.0));
  sem = sd / std::sqrt(n);
}

// ----------------------------------------------------------------------------
// MakeNoisyWindow
//   nx x ny window with truth I*Shape(ix, iy) plus Gaussian noise of width sigma_pix_e.
// ----------------------------------------------------------------------------
std::vector<double> MakeNoisyWindow(int nx, int ny, double mux, double muy, double sigma,
                                     double I, double sigma_pix_e, std::mt19937_64& rng) {
  std::normal_distribution<double> gaus(0.0, 1.0);
  std::vector<double> pixels(static_cast<std::size_t>(nx * ny));
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {
      const double truth = I * Shape(ix, iy, mux, muy, sigma);
      pixels[static_cast<std::size_t>(iy * nx + ix)] = truth + sigma_pix_e * gaus(rng);
    }
  }
  return pixels;
}

// MakeNoisyWindow1D: the 1D window for the collapsed-row mode (I*Shape1D(ix) plus Gaussian noise).
// 1x100-mode analog: truth from Shape1D (the y-marginal), one collapsed row.
std::vector<double> MakeNoisyWindow1D(int nx, double mux, double sigma, double I,
                                       double sigma_pix_e, std::mt19937_64& rng) {
  std::normal_distribution<double> gaus(0.0, 1.0);
  std::vector<double> pixels(static_cast<std::size_t>(nx));
  for (int ix = 0; ix < nx; ++ix) {
    const double truth = I * Shape1D(ix, mux, sigma);
    pixels[static_cast<std::size_t>(ix)] = truth + sigma_pix_e * gaus(rng);
  }
  return pixels;
}

// One 1D injection point: true centre, width [pixels] and total charge [e-].
struct TestPoint1D {
  const char* label;
  double mux_px, sigma_xy_px, I_e;
};

// Fit minus truth for every repeat of one 1D test point.
struct Residuals1D {
  std::vector<double> d_mux, d_sigma, d_I;
};

}  // namespace

// ----------------------------------------------------------------------------
// main
//   Validation of ClusterFitEngine against known truth: I inject noisy windows with a
//   known centre, width and charge at five 2D test points, refit them with Nelder-Mead
//   and/or Minuit2, and require each parameter's mean bias to be below 5 standard errors.
//   When both minimizers run, they must agree to 0.02 px in mu and sigma and 1 e- in charge.
//   The same recovery test is repeated for the 1D (1x100) mode. Options: --fit_method
//   nelder_mead|minuit2|both and --n_repeat N. Returns 0 if all checks pass.
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
  std::string fit_method_arg = "both";
  int n_repeat = 200;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--fit_method" && i + 1 < argc) { fit_method_arg = argv[++i]; continue; }
    if (arg == "--n_repeat" && i + 1 < argc) { n_repeat = std::stoi(argv[++i]); continue; }
  }
  const bool run_nm = (fit_method_arg == "nelder_mead" || fit_method_arg == "both");
  const bool run_mn = (fit_method_arg == "minuit2" || fit_method_arg == "both");
  if (!run_nm && !run_mn) {
    std::fprintf(stderr, "unknown --fit_method %s (expected nelder_mead|minuit2|both)\n",
                 fit_method_arg.c_str());
    return 2;
  }

  const int nx = 15, ny = 15;  // fit window [pixels]
  const double sigma_pix_e = 0.16;

  std::vector<TestPoint> points = {
      {"bright, mid-width",      7.30, 6.80, 1.00, 500.0},
      {"medium, mid-width",      7.30, 6.80, 1.00, 100.0},
      {"faint, mid-width",       7.30, 6.80, 1.00,  20.0},
      {"medium, narrow (near fiducial floor)", 7.30, 6.80, 0.40, 100.0},
      {"medium, wide (near fiducial ceiling)", 7.30, 6.80, 1.80, 100.0},
  };

  std::mt19937_64 rng(20260825ULL);  // fixed seed -- reproducible regression baseline

  ClusterFitConfig cfg_nm; cfg_nm.sigma_pix_e = sigma_pix_e; cfg_nm.method = ClusterFitConfig::Method::kNelderMead;
  ClusterFitConfig cfg_mn; cfg_mn.sigma_pix_e = sigma_pix_e; cfg_mn.method = ClusterFitConfig::Method::kMinuit2;
  ClusterFitEngine engine_nm(cfg_nm);
  ClusterFitEngine engine_mn(cfg_mn);

  bool ok = true;  // stays true while every check passes
  constexpr double kBiasSemThreshold = 5.0;      // |mean bias| / SEM
  constexpr double kCrossCheckTolMu = 0.02;      // px
  constexpr double kCrossCheckTolSigma = 0.02;   // px
  constexpr double kCrossCheckTolI = 1.0;        // electrons
  double max_cross_mu = 0.0, max_cross_sigma = 0.0, max_cross_I = 0.0;

  for (const auto& tp : points) {
    Residuals res_nm, res_mn;
    for (int rep = 0; rep < n_repeat; ++rep) {
      auto pixels = MakeNoisyWindow(nx, ny, tp.mux_px, tp.muy_px, tp.sigma_xy_px, tp.I_e,
                                     sigma_pix_e, rng);
      ClusterFitResult r_nm, r_mn;
      if (run_nm) r_nm = engine_nm.Fit(pixels, nx, ny);
      if (run_mn) r_mn = engine_mn.Fit(pixels, nx, ny);

      if (run_nm) {
        res_nm.d_mux.push_back(r_nm.mux_px - tp.mux_px);
        res_nm.d_muy.push_back(r_nm.muy_px - tp.muy_px);
        res_nm.d_sigma.push_back(r_nm.sigma_xy_px - tp.sigma_xy_px);
        res_nm.d_I.push_back(r_nm.I_hat_e - tp.I_e);
      }
      if (run_mn) {
        res_mn.d_mux.push_back(r_mn.mux_px - tp.mux_px);
        res_mn.d_muy.push_back(r_mn.muy_px - tp.muy_px);
        res_mn.d_sigma.push_back(r_mn.sigma_xy_px - tp.sigma_xy_px);
        res_mn.d_I.push_back(r_mn.I_hat_e - tp.I_e);
      }
      if (run_nm && run_mn) {
        max_cross_mu = std::max({max_cross_mu, std::abs(r_nm.mux_px - r_mn.mux_px),
                                  std::abs(r_nm.muy_px - r_mn.muy_px)});
        max_cross_sigma = std::max(max_cross_sigma, std::abs(r_nm.sigma_xy_px - r_mn.sigma_xy_px));
        max_cross_I = std::max(max_cross_I, std::abs(r_nm.I_hat_e - r_mn.I_hat_e));
      }
    }

    auto report = [&](const char* method_name, const Residuals& res) {
      double mean, sem;
      bool pt_ok = true;
      std::printf("  [%s] %s:\n", method_name, tp.label);
      const char* names[4] = {"mux", "muy", "sigma", "I"};
      const std::vector<double>* vecs[4] = {&res.d_mux, &res.d_muy, &res.d_sigma, &res.d_I};
      for (int k = 0; k < 4; ++k) {
        MeanAndSem(*vecs[k], mean, sem);
        const double z = (sem > 0.0) ? std::abs(mean) / sem : 0.0;
        std::printf("    %-6s bias=%+.4f  sem=%.4f  |bias|/sem=%.2f  %s\n",
                    names[k], mean, sem, z, (z < kBiasSemThreshold) ? "ok" : "FAIL");
        if (z >= kBiasSemThreshold) pt_ok = false;
      }
      return pt_ok;
    };

    if (run_nm) ok &= report("nelder_mead", res_nm);
    if (run_mn) ok &= report("minuit2    ", res_mn);
  }

  if (run_nm && run_mn) {
    std::printf("\n[cross-check] max |NM - Minuit2| over all trials: "
                "mu=%.4f px  sigma=%.4f px  I=%.4f e-\n",
                max_cross_mu, max_cross_sigma, max_cross_I);
    if (max_cross_mu > kCrossCheckTolMu || max_cross_sigma > kCrossCheckTolSigma ||
        max_cross_I > kCrossCheckTolI) {
      std::printf("[cross-check] FAIL (tolerance: mu=%.3f sigma=%.3f I=%.3f)\n",
                  kCrossCheckTolMu, kCrossCheckTolSigma, kCrossCheckTolI);
      ok = false;
    } else {
      std::printf("[cross-check] PASS\n");
    }
  }

  // ---- 1D (DAMIC 1x100 mode) recovery check ----
  // Same idea as the 2D loop above, but on a single collapsed row using
  // Shape1D truth and ClusterFitConfig::one_dimensional=true -- validates
  // Slices 1-3 of the joint-likelihood plan (ClusterFitModel's 1D
  // primitives, PixelSimulator's collapse_y mode is exercised separately by
  // ccdarksens_validate_wimp_nucleon_paper_repro against a real 1x100
  // config; this checks the fit engine itself against known truth).
  {
    const int nx1d = 41;
    ClusterFitConfig cfg_1d = cfg_nm;
    cfg_1d.one_dimensional = true;
    ClusterFitEngine engine_1d(cfg_1d);

    std::vector<TestPoint1D> points_1d = {
        {"1x100 bright, mid-width", 20.3, 1.00, 500.0},
        {"1x100 medium, mid-width", 20.3, 1.00, 100.0},
        {"1x100 faint, mid-width",  20.3, 1.00,  20.0},
        {"1x100 medium, narrow",    20.3, 0.40, 100.0},
        {"1x100 medium, wide",      20.3, 1.80, 100.0},
    };

    for (const auto& tp : points_1d) {
      Residuals1D res;
      for (int rep = 0; rep < n_repeat; ++rep) {
        auto pixels = MakeNoisyWindow1D(nx1d, tp.mux_px, tp.sigma_xy_px, tp.I_e, sigma_pix_e, rng);
        ClusterFitResult r = engine_1d.Fit(pixels, nx1d, 1);
        res.d_mux.push_back(r.mux_px - tp.mux_px);
        res.d_sigma.push_back(r.sigma_xy_px - tp.sigma_xy_px);
        res.d_I.push_back(r.I_hat_e - tp.I_e);
      }
      double mean, sem;
      bool pt_ok = true;
      std::printf("  [1D nelder_mead] %s:\n", tp.label);
      const char* names[3] = {"mux", "sigma", "I"};
      const std::vector<double>* vecs[3] = {&res.d_mux, &res.d_sigma, &res.d_I};
      for (int k = 0; k < 3; ++k) {
        MeanAndSem(*vecs[k], mean, sem);
        const double z = (sem > 0.0) ? std::abs(mean) / sem : 0.0;
        std::printf("    %-6s bias=%+.4f  sem=%.4f  |bias|/sem=%.2f  %s\n",
                    names[k], mean, sem, z, (z < kBiasSemThreshold) ? "ok" : "FAIL");
        if (z >= kBiasSemThreshold) pt_ok = false;
      }
      ok &= pt_ok;
    }
  }

  std::printf("\n%s\n", ok ? "ALL CHECKS PASS" : "SOME CHECKS FAILED");
  return ok ? 0 : 1;
}
