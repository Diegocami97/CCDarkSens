// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ClusterFitEngine.hh -- 3-parameter (mux, muy, sigma_xy) search that
//  minimizes ClusterFitModel's ObjectiveG, recovering the best-fit pixel-
//  cloud position/width (and, via the closed-form IHat, amplitude) and the
//  resulting ΔLL test statistic. Two interchangeable methods: a dependency-
//  free Nelder-Mead (default) and an optional ROOT::Math::Minuit2 path,
//  selectable per call.
// ===========================================================================

#pragma once

#include <vector>
#include "ccdarksens/response/ClusterFitModel.hh"

namespace ccdarksens {

// ----------------------------------------------------------------------------
// ClusterFitConfig
//   Settings of the cluster fit: assumed noise, sigma bounds, minimizer and its limits.
// ----------------------------------------------------------------------------
struct ClusterFitConfig {
  /// Per-pixel Gaussian readout-noise sigma, electrons. Must match the
  /// sigma_pix_e the pixel window under fit was actually generated/read
  /// with (PixelSimulatorConfig::sigma_readout_e in the caller).
  double sigma_pix_e = 0.16;

  /// Search bounds on sigma_xy, in pixel units. Values tried outside this
  /// range are clamped for evaluation purposes (soft bound, not a hard
  /// error) so the search cannot silently wander to a degenerate width.
  double sigma_xy_lo_px = 0.1;
  double sigma_xy_hi_px = 2.0;

  enum class Method { kNelderMead, kMinuit2 };  // minimizer choice
  Method method = Method::kNelderMead;  // minimizer actually used

  /// When true, fit a 1D Gaussian (mux, sigma only -- muy fixed, unused) over
  /// a single collapsed row instead of the normal 2D (mux, muy, sigma) fit --
  /// the DAMIC 1x100 readout mode (PhysRevD.94.082006 Sec. VI). The window
  /// passed to Fit() must have ny == 1 in this mode (see
  /// PixelSimulatorConfig::collapse_y, which produces such a window).
  /// Only Method::kNelderMead is implemented for one_dimensional; kMinuit2
  /// falls back to it (matches the existing Minuit2-unavailable fallback
  /// pattern below).
  bool one_dimensional = false;

  // Nelder-Mead iteration cap (per start; kMinuit2 ignores this). 300 is
  // the accuracy-first default for fitting real (or synthetic-but-precise)
  // events -- see docs/ClusterFitMC_Design.md §3 for why pure-noise toy
  // ensembles (NoiseTailCalibrator) deliberately override this lower: a
  // flat noise landscape has no real minimum to converge toward, so it
  // burns the full iteration budget on nearly every toy, and per-toy
  // precision matters far less there than it does for a real candidate's
  // reported energy (the calibration is a statistical quantile over many
  // toys, not a single precise measurement). Do not lower this default
  // without re-running ccdarksens_validate_cluster_fit_engine -- 60 is
  // measurably unsafe on real signal (bias significance > 15 sigma on the
  // near-fiducial-ceiling test case).
  int max_iterations = 300;
  double tol = 1e-8;         // convergence tolerance (spread of objective across the simplex)
};

// ----------------------------------------------------------------------------
// ClusterFitResult
//   Best-fit parameters of one window and the resulting test statistic.
// ----------------------------------------------------------------------------
struct ClusterFitResult {
  bool ok = false;  // true if the fit ran
  double mux_px = 0.0;  // fitted centre, x [pixels]
  double muy_px = 0.0;  // fitted centre, y [pixels] (0 in the 1D mode)
  double sigma_xy_px = 0.0;  // fitted Gaussian width [pixels]
  double I_hat_e = 0.0;  // fitted total charge [e-]
  double delta_ll = 0.0;  // <= 0; more negative = more signal-like
  int n_eval = 0;         // objective-function calls, for perf diagnostics
};

// ----------------------------------------------------------------------------
// ClusterFitEngine
//   I fit an isotropic 2D pixel-integrated Gaussian (or its 1D collapsed-row
//   version) to a pixel window by minimizing ObjectiveG, and return the best
//   parameters and Delta LL. The minimizer is my own Nelder-Mead (default) or
//   ROOT's Minuit2 when it was compiled in.
// ----------------------------------------------------------------------------
class ClusterFitEngine {
 public:
  // Constructor: keep the fit settings.
  explicit ClusterFitEngine(const ClusterFitConfig& cfg);

  /**
   * Fit a pixel window. pixels is a flat, row-major nx*ny array (pixel
   * (ix,iy) at index iy*nx+ix), same layout PixelSimulator::PixelCharges()
   * returns.
   */
  ClusterFitResult Fit(const std::vector<double>& pixels, int nx, int ny) const;

 private:
  ClusterFitConfig cfg_;  // fit settings

  /// Cheap starting guess: brightness-weighted centroid for (mux,muy), then
  /// a short grid scan over sigma at that centroid. Falls back to the
  /// window center / bound midpoint if the window has no positive signal
  /// (e.g. a pure-noise toy with net-negative fluctuation).
  void Seed_(const FitWindow& w, double& mux0, double& muy0, double& sigma0) const;
  /// 1D analog: weighted centroid over x only, sigma grid scan using
  /// ObjectiveG1D.
  void Seed1D_(const FitWindow& w, double& mux0, double& sigma0) const;

  // 2D Nelder-Mead fit from the given seed (multi-start, sigmoid-bounded width).
  ClusterFitResult FitNelderMead_(const FitWindow& w, double mux0, double muy0, double sigma0) const;
  ClusterFitResult FitMinuit2_(const FitWindow& w, double mux0, double muy0, double sigma0) const;
  ClusterFitResult FitNelderMead1D_(const FitWindow& w, double mux0, double sigma0) const;
};

}  // namespace ccdarksens
