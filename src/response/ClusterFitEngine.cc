// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ClusterFitEngine.cc -- 3-parameter (mux, muy, sigma_xy) search minimizing
//  ClusterFitModel's ObjectiveG: a dependency-free Nelder-Mead default, and
//  an optional Minuit2 path mirroring src/stats/ProfileLikelihood.cc's
//  pattern.
// ===========================================================================

#include "ccdarksens/response/ClusterFitEngine.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>

#ifdef CCDARKSENS_USE_MINUIT2
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/IFunction.h>
#include <Math/Minimizer.h>
#endif

namespace ccdarksens {

namespace {

// Clamp sigma into [lo, hi].
double ClampSigma(double sigma, double lo, double hi) {
  return std::max(lo, std::min(hi, sigma));
}

// Objective evaluated with sigma soft-clamped into [lo,hi] -- the search is
// free to propose sigma outside the bound, but the objective it sees there
// is flat (evaluated at the clamped value), so there's no gradient pulling
// it further out and no discontinuity/NaN at the boundary.
// ----------------------------------------------------------------------------
// EvalObjective
//   ObjectiveG at a soft-clamped sigma (used by the Minuit2 functor).
// ----------------------------------------------------------------------------
double EvalObjective(const FitWindow& w, double mux, double muy, double sigma,
                      double sigma_lo, double sigma_hi) {
  return ObjectiveG(w, mux, muy, ClampSigma(sigma, sigma_lo, sigma_hi));
}

}  // namespace

// Constructor: keep the fit settings.
ClusterFitEngine::ClusterFitEngine(const ClusterFitConfig& cfg) : cfg_(cfg) {}

// ----------------------------------------------------------------------------
// ClusterFitEngine::Seed_
//   Starting point for the 2D fit: the charge-weighted centroid (negative
//   pixels ignored; the window centre if nothing is positive) and, at that
//   centroid, the best of ten log-spaced sigma values in the allowed range.
// ----------------------------------------------------------------------------
void ClusterFitEngine::Seed_(const FitWindow& w, double& mux0, double& muy0, double& sigma0) const {
  double sum_w = 0.0, sum_wx = 0.0, sum_wy = 0.0;
  for (int iy = 0; iy < w.ny; ++iy) {
    for (int ix = 0; ix < w.nx; ++ix) {
      const double d = (*w.pixels)[static_cast<std::size_t>(iy * w.nx + ix)];
      const double wgt = std::max(0.0, d);  // ignore negative noise fluctuations
      sum_w += wgt;
      sum_wx += wgt * ix;
      sum_wy += wgt * iy;
    }
  }
  if (sum_w > 0.0) {
    mux0 = sum_wx / sum_w;
    muy0 = sum_wy / sum_w;
  } else {
    mux0 = 0.5 * (w.nx - 1);
    muy0 = 0.5 * (w.ny - 1);
  }

  // Short grid scan over sigma at the centroid to pick a sensible starting width.
  double best_g = 0.0;
  double best_sigma = 0.5 * (cfg_.sigma_xy_lo_px + cfg_.sigma_xy_hi_px);
  constexpr int kScanPoints = 10;
  for (int i = 0; i < kScanPoints; ++i) {
    const double t = static_cast<double>(i) / (kScanPoints - 1);
    // log-spaced scan across the configured sigma range
    const double sigma = cfg_.sigma_xy_lo_px *
        std::pow(cfg_.sigma_xy_hi_px / cfg_.sigma_xy_lo_px, t);
    const double g = ObjectiveG(w, mux0, muy0, sigma);
    if (g < best_g) {
      best_g = g;
      best_sigma = sigma;
    }
  }
  sigma0 = best_sigma;
}

// ----------------------------------------------------------------------------
// ClusterFitEngine::Seed1D_
//   1D version of Seed_: x centroid plus a ten-point sigma scan with ObjectiveG1D.
// ----------------------------------------------------------------------------
void ClusterFitEngine::Seed1D_(const FitWindow& w, double& mux0, double& sigma0) const {
  double sum_w = 0.0, sum_wx = 0.0;
  for (int ix = 0; ix < w.nx; ++ix) {
    const double d = (*w.pixels)[static_cast<std::size_t>(ix)];
    const double wgt = std::max(0.0, d);
    sum_w += wgt;
    sum_wx += wgt * ix;
  }
  mux0 = (sum_w > 0.0) ? (sum_wx / sum_w) : 0.5 * (w.nx - 1);

  double best_g = 0.0;
  double best_sigma = 0.5 * (cfg_.sigma_xy_lo_px + cfg_.sigma_xy_hi_px);
  constexpr int kScanPoints = 10;
  for (int i = 0; i < kScanPoints; ++i) {
    const double t = static_cast<double>(i) / (kScanPoints - 1);
    const double sigma = cfg_.sigma_xy_lo_px *
        std::pow(cfg_.sigma_xy_hi_px / cfg_.sigma_xy_lo_px, t);
    const double g = ObjectiveG1D(w, mux0, sigma);
    if (g < best_g) {
      best_g = g;
      best_sigma = sigma;
    }
  }
  sigma0 = best_sigma;
}

namespace {

// Parameter vector of the 2D fit.
using Point3 = std::array<double, 3>;  // {mux, muy, u} -- u is an UNCONSTRAINED sigma surrogate

// One downhill-simplex run from a single starting point. Returns the best
// vertex found, its objective value, and the evaluation count.
// Result of one Nelder-Mead run.
struct SimplexRun {
  Point3 best;
  double f_best;
  int n_eval;
};

// ----------------------------------------------------------------------------
// RunSimplex
//   One downhill-simplex (Nelder-Mead) run in {mux, muy, u}, where u is an
//   unconstrained stand-in for sigma. It stops when both the objective spread
//   and the physical extent of the simplex are small, or after max_iterations.
//   Returns the best vertex, its objective value and the number of evaluations.
// ----------------------------------------------------------------------------
SimplexRun RunSimplex(const std::function<double(const Point3&)>& eval,
                       double mux0, double muy0, double u0,
                       int max_iterations, double tol) {
  // Standard Nelder-Mead coefficients.
  constexpr double kAlpha = 1.0;   // reflection
  constexpr double kGamma = 2.0;   // expansion
  constexpr double kRho = 0.5;     // contraction
  constexpr double kSigmaShrink = 0.5;  // shrink

  std::array<Point3, 4> simplex = {{
      {mux0, muy0, u0},
      {mux0 + 0.6, muy0, u0},
      {mux0, muy0 + 0.6, u0},
      {mux0, muy0, u0 + 0.5},
  }};
  std::array<double, 4> fval;
  int n_eval = 0;
  for (int i = 0; i < 4; ++i) { fval[i] = eval(simplex[i]); ++n_eval; }

  auto sort_simplex = [&]() {
    std::array<int, 4> idx = {0, 1, 2, 3};
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return fval[a] < fval[b]; });
    std::array<Point3, 4> s2;
    std::array<double, 4> f2;
    for (int i = 0; i < 4; ++i) { s2[i] = simplex[idx[i]]; f2[i] = fval[idx[i]]; }
    simplex = s2;
    fval = f2;
  };

  for (int iter = 0; iter < max_iterations; ++iter) {
    sort_simplex();
    // Convergence: spread of BOTH the objective values AND the simplex's
    // own physical extent (mux,muy,u) must be tiny. Objective-only spread
    // is not sufficient near the sigmoid's saturating tails (u -> +-inf,
    // sigma -> lo/hi): the u-space gradient vanishes there even when the
    // true sigma-space objective has not actually converged, which let the
    // simplex falsely "converge" pinned near a bound -- caught by the
    // Nelder-Mead-vs-Minuit2 cross-check (see docs/ClusterFitMC_Design.md).
    double max_extent = 0.0;
    for (int k = 0; k < 3; ++k) {
      double lo_k = simplex[0][k], hi_k = simplex[0][k];
      for (int i = 1; i < 4; ++i) { lo_k = std::min(lo_k, simplex[i][k]); hi_k = std::max(hi_k, simplex[i][k]); }
      max_extent = std::max(max_extent, hi_k - lo_k);
    }
    if (std::abs(fval[3] - fval[0]) < tol && max_extent < 1e-5) break;

    Point3 centroid = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i)  // all but the worst (index 3)
      for (int k = 0; k < 3; ++k) centroid[k] += simplex[i][k] / 3.0;

    Point3 worst = simplex[3];
    Point3 reflected;
    for (int k = 0; k < 3; ++k) reflected[k] = centroid[k] + kAlpha * (centroid[k] - worst[k]);
    const double f_reflected = eval(reflected); ++n_eval;

    if (f_reflected < fval[0]) {
      Point3 expanded;
      for (int k = 0; k < 3; ++k) expanded[k] = centroid[k] + kGamma * (reflected[k] - centroid[k]);
      const double f_expanded = eval(expanded); ++n_eval;
      if (f_expanded < f_reflected) { simplex[3] = expanded; fval[3] = f_expanded; }
      else { simplex[3] = reflected; fval[3] = f_reflected; }
    } else if (f_reflected < fval[2]) {
      simplex[3] = reflected; fval[3] = f_reflected;
    } else {
      Point3 contracted;
      for (int k = 0; k < 3; ++k) contracted[k] = centroid[k] + kRho * (worst[k] - centroid[k]);
      const double f_contracted = eval(contracted); ++n_eval;
      if (f_contracted < fval[3]) {
        simplex[3] = contracted; fval[3] = f_contracted;
      } else {
        for (int i = 1; i < 4; ++i) {
          for (int k = 0; k < 3; ++k)
            simplex[i][k] = simplex[0][k] + kSigmaShrink * (simplex[i][k] - simplex[0][k]);
          fval[i] = eval(simplex[i]); ++n_eval;
        }
      }
    }
  }
  sort_simplex();
  return {simplex[0], fval[0], n_eval};
}

// Parameter vector of the 1D fit.
using Point2 = std::array<double, 2>;  // {mux, u} -- 1D fit's unconstrained sigma surrogate

// Result of one 2-parameter Nelder-Mead run.
struct SimplexRun2 {
  Point2 best;
  double f_best;
  int n_eval;
};

// 2-parameter Nelder-Mead, same coefficients/convergence logic as RunSimplex
// above, just one fewer dimension (a 3-vertex simplex instead of 4).
// Nelder-Mead for the 1D fit (a 3-vertex simplex in {mux, u}); same logic and stopping rule as RunSimplex.
SimplexRun2 RunSimplex2(const std::function<double(const Point2&)>& eval,
                        double mux0, double u0, int max_iterations, double tol) {
  constexpr double kAlpha = 1.0, kGamma = 2.0, kRho = 0.5, kSigmaShrink = 0.5;

  std::array<Point2, 3> simplex = {{
      {mux0, u0},
      {mux0 + 0.6, u0},
      {mux0, u0 + 0.5},
  }};
  std::array<double, 3> fval;
  int n_eval = 0;
  for (int i = 0; i < 3; ++i) { fval[i] = eval(simplex[i]); ++n_eval; }

  auto sort_simplex = [&]() {
    std::array<int, 3> idx = {0, 1, 2};
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return fval[a] < fval[b]; });
    std::array<Point2, 3> s2;
    std::array<double, 3> f2;
    for (int i = 0; i < 3; ++i) { s2[i] = simplex[idx[i]]; f2[i] = fval[idx[i]]; }
    simplex = s2;
    fval = f2;
  };

  for (int iter = 0; iter < max_iterations; ++iter) {
    sort_simplex();
    double max_extent = 0.0;
    for (int k = 0; k < 2; ++k) {
      double lo_k = simplex[0][k], hi_k = simplex[0][k];
      for (int i = 1; i < 3; ++i) { lo_k = std::min(lo_k, simplex[i][k]); hi_k = std::max(hi_k, simplex[i][k]); }
      max_extent = std::max(max_extent, hi_k - lo_k);
    }
    if (std::abs(fval[2] - fval[0]) < tol && max_extent < 1e-5) break;

    Point2 centroid = {0.0, 0.0};
    for (int i = 0; i < 2; ++i)  // all but the worst (index 2)
      for (int k = 0; k < 2; ++k) centroid[k] += simplex[i][k] / 2.0;

    Point2 worst = simplex[2];
    Point2 reflected;
    for (int k = 0; k < 2; ++k) reflected[k] = centroid[k] + kAlpha * (centroid[k] - worst[k]);
    const double f_reflected = eval(reflected); ++n_eval;

    if (f_reflected < fval[0]) {
      Point2 expanded;
      for (int k = 0; k < 2; ++k) expanded[k] = centroid[k] + kGamma * (reflected[k] - centroid[k]);
      const double f_expanded = eval(expanded); ++n_eval;
      if (f_expanded < f_reflected) { simplex[2] = expanded; fval[2] = f_expanded; }
      else { simplex[2] = reflected; fval[2] = f_reflected; }
    } else if (f_reflected < fval[1]) {
      simplex[2] = reflected; fval[2] = f_reflected;
    } else {
      Point2 contracted;
      for (int k = 0; k < 2; ++k) contracted[k] = centroid[k] + kRho * (worst[k] - centroid[k]);
      const double f_contracted = eval(contracted); ++n_eval;
      if (f_contracted < fval[2]) {
        simplex[2] = contracted; fval[2] = f_contracted;
      } else {
        for (int i = 1; i < 3; ++i) {
          for (int k = 0; k < 2; ++k)
            simplex[i][k] = simplex[0][k] + kSigmaShrink * (simplex[i][k] - simplex[0][k]);
          fval[i] = eval(simplex[i]); ++n_eval;
        }
      }
    }
  }
  sort_simplex();
  return {simplex[0], fval[0], n_eval};
}

}  // namespace

// ----------------------------------------------------------------------------
// ClusterFitEngine::FitNelderMead_
//   2D fit with my Nelder-Mead. I map the width through a sigmoid,
//   sigma = lo + (hi-lo)/(1+exp(-u)), so it can never leave the allowed range,
//   and start from the seed plus two perturbed widths (u -/+ 1.5), keeping the
//   best result. Fills the fit result and Delta LL = G_min/(2*sigma_pix^2).
// ----------------------------------------------------------------------------
ClusterFitResult ClusterFitEngine::FitNelderMead_(const FitWindow& w, double mux0, double muy0,
                                                    double sigma0) const {
  // sigma = lo + (hi-lo)*sigmoid(u) maps any real u onto the open interval
  // (lo,hi), so sigma can never leave the valid range no matter what the
  // simplex tries -- replaces an earlier soft-clamp-at-evaluation approach
  // that made the objective perfectly FLAT beyond the bound (no gradient to
  // climb back down with, confirmed by the Minuit2 cross-check finding a
  // strictly better minimum in exactly that situation).
  const double lo = cfg_.sigma_xy_lo_px, hi = cfg_.sigma_xy_hi_px;
  auto u_to_sigma = [lo, hi](double u) { return lo + (hi - lo) / (1.0 + std::exp(-u)); };
  auto sigma_to_u = [lo, hi](double sigma) {
    const double frac = std::clamp((sigma - lo) / (hi - lo), 1e-6, 1.0 - 1e-6);
    return std::log(frac / (1.0 - frac));
  };
  std::function<double(const Point3&)> eval = [&](const Point3& p) {
    return ObjectiveG(w, p[0], p[1], u_to_sigma(p[2]));
  };
  const double u0 = sigma_to_u(sigma0);

  // Multi-start: run from the seed plus two perturbed sigma starts (one
  // narrower, one wider) and keep whichever converges to the best
  // objective. Guards against any single run settling in a locally flat
  // region near the sigmoid's saturating tails -- cheap here since each
  // run is a few hundred evaluations of a closed-form objective.
  SimplexRun best = RunSimplex(eval, mux0, muy0, u0, cfg_.max_iterations, cfg_.tol);
  for (double du : {-1.5, 1.5}) {
    SimplexRun alt = RunSimplex(eval, mux0, muy0, u0 + du, cfg_.max_iterations, cfg_.tol);
    if (alt.f_best < best.f_best) best = alt;
  }

  ClusterFitResult r;
  r.ok = true;
  r.mux_px = best.best[0];
  r.muy_px = best.best[1];
  r.sigma_xy_px = u_to_sigma(best.best[2]);
  r.I_hat_e = IHat(w, r.mux_px, r.muy_px, r.sigma_xy_px);
  r.delta_ll = DeltaLLFromG(best.f_best, cfg_.sigma_pix_e);
  r.n_eval = best.n_eval;
  return r;
}

// ----------------------------------------------------------------------------
// ClusterFitEngine::FitNelderMead1D_
//   1D fit (single collapsed row) with the same sigmoid mapping and multi-start strategy. muy is reported as 0.
// ----------------------------------------------------------------------------
ClusterFitResult ClusterFitEngine::FitNelderMead1D_(const FitWindow& w, double mux0, double sigma0) const {
  const double lo = cfg_.sigma_xy_lo_px, hi = cfg_.sigma_xy_hi_px;
  auto u_to_sigma = [lo, hi](double u) { return lo + (hi - lo) / (1.0 + std::exp(-u)); };
  auto sigma_to_u = [lo, hi](double sigma) {
    const double frac = std::clamp((sigma - lo) / (hi - lo), 1e-6, 1.0 - 1e-6);
    return std::log(frac / (1.0 - frac));
  };
  std::function<double(const Point2&)> eval = [&](const Point2& p) {
    return ObjectiveG1D(w, p[0], u_to_sigma(p[1]));
  };
  const double u0 = sigma_to_u(sigma0);

  SimplexRun2 best = RunSimplex2(eval, mux0, u0, cfg_.max_iterations, cfg_.tol);
  for (double du : {-1.5, 1.5}) {
    SimplexRun2 alt = RunSimplex2(eval, mux0, u0 + du, cfg_.max_iterations, cfg_.tol);
    if (alt.f_best < best.f_best) best = alt;
  }

  ClusterFitResult r;
  r.ok = true;
  r.mux_px = best.best[0];
  r.muy_px = 0.0;
  r.sigma_xy_px = u_to_sigma(best.best[1]);
  r.I_hat_e = IHat1D(w, r.mux_px, r.sigma_xy_px);
  r.delta_ll = DeltaLLFromG(best.f_best, cfg_.sigma_pix_e);
  r.n_eval = best.n_eval;
  return r;
}

#ifdef CCDARKSENS_USE_MINUIT2
namespace {

// Adapter that exposes the 3-parameter objective to ROOT's Minuit2.
class ClusterObjectiveFunctor : public ROOT::Math::IBaseFunctionMultiDim {
 public:
  ClusterObjectiveFunctor(const FitWindow& w, double sigma_lo, double sigma_hi)
      : w_(w), sigma_lo_(sigma_lo), sigma_hi_(sigma_hi) {}
  double DoEval(const double* x) const override {
    return EvalObjective(w_, x[0], x[1], x[2], sigma_lo_, sigma_hi_);
  }
  unsigned int NDim() const override { return 3; }
  ROOT::Math::IBaseFunctionMultiDim* Clone() const override {
    return new ClusterObjectiveFunctor(w_, sigma_lo_, sigma_hi_);
  }

 private:
  FitWindow w_;
  double sigma_lo_, sigma_hi_;
};

}  // namespace
#endif

// ----------------------------------------------------------------------------
// ClusterFitEngine::FitMinuit2_
//   2D fit with Minuit2 (Simplex) with sigma limited to the allowed range. I fall
//   back to Nelder-Mead when Minuit2 is not compiled in, cannot be created, or
//   fails to converge.
// ----------------------------------------------------------------------------
ClusterFitResult ClusterFitEngine::FitMinuit2_(const FitWindow& w, double mux0, double muy0,
                                                 double sigma0) const {
#ifdef CCDARKSENS_USE_MINUIT2
  std::unique_ptr<ROOT::Math::Minimizer> min(
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex"));
  if (!min) return FitNelderMead_(w, mux0, muy0, sigma0);  // fallback, same pattern as ProfileLikelihood.cc

  ClusterObjectiveFunctor obj(w, cfg_.sigma_xy_lo_px, cfg_.sigma_xy_hi_px);
  ROOT::Math::Functor f(obj, 3);
  min->SetFunction(f);
  min->SetVariable(0, "mux", mux0, 0.1);
  min->SetVariable(1, "muy", muy0, 0.1);
  min->SetVariable(2, "sigma", sigma0, 0.05);
  min->SetVariableLimits(2, cfg_.sigma_xy_lo_px, cfg_.sigma_xy_hi_px);
  min->SetPrintLevel(0);
  min->SetTolerance(cfg_.tol);
  if (!min->Minimize()) return FitNelderMead_(w, mux0, muy0, sigma0);

  const double* x = min->X();
  ClusterFitResult r;
  r.ok = true;
  r.mux_px = x[0];
  r.muy_px = x[1];
  r.sigma_xy_px = ClampSigma(x[2], cfg_.sigma_xy_lo_px, cfg_.sigma_xy_hi_px);
  r.I_hat_e = IHat(w, r.mux_px, r.muy_px, r.sigma_xy_px);
  r.delta_ll = DeltaLLFromG(min->MinValue(), cfg_.sigma_pix_e);
  r.n_eval = static_cast<int>(min->NCalls());
  return r;
#else
  return FitNelderMead_(w, mux0, muy0, sigma0);
#endif
}

// ----------------------------------------------------------------------------
// ClusterFitEngine::Fit
//   Public entry point. pixels is the row-major nx*ny window (as returned by
//   PixelSimulator::PixelCharges()). I seed, then run the 1D Nelder-Mead in
//   one-dimensional mode, otherwise Minuit2 or the 2D Nelder-Mead as configured.
// ----------------------------------------------------------------------------
ClusterFitResult ClusterFitEngine::Fit(const std::vector<double>& pixels, int nx, int ny) const {
  FitWindow w{&pixels, nx, ny, cfg_.sigma_pix_e};

  if (cfg_.one_dimensional) {
    // kMinuit2 has no 1D implementation -- falls back to Nelder-Mead, same
    // pattern as the 2D path's Minuit2-unavailable fallback below.
    double mux0, sigma0;
    Seed1D_(w, mux0, sigma0);
    return FitNelderMead1D_(w, mux0, sigma0);
  }

  double mux0, muy0, sigma0;
  Seed_(w, mux0, muy0, sigma0);

  if (cfg_.method == ClusterFitConfig::Method::kMinuit2) {
    return FitMinuit2_(w, mux0, muy0, sigma0);
  }
  return FitNelderMead_(w, mux0, muy0, sigma0);
}

}  // namespace ccdarksens
