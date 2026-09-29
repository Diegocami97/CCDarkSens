// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ProfileLikelihood.cc -- Implements pydme-compatible Poisson profile
//  likelihood with Bp+θ·Br backgrounds, optional priors, and Brent/Minuit
//  profiling over θ and (log σ, θ).
// ===========================================================================

#include "ccdarksens/stats/ProfileLikelihood.hh"
#include "ccdarksens/stats/StatsUtils.hh"
#include <cmath>
#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>

#ifdef CCDARKSENS_USE_MINUIT2
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/IFunction.h>
#endif

namespace ccdarksens::stats {

// Set the observed counts per bin (use the background for an Asimov data set).
void ProfileLikelihood::SetData(const std::vector<double>& data) {
  data_ = data;
}

// Switch to scale mode: B = scale * B_template (this clears any Bp/Br).
void ProfileLikelihood::SetBTemplate(const std::vector<double>& B_template) {
  B_template_ = B_template;
  Bp_.clear();
  Br_.clear();
}

// Switch to theta mode: B = Bp + theta * Br. Throws std::invalid_argument if the sizes differ.
void ProfileLikelihood::SetBpBr(const std::vector<double>& Bp, const std::vector<double>& Br) {
  if (Bp.size() != Br.size())
    throw std::invalid_argument("ProfileLikelihood::SetBpBr: Bp and Br size mismatch");
  Bp_ = Bp;
  Br_ = Br;
}

// Set an extra constraint term added to the NLL (for example a Gaussian prior).
void ProfileLikelihood::SetConstrain(std::function<double(double param)> constrain) {
  constrain_ = std::move(constrain);
}

// ----------------------------------------------------------------------------
// ProfileLikelihood::NLL
//   NLL = sum_i [ mu_i - n_i ln(mu_i) ] with mu_i = S_i + B_i(param), where B is
//   scale*B_template (scale mode) or Bp + theta*Br (theta mode). A bin with
//   mu <= 0 and data > 0 costs a 1e9 penalty. In theta mode the pydme constraint
//   on theta is added (plain, Gamma-sign, tau-weighted or multi-bin form, as
//   configured), then the optional user constraint. Throws
//   std::invalid_argument if the vector sizes do not match.
// ----------------------------------------------------------------------------
double ProfileLikelihood::NLL(const std::vector<double>& S, double param) const {
  const std::size_t n = data_.size();
  if (n != S.size())
    throw std::invalid_argument("ProfileLikelihood::NLL: data/S size mismatch");

  if (UseBpBr()) {
    if (n != Bp_.size() || n != Br_.size())
      throw std::invalid_argument("ProfileLikelihood::NLL: Bp/Br size mismatch");
    double nll = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double B_i = Bp_[i] + param * Br_[i];
      const double mu = S[i] + B_i;
      const double n_i = data_[i];
      if (mu <= 0.0) {
        if (n_i > 0.0) nll += 1e9;
        continue;
      }
      nll += mu - n_i * safe_log(mu);
    }
    // Constraint (prior on theta): see pydme constraint_pattern (tau-weighted) vs background_pattern L_array (plain).
    // Tau-weighted: tau = N/sum(Br), constraint_i = -theta*tau*Br_i + (N*Br_i/sum(Br))*ln(theta*tau*Br_i); minimum at theta=1.
    if (constrain_prior_strength_ > 0.0) {
      if (constrain_use_tau_weighted_) {
        double sum_Br = 0.0;
        for (std::size_t i = 0; i < n; ++i) sum_Br += Br_[i];
        if (sum_Br <= 0.0) nll += 1e9;
        else {
          const double tau = constrain_prior_strength_ / sum_Br;
          for (std::size_t i = 0; i < n; ++i) {
            const double x = param * tau * Br_[i];
            if (x <= 0.0) nll += 1e9;
            else {
              const double Nrc_weighted = constrain_prior_strength_ * Br_[i] / sum_Br;
              // pydme: nLnL += -sum(poisson + constraint), so we add -constraint = +x - Nrc_weighted*ln(x) to NLL
              nll += x - Nrc_weighted * safe_log(x);
            }
          }
        }
      } else {
        const double n_bins = static_cast<double>(constrain_n_bins_);
        for (std::size_t i = 0; i < n; ++i) {
          const double tBr = param * Br_[i];
          if (tBr <= 0.0) nll += 1e9;
          else if (constrain_gamma_sign_) {
            if (constrain_n_bins_ <= 1)
              nll += tBr - constrain_prior_strength_ * safe_log(tBr);
            else
              nll += tBr - constrain_prior_strength_ * n_bins * safe_log(tBr / n_bins);
          } else {
            if (constrain_n_bins_ <= 1)
              nll += -tBr + constrain_prior_strength_ * safe_log(tBr);
            else
              nll += -tBr + constrain_prior_strength_ * n_bins * safe_log(tBr / n_bins);
          }
        }
      }
    }
    if (constrain_) nll += constrain_(param);
    return nll;
  }

  if (n != B_template_.size())
    throw std::invalid_argument("ProfileLikelihood::NLL: B_template size mismatch");
  double nll = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double mu = S[i] + param * B_template_[i];
    const double n_i = data_[i];
    if (mu <= 0.0) {
      if (n_i > 0.0) nll += 1e9;
      continue;
    }
    nll += mu - n_i * safe_log(mu);
  }
  if (constrain_) nll += constrain_(param);
  return nll;
}

// Brent's method (1D minimization) for scale in [scale_lo, scale_hi]
// ----------------------------------------------------------------------------
// ProfileLikelihood::MinimizeOverScale
//   Brent's 1D minimization of the NLL over the nuisance parameter in
//   [scale_lo, scale_hi] (200 iterations at most). Returns {param_hat, nll_min}.
// ----------------------------------------------------------------------------
std::pair<double, double> ProfileLikelihood::MinimizeOverScale(
    const std::vector<double>& S,
    double scale_lo, double scale_hi, double tol) const {
  constexpr double cgold = 0.381966;
  double a = scale_lo, b = scale_hi;
  double x = a + cgold * (b - a), w = x, v = x;
  double fx = NLL(S, x), fw = fx, fv = fx;
  double d = 0.0, e = 0.0;
  const int maxiter = 200;
  for (int iter = 0; iter < maxiter; ++iter) {
    const double xm = 0.5 * (a + b);
    if (std::abs(x - xm) <= 2.0 * tol - 0.5 * (b - a))
      return {x, fx};
    if (std::abs(e) > tol) {
      const double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      p = 2 * (q - r) >= 0 ? p : -p;
      q = std::abs(q - r);
      if (p > 0) q = -q;
      p = std::abs(p);
      const double etemp = e;
      e = d;
      if (std::abs(p) >= std::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
        d = cgold * (e = (x >= xm ? a - x : b - x));
      else {
        d = p / q;
        const double u = x + d;
        if (u - a < 2 * tol || b - u < 2 * tol) d = (xm - x > 0 ? tol : -tol);
      }
    } else {
      d = cgold * (e = (x >= xm ? a - x : b - x));
    }
    const double u = (std::abs(d) >= tol ? x + d : x + (d > 0 ? tol : -tol));
    const double fu = NLL(S, u);
    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      v = w; w = x; x = u;
      fv = fw; fw = fx; fx = fu;
    } else {
      if (u >= x) b = u; else a = u;
      if (fu <= fw || w == x) { v = w; w = u; fv = fw; fw = fu; }
      else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
    }
  }
  return {x, fx};
}

// ----------------------------------------------------------------------------
// ProfileLikelihood::EvaluateRatio
//   q = 2*(NLL_test - NLL_null), each profiled over the nuisance parameter; floored at 0.
// ----------------------------------------------------------------------------
double ProfileLikelihood::EvaluateRatio(const std::vector<double>& S_null,
                                       const std::vector<double>& S_test,
                                       double scale_lo, double scale_hi) const {
  const auto [scale_hat_null, nll_null] = MinimizeOverScale(S_null, scale_lo, scale_hi);
  const auto [scale_hat_test, nll_test] = MinimizeOverScale(S_test, scale_lo, scale_hi);
  (void)scale_hat_null;
  (void)scale_hat_test;
  const double q = 2.0 * (nll_test - nll_null);
  return q < 0.0 ? 0.0 : q;
}

#ifdef CCDARKSENS_USE_MINUIT2
namespace {
// Adapter exposing the 1D NLL(theta) at fixed signal to Minuit2.
class NLLFunctor : public ROOT::Math::IBaseFunctionMultiDim {
public:
  NLLFunctor(const ProfileLikelihood* pl, const std::vector<double>& S) : pl_(pl), S_(S) {}
  double DoEval(const double* x) const override { return pl_->NLL(S_, x[0]); }
  unsigned int NDim() const override { return 1; }
  ROOT::Math::IBaseFunctionMultiDim* Clone() const override { return new NLLFunctor(pl_, S_); }

private:
  const ProfileLikelihood* pl_;
  std::vector<double> S_;
};

// Adapter exposing NLL(log10 sigma, theta) to Minuit2 (the signal is rebuilt at every call).
class NLL2DFunctor : public ROOT::Math::IBaseFunctionMultiDim {
public:
  NLL2DFunctor(const ProfileLikelihood* pl,
               ProfileLikelihood::S_from_log10_sigma_t S_from_log10_sigma)
      : pl_(pl), S_from_log10_sigma_(std::move(S_from_log10_sigma)) {}
  double DoEval(const double* x) const override {
    S_buf_ = S_from_log10_sigma_(x[0]);
    return pl_->NLL(S_buf_, x[1]);
  }
  unsigned int NDim() const override { return 2; }
  ROOT::Math::IBaseFunctionMultiDim* Clone() const override {
    return new NLL2DFunctor(pl_, S_from_log10_sigma_);
  }

private:
  const ProfileLikelihood* pl_;
  ProfileLikelihood::S_from_log10_sigma_t S_from_log10_sigma_;
  mutable std::vector<double> S_buf_;
};
}  // namespace
#endif

// ----------------------------------------------------------------------------
// ProfileLikelihood::MinimizeOverScaleMinuit
//   Minuit2 (Simplex) version of MinimizeOverScale. Boundary solutions are
//   re-done with Brent unless accept_boundary is set; with accept_boundary
//   I also compare against both bounds explicitly and keep the lowest NLL, so
//   the 1D and 2D fits treat boundaries identically. Falls back to Brent when
//   Minuit2 is not compiled in or fails.
// ----------------------------------------------------------------------------
std::pair<double, double> ProfileLikelihood::MinimizeOverScaleMinuit(
    const std::vector<double>& S,
    double scale_lo, double scale_hi) const {
#ifdef CCDARKSENS_USE_MINUIT2
  // Use Simplex to avoid "VariableMetricBuilder Initial matrix not pos.def" (Migrad needs pos.def Hessian).
  std::unique_ptr<ROOT::Math::Minimizer> min(
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex"));
  if (!min) return MinimizeOverScale(S, scale_lo, scale_hi);

  NLLFunctor nll(this, S);
  ROOT::Math::Functor f(nll, 1);
  min->SetFunction(f);
  const double start = 0.5 * (scale_lo + scale_hi);
  const double step = 0.01 * (scale_hi - scale_lo);
  if (step <= 0.0) return MinimizeOverScale(S, scale_lo, scale_hi);
  min->SetVariable(0, "theta", start, step);
  min->SetVariableLimits(0, scale_lo, scale_hi);
  min->SetPrintLevel(0);
  min->SetTolerance(1e-8);  // avoid stopping at a rough minimum
  if (!min->Minimize()) return MinimizeOverScale(S, scale_lo, scale_hi);
  const double* x = min->X();
  double theta_hat = x[0];
  double nll_hat = min->MinValue();
  // If Minuit converged to a boundary, the minimum may be spurious -> fall
  // back to Brent. However, when accept_boundary_ is enabled (pydme-style),
  // boundary minima are legitimate (the constraint prior can drive theta to
  // the lower bound), so we keep them and skip the fallback. This is
  // critical: if 1D and 2D fits use *different* policies for boundaries,
  // the 1D result may sit slightly inside the bound while the 2D result
  // sits exactly on the bound; with a steep constraint gradient at the
  // boundary that mismatch produces a spurious finite NLL gap and inflates
  // q = 2*(NLL_top - NLL_glob).
  const double margin = 0.01 * (scale_hi - scale_lo);
  const bool at_boundary =
      (theta_hat <= scale_lo + margin || theta_hat >= scale_hi - margin);
  if (at_boundary && !accept_boundary_) {
    return MinimizeOverScale(S, scale_lo, scale_hi);
  }
  // When accept_boundary_ is enabled, also explicitly probe the lower and
  // upper θ boundaries and keep whichever NLL is lowest. Simplex is a local
  // search and may converge near (but not exactly on) the boundary even when
  // the true minimum is on the boundary. Without this snap-to-boundary step
  // the 1D θ-only fit can disagree with the 2D (σ, θ) fit by a few floating
  // units of θ; with a steep constraint gradient at the bound that maps to
  // an O(10) artificial NLL gap and an inflated test statistic.
  if (accept_boundary_) {
    const double nll_lo = NLL(S, scale_lo);
    if (nll_lo < nll_hat) {
      theta_hat = scale_lo;
      nll_hat = nll_lo;
    }
    const double nll_hi = NLL(S, scale_hi);
    if (nll_hi < nll_hat) {
      theta_hat = scale_hi;
      nll_hat = nll_hi;
    }
  }
  return {theta_hat, nll_hat};
#else
  return MinimizeOverScale(S, scale_lo, scale_hi);
#endif
}

// ----------------------------------------------------------------------------
// ProfileLikelihood::MinimizeOverSigmaAndTheta (default start)
//   Joint minimization with the seed at the middle of the box; forwards to the version below.
// ----------------------------------------------------------------------------
ProfileLikelihood::Minimize2DResult ProfileLikelihood::MinimizeOverSigmaAndTheta(
    S_from_log10_sigma_t S_from_log10_sigma,
    double log10_sigma_lo, double log10_sigma_hi,
    double theta_lo, double theta_hi) const {
  return MinimizeOverSigmaAndTheta(std::move(S_from_log10_sigma),
                                   log10_sigma_lo, log10_sigma_hi,
                                   theta_lo, theta_hi,
                                   std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN());
}

// ----------------------------------------------------------------------------
// ProfileLikelihood::MinimizeOverSigmaAndTheta
//   Minuit2 (Simplex) minimization over (log10 sigma, theta) inside the given
//   box, starting at the supplied point (NaN = box centre, clamped into the
//   box). Unless accept_boundary is set, a minimum on the box edge is rejected
//   as spurious. result.ok is false if Minuit2 is unavailable or fails.
// ----------------------------------------------------------------------------
ProfileLikelihood::Minimize2DResult ProfileLikelihood::MinimizeOverSigmaAndTheta(
    S_from_log10_sigma_t S_from_log10_sigma,
    double log10_sigma_lo, double log10_sigma_hi,
    double theta_lo, double theta_hi,
    double log10_sigma_start, double theta_start) const {
  Minimize2DResult out;
  out.ok = false;
#ifdef CCDARKSENS_USE_MINUIT2
  std::unique_ptr<ROOT::Math::Minimizer> min(
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex"));
  if (!min) return out;

  S_from_log10_sigma_t S_for_check = S_from_log10_sigma;  // copy for validation after move
  NLL2DFunctor nll2(this, std::move(S_from_log10_sigma));
  ROOT::Math::Functor f(nll2, 2);
  min->SetFunction(f);
  if (!std::isfinite(log10_sigma_start))
    log10_sigma_start = 0.5 * (log10_sigma_lo + log10_sigma_hi);
  else
    log10_sigma_start =
        std::max(log10_sigma_lo, std::min(log10_sigma_hi, log10_sigma_start));
  if (!std::isfinite(theta_start))
    theta_start = 0.5 * (theta_lo + theta_hi);
  else
    theta_start = std::max(theta_lo, std::min(theta_hi, theta_start));
  const double step_log10 = 0.01 * (log10_sigma_hi - log10_sigma_lo);
  const double step_theta = 0.01 * (theta_hi - theta_lo);
  if (step_log10 <= 0.0 || step_theta <= 0.0) return out;
  min->SetVariable(0, "log10_sigma", log10_sigma_start, step_log10);
  min->SetVariableLimits(0, log10_sigma_lo, log10_sigma_hi);
  min->SetVariable(1, "theta", theta_start, step_theta);
  min->SetVariableLimits(1, theta_lo, theta_hi);
  min->SetPrintLevel(0);
  min->SetTolerance(1e-8);
  if (!min->Minimize()) return out;
  const double* x = min->X();
  out.log10_sigma_hat = x[0];
  out.theta_hat = x[1];
  out.nll_min = min->MinValue();

  // Reject if minimum is at a boundary (spurious for Simplex), unless accept_boundary_ (pydme-style).
  if (!accept_boundary_) {
    const double margin_log10 = 0.02 * (log10_sigma_hi - log10_sigma_lo);
    const double margin_theta = 0.02 * (theta_hi - theta_lo);
    if (out.log10_sigma_hat <= log10_sigma_lo + margin_log10 ||
        out.log10_sigma_hat >= log10_sigma_hi - margin_log10 ||
        out.theta_hat <= theta_lo + margin_theta ||
        out.theta_hat >= theta_hi - margin_theta)
      return out;  // ok stays false
  }

  // Sanity check: NLL at found point should match MinValue.
  std::vector<double> S_check = S_for_check(out.log10_sigma_hat);
  if (S_check.size() != data_.size()) return out;
  const double nll_check = NLL(S_check, out.theta_hat);
  if (std::abs(nll_check - out.nll_min) > 0.01 * (1.0 + std::abs(out.nll_min)))
    return out;  // ok stays false

  out.ok = true;
#endif
  return out;
}

}  // namespace ccdarksens::stats
