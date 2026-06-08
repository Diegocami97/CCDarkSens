// ============================================================================
//  CCDarkSens — ProfileLikelihood
//  Header declaring the pattern-space profile likelihood API (data/B templates, constraints, and 1D/2D minimizers).
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#pragma once

#include <vector>
#include <functional>

namespace ccdarksens::stats {

/**
 * Profile likelihood in pattern space, matching pydme:
 * - Data: either Asimov (data = B) or real observed counts D.
 * - Background: either B(scale) = scale * B_template, or B(theta) = Bp + theta * Br (pydme).
 * - NLL = sum_i [ mu_i - data_i * ln(mu_i) ] + optional L(theta/scale).
 *
 * Use MinimizeOverScale (or over theta in Bp_Br mode) to get (param_hat, nll) for a given S.
 * Then q = 2*(nll_test - nll_null) with null = (S=0), test = (S=S_pat).
 */
class ProfileLikelihood {
public:
  ProfileLikelihood() = default;

  /// Set observed counts (length = n_bins). For Asimov use background as data.
  void SetData(const std::vector<double>& data);

  /// Scale mode: B = scale * B_template. Model is S + scale * B_template.
  void SetBTemplate(const std::vector<double>& B_template);

  /// Pydme mode: B_i = Bp_i + theta * Br_i. One global theta. Call after SetData (same size).
  void SetBpBr(const std::vector<double>& Bp, const std::vector<double>& Br);

  /// Pydme constrain: L(theta) = sum_i ( ±(-theta*Br_i + prior_strength*ln(theta*Br_i)) ). Use 0 to disable.
  void SetConstrainPriorStrength(double prior_strength) { constrain_prior_strength_ = prior_strength; }
  /// When true, use Gamma-prior sign (+theta*Br - N*ln(theta*Br)) so prior pulls theta toward mode N/Br; when false, match pydme (-theta*Br + N*ln(theta*Br)).
  void SetConstrainGammaSign(bool use_gamma_sign) { constrain_gamma_sign_ = use_gamma_sign; }
  /// When true, use pydme tau-weighted form: tau=N/sum(Br), constraint_i = -theta*tau*Br_i + (N*Br_i/sum(Br))*ln(theta*tau*Br_i); prior mode at theta=1 (nominal).
  void SetConstrainUseTauWeighted(bool use_tau) { constrain_use_tau_weighted_ = use_tau; }
  /// When > 1, use pydme multi-bin form: -theta*Br_i + prior_strength*n_bins*ln(theta*Br_i/n_bins) per pattern (matches pydme when len(gamma)=n_bins). Default 1 = single-bin form.
  void SetConstrainNBins(int n_bins) { constrain_n_bins_ = n_bins <= 0 ? 1 : n_bins; }
  /// When true, 2D (log10_sigma, theta) fit is not rejected when minimum is at a boundary; matches pydme which does not reject. Default false.
  // When true, both the 1D θ-minimization and the 2D (σ, θ) minimization may
  // return solutions on the parameter boundary instead of being rejected as
  // spurious. This matches pydme's behavior, where the constraint's gradient
  // at the lower θ boundary can drive the global NLL minimum to the boundary
  // legitimately. The flag affects MinimizeOverScaleMinuit() and
  // MinimizeOverSigmaAndTheta(); the legacy method name is kept for API
  // compatibility with existing call sites.
  void SetAcceptBoundary2DMinimum(bool accept) { accept_boundary_ = accept; }

  /// Optional extra constrain term (e.g. Gaussian prior). Applied on top of pydme L if in Bp_Br mode.
  void SetConstrain(std::function<double(double param)> constrain);

  /// NLL(data | S + B(param)) + optional constrain(param). param = scale or theta.
  double NLL(const std::vector<double>& S, double param) const;
  /// Minimize NLL over scale/theta for given S. Returns (param_hat, nll_min). Bounds [lo, hi].
  std::pair<double, double> MinimizeOverScale(const std::vector<double>& S,
                                               double scale_lo = 0.01,
                                               double scale_hi = 10.0,
                                               double tol = 1e-6) const;
  /// Same as MinimizeOverScale but using ROOT Minuit2 (matches pydme minimizer). Falls back to Brent if Minuit2 unavailable.
  std::pair<double, double> MinimizeOverScaleMinuit(const std::vector<double>& S,
                                                    double scale_lo = 0.01,
                                                    double scale_hi = 10.0) const;

  /// Pydme-style: minimize NLL over (log10(sigma_e), theta) in one 2D Minuit run.
  /// S_from_log10_sigma(log10_sigma) returns signal vector for that cross-section (e.g. from interpolation).
  /// Returns (log10_sigma_hat, theta_hat, nll_min). ok is false if Minuit2 unavailable or minimization failed.
  using S_from_log10_sigma_t = std::function<std::vector<double>(double log10_sigma)>;
  struct Minimize2DResult {
    double log10_sigma_hat = 0.0;
    double theta_hat = 0.0;
    double nll_min = 0.0;
    bool ok = false;
  };
  Minimize2DResult MinimizeOverSigmaAndTheta(S_from_log10_sigma_t S_from_log10_sigma,
                                             double log10_sigma_lo,
                                             double log10_sigma_hi,
                                             double theta_lo,
                                             double theta_hi) const;

  /// Same as above, but with explicit starting points. Pass NaN to use the
  /// default box-midpoint seed (current behavior). Useful for sub-toy loops
  /// where the global NLL minimum is known a priori to lie near a specific
  /// (log10_sigma, theta) — e.g. threshold toys generated at a fixed σ_thr.
  Minimize2DResult MinimizeOverSigmaAndTheta(S_from_log10_sigma_t S_from_log10_sigma,
                                             double log10_sigma_lo,
                                             double log10_sigma_hi,
                                             double theta_lo,
                                             double theta_hi,
                                             double log10_sigma_start,
                                             double theta_start) const;

  /// Profile likelihood ratio: q = 2*(NLL_test - NLL_null). Null and test both minimize over param.
  double EvaluateRatio(const std::vector<double>& S_null,
                      const std::vector<double>& S_test,
                      double scale_lo = 0.01,
                      double scale_hi = 10.0) const;

  std::size_t Nbins() const { return data_.size(); }
  const std::vector<double>& Data() const { return data_; }
  const std::vector<double>& BTemplate() const { return B_template_; }
  /// True if Bp+theta*Br mode is active.
  bool UseBpBr() const { return !Bp_.empty() && !Br_.empty(); }

private:
  std::vector<double> data_;
  std::vector<double> B_template_;
  std::vector<double> Bp_;
  std::vector<double> Br_;
  double constrain_prior_strength_ = 0.0;
  bool constrain_gamma_sign_ = false;  // false = pydme sign; true = Gamma prior sign (+tBr - N*ln(tBr))
  bool constrain_use_tau_weighted_ = false;  // true = tau-weighted form so prior mode at theta=1
  int constrain_n_bins_ = 1;  // pydme len(gamma); 1 = single-bin constraint form
  bool accept_boundary_ = false;  // if true, both 1D and 2D fits may return boundary solutions (pydme-style)
  std::function<double(double)> constrain_{};
};

}  // namespace ccdarksens::stats
