// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ClusterFitModel.cc -- Closed-form math for the pixel-window signal-vs-
//  noise ΔLL test: pixel-integrated 2D Gaussian shape, analytic amplitude
//  profiling, and the resulting ΔLL objective.
// ===========================================================================

#include "ccdarksens/response/ClusterFitModel.hh"
#include <cmath>

namespace ccdarksens {

namespace {

// Standard normal CDF, Phi(x) = 0.5*(1 + erf(x/sqrt(2))).
double NormalCdf(double x) {
  // Phi(x) = 0.5 * [1 + erf(x / sqrt(2))] -- same formula as Diffusion::normal_cdf_
  return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

}  // namespace

// ----------------------------------------------------------------------------
// PixelIntegral1D
//   Probability mass of a 1D Gaussian(mu, sigma) inside the pixel
//   [center-0.5, center+0.5): the difference of two normal CDFs. For a
//   non-positive sigma the Gaussian is a point at mu, so I return 1 or 0.
// ----------------------------------------------------------------------------
double PixelIntegral1D(double center_px, double mu_px, double sigma_px) {
  if (!(sigma_px > 0.0)) {
    return (center_px - 0.5 <= mu_px && mu_px < center_px + 0.5) ? 1.0 : 0.0;
  }
  const double z_hi = (center_px + 0.5 - mu_px) / sigma_px;
  const double z_lo = (center_px - 0.5 - mu_px) / sigma_px;
  return std::max(0.0, NormalCdf(z_hi) - NormalCdf(z_lo));
}

// ----------------------------------------------------------------------------
// Shape
//   2D pixel-integrated Gaussian: PixelIntegral1D in x times PixelIntegral1D in y.
// ----------------------------------------------------------------------------
double Shape(int ix, int iy, double mux_px, double muy_px, double sigma_px) {
  return PixelIntegral1D(static_cast<double>(ix), mux_px, sigma_px) *
         PixelIntegral1D(static_cast<double>(iy), muy_px, sigma_px);
}

// ----------------------------------------------------------------------------
// IHat
//   Best-fit amplitude at fixed shape parameters: I_hat = sum(data*shape) /
//   sum(shape^2), or 0 if the shape has no weight in the window.
// ----------------------------------------------------------------------------
double IHat(const FitWindow& w, double mux_px, double muy_px, double sigma_px) {
  double num = 0.0;
  double den = 0.0;
  for (int iy = 0; iy < w.ny; ++iy) {
    for (int ix = 0; ix < w.nx; ++ix) {
      const double s = Shape(ix, iy, mux_px, muy_px, sigma_px);
      const double d = (*w.pixels)[static_cast<std::size_t>(iy * w.nx + ix)];
      num += d * s;
      den += s * s;
    }
  }
  return (den > 0.0) ? (num / den) : 0.0;
}

// ----------------------------------------------------------------------------
// ObjectiveG
//   Fit objective G = -I_hat^2 * sum(shape^2) that the fitter minimizes over
//   (mux, muy, sigma); I compute I_hat inline to avoid a second pass over the
//   window. Returns 0 if the shape has no weight in the window.
// ----------------------------------------------------------------------------
double ObjectiveG(const FitWindow& w, double mux_px, double muy_px, double sigma_px) {
  double sum_shape2 = 0.0;
  double sum_data_shape = 0.0;
  for (int iy = 0; iy < w.ny; ++iy) {
    for (int ix = 0; ix < w.nx; ++ix) {
      const double s = Shape(ix, iy, mux_px, muy_px, sigma_px);
      const double d = (*w.pixels)[static_cast<std::size_t>(iy * w.nx + ix)];
      sum_shape2 += s * s;
      sum_data_shape += d * s;
    }
  }
  if (sum_shape2 <= 0.0) return 0.0;
  const double i_hat = sum_data_shape / sum_shape2;
  return -(i_hat * i_hat) * sum_shape2;
}

// ----------------------------------------------------------------------------
// DeltaLLFromG
//   Delta LL = G_min / (2*sigma_pix^2), which is <= 0. Returns 0 for sigma_pix <= 0.
// ----------------------------------------------------------------------------
double DeltaLLFromG(double g_min, double sigma_pix_e) {
  if (!(sigma_pix_e > 0.0)) return 0.0;
  return g_min / (2.0 * sigma_pix_e * sigma_pix_e);
}

// ----------------------------------------------------------------------------
// Shape1D
//   1D shape for the collapsed-row (1x100) mode: the pixel-integrated Gaussian in x only.
// ----------------------------------------------------------------------------
double Shape1D(int ix, double mux_px, double sigma_px) {
  return PixelIntegral1D(static_cast<double>(ix), mux_px, sigma_px);
}

// ----------------------------------------------------------------------------
// IHat1D
//   1D counterpart of IHat over the single collapsed row.
// ----------------------------------------------------------------------------
double IHat1D(const FitWindow& w, double mux_px, double sigma_px) {
  double num = 0.0;
  double den = 0.0;
  for (int ix = 0; ix < w.nx; ++ix) {
    const double s = Shape1D(ix, mux_px, sigma_px);
    const double d = (*w.pixels)[static_cast<std::size_t>(ix)];
    num += d * s;
    den += s * s;
  }
  return (den > 0.0) ? (num / den) : 0.0;
}

// ----------------------------------------------------------------------------
// ObjectiveG1D
//   1D counterpart of ObjectiveG over the single collapsed row.
// ----------------------------------------------------------------------------
double ObjectiveG1D(const FitWindow& w, double mux_px, double sigma_px) {
  double sum_shape2 = 0.0;
  double sum_data_shape = 0.0;
  for (int ix = 0; ix < w.nx; ++ix) {
    const double s = Shape1D(ix, mux_px, sigma_px);
    const double d = (*w.pixels)[static_cast<std::size_t>(ix)];
    sum_shape2 += s * s;
    sum_data_shape += d * s;
  }
  if (sum_shape2 <= 0.0) return 0.0;
  const double i_hat = sum_data_shape / sum_shape2;
  return -(i_hat * i_hat) * sum_shape2;
}

}  // namespace ccdarksens
