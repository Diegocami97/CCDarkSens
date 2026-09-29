// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ClusterFitModel.hh -- Closed-form math for the pixel-window signal-vs-
//  noise ΔLL test: pixel-integrated 2D Gaussian shape, analytic amplitude
//  profiling, and the resulting ΔLL objective. No ROOT dependency — pure
//  arithmetic, reused by both the toy-noise calibration
//  (NoiseTailCalibrator) and the real-event fit (ClusterFitMC) via
//  ClusterFitEngine.
// ===========================================================================

#pragma once

#include <vector>

namespace ccdarksens {

/**
 * A pixel-charge window under test: a flat, row-major nx*ny array of
 * per-pixel charge (electrons), plus the per-pixel Gaussian readout-noise
 * sigma assumed for every pixel. Pixel (0,0) is index 0; pixel (ix,iy) is
 * index iy*nx + ix. Coordinates in Shape()/PixelIntegral1D() are in pixel
 * units (pitch = 1), with pixel ix spanning [ix-0.5, ix+0.5).
 */
// ----------------------------------------------------------------------------
// FitWindow
//   A view of one pixel window handed to the fit functions (documented above).
// ----------------------------------------------------------------------------
struct FitWindow {
  const std::vector<double>* pixels = nullptr;  // row-major pixel charges [e-] (not owned)
  int nx = 0;  // window width [pixels]
  int ny = 0;  // window height [pixels]
  double sigma_pix_e = 0.0;  // per-pixel Gaussian readout noise [e-]
};

/**
 * Fraction of a 1D Gaussian(mu, sigma)'s probability mass landing inside
 * the pixel centered at center_px (i.e. the interval [center-0.5, center+0.5)).
 * Exact erf-based CDF difference -- same idiom as Diffusion::normal_cdf_ /
 * interval_prob_ (src/response/Diffusion.cc), not a point-sampled density.
 */
double PixelIntegral1D(double center_px, double mu_px, double sigma_px);

/**
 * 2D pixel-integrated Gaussian shape at pixel (ix,iy): the fraction of a
 * unit-charge, isotropic 2D Gaussian cloud centered at (mux_px,muy_px) with
 * width sigma_px landing in that one pixel. Product of two PixelIntegral1D
 * calls (x and y separate for an isotropic Gaussian).
 */
double Shape(int ix, int iy, double mux_px, double muy_px, double sigma_px);

/**
 * Closed-form amplitude MLE at fixed (mux,muy,sigma): the pixel model is
 * data_ij = I*shape_ij + noise_ij, noise ~ Gaussian(0,sigma_pix) i.i.d., so
 * this is the weighted-least-squares solution I_hat = sum(data*shape) /
 * sum(shape^2). Returns 0 if sum(shape^2) is degenerate (e.g. window has
 * zero overlap with the signal shape).
 */
double IHat(const FitWindow& w, double mux_px, double muy_px, double sigma_px);

/**
 * The fit objective actually minimized by ClusterFitEngine:
 *   G(theta) = -IHat(theta)^2 * sum_ij shape_ij(theta)^2
 * Minimizing G over (mux,muy,sigma) is exactly equivalent to maximizing the
 * signal likelihood L_G(theta, IHat(theta)) -- see docs/ClusterFitMC_Design.md
 * for the full derivation of why this is correct and why L_n (the white-
 * noise null) never needs separate evaluation.
 */
double ObjectiveG(const FitWindow& w, double mux_px, double muy_px, double sigma_px);

/**
 * Converts a minimized objective value G(theta_hat) into the actual ΔLL
 * test statistic: ΔLL = G(theta_hat) / (2*sigma_pix^2). ΔLL <= 0; more
 * negative means more signal-like (better fit than the white-noise null).
 * Dimension-agnostic (used by both the 2D and 1D fit paths below).
 */
double DeltaLLFromG(double g_min, double sigma_pix_e);

/**
 * 1D analog of Shape()/IHat()/ObjectiveG(), for the DAMIC 1x100 readout mode:
 * 100 pixel rows are summed together in hardware before the single readout
 * (PhysRevD.94.082006 Sec. VI), so the truth landing in one such row is the
 * y-marginal of the 2D charge cloud -- and because the 2D shape is separable
 * (Shape = PixelIntegral1D(x) * PixelIntegral1D(y)), that marginal is exactly
 * PixelIntegral1D(x) alone, not a fixed-row slice of the 2D shape. `w` here
 * must have w.ny == 1 (one collapsed row); pixel ix is at flat index ix.
 */
double Shape1D(int ix, double mux_px, double sigma_px);
double IHat1D(const FitWindow& w, double mux_px, double sigma_px);
double ObjectiveG1D(const FitWindow& w, double mux_px, double sigma_px);

}  // namespace ccdarksens
