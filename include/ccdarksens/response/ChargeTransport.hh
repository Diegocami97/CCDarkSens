// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ChargeTransport.hh -- Header for depth sampling and 2D electron-cloud
//  generation with configurable diffusion parameters.
// ===========================================================================

#pragma once

#include <vector>
#include <random>

namespace ccdarksens {

// ----------------------------------------------------------------------------
// ChargeTransportConfig
//   Sensor thickness and diffusion parameters for the charge-transport model.
// ----------------------------------------------------------------------------
struct ChargeTransportConfig {
  double thickness_um   = 670.0;  // sensor thickness [um]
  double A_um2          = 803.25;   ///< match Diffusion config
  double b_umInv        = 6.5e-4;  // diffusion constant b [1/um]
  double alpha          = 1.0;  // energy-dependence offset of the width
  double beta_per_keV   = 0.0;  // energy-dependence slope of the width [1/keV]
  std::uint64_t rng_seed= 12345ULL;  // RNG seed
};

/// ChargeTransport: convert (E, depth) → 2D Gaussian cloud in (x,y).
/// Intent is to reuse the same functional form as Diffusion::sigma_xy_um_.
class ChargeTransport {
public:
  explicit ChargeTransport(const ChargeTransportConfig& cfg);

  /// Sample a depth uniformly in [0, thickness].
  double SampleDepthUm();

  /// Compute σ_xy(z, E) in μm using my standard parametrization.
  double SigmaXYUm(double z_um, double Ee_eV) const;

  /// Sample n_e electrons around (x0,y0) with Gaussian σ_xy.
  void SampleCloudXY(double x0_um, double y0_um,
                     double sigma_xy_um,
                     int n_e,
                     std::vector<double>& xs_um,
                     std::vector<double>& ys_um);

private:
  ChargeTransportConfig cfg_;  // configuration
  std::mt19937_64 rng_;  // random-number generator
  std::uniform_real_distribution<double> uni_{0.0, 1.0};  // uniform [0,1) draws (depth sampling)
  std::normal_distribution<double> gaus_{0.0, 1.0};  // unit Gaussian draws (cloud sampling)
};

} // namespace ccdarksens
