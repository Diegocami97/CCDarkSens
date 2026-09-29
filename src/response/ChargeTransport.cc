// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ChargeTransport.cc -- Samples electron depth and 2D Gaussian charge
//  clouds with the same σ_xy(z,E) parametrization used by diffusion
//  modeling.
// ===========================================================================

#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/DiffusionPhysics.hh"

#include <cmath>

namespace ccdarksens {

// Constructor: store the configuration and seed the RNG.
ChargeTransport::ChargeTransport(const ChargeTransportConfig& cfg)
  : cfg_(cfg),
    rng_(cfg.rng_seed)
{
}

// ----------------------------------------------------------------------------
// ChargeTransport::SampleDepthUm
//   Depth of an interaction, uniform between 0 and the sensor thickness [um].
// ----------------------------------------------------------------------------
double ChargeTransport::SampleDepthUm() {
  if (cfg_.thickness_um <= 0.0) return 0.0;
  const double u = uni_(rng_);
  return u * cfg_.thickness_um;
}

// ----------------------------------------------------------------------------
// ChargeTransport::SigmaXYUm
//   Lateral diffusion width sigma_xy(z, E) [um]; I use the shared ComputeSigmaXYUm().
// ----------------------------------------------------------------------------
double ChargeTransport::SigmaXYUm(double z_um, double Ee_eV) const {
  return ComputeSigmaXYUm(z_um, Ee_eV, cfg_.A_um2, cfg_.b_umInv, cfg_.alpha, cfg_.beta_per_keV);
}

// ----------------------------------------------------------------------------
// ChargeTransport::SampleCloudXY
//   I append n_e electron positions to xs_um / ys_um, each a Gaussian
//   displacement of width sigma_xy_um around (x0_um, y0_um).
// ----------------------------------------------------------------------------
void ChargeTransport::SampleCloudXY(double x0_um, double y0_um,
                                    double sigma_xy_um,
                                    int n_e,
                                    std::vector<double>& xs_um,
                                    std::vector<double>& ys_um)
{
  if (n_e <= 0) return;
  xs_um.reserve(xs_um.size() + static_cast<std::size_t>(n_e));
  ys_um.reserve(ys_um.size() + static_cast<std::size_t>(n_e));

  for (int i = 0; i < n_e; ++i) {
    const double dx = sigma_xy_um * gaus_(rng_);
    const double dy = sigma_xy_um * gaus_(rng_);
    xs_um.push_back(x0_um + dx);
    ys_um.push_back(y0_um + dy);
  }
}

} // namespace ccdarksens
