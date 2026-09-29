// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PixelSimulator.cc -- Maps continuous (x,y) electron positions to a
//  discrete pixel charge array with optional Poisson DC and Gaussian readout
//  noise.
// ===========================================================================

#include "ccdarksens/response/PixelSimulator.hh"

#include <cmath>

namespace ccdarksens {

// Constructor: allocate a zeroed nx*ny map and seed the RNG and the dark-current Poisson.
PixelSimulator::PixelSimulator(const PixelSimulatorConfig& cfg)
  : cfg_(cfg),
    qpix_(static_cast<std::size_t>(cfg.nx * cfg.ny), 0.0),
    rng_(cfg.rng_seed),
    pois_dc_( (cfg.lambda_dc > 0.0) ? cfg.lambda_dc : 0.0 )
{
}

// Set every pixel charge to zero.
void PixelSimulator::Reset() {
  std::fill(qpix_.begin(), qpix_.end(), 0.0);
}

// True if pixel (ix, iy) lies inside the patch.
bool PixelSimulator::InBounds_(int ix, int iy) const {
  return (ix >= 0 && ix < cfg_.nx && iy >= 0 && iy < cfg_.ny);
}

// ----------------------------------------------------------------------------
// PixelSimulator::DepositElectron
//   Add one electron at (x_um, y_um), with (0,0) at the patch centre, to the
//   nearest pixel. Electrons outside the patch are dropped. In collapse_y (1x100)
//   mode the y position is ignored and the charge always goes to row 0.
// ----------------------------------------------------------------------------
void PixelSimulator::DepositElectron(double x_um, double y_um) {
  // For now we treat (x_um, y_um) such that (0,0) is the center of the patch
  // (or sensor), and map to nearest pixel.
  const double pitch = cfg_.pixel_size_um;
  const int ix0 = static_cast<int>(std::lround(x_um / pitch)) + cfg_.nx / 2;

  if (cfg_.collapse_y) {
    // 1x100 readout: 100 physical rows are summed in hardware before the
    // single readout, so every electron's charge lands in row 0 regardless
    // of its y offset -- no y bounds-check/rejection.
    if (ix0 < 0 || ix0 >= cfg_.nx) return;
    qpix_[static_cast<std::size_t>(ix0)] += 1.0;
    return;
  }

  const int iy0 = static_cast<int>(std::lround(y_um / pitch)) + cfg_.ny / 2;
  if (!InBounds_(ix0, iy0)) return;
  const std::size_t idx = static_cast<std::size_t>(iy0 * cfg_.nx + ix0);
  qpix_[idx] += 1.0;
}

// ----------------------------------------------------------------------------
// PixelSimulator::AddDarkCurrent
//   Add Poisson(lambda_dc) electrons to every pixel. Does nothing (and draws no random numbers) if lambda_dc <= 0.
// ----------------------------------------------------------------------------
void PixelSimulator::AddDarkCurrent() {
  if (cfg_.lambda_dc <= 0.0) return;
  const int N = cfg_.nx * cfg_.ny;
  for (int i = 0; i < N; ++i) {
    const int k = pois_dc_(rng_);
    if (k > 0) qpix_[static_cast<std::size_t>(i)] += static_cast<double>(k);
  }
}

// ----------------------------------------------------------------------------
// PixelSimulator::AddReadoutNoise
//   Add Gaussian noise of width sigma_readout_e to every pixel. Does nothing if sigma_readout_e <= 0.
// ----------------------------------------------------------------------------
void PixelSimulator::AddReadoutNoise() {
  if (cfg_.sigma_readout_e <= 0.0) return;
  const int N = cfg_.nx * cfg_.ny;
  for (int i = 0; i < N; ++i) {
    qpix_[static_cast<std::size_t>(i)] += cfg_.sigma_readout_e * gaus_(rng_);
  }
}

// Change the dark-current mean and rebuild the Poisson draw.
void PixelSimulator::SetDarkCurrent(double lambda_dc) {
  cfg_.lambda_dc = lambda_dc;
  pois_dc_ = std::poisson_distribution<int>( (lambda_dc > 0.0) ? lambda_dc : 0.0 );
}

// Change the readout-noise width.
void PixelSimulator::SetReadoutNoise(double sigma_readout_e) {
  cfg_.sigma_readout_e = sigma_readout_e;
}

} // namespace ccdarksens
