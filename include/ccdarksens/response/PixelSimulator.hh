// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PixelSimulator.hh -- Header for pixel-grid charge deposition, dark
//  current, and readout-noise simulation settings.
// ===========================================================================

#pragma once

#include <vector>
#include <random>

namespace ccdarksens {

/// Simulation mode:
///  - LocalPatch: small pixel patch around interaction (EfficiencyMC / efficiencies)
///  - FullCCD: full sensor for validation / global studies
enum class PixelSimMode {
  LocalPatch,
  FullCCD,
  RowSegment
};

// ----------------------------------------------------------------------------
// PixelSimulatorConfig
//   Pixel-patch geometry, noise settings and RNG seed for PixelSimulator.
// ----------------------------------------------------------------------------
struct PixelSimulatorConfig {
  PixelSimMode mode         = PixelSimMode::LocalPatch;  // which kind of patch to simulate
  int          nx           = 7;       ///< pixels in x (patch or full sensor)
  int          ny           = 7;       ///< pixels in y
  double       pixel_size_um= 15.0;    ///< pixel pitch (μm)
  double       lambda_dc    = 0.0;     ///< DC mean e-/pixel/exposure
  double       sigma_readout_e = 0.0;  ///< Gaussian readout noise [e-]
  std::uint64_t rng_seed    = 1234567;  // RNG seed

  /// When true, DepositElectron() forces every electron into row 0
  /// unconditionally (no y bounds-check) instead of binning/rejecting by y --
  /// the DAMIC 1x100 readout mode's hardware row-summing (100 physical rows
  /// combined before the single readout; PhysRevD.94.082006 Sec. VI), not a
  /// literal ny=1 window (which would silently drop most charge instead of
  /// summing it). ny should be 1 in this mode; x-binning is unaffected.
  bool         collapse_y = false;
};

/// PixelSimulator: take continuous electron positions in μm and convert them
/// into a discrete pixel charge map, with optional DC and readout noise.
// ----------------------------------------------------------------------------
// PixelSimulator
//   Pixel charge map: deposit electrons, then optionally add Poisson dark current and Gaussian readout noise.
// ----------------------------------------------------------------------------
class PixelSimulator {
public:
  explicit PixelSimulator(const PixelSimulatorConfig& cfg);

  /// Reset all pixel charges to zero.
  void Reset();

  /// Deposit a single electron at (x_um, y_um).
  /// Interpretation of coordinates depends on mode (local vs full).
  void DepositElectron(double x_um, double y_um);

  /// Optionally add Poisson dark current to each pixel.
  void AddDarkCurrent();

  /// Optionally add Gaussian readout noise to each pixel.
  void AddReadoutNoise();

  /// Accessors
  int NX() const noexcept { return cfg_.nx; }
  int NY() const noexcept { return cfg_.ny; }

  /// Flat, row-major pixel array: size = NX()*NY().
  const std::vector<double>& PixelCharges() const noexcept { return qpix_; }

  /// Mutable access if needed by higher-level logic.
  std::vector<double>& PixelCharges() noexcept { return qpix_; }

  /// Dynamic adjustment of DC and noise.
  void SetDarkCurrent(double lambda_dc);
  void SetReadoutNoise(double sigma_readout_e);

private:
  PixelSimulatorConfig cfg_;  // settings
  std::vector<double> qpix_;  // row-major pixel charges [e-]
  

  std::mt19937_64 rng_;  // random-number generator
  std::poisson_distribution<int> pois_dc_;  // Poisson draw for the dark current
  std::normal_distribution<double> gaus_{0.0, 1.0};  // unit Gaussian draw for the readout noise

  // True if pixel (ix, iy) lies inside the patch.
  bool InBounds_(int ix, int iy) const;
};

} // namespace ccdarksens
