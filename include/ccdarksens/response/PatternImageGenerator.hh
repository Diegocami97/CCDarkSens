// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PatternImageGenerator.hh -- Header for 2D binned image generation config
//  and ChargeTransport-backed image simulation.
// ===========================================================================

#pragma once

#include <memory>
#include <random>
#include <vector>

class TH2D;

namespace ccdarksens {

class ChargeTransport;

/// Config for 2D binned image (notebook-style generate_image_E).
/// Image size: either from detector (raw_rows, raw_cols) or from binned dimensions (nrows_binned, ncols).
// ----------------------------------------------------------------------------
// PatternImageConfig
//   Geometry, noise and randomization settings of the 2D image generator (see the comment above).
// ----------------------------------------------------------------------------
struct PatternImageConfig {
  /// When both > 0: raw image size from detector (rows, cols). Binned size = raw/row_binning, raw/col_binning.
  int raw_rows = 0;  // raw image height from the detector [pixels] (0 = not set)
  int raw_cols = 0;  // raw image width from the detector [pixels] (0 = not set)
  /// When raw_rows/raw_cols not set: binned dimensions (raw then = nrows_binned*row_binning, ncols*col_binning).
  int nrows_binned  = 3;     ///< rows after row binning (ignored if raw_rows > 0)
  int ncols         = 50;    ///< columns after column binning (ignored if raw_cols > 0)
  int row_binning   = 100;   ///< raw rows per binned row
  int col_binning   = 1;     ///< raw columns per binned column (1 = no col bin)
  double pixel_size_um   = 15.0;  // pixel pitch [um]
  double sigma_readout_e = 0.21;  // readout noise per binned pixel [e-]
  double lambda_dc      = 0.0;   ///< dark current (0 = off)
  uint64_t rng_seed     = 987654321ULL;  // RNG seed
  bool include_dark_current = false;  // add Poisson dark current (needs lambda_dc > 0)

  /// If true, the cloud center (cx, cy) is sampled from a uniform distribution
  /// within the image interior on each GenerateImage() call, matching the Python
  /// reference:  x0 = uniform(x_min+15, x_max-15),  y0 = uniform(y_min+15, y_max-15).
  /// The margin in raw-pixel units is given by center_margin_pix.
  /// Default false (backward-compatible: fixed center at image midpoint).
  bool randomize_center    = false;
  double center_margin_pix = 15.0; ///< minimum raw-pixel distance from edge when randomising
};

/**
 * Generates a single 2D binned image (e.g. 3×50) matching the notebook's
 * generate_image_E(ne, Ee): diffusion from random z and σ_xy(z,E), row binning,
 * readout noise, optional dark current.
 */
class PatternImageGenerator {
public:
  // Constructor: settings and the charge-transport model.
  PatternImageGenerator(const PatternImageConfig& cfg,
                        std::shared_ptr<ChargeTransport> ct);

  /// Generate one 2D image: (nrows_binned × ncols). Row-major: [row][col].
  std::vector<std::vector<double>> GenerateImage(int n_e, double Ee_eV);

  /// Same as GenerateImage but fill a ROOT TH2D (bins 1..ncols, 1..nrows_binned).
  void GenerateImageToTH2D(int n_e, double Ee_eV, TH2D* h2);

  /// Simulate ideal 3×5 cluster with middle row [1,1]=b, [1,2]=c, [1,3]=d; add readout noise.
  /// Returns 3×5 row-major (rows then cols). Notebook: simulate_cluster(b,c,d).
  std::vector<std::vector<double>> SimulateCluster(double b, double c, double d);

  /// Effective binned dimensions (from raw/row_binning, raw/col_binning when raw set; else from config).
  int NrowsBinned() const;
  int Ncols() const;
  int RowBinning() const { return cfg_.row_binning; }
  int ColBinning() const { return cfg_.col_binning; }

private:
  PatternImageConfig cfg_;  // settings
  std::shared_ptr<ChargeTransport> ct_;  // depth / diffusion sampler
  std::mt19937_64 rng_;  // random-number generator
  std::normal_distribution<double> gaus_{0.0, 1.0};  // unit Gaussian for the readout noise
};

} // namespace ccdarksens
