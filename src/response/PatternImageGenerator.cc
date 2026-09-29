// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PatternImageGenerator.cc -- Generates binned 2D charge images (and ideal
//  3×5 clusters) from deposited electrons, diffusion, readout noise, and
//  optional dark current.
// ===========================================================================

#include "ccdarksens/response/PatternImageGenerator.hh"
#include "ccdarksens/response/ChargeTransport.hh"

#include <TH2D.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ccdarksens {

// ----------------------------------------------------------------------------
// PatternImageGenerator::PatternImageGenerator
//   I check the configuration: with a raw size, the binning factors must be
//   positive and divide the raw size; otherwise every dimension must be positive.
//   The charge-transport pointer must not be null. Throws std::runtime_error.
// ----------------------------------------------------------------------------
PatternImageGenerator::PatternImageGenerator(const PatternImageConfig& cfg,
                                             std::shared_ptr<ChargeTransport> ct)
  : cfg_(cfg),
    ct_(std::move(ct)),
    rng_(cfg_.rng_seed)
{
  const bool use_raw = (cfg_.raw_rows > 0 && cfg_.raw_cols > 0);
  if (use_raw) {
    if (cfg_.row_binning <= 0 || cfg_.col_binning <= 0)
      throw std::runtime_error("PatternImageGenerator: row_binning and col_binning must be > 0 when using raw size");
    if (cfg_.raw_rows % cfg_.row_binning != 0 || cfg_.raw_cols % cfg_.col_binning != 0)
      throw std::runtime_error("PatternImageGenerator: raw_rows/raw_cols must be divisible by row_binning/col_binning");
  } else {
    if (cfg_.nrows_binned <= 0 || cfg_.ncols <= 0 || cfg_.row_binning <= 0 || cfg_.col_binning <= 0)
      throw std::runtime_error("PatternImageGenerator: nrows_binned, ncols, row_binning, col_binning must be > 0");
  }
  if (!ct_)
    throw std::runtime_error("PatternImageGenerator: ChargeTransport is null");
}

// ----------------------------------------------------------------------------
// PatternImageGenerator::NrowsBinned
//   Number of binned rows: raw_rows / row_binning if the raw size is set, otherwise nrows_binned.
// ----------------------------------------------------------------------------
int PatternImageGenerator::NrowsBinned() const
{
  if (cfg_.raw_rows > 0 && cfg_.raw_cols > 0)
    return cfg_.raw_rows / cfg_.row_binning;
  return cfg_.nrows_binned;
}

// ----------------------------------------------------------------------------
// PatternImageGenerator::Ncols
//   Number of binned columns: raw_cols / col_binning if the raw size is set, otherwise ncols.
// ----------------------------------------------------------------------------
int PatternImageGenerator::Ncols() const
{
  if (cfg_.raw_rows > 0 && cfg_.raw_cols > 0)
    return cfg_.raw_cols / cfg_.col_binning;
  return cfg_.ncols;
}

// ----------------------------------------------------------------------------
// PatternImageGenerator::GenerateImage
//   One simulated event as a binned image [row][col]. I sample a depth, get
//   sigma_xy(z, E), place a Gaussian cloud of n_e electrons at the image centre (or a
//   random position if randomize_center is set), drop each electron in its raw
//   pixel, sum row_binning raw rows and col_binning raw columns, add Gaussian
//   readout noise to every binned pixel, and optionally Poisson dark current.
// ----------------------------------------------------------------------------
std::vector<std::vector<double>> PatternImageGenerator::GenerateImage(int n_e, double Ee_eV)
{
  const int ny_raw = (cfg_.raw_rows > 0 && cfg_.raw_cols > 0)
    ? cfg_.raw_rows
    : (cfg_.nrows_binned * cfg_.row_binning);
  const int nx_raw = (cfg_.raw_rows > 0 && cfg_.raw_cols > 0)
    ? cfg_.raw_cols
    : (cfg_.ncols * cfg_.col_binning);
  const int nrows_binned = NrowsBinned();
  const int ncols = Ncols();
  const double pitch = cfg_.pixel_size_um;

  // Raw grid (ny_raw × nx_raw), row-major as [iy][ix]
  std::vector<std::vector<double>> raw(static_cast<std::size_t>(ny_raw),
                                       std::vector<double>(static_cast<std::size_t>(nx_raw), 0.0));

  if (n_e > 0) {
    const double z_um = ct_->SampleDepthUm();
    const double sigma_xy_um = ct_->SigmaXYUm(z_um, Ee_eV);

    // Cloud centre: either randomised per-trial (matching Python x0/y0 = uniform)
    // or fixed at the image midpoint.
    double cx_um, cy_um;
    if (cfg_.randomize_center) {
      // Match Python generate_image_E: x0=uniform(x_min+margin, x_max-margin),
      // y0=uniform(y_min+margin, y_max-margin) where y is within the middle
      // binned row (raw rows [row_binning, 2*row_binning)).
      const double m  = cfg_.center_margin_pix;
      const double lo_x = m * pitch;
      const double hi_x = (nx_raw - m) * pitch;
      // y: restrict to middle binned row (Python: bining+y, y in [margin, row_binning-margin])
      const double lo_y = (cfg_.row_binning + m) * pitch;
      const double hi_y = (2.0 * cfg_.row_binning - m) * pitch;
      cx_um = std::uniform_real_distribution<double>(lo_x, hi_x)(rng_);
      cy_um = std::uniform_real_distribution<double>(lo_y, hi_y)(rng_);
    } else {
      cx_um = (nx_raw / 2.0) * pitch;
      cy_um = (ny_raw / 2.0) * pitch;
    }

    std::vector<double> xs_um, ys_um;
    ct_->SampleCloudXY(cx_um, cy_um, sigma_xy_um, n_e, xs_um, ys_um);

    for (int i = 0; i < n_e && i < (int)xs_um.size(); ++i) {
      int ix = static_cast<int>(xs_um[static_cast<std::size_t>(i)] / pitch);
      int iy = static_cast<int>(ys_um[static_cast<std::size_t>(i)] / pitch);
      ix = std::max(0, std::min(nx_raw - 1, ix));
      iy = std::max(0, std::min(ny_raw - 1, iy));
      raw[static_cast<std::size_t>(iy)][static_cast<std::size_t>(ix)] += 1.0;
    }
  }

  // Row binning: sum every row_binning rows → (nrows_binned × nx_raw)
  std::vector<std::vector<double>> after_row(static_cast<std::size_t>(nrows_binned),
                                             std::vector<double>(static_cast<std::size_t>(nx_raw), 0.0));
  for (int r = 0; r < nrows_binned; ++r) {
    for (int j = 0; j < cfg_.row_binning; ++j) {
      const int iy = r * cfg_.row_binning + j;
      for (int c = 0; c < nx_raw; ++c)
        after_row[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] += raw[static_cast<std::size_t>(iy)][static_cast<std::size_t>(c)];
    }
  }

  // Column binning: sum every col_binning columns → (nrows_binned × ncols)
  std::vector<std::vector<double>> out(static_cast<std::size_t>(nrows_binned),
                                        std::vector<double>(static_cast<std::size_t>(ncols), 0.0));
  for (int r = 0; r < nrows_binned; ++r) {
    for (int c = 0; c < ncols; ++c) {
      for (int k = 0; k < cfg_.col_binning; ++k)
        out[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] += after_row[static_cast<std::size_t>(r)][static_cast<std::size_t>(c * cfg_.col_binning + k)];
    }
  }

  // Readout noise
  for (int r = 0; r < nrows_binned; ++r) {
    for (int c = 0; c < ncols; ++c) {
      out[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] += cfg_.sigma_readout_e * static_cast<double>(gaus_(rng_));
    }
  }
  // Optional dark current
  if (cfg_.include_dark_current && cfg_.lambda_dc > 0.0) {
    std::poisson_distribution<int> pois(cfg_.lambda_dc);
    for (int r = 0; r < nrows_binned; ++r) {
      for (int col = 0; col < ncols; ++col) {
        const int k = pois(rng_);
        if (k > 0) out[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)] += static_cast<double>(k);
      }
    }
  }

  return out;
}

// ----------------------------------------------------------------------------
// PatternImageGenerator::GenerateImageToTH2D
//   Same as GenerateImage, written into a ROOT TH2D with the first row at the top (a null pointer is ignored).
// ----------------------------------------------------------------------------
void PatternImageGenerator::GenerateImageToTH2D(int n_e, double Ee_eV, TH2D* h2)
{
  if (!h2) return;
  auto img = GenerateImage(n_e, Ee_eV);
  const int ncols = Ncols();
  const int ny = NrowsBinned();
  for (int r = 0; r < ny; ++r) {
    for (int c = 0; c < ncols; ++c) {
      h2->SetBinContent(c + 1, ny - r, img[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]);
    }
  }
}

// ----------------------------------------------------------------------------
// PatternImageGenerator::SimulateCluster
//   An ideal 3x5 cluster whose middle row is [0, b, c, d, 0], with Gaussian readout
//   noise added to every pixel and rounded to five decimals.
// ----------------------------------------------------------------------------
std::vector<std::vector<double>> PatternImageGenerator::SimulateCluster(double b, double c, double d)
{
  // Match notebook: 3x5 array, middle row [0, b, c, d, 0], add N(0, sigma) then round to 5 decimals
  std::vector<std::vector<double>> cl(3, std::vector<double>(5, 0.0));
  cl[1][1] = b;
  cl[1][2] = c;
  cl[1][3] = d;
  const double sigma = cfg_.sigma_readout_e;
  for (int ri = 0; ri < 3; ++ri) {
    for (int ci = 0; ci < 5; ++ci) {
      double v = cl[static_cast<std::size_t>(ri)][static_cast<std::size_t>(ci)] + sigma * static_cast<double>(gaus_(rng_));
      cl[static_cast<std::size_t>(ri)][static_cast<std::size_t>(ci)] = std::round(v * 100000.0) / 100000.0;  // 5 decimals like notebook
    }
  }
  return cl;
}

} // namespace ccdarksens
