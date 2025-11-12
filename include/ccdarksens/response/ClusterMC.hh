#pragma once
#include <memory>
#include <random>
#include <string>
#include <vector>

class TH1D;

namespace ccdarksens {

// Lightweight geometry view (keep in sync with Detector)
struct GeometryLite {
  int    rows = 0;
  int    cols = 0;
  double pixel_size_um = 15.0; // square pixels
  double thickness_um  = 0.0;
  double active_fraction = 1.0;
};

struct DiffusionPars {
  double A_um2 = 803.25;
  double b_umInv = 6.5e-4;
  double alpha = 1.0;
  double beta_per_keV = 0.0; // MVP kept 0 for DM-e low E
};

struct ClusterMCConfig {
  // Sampling
  int     n_events_per_ne = 20000;
  uint64_t rng_seed = 987654321ULL;

  // Readout / thresholds
  double  sigma_readout_e = 0.16; // per pixel, RMS in e-
  double  Qmin_e = 0.3;           // detection threshold (per pixel)
  double  Qmax_e = 4.0;           // upper guard (reject bright tails if needed)

  // Diffusion / geometry
  DiffusionPars diffusion;
  int rows_bin = 1;
  int cols_bin = 1;

  // Dark current pileup (optional)
  bool   pileup_with_dc = false;
  double lambda_dc_per_pix_per_exposure = 0.0; // mean e- / pixel / exposure

  // Pixel window for collecting a cluster "strip" around the impact point
  int    half_window_pix = 2; // simulate +/- this many pixels in each axis

  // If true, accept MN or MNL style patterns (1-pix over Qmin, or 2-3 adjacents)
  bool   enable_MN   = true;
  bool   enable_MNL  = true;
};

class ClusterMC {
public:
  ClusterMC(const ClusterMCConfig& cfg, const GeometryLite& geo);

  // Precompute ε_MC(n_e) on integer bins [ne_min, ne_max]
  // Name: "eps_mc"; Title: "Cluster-MC efficiency; n_{e}; #epsilon"
  std::unique_ptr<TH1D> PrecomputeEpsilon(int ne_min, int ne_max);

private:
  ClusterMCConfig cfg_;
  GeometryLite    geo_;

  std::mt19937_64 rng_;
  std::uniform_real_distribution<double> uni_{0.0,1.0};
  std::normal_distribution<double> gaus_{0.0, 1.0};

  // Depth sampling (uniform 0..thickness for MVP)
  double sample_depth_um();

  // Lateral diffusion sigma at depth z (and weak E dependence via beta if used)
  double sigma_xy_um(double z_um, double E_eV) const;

  // Simulate one event with true n_e and reference energy (only for beta ≠ 0)
  // Returns true if the pattern passes MN/MNL after DC pileup and readout noise.
  bool simulate_event(int n_true, double E_eV);

  // Helpers
  void add_dark_current(std::vector<double>& qpix) ; // Poisson DC per pixel
  void add_readout_noise(std::vector<double>& qpix);      // Gaussian per pixel
  bool pass_MN_MNL(const std::vector<double>& qpix) const;

  // Deposit n_true electrons around (0,0) with Gaussian σ, pixelize on a local grid
  void deposit_electrons(int n_true, double sigma_um, std::vector<double>& qpix);

  // Build a local pixel grid with binning; returns total simulated pixels
  int local_grid_size_1d() const; // = (2*half_window_pix+1)/bin
};

} // namespace ccdarksens
