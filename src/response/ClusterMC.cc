#include "ccdarksens/response/ClusterMC.hh"
#include <TH1D.h>
#include <algorithm>
#include <cmath>

namespace ccdarksens {

ClusterMC::ClusterMC(const ClusterMCConfig& cfg, const GeometryLite& geo)
: cfg_(cfg), geo_(geo), rng_(cfg.rng_seed)
{
  if (cfg_.rows_bin < 1) cfg_.rows_bin = 1;
  if (cfg_.cols_bin < 1) cfg_.cols_bin = 1;
  if (cfg_.half_window_pix < 1) cfg_.half_window_pix = 1;
}

double ClusterMC::sample_depth_um() {
  // Uniform in [0, thickness] for MVP
  return uni_(rng_) * std::max(0.0, geo_.thickness_um);
}

double ClusterMC::sigma_xy_um(double z_um, double E_eV) const {
  // E_ref in keV for beta scaling (if used)
  const double E_keV = std::max(0.0, E_eV) * 1e-3;
  const double term  = std::max(1.0 - cfg_.diffusion.b_umInv * z_um, 1e-9);
  const double geom  = std::sqrt(-cfg_.diffusion.A_um2 * std::log(term));
  const double scale = (cfg_.diffusion.alpha + cfg_.diffusion.beta_per_keV * E_keV);
  return std::max(0.0, geom * scale);
}

int ClusterMC::local_grid_size_1d() const {
  const int full = 2*cfg_.half_window_pix + 1;
  const int rb = std::max(1, cfg_.rows_bin);
  // Treat square binning (rows_bin = cols_bin) for the local patch
  const int binned = (full + rb - 1) / rb; // ceil
  return binned;
}

void ClusterMC::deposit_electrons(int n_true, double sigma_um, std::vector<double>& qpix)
{
  // Local grid centered at (0,0). Indexing: [iy * Nx + ix].
  const int Nx = local_grid_size_1d();
  const int Ny = Nx; // square patch
  qpix.assign(Nx*Ny, 0.0);

  // Pixel pitch after binning
  const double pitch_um = geo_.pixel_size_um * std::max(1, cfg_.rows_bin);

  // For each electron: sample (x,y) ~ N(0, sigma), map to nearest pixel (MVP)
  // (You can replace nearest-pixel with Gaussian integration over pixels later.)
  for (int i = 0; i < n_true; ++i) {
    const double x_um = sigma_um * gaus_(rng_);
    const double y_um = sigma_um * gaus_(rng_);

    int ix0 = static_cast<int>(std::lround(x_um / pitch_um));
    int iy0 = static_cast<int>(std::lround(y_um / pitch_um));

    // Shift to local grid coordinates centered at 0
    ix0 += Nx/2;  // center pixel index
    iy0 += Ny/2;

    if (ix0 >= 0 && ix0 < Nx && iy0 >= 0 && iy0 < Ny) {
      qpix[iy0 * Nx + ix0] += 1.0; // 1 electron unit
    }
  }
}

void ClusterMC::add_dark_current(std::vector<double>& qpix) {
  if (!cfg_.pileup_with_dc) return;
  const int N = static_cast<int>(qpix.size());
  if (cfg_.lambda_dc_per_pix_per_exposure <= 0.0) return;

  // Poisson(dc) electrons per pixel; approximate by lambda if tiny? Use Poisson exactly.
  std::poisson_distribution<int> pois(cfg_.lambda_dc_per_pix_per_exposure);
  for (int i = 0; i < N; ++i) {
    const int k = pois(rng_);
    if (k > 0) qpix[i] += static_cast<double>(k);
  }
}

void ClusterMC::add_readout_noise(std::vector<double>& qpix) {
  if (cfg_.sigma_readout_e <= 0.0) return;
  const int N = static_cast<int>(qpix.size());
  for (int i = 0; i < N; ++i) {
    qpix[i] += cfg_.sigma_readout_e * gaus_(rng_);
  }
}

bool ClusterMC::pass_MN_MNL(const std::vector<double>& qpix) const {
  // Apply simple MN / MNL classification:
  // - MN: at least one pixel with Q in [Qmin,Qmax]
  // - MNL: at least two adjacent pixels with Q in [Qmin,Qmax]
  const double qmin = cfg_.Qmin_e;
  const double qmax = cfg_.Qmax_e;
  const int Nx = local_grid_size_1d();
  const int Ny = Nx;

  auto in_window = [&](double q){ return (q >= qmin && q <= qmax); };

  bool mn_ok = false;
  bool mnl_ok = false;

  for (int iy = 0; iy < Ny; ++iy) {
    for (int ix = 0; ix < Nx; ++ix) {
      const double q0 = qpix[iy * Nx + ix];
      if (cfg_.enable_MN && in_window(q0)) mn_ok = true;

      if (cfg_.enable_MNL && in_window(q0)) {
        // check 4-neighbors for adjacency
        const int dx[4] = {+1,-1,0,0};
        const int dy[4] = {0,0,+1,-1};
        for (int k=0;k<4;++k){
          int jx = ix + dx[k];
          int jy = iy + dy[k];
          if (jx<0 || jx>=Nx || jy<0 || jy>=Ny) continue;
          const double q1 = qpix[jy * Nx + jx];
          if (in_window(q1)) { mnl_ok = true; break; }
        }
      }
    }
  }

  return ( (cfg_.enable_MN  && mn_ok) ||
           (cfg_.enable_MNL && mnl_ok) );
}

bool ClusterMC::simulate_event(int n_true, double E_eV)
{
  // 1) sample mean sigma from depth (single z per event for MVP)
  const double z_um = sample_depth_um();
  const double sig_um = sigma_xy_um(z_um, E_eV);

  // 2) deposit electrons in local pixel patch
  std::vector<double> qpix;
  deposit_electrons(n_true, sig_um, qpix);

  // 3) optional DC pileup
  add_dark_current(qpix);

  // 4) readout noise
  add_readout_noise(qpix);

  // 5) MN/MNL decision
  return pass_MN_MNL(qpix);
}

std::unique_ptr<TH1D> ClusterMC::PrecomputeEpsilon(int ne_min, int ne_max)
{
  const int nbin = ne_max - ne_min + 1;
  std::vector<double> edges(nbin + 1);
  for (int i = 0; i <= nbin; ++i) edges[i] = (ne_min - 0.5) + i;

  auto h = std::make_unique<TH1D>("eps_mc",
                                  "Cluster-MC efficiency; n_{e}; #epsilon",
                                  nbin, edges.data());
  h->Sumw2(false);

  // MVP reference energy for sigma(beta) if ever used (kept 0 by default)
  const double E_ref_eV = 50.0;

  for (int ne = ne_min; ne <= ne_max; ++ne) {
    const int trials = std::max(1, cfg_.n_events_per_ne);
    int accepted = 0;

    for (int i = 0; i < trials; ++i) {
      if (simulate_event(ne, E_ref_eV)) ++accepted;
    }

    const double eps = (trials > 0) ? (static_cast<double>(accepted) / trials) : 0.0;
    const int bin = h->FindBin(ne);
    h->SetBinContent(bin, std::clamp(eps, 0.0, 1.0));
  }

  return h;
}

} // namespace ccdarksens
