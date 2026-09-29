// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ClusterFitMC.cc -- Forward-sim + fit + kernel builder. See
//  docs/ClusterFitMC_Design.md.
// ===========================================================================

#include "ccdarksens/response/ClusterFitMC.hh"

#include <algorithm>
#include <cmath>

namespace ccdarksens {

// Constructor: keep the settings and the shared charge-transport model, seed the RNG.
ClusterFitMC::ClusterFitMC(const ClusterFitMCConfig& cfg, std::shared_ptr<ChargeTransport> ct)
    : cfg_(cfg), ct_(std::move(ct)), rng_(cfg.rng_seed) {}

// ----------------------------------------------------------------------------
// ClusterFitMC::BuildKernel
//   For every true energy E: draw n_e from a Gaussian with mean E/eh_pair and
//   variance fano*mean, sample a depth and the diffusion width, deposit the
//   electrons, add dark current (if configured) and readout noise, and fit.
//   An event is accepted if Delta LL beats the calibrated cut and the fitted width
//   is inside the fiducial window; it is then filled into the E_reco bin of
//   I_hat * eh_pair (events reconstructed outside the requested range count as
//   missed). Each row is divided by ne_trials_per_point, so K[i][j] is the
//   probability of E_true[i] -> E_reco bin j and the row sum is the efficiency.
// ----------------------------------------------------------------------------
KernelMatrix ClusterFitMC::BuildKernel(const std::vector<double>& Etrue_grid_eV,
                                        const std::vector<double>& Ereco_edges_eV) const {
  KernelMatrix result;
  result.Etrue_grid_eV = Etrue_grid_eV;
  result.Ereco_edges_eV = Ereco_edges_eV;
  const int n_ereco_bins = static_cast<int>(Ereco_edges_eV.size()) - 1;
  result.K.assign(Etrue_grid_eV.size(), std::vector<double>(static_cast<std::size_t>(n_ereco_bins), 0.0));
  result.efficiency.assign(Etrue_grid_eV.size(), 0.0);

  PixelSimulator pixSim(cfg_.pix_cfg);
  ClusterFitEngine engine(cfg_.fit_cfg);

  for (std::size_t iE = 0; iE < Etrue_grid_eV.size(); ++iE) {
    const double E_true = Etrue_grid_eV[iE];
    const double mean_ne = E_true / cfg_.eh_pair_eV;
    const double sigma_ne = std::sqrt(std::max(0.0, cfg_.fano_factor * mean_ne));

    int n_accept = 0;
    for (int trial = 0; trial < cfg_.ne_trials_per_point; ++trial) {
      const int n_e = std::max(0, static_cast<int>(std::lround(mean_ne + sigma_ne * gaus_(rng_))));
      if (n_e == 0) continue;  // nothing deposited -- automatically a miss

      pixSim.Reset();
      const double z_um = ct_->SampleDepthUm();
      const double sigma_xy_um = ct_->SigmaXYUm(z_um, E_true);
      std::vector<double> xs_um, ys_um;
      ct_->SampleCloudXY(0.0, 0.0, sigma_xy_um, n_e, xs_um, ys_um);
      for (std::size_t i = 0; i < xs_um.size(); ++i) pixSim.DepositElectron(xs_um[i], ys_um[i]);
      pixSim.AddDarkCurrent();   // no-op (no RNG draws) unless pix_cfg.lambda_dc > 0
      pixSim.AddReadoutNoise();

      const auto r = engine.Fit(pixSim.PixelCharges(), pixSim.NX(), pixSim.NY());
      if (r.delta_ll > cfg_.delta_ll_cut) continue;  // didn't beat the noise-tail detection cut
      if (r.sigma_xy_px < cfg_.sigma_xy_fid_min_px || r.sigma_xy_px > cfg_.sigma_xy_fid_max_px) continue;  // fiducial cut

      const double E_reco = r.I_hat_e * cfg_.eh_pair_eV;
      const auto it = std::upper_bound(Ereco_edges_eV.begin(), Ereco_edges_eV.end(), E_reco);
      if (it == Ereco_edges_eV.begin() || it == Ereco_edges_eV.end()) continue;  // reconstructed outside the requested Ereco range
      const int bin = static_cast<int>(it - Ereco_edges_eV.begin()) - 1;

      result.K[iE][static_cast<std::size_t>(bin)] += 1.0;
      ++n_accept;
    }

    const double norm = 1.0 / static_cast<double>(cfg_.ne_trials_per_point);
    for (double& v : result.K[iE]) v *= norm;
    result.efficiency[iE] = static_cast<double>(n_accept) * norm;
  }

  return result;
}

}  // namespace ccdarksens
