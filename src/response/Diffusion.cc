#include "ccdarksens/response/Diffusion.hh"
#include <TH1D.h>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace ccdarksens {

double Diffusion::sigma_xy_um_(double z_um, double Ee_eV) const {
  // Ee in keV for the linear term
  const double Ee_keV = Ee_eV * 1e-3;
  const double inside = 1.0 - b_umInv_ * z_um;
  if (inside <= 0.0) return std::numeric_limits<double>::quiet_NaN();
  const double geom = std::sqrt(-A_um2_ * std::log(inside));
  return geom * (alpha_ + beta_per_keV_ * Ee_keV);
}

double Diffusion::ComputeSigmaElectrons(double Ee_eV, std::size_t nz_steps) const {
  if (thickness_um_ <= 0.0 || kappa_e_per_um_ < 0.0)
    return std::sqrt(std::max(0.0, sigma_readout_e_ * sigma_readout_e_));

  // z-average sigma_xy over [0, thickness], simple trapezoidal integration
  const double dz = thickness_um_ / static_cast<double>(nz_steps);
  double acc = 0.0;
  for (std::size_t i=0; i<=nz_steps; ++i) {
    const double z = dz * static_cast<double>(i);
    double w = (i==0 || i==nz_steps) ? 0.5 : 1.0;
    const double sxy = sigma_xy_um_(z, Ee_eV);
    if (std::isfinite(sxy)) acc += w * sxy;
  }
  const double mean_sigma_xy_um = acc * dz / thickness_um_;

  // map to electrons and add readout noise in quadrature
  const double sigma_geom_e = kappa_e_per_um_ * mean_sigma_xy_um;
  const double sigma_e2 = sigma_geom_e * sigma_geom_e + sigma_readout_e_ * sigma_readout_e_;
  return std::sqrt(std::max(0.0, sigma_e2));
}

double Diffusion::normal_cdf_(double x) {
  // Phi(x) = 0.5 * [1 + erf(x / sqrt(2))]
  return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

double Diffusion::interval_prob_(double mu, double sigma, double a, double b) {
  if (!(sigma > 0.0)) return (a <= mu && mu < b) ? 1.0 : 0.0;
  const double z1 = (a - mu) / sigma;
  const double z2 = (b - mu) / sigma;
  return std::max(0.0, normal_cdf_(z2) - normal_cdf_(z1));
}

void Diffusion::Apply(TH1D& h_ne, double Ee_eV) const {
  const int nb = h_ne.GetNbinsX();
  if (nb <= 0) return;

  const double sigma_e = ComputeSigmaElectrons(Ee_eV);
  if (!(sigma_e > 0.0)) return; // no-op if sigma is zero

  // integer binning assumed: bins centered at ..., n-0.5 .. n+0.5 ...
  const double x_min = h_ne.GetXaxis()->GetXmin();
  const double x_max = h_ne.GetXaxis()->GetXmax();

  // Copy source contents
  std::vector<double> src(nb+1, 0.0); // 1-indexed for ROOT bins
  for (int i=1; i<=nb; ++i) src[i] = h_ne.GetBinContent(i);

  // Destination accumulator
  std::vector<double> dst(nb+1, 0.0);

  // For each source bin centered at n (integer), distribute to neighbors
  for (int i=1; i<=nb; ++i) {
    const double n = h_ne.GetBinCenter(i);        // should be integer
    const double val = src[i];
    if (val == 0.0) continue;

    // truncate kernel at ±5σ for speed
    const int q_min = std::max(1, h_ne.FindBin(n - 5.0*sigma_e));
    const int q_max = std::min(nb, h_ne.FindBin(n + 5.0*sigma_e));
    for (int q=q_min; q<=q_max; ++q) {
      const double q_center = h_ne.GetBinCenter(q);
      const double a = std::max(x_min, q_center - 0.5);
      const double b = std::min(x_max, q_center + 0.5);
      const double p = interval_prob_(/*mu=*/n, /*sigma=*/sigma_e, a, b);
      dst[q] += val * p;
    }
  }

  // Write back
  for (int i=1; i<=nb; ++i) h_ne.SetBinContent(i, dst[i]);
}

} // namespace ccdarksens
