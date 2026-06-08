// ============================================================================
//  CCDarkSens — PCDCalculator
//  Folds signal/background n_e spectra through P(q|n_e) into q-space rates and builds P(n_obs|n_true) reconstruction kernels from PCD tables.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

// #include "ccdarksens/response/PCDCalculator.hh"

// #include <TH1D.h>
// #include <cmath>
// #include <stdexcept>
// #include <algorithm>

// namespace ccdarksens {

// namespace {

// // Gaussian CDF
// inline double normal_cdf(double x) {
//   return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
// }

// // Probability that a Gaussian N(mu, sigma) falls in [a,b]
// inline double gaussian_interval_prob(double mu, double sigma,
//                                      double a, double b)
// {
//   if (!(sigma > 0.0)) {
//     return (a <= mu && mu < b) ? 1.0 : 0.0;
//   }
//   const double z1 = (a - mu) / sigma;
//   const double z2 = (b - mu) / sigma;
//   return std::max(0.0, normal_cdf(z2) - normal_cdf(z1));
// }

// } // namespace

// PCDCalculator::PCDCalculator(const PCDCalculatorConfig& cfg)
//   : cfg_(cfg)
// {
// }

// void PCDCalculator::SetConfig(const PCDCalculatorConfig& cfg) {
//   cfg_ = cfg;
// }

// const PCDCalculatorConfig& PCDCalculator::GetConfig() const noexcept {
//   return cfg_;
// }

// void PCDCalculator::SetSignalSpectrum(const TH1D& h_ne) {
//   // Interpret h_ne as expected *signal* counts per n_e bin over the whole detector.
//   // We convert to "per pixel" expectation by dividing by n_pixels elsewhere,
//   // so here we just store shape proportional to counts.
//   const int nb = h_ne.GetNbinsX();
//   if (nb <= 0) {
//     p_signal_ne_.clear();
//     return;
//   }

//   // Determine max n_e from axis (assumes integer bin centers).
//   int max_ne = static_cast<int>(std::round(h_ne.GetXaxis()->GetXmax()));
//   if (max_ne < 0) max_ne = 0;
//   p_signal_ne_.assign(static_cast<std::size_t>(max_ne + 1), 0.0);

//   for (int i = 1; i <= nb; ++i) {
//     const double n_center = h_ne.GetBinCenter(i);
//     const int    n        = static_cast<int>(std::lround(n_center));
//     if (n < 0 || n > max_ne) continue;
//     p_signal_ne_[static_cast<std::size_t>(n)] += h_ne.GetBinContent(i);
//   }
// }

// void PCDCalculator::SetSignalSpectrum(const std::vector<double>& p_ne) {
//   p_signal_ne_ = p_ne;
// }

// void PCDCalculator::SetDarkCurrent(double lambda_dc) {
//   cfg_.lambda_dc = lambda_dc;
// }

// void PCDCalculator::SetReadoutNoise(double sigma_readout_e) {
//   cfg_.sigma_readout_e = sigma_readout_e;
// }

// std::unique_ptr<TH1D> PCDCalculator::ComputePCD() const {
//   if (cfg_.q_max < cfg_.q_min) {
//     throw std::invalid_argument("PCDCalculator: q_max < q_min");
//   }

//   const int q_min = cfg_.q_min;
//   const int q_max = cfg_.q_max;
//   const int nb_q  = q_max - q_min + 1;

//   std::vector<double> edges(nb_q + 1);
//   for (int i = 0; i <= nb_q; ++i) {
//     edges[i] = (q_min - 0.5) + i; // bins [q-0.5, q+0.5]
//   }

//   auto h = std::make_unique<TH1D>(cfg_.hist_name.c_str(),
//                                   cfg_.hist_name.c_str(),
//                                   nb_q, edges.data());
//   h->Sumw2();

//   if (p_signal_ne_.empty()) {
//     // No signal => only DC + noise. For now we treat signal-only PCD;
//     // user can extend this for pure-DC backgrounds later.
//     return h;
//   }

//   // We interpret p_signal_ne_[n] as the *expected number of signal events per pixel*
//   // yielding n electrons, so sum_n p_signal_ne_[n] = mu_sig (per pixel).
//   const std::size_t n_max = p_signal_ne_.size() - 1;

//   const double lambda_dc = std::max(0.0, cfg_.lambda_dc);
//   const double sigma     = std::max(0.0, cfg_.sigma_readout_e);

//   // Precompute Poisson DC pmf up to a reasonable k_max.
//   // For small occupancies this can be modest; we choose k_max ~ q_max + 5.
//   const int k_max = std::max(0, q_max + 5);
//   std::vector<double> p_dc(static_cast<std::size_t>(k_max + 1), 0.0);
//   if (lambda_dc > 0.0) {
//     p_dc[0] = std::exp(-lambda_dc);
//     for (int k = 1; k <= k_max; ++k) {
//       p_dc[static_cast<std::size_t>(k)] =
//         p_dc[static_cast<std::size_t>(k - 1)] * lambda_dc / static_cast<double>(k);
//     }
//   } else {
//     p_dc[0] = 1.0;
//   }

//   // Accumulator for PCD counts per pixel.
//   std::vector<double> pcd_per_pixel(static_cast<std::size_t>(nb_q), 0.0);

//   // Triple sum over n_e signal, DC electrons, and q bins with Gaussian readout noise.
//   for (std::size_t n = 0; n <= n_max; ++n) {
//     const double w_sig = p_signal_ne_[n];
//     if (w_sig <= 0.0) continue;

//     for (int k = 0; k <= k_max; ++k) {
//       const double w_dc = p_dc[static_cast<std::size_t>(k)];
//       if (w_dc <= 0.0) continue;

//       const double mean_e = static_cast<double>(n + k); // mean before noise

//       for (int ibin = 0; ibin < nb_q; ++ibin) {
//         const double q_center = static_cast<double>(q_min + ibin);
//         const double a = q_center - 0.5;
//         const double b = q_center + 0.5;

//         const double p_gauss = gaussian_interval_prob(mean_e, sigma, a, b);
//         if (p_gauss <= 0.0) continue;

//         pcd_per_pixel[static_cast<std::size_t>(ibin)] += w_sig * w_dc * p_gauss;
//       }
//     }
//   }

//   // Scale by number of pixels if requested: we then obtain expected counts per bin.
//   const double scale = (cfg_.n_pixels > 0) ? static_cast<double>(cfg_.n_pixels) : 1.0;
//   for (int ibin = 0; ibin < nb_q; ++ibin) {
//     const double val = pcd_per_pixel[static_cast<std::size_t>(ibin)] * scale;
//     h->SetBinContent(ibin + 1, val);
//   }

//   return h;
// }

// } // namespace ccdarksens

#include "ccdarksens/response/PCDCalculator.hh"

#include <stdexcept>
#include <TH1D.h>

namespace ccdarksens {

std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>>
PCDCalculator::FoldSpectra(
    const std::map<int, std::unique_ptr<TH1D>>& pcd_table,
    const std::map<int, double>& Sn,
    const std::map<int, double>& Bn) const
{
  if (pcd_table.empty()) {
    throw std::runtime_error("PCDCalculator::FoldSpectra: empty P(q | n_e) table.");
  }

  // Use the first histogram as a template for binning.
  const TH1D* htemplate = pcd_table.begin()->second.get();
  if (!htemplate) {
    throw std::runtime_error("PCDCalculator::FoldSpectra: null histogram template.");
  }

  int nbins = htemplate->GetNbinsX();
  double qmin = htemplate->GetXaxis()->GetXmin();
  double qmax = htemplate->GetXaxis()->GetXmax();

  auto h_sig = std::make_unique<TH1D>("pcd_signal_q",
                                      "Signal in q-space",
                                      nbins, qmin, qmax);
  auto h_bkg = std::make_unique<TH1D>("pcd_background_q",
                                      "Background in q-space",
                                      nbins, qmin, qmax);

  h_sig->Sumw2();
  h_bkg->Sumw2();

  // Loop over each n_e entry in the PCD table
  for (const auto& kv : pcd_table) {
    int ne = kv.first;
    const TH1D* h_pcd = kv.second.get();

    if (!h_pcd) continue;

    // Find S(n_e), B(n_e)
    double S_ne = 0.0;
    double B_ne = 0.0;

    auto itS = Sn.find(ne);
    if (itS != Sn.end()) S_ne = itS->second;

    auto itB = Bn.find(ne);
    if (itB != Bn.end()) B_ne = itB->second;

    // Fold P(q | n_e) with S(n_e) and B(n_e)
    for (int ibin = 1; ibin <= nbins; ++ibin) {
      double Pq = h_pcd->GetBinContent(ibin);
      h_sig->AddBinContent(ibin, S_ne * Pq);
      h_bkg->AddBinContent(ibin, B_ne * Pq);
    }
  }

  return { std::move(h_sig), std::move(h_bkg) };
}

// =====================================================================
// NEW: Build P(n_obs | n_true) from P(q | n_true)
// =====================================================================

PCDCalculator::NeKernel
PCDCalculator::BuildNeKernelFromPCD(const std::map<int, std::unique_ptr<TH1D>>& pcd_table,
                                    int ne_min, int ne_max,
                                    double sigma_res,
                                    double Dqmin,
                                    double Dqmax) const
{
  if (pcd_table.empty()) {
    throw std::runtime_error("PCDCalculator::BuildNeKernelFromPCD: empty P(q | n_e) table.");
  }
  if (ne_max < ne_min) {
    throw std::runtime_error("PCDCalculator::BuildNeKernelFromPCD: ne_max < ne_min.");
  }
  if (!(sigma_res > 0.0)) {
    throw std::runtime_error("PCDCalculator::BuildNeKernelFromPCD: sigma_res must be > 0.");
  }

  const int n_true_bins = ne_max - ne_min + 1;
  const int n_obs_bins  = n_true_bins;

  NeKernel kernel(static_cast<std::size_t>(n_true_bins),
                  std::vector<double>(static_cast<std::size_t>(n_obs_bins), 0.0));

  const TH1D* htemplate = pcd_table.begin()->second.get();
  if (!htemplate) {
    throw std::runtime_error("PCDCalculator::BuildNeKernelFromPCD: null histogram template.");
  }

  const int nbins_q = htemplate->GetNbinsX();

  // Loop over true n_e values present in the PCD table
  for (const auto& kv : pcd_table) {
    const int ne_true = kv.first;
    if (ne_true < ne_min || ne_true > ne_max) continue;

    const TH1D* h_pcd = kv.second.get();
    if (!h_pcd) continue;

    const int i_true = ne_true - ne_min;

    // For each possible observed n_obs (our "qe")
    for (int n_obs = ne_min; n_obs <= ne_max; ++n_obs) {
      const double qe = static_cast<double>(n_obs);

      // Define window [qmin, qmax] around qe, in analogy to Python
      const double qmin = qe - Dqmin * sigma_res;
      const double qmax = qe + Dqmax * sigma_res;

      double prob = 0.0;

      // Integrate P(q | n_true) over [qmin, qmax] by summing bins
      for (int ibin = 1; ibin <= nbins_q; ++ibin) {
        const double q_center = h_pcd->GetXaxis()->GetBinCenter(ibin);
        if (q_center < qmin || q_center > qmax) continue;

        const double Pq = h_pcd->GetBinContent(ibin);
        if (Pq <= 0.0) continue;

        prob += Pq;
      }

      const int j_obs = n_obs - ne_min;
      kernel[static_cast<std::size_t>(i_true)][static_cast<std::size_t>(j_obs)] = prob;
    }

    // Renormalize this row so Σ_{n_obs} P(n_obs | n_true) = 1 (if > 0)
    double row_sum = 0.0;
    for (double v : kernel[static_cast<std::size_t>(i_true)]) row_sum += v;

    if (row_sum > 0.0) {
      for (double& v : kernel[static_cast<std::size_t>(i_true)]) {
        v /= row_sum;
      }
    }
  }

  return kernel;
}


// =====================================================================
// NEW: Fold S_true(n_e) through P(n_obs | n_true)
// =====================================================================

std::unique_ptr<TH1D>
PCDCalculator::FoldNeSpectrum(TH1D& h_true,
                              const NeKernel& kernel,
                              int ne_min, int ne_max,
                              const std::string& name) const
{
  if (ne_max < ne_min) {
    throw std::runtime_error("PCDCalculator::FoldNeSpectrum: ne_max < ne_min.");
  }

  const int n_true_bins = ne_max - ne_min + 1;
  if (static_cast<int>(kernel.size()) != n_true_bins) {
    throw std::runtime_error("PCDCalculator::FoldNeSpectrum: kernel has wrong number of true-n_e rows.");
  }

  // Clone the input histogram to preserve binning, then reset contents
  auto h_rec = std::unique_ptr<TH1D>(
      static_cast<TH1D*>(h_true.Clone(name.c_str())));
  h_rec->Reset("ICES"); // reset contents, keep axis, errors

  // For convenience, we assume obs range is also [ne_min, ne_max]
  for (int ne_true = ne_min; ne_true <= ne_max; ++ne_true) {
    const int i_true = ne_true - ne_min;
    const auto& row  = kernel[static_cast<std::size_t>(i_true)];
    if (static_cast<int>(row.size()) != n_true_bins) {
      throw std::runtime_error("PCDCalculator::FoldNeSpectrum: kernel row has wrong length.");
    }

    const int bin_true = h_true.FindBin(ne_true);
    const double S_true = h_true.GetBinContent(bin_true);
    if (S_true <= 0.0) continue;

    for (int ne_obs = ne_min; ne_obs <= ne_max; ++ne_obs) {
      const int j_obs = ne_obs - ne_min;
      const double P = row[static_cast<std::size_t>(j_obs)];
      if (P <= 0.0) continue;

      const int bin_obs = h_rec->FindBin(ne_obs);
      const double delta = S_true * P;
      h_rec->AddBinContent(bin_obs, delta);
    }
  }

  return h_rec;
}

} // namespace ccdarksens

