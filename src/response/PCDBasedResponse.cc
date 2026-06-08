// ============================================================================
//  CCDarkSens — PCDBasedResponse
//  Monte Carlo simulation of normalized P(q|n_e) pixel-charge distributions using ChargeTransport and PixelSimulator.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#include "ccdarksens/response/PCDBasedResponse.hh"

#include <TH1D.h>
#include <TRandom3.h>
#include <stdexcept>
#include <numeric>   // for std::accumulate
#include <cmath>

namespace ccdarksens {

PCDBasedResponse::PCDBasedResponse(const PCDResponseConfig& cfg,
                                   std::shared_ptr<ChargeTransport> ct)
: cfg_(cfg),
  ct_(std::move(ct))
{
  if (!ct_) {
    throw std::runtime_error("PCDBasedResponse: ChargeTransport pointer is null.");
  }
}

/**
 * Build P(q | n_e) for ne_min <= n_e <= ne_max.
 *
 * For each n_e:
 *   - Create a TH1D from cfg_.q_min to cfg_.q_max with cfg_.nbins
 *   - Run cfg_.mc_trials events
 *       * For each event: deposit n_e electrons via ChargeTransport,
 *         accumulate pixel charges in PixelSimulator
 *         total charge q = sum of all pixel charges
 *         Fill the histogram
 *   - Normalize the histogram to unity
 */
const std::map<int, std::unique_ptr<TH1D>>&
PCDBasedResponse::BuildPCDTable(int ne_min, int ne_max, double Ee_ref_eV)
{
  pcd_table_.clear();

  // Prepare PixelSimulator (reused for each event)
  PixelSimulator pixSim(cfg_.pix_cfg);

  // RNG for reproducibility
  TRandom3 rng(cfg_.seed);

  // Loop over electron counts
  for (int ne_true = ne_min; ne_true <= ne_max; ++ne_true) {

    // Create histogram for this n_e
    auto h = std::make_unique<TH1D>(
      Form("pcd_qdist_ne%d", ne_true),
      Form("P(q | n_e = %d)", ne_true),
      cfg_.nbins, cfg_.q_min, cfg_.q_max
    );
    h->Sumw2();

    // Monte Carlo
    for (int it = 0; it < cfg_.mc_trials; ++it) { 

      pixSim.Reset();

      // Sample electron cloud
      std::vector<double> xs(ne_true), ys(ne_true);
      ct_->SampleCloudXY(0.0, 0.0, Ee_ref_eV, ne_true, xs, ys);

      // Deposit electrons
      for (int i = 0; i < ne_true; ++i) {
        pixSim.DepositElectron(xs[i], ys[i]);
      }

      // Add DC + readout noise
      // pixSim.AddDarkCurrent();
      pixSim.AddReadoutNoise();

      // Get pixel charges
      const std::vector<double>& qpix = pixSim.PixelCharges();

      // Total charge q = sum of all pixels
      double q_total = std::accumulate(qpix.begin(), qpix.end(), 0.0);

      // Fill histogram
      h->Fill(q_total);
    }

    // Normalize histogram to unity
    double integral = h->Integral();
    if (integral > 0) {
      h->Scale(1.0 / integral);
    }

    // Store histogram
    pcd_table_[ne_true] = std::move(h);
  }

  return pcd_table_;
}

} // namespace ccdarksens
