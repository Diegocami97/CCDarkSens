// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PCDBasedResponse.hh -- Header for Monte Carlo construction of P(q|n_e)
//  via ChargeTransport and PixelSimulator.
// ===========================================================================

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <TH1D.h>

#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/PixelSimulator.hh"

namespace ccdarksens {

/**
 * Configuration for PCDBasedResponse: how to build P(q | n_e).
 *
 * q_min, q_max, nbins define the q-axis histogram.
 * mc_trials defines how many Monte Carlo events to run for each n_e.
 */
struct PCDResponseConfig {
  double q_min  = 0.0;  // lowest total charge [e-]
  double q_max  = 20.0;  // highest total charge [e-]
  int    nbins  = 200;  // charge bins
  int    mc_trials = 50000;  // Monte Carlo events per n_e

  // Pixel simulator settings (geometry, noise, DC, etc.)
  PixelSimulatorConfig pix_cfg;

  // RNG seed for reproducibility
  unsigned long seed = 12345;
};

/**
 * PCDBasedResponse:
 *
 *  Builds P(q | n_e) from full PCD simulation:
 *
 *    - ChargeTransport → pixel (x,y) diffusion
 *    - PixelSimulator  → pixel charges (DC + readout noise)
 *    - total q = sum of all pixel charges
 *
 *  For each n_e in a specified [ne_min, ne_max] range, we generate 'mc_trials'
 *  events, fill a TH1D histogram, and normalize it to obtain P(q | n_e).
 *
 *  Output is stored as:
 *
 *    pcd_table_[n_e] = unique_ptr<TH1D> (owned, normalized)
 */
class PCDBasedResponse {
public:
  // Constructor: settings and the charge-transport model (must not be null).
  PCDBasedResponse(const PCDResponseConfig& cfg,
                   std::shared_ptr<ChargeTransport> ct);

  /**
   * Build P(q | n_e) for ne_min ≤ n_e ≤ ne_max.
   *
   * Returns a map<int, TH1D*> referencing the internally owned histograms.
   */
  const std::map<int, std::unique_ptr<TH1D>>&
  BuildPCDTable(int ne_min, int ne_max, double Ee_ref_eV = 0.0);

  /// Access the PCD table after building.
  const std::map<int, std::unique_ptr<TH1D>>& GetPCDTable() const {
    return pcd_table_;
  }

private:
  PCDResponseConfig cfg_;  // settings
  std::shared_ptr<ChargeTransport> ct_;  // charge-transport model

  /// pcd_table_[n_e] holds the normalized TH1D for that n_e.
  std::map<int, std::unique_ptr<TH1D>> pcd_table_;
};

} // namespace ccdarksens
