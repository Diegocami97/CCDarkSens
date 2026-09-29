// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ClusterEnergyResponseFold.hh -- Header for the cluster-energy
//  ResponseFold implementation, wrapping FoldEtrueToErecoRates against a
//  precomputed ClusterFitMC kernel.
// ===========================================================================

#pragma once

#include "ccdarksens/response/ClusterEnergyRates.hh"
#include "ccdarksens/response/ClusterFitMC.hh"
#include "ccdarksens/response/ResponseFold.hh"

namespace ccdarksens {

/**
 * Wraps FoldEtrueToErecoRates against a KernelMatrix built ONCE (expensive,
 * ~minutes, via ClusterFitMC::BuildKernel) before the grid loop -- the
 * cluster-energy analogue of PatternResponseFold/NeSpaceResponseFold. No new
 * physics: the kernel and the fold function are both already validated
 * (docs/ClusterFitMC_Design.md); this is exclusively the ResponseFold
 * wrapping so the scan app can treat all three analysis spaces identically.
 */
// ResponseFold for the cluster_energy analysis space (see the comment above).
class ClusterEnergyResponseFold : public ResponseFold {
 public:
  // Takes the already-built kernel by value.
  explicit ClusterEnergyResponseFold(KernelMatrix kernel);

  std::vector<double> Fold(const TH1D& dRdE_density, double exposure_kg_year) const override;
  std::size_t NumBins() const override;
  std::vector<std::string> BinLabels() const override;
  std::vector<double> ErecoEdgesEV() const override { return kernel_.Ereco_edges_eV; }

 private:
  KernelMatrix kernel_;  // K[E_true][E_reco] and its grids
};

}  // namespace ccdarksens
