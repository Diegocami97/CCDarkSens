// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ClusterEnergyResponseFold.cc -- Cluster-energy ResponseFold:
//  FoldEtrueToErecoRates against a precomputed ClusterFitMC kernel.
// ===========================================================================

#include "ccdarksens/response/ClusterEnergyResponseFold.hh"

namespace ccdarksens {

// Constructor: keep the kernel.
ClusterEnergyResponseFold::ClusterEnergyResponseFold(KernelMatrix kernel)
    : kernel_(std::move(kernel)) {}

// ----------------------------------------------------------------------------
// ClusterEnergyResponseFold::Fold
//   Expected counts per E_reco bin for dR/dE_true and the exposure [kg*year].
// ----------------------------------------------------------------------------
std::vector<double> ClusterEnergyResponseFold::Fold(const TH1D& dRdE_density,
                                                      double exposure_kg_year) const {
  return FoldEtrueToErecoRates(dRdE_density, exposure_kg_year, kernel_);
}

// ----------------------------------------------------------------------------
// ClusterEnergyResponseFold::NumBins
//   Number of E_reco bins.
// ----------------------------------------------------------------------------
std::size_t ClusterEnergyResponseFold::NumBins() const {
  return kernel_.Ereco_edges_eV.empty() ? 0 : kernel_.Ereco_edges_eV.size() - 1;
}

// ----------------------------------------------------------------------------
// ClusterEnergyResponseFold::BinLabels
//   Labels "Ereco_<lo>_<hi>eV", one per E_reco bin.
// ----------------------------------------------------------------------------
std::vector<std::string> ClusterEnergyResponseFold::BinLabels() const {
  std::vector<std::string> labels;
  const std::size_t n = NumBins();
  labels.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    labels.push_back("Ereco_" + std::to_string(static_cast<int>(kernel_.Ereco_edges_eV[i])) +
                      "_" + std::to_string(static_cast<int>(kernel_.Ereco_edges_eV[i + 1])) + "eV");
  }
  return labels;
}

}  // namespace ccdarksens
