// ============================================================================
//  CCDarkSens — DarkPhotonModel
//  Header for hidden-photon absorption model configuration and dR/dE rate-table spectrum loading.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#pragma once
#include <memory>
#include <string>
#include "ccdarksens/io/RateTable.hh"

class TH1D;

namespace ccdarksens {

struct DarkPhotonConfig {
  std::string material;           // "Si"
  std::string mediator;           // fixed tag, e.g. "absorption"
  std::string rates_dir;          // e.g. data/darkphoton_rates/Si
  std::string filename_template;  // e.g. dRdE_{material}_{mediator}_m{mA_eV}_e{epsilon}.csv
  double      mA_eV = 0.0;        // hidden-photon mass, used in filename
  std::string epsilon;            // kinetic-mixing parameter, kept as string to match filename exactly
  std::string epsilon_ref = "";   // if set, load file at epsilon_ref and scale rate by (epsilon/epsilon_ref)^2

  // output spectrum binning
  double Emin_eV = 0.0;
  double Emax_eV = 20.0;
  int    nbins   = 200;
};

class DarkPhotonModel {
public:
  DarkPhotonModel() = default;
  ~DarkPhotonModel();
  bool Configure(const DarkPhotonConfig& cfg);

  /// dR/dE in events / kg / day / eV
  std::unique_ptr<TH1D> MakeSpectrum_E() const;

  const DarkPhotonConfig& cfg() const noexcept { return cfg_; }

private:
  std::string ResolvePath_() const;
  DarkPhotonConfig              cfg_{};
  std::unique_ptr<RateTable>    table_;
};

} // namespace ccdarksens
