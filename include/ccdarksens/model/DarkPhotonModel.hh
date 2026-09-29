// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  DarkPhotonModel.hh -- I declare DarkPhotonModel: it loads the pre-
//  computed rate table for one grid point (hidden-photon absorption) and
//  returns the dR/dE_e spectrum as a ROOT histogram.
// ===========================================================================

#pragma once
#include <memory>
#include <string>
#include "ccdarksens/io/RateTable.hh"

class TH1D;

namespace ccdarksens {

// ----------------------------------------------------------------------------
// DarkPhotonConfig
//   Where the dark-photon absorption rate files live and which (m_A, epsilon)
//   grid point to load.
// ----------------------------------------------------------------------------
struct DarkPhotonConfig {
  std::string material;           // "Si"
  std::string mediator;           // fixed tag, e.g. "absorption"
  std::string rates_dir;          // e.g. data/darkphoton_rates/Si
  std::string filename_template;  // e.g. dRdE_{material}_{mediator}_m{mA_eV}_e{epsilon}.csv
  double      mA_eV = 0.0;        // hidden-photon mass, used in filename
  std::string epsilon;            // kinetic-mixing parameter, kept as string to match filename exactly
  std::string epsilon_ref = "";   // if set, load file at epsilon_ref and scale rate by (epsilon/epsilon_ref)^2

  // output spectrum binning
  double Emin_eV = 0.0;  // lower edge of the output spectrum [eV]
  double Emax_eV = 20.0;  // upper edge of the output spectrum [eV]
  int    nbins   = 200;  // number of output bins
};

// ----------------------------------------------------------------------------
// DarkPhotonModel
//   Signal-rate model for hidden-photon absorption. I resolve the rate-table CSV
//   for one (mass, coupling) grid point from the config's directory and
//   filename template, load it with RateTable, and return dR/dE_e as a ROOT
//   histogram in events/(kg*year*eV) on the requested energy grid.
// ----------------------------------------------------------------------------
class DarkPhotonModel {
public:
  DarkPhotonModel() = default;
  ~DarkPhotonModel();
  // Resolve the CSV path from the config and load the table.
  // Returns false if the file is missing or unreadable (the spectrum is then null).
  bool Configure(const DarkPhotonConfig& cfg);

  /// dR/dE in events / kg / year / eV
  std::unique_ptr<TH1D> MakeSpectrum_E() const;

  const DarkPhotonConfig& cfg() const noexcept { return cfg_; }

private:
  // Build the CSV path by substituting the grid-point values into the filename template.
  std::string ResolvePath_() const;
  DarkPhotonConfig              cfg_{};
  std::unique_ptr<RateTable>    table_;  // the loaded rate table (null until Configure() succeeds)
};

} // namespace ccdarksens
