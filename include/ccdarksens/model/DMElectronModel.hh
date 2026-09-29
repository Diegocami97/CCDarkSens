// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  DMElectronModel.hh -- I declare DMElectronModel: it loads the pre-
//  computed rate table for one grid point (DM-electron scattering) and
//  returns the dR/dE_e spectrum as a ROOT histogram.
// ===========================================================================

#pragma once
#include <memory>
#include <string>
#include "ccdarksens/io/RateTable.hh" 

class TH1D;

namespace ccdarksens {

// ----------------------------------------------------------------------------
// DMElectronConfig
//   Where the DM-electron rate files live and which grid point to load.
// ----------------------------------------------------------------------------
struct DMElectronConfig {
  std::string material;           // "Si"
  std::string mediator;           // "heavy" | "massless"
  std::string rates_dir;          // e.g. data/qedark_rates/Si/heavy
  std::string filename_template;  // e.g. dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv
  double      mchi_MeV = 0.0;     // used in filename
  std::string sigma_e_cm2;        // keep as string to match filename exactly

  // output spectrum binning
  double Emin_eV = 0.0;  // lower edge of the output spectrum [eV]
  double Emax_eV = 20.0;  // upper edge of the output spectrum [eV]
  int    nbins   = 200;  // number of output bins
};

// class RateTable;

// ----------------------------------------------------------------------------
// DMElectronModel
//   Signal-rate model for DM-electron scattering. I resolve the rate-table CSV
//   for one (mass, coupling) grid point from the config's directory and
//   filename template, load it with RateTable, and return dR/dE_e as a ROOT
//   histogram in events/(kg*year*eV) on the requested energy grid.
// ----------------------------------------------------------------------------
class DMElectronModel {
public:
  DMElectronModel() = default;
  ~DMElectronModel();                 // <-- declare only (no = default here)
  // Resolve the CSV path from the config and load the table.
  // Returns false if the file is missing or unreadable (the spectrum is then null).
  bool Configure(const DMElectronConfig& cfg);

  /// dR/dE in events / kg / year / eV
  std::unique_ptr<TH1D> MakeSpectrum_E() const;

  const DMElectronConfig& cfg() const noexcept { return cfg_; }

private:
  // Build the CSV path by substituting the grid-point values into the filename template.
  std::string ResolvePath_() const;
  DMElectronConfig              cfg_{};
  std::unique_ptr<RateTable>    table_;  // the loaded rate table (null until Configure() succeeds)
};

} // namespace ccdarksens
