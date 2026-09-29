// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  MigdalModel.hh -- I declare MigdalModel: it loads the pre-computed rate
//  table for one grid point (Migdal-effect ionization) and returns the
//  dR/dE_e spectrum as a ROOT histogram.
// ===========================================================================

#pragma once
#include <memory>
#include <string>
#include "ccdarksens/io/RateTable.hh"
#include "ccdarksens/model/DMNucleonConfig.hh"

class TH1D;

namespace ccdarksens {

// ----------------------------------------------------------------------------
// MigdalModel
//   Signal-rate model for Migdal-effect ionization. I resolve the rate-table CSV
//   for one (mass, coupling) grid point from the config's directory and
//   filename template, load it with RateTable, and return dR/dE_e as a ROOT
//   histogram in events/(kg*year*eV) on the requested energy grid.
// ----------------------------------------------------------------------------
class MigdalModel {
public:
  MigdalModel() = default;
  ~MigdalModel();
  // Resolve the CSV path from the config and load the table.
  // Returns false if the file is missing or unreadable (the spectrum is then null).
  bool Configure(const DMNucleonConfig& cfg);

  /// dR/dE_e in events / kg / year / eV
  std::unique_ptr<TH1D> MakeSpectrum_E() const;

  const DMNucleonConfig& cfg() const noexcept { return cfg_; }

private:
  // Build the CSV path by substituting the grid-point values into the filename template.
  std::string ResolvePath_() const;
  DMNucleonConfig                cfg_{};
  std::unique_ptr<RateTable>     table_;  // the loaded rate table (null until Configure() succeeds)
};

} // namespace ccdarksens
