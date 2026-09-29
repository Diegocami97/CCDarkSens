// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  WimpNucleonModel.hh -- I declare WimpNucleonModel: it loads the pre-
//  computed rate table for one grid point (WIMP-nucleus spin-independent
//  elastic scattering) and returns the dR/dE_ee spectrum as a ROOT
//  histogram.
// ===========================================================================

#pragma once
#include <memory>
#include <string>
#include "ccdarksens/io/RateTable.hh"
#include "ccdarksens/model/DMNucleonConfig.hh"

class TH1D;

namespace ccdarksens {

// ----------------------------------------------------------------------------
// WimpNucleonModel
//   Signal-rate model for WIMP-nucleus spin-independent elastic scattering. I resolve the rate-table CSV
//   for one (mass, coupling) grid point from the config's directory and
//   filename template, load it with RateTable, and return dR/dE_ee as a ROOT
//   histogram in events/(kg*year*eV) on the requested energy grid.
// ----------------------------------------------------------------------------
class WimpNucleonModel {
public:
  WimpNucleonModel() = default;
  ~WimpNucleonModel();
  // Resolve the CSV path from the config and load the table.
  // Returns false if the file is missing or unreadable (the spectrum is then null).
  bool Configure(const DMNucleonConfig& cfg);

  /// dR/dE_ee in events / kg / year / eV
  std::unique_ptr<TH1D> MakeSpectrum_E() const;

  const DMNucleonConfig& cfg() const noexcept { return cfg_; }

private:
  // Build the CSV path by substituting the grid-point values into the filename template.
  std::string ResolvePath_() const;
  DMNucleonConfig                cfg_{};
  std::unique_ptr<RateTable>     table_;  // the loaded rate table (null until Configure() succeeds)
};

} // namespace ccdarksens
