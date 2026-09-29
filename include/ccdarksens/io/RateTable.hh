// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  RateTable.hh -- I declare RateTable: a two-column rate table (E [eV],
//  dR/dE [events/(kg*year*eV)]) loaded from CSV, and its conversion into a
//  ROOT spectrum on a uniform energy grid.
// ===========================================================================

#pragma once
#include <memory>
#include <string>
#include <vector>

class TH1D;

namespace ccdarksens {

// Descriptive metadata of a rate table (which model/point it belongs to).
struct RateMeta {
  std::string material;  // target material, e.g. "Si"
  std::string mediator;  // mediator type: heavy | light | ultralight
  double      mchi_MeV     = 0.0;  // dark-matter mass [MeV]
  double      sigma_e_cm2  = 0.0;  // reference DM-electron cross section [cm^2]
  double      binsize_eV   = 0.1;  // native energy step of the table [eV]
};

/// Holds a QEDark rate table (E [eV], dR/dE [events/(kg·year·eV)]). CSV must use events/(kg·year·eV).
class RateTable {
public:
  RateTable() = default;

  /// Load a 2-column CSV with a single header row; comments (#) are ignored.
  /// Columns: E [eV], dRdE [events/(kg·year·eV)].
  /// Returns true on success.
  bool LoadCSV(const std::string& path);

  const std::vector<double>& E() const noexcept { return E_eV_; }
  const std::vector<double>& R_kg_year_eV() const noexcept { return R_kg_year_eV_; }
  const RateMeta& meta() const noexcept { return meta_; }

  /// Build a ROOT spectrum: bin content = dR/dE [events/(kg·year·eV)]. Multiply by exposure_kg_year and dE to get counts.
  /// Performs a simple linear interpolation onto [Emin, Emax] with nbins.
  std::unique_ptr<TH1D> MakeTH1D(const std::string& name,
                                 double Emin_eV, double Emax_eV,
                                 int nbins) const;

private:
  std::vector<double> E_eV_;  // energy column [eV]
  std::vector<double> R_kg_year_eV_;  // dR/dE column [events/(kg*year*eV)]
  RateMeta            meta_;  // metadata of the loaded table
};

} // namespace ccdarksens
