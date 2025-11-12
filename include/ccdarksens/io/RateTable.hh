#pragma once
#include <memory>
#include <string>
#include <vector>

class TH1D;

namespace ccdarksens {

struct RateMeta {
  std::string material;
  std::string mediator;
  double      mchi_MeV     = 0.0;
  double      sigma_e_cm2  = 0.0;
  double      binsize_eV   = 0.1;
};

/// Holds a QEDark rate table (E [eV], dR/dE [events / g / day / eV])
class RateTable {
public:
  RateTable() = default;

  /// Load a 2-column CSV with a single header row; comments (#) are ignored.
  /// Columns: E_eV, dRdE_g_day_eV
  /// Returns true on success.
  bool LoadCSV(const std::string& path);

  const std::vector<double>& E() const noexcept { return E_eV_; }
  const std::vector<double>& R_kg_year_eV() const noexcept { return R_kg_year_eV_; }
  const RateMeta& meta() const noexcept { return meta_; }

  /// Build a ROOT spectrum in units of events / kg / day / eV.
  /// Performs a simple linear interpolation onto [Emin, Emax] with nbins.
  std::unique_ptr<TH1D> MakeTH1D(const std::string& name,
                                 double Emin_eV, double Emax_eV,
                                 int nbins) const;

private:
  std::vector<double> E_eV_;
  std::vector<double> R_kg_year_eV_;
  RateMeta            meta_;
};

} // namespace ccdarksens
