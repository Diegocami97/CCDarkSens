// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ChargeIonization.hh -- Header for table-driven P(n_e|E) ionization and
//  dR/dE→n_e folding.
// ===========================================================================

#pragma once
#include <memory>
#include <string>
#include <utility>
#include <vector>

class TH1D;

namespace ccdarksens {

// Table-driven P(n_e | E): CSV columns = Er_eV, P1, P2, ..., Pk
// ----------------------------------------------------------------------------
// ChargeIonization
//   Table-driven ionization response P(n_e | E). The CSV columns are
//   Er_eV, P1, P2, ..., Pk (probability of producing n = 1..k electrons at
//   energy Er). I use it to fold a dR/dE spectrum into an n_e spectrum.
// ----------------------------------------------------------------------------
class ChargeIonization {
public:
  // Load the P(n_e | E) table from a CSV file (throws if it cannot be read or is malformed).
  explicit ChargeIonization(std::string table_csv);

  // Fold dR/dE [events/(kg*year*eV)] * exposure into an n_e histogram over [ne_min, ne_max].
  std::unique_ptr<TH1D> FoldToNe(const TH1D& dRdE,
                                 double exposure_kg_year,
                                 int ne_min, int ne_max) const;

  // Largest n_e the table has a column for (k).
  int MaxNeFromTable() const { return static_cast<int>(pn_given_E_.size()); }

  // NEW: return P(n_e | E) for ne_min ≤ n_e ≤ ne_max
  // P(n_e | E_eV) for every n_e in [ne_min, ne_max] (zero outside the table's n_e columns).
  std::vector<double> ProbNeGivenE(double E_eV,
                                  int ne_min,
                                  int ne_max) const;

private:
  std::vector<std::pair<std::vector<double>, std::vector<double>>> pn_given_E_;  // one (E grid, P) pair per n_e = 1..k
  // Linear interpolation of y(x) at xq, clamped to the end values and to [0,1].
  static double InterpLinearClamped(const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    double xq);
  // Read the table from the CSV file.
  void LoadCSV_(const std::string& path);
};

} // namespace ccdarksens
