// ============================================================================
//  CCDarkSens — ChargeIonization
//  Header for table-driven P(n_e|E) ionization and dR/dE→n_e folding.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#pragma once
#include <memory>
#include <string>
#include <utility>
#include <vector>

class TH1D;

namespace ccdarksens {

// Table-driven P(n_e | E): CSV columns = Er_eV, P1, P2, ..., Pk
class ChargeIonization {
public:
  explicit ChargeIonization(std::string table_csv);

  std::unique_ptr<TH1D> FoldToNe(const TH1D& dRdE,
                                 double exposure_kg_year,
                                 int ne_min, int ne_max) const;

  int MaxNeFromTable() const { return static_cast<int>(pn_given_E_.size()); }

  // NEW: return P(n_e | E) for ne_min ≤ n_e ≤ ne_max
  std::vector<double> ProbNeGivenE(double E_eV,
                                  int ne_min,
                                  int ne_max) const;

private:
  std::vector<std::pair<std::vector<double>, std::vector<double>>> pn_given_E_;
  static double InterpLinearClamped(const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    double xq);
  void LoadCSV_(const std::string& path);
};

} // namespace ccdarksens
