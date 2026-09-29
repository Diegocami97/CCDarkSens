// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PatternEfficiency.hh -- Header for applying a stored ε(n_e) efficiency
//  histogram to spectra.
// ===========================================================================

#pragma once
#include <memory>

class TH1D;

namespace ccdarksens {

// ----------------------------------------------------------------------------
// PatternEfficiency
//   Holds a detection-efficiency histogram epsilon(n_e) and multiplies a
//   spectrum by it bin by bin.
// ----------------------------------------------------------------------------
class PatternEfficiency {
public:
  void SetEfficiencyHist(const TH1D& epsilon_ne); // clones input
  void Apply(TH1D& target_ne) const;              // multiplies bin-by-bin

private:
  std::unique_ptr<TH1D> eps_;  // private copy of epsilon(n_e)
};

} // namespace ccdarksens
