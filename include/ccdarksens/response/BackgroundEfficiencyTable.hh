// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  BackgroundEfficiencyTable.hh -- Header for loading and evaluating a
//  (E_keV, efficiency) background detection-efficiency curve -- distinct
//  from the signal's kernel-based efficiency (see
//  docs/ClusterFitMC_Design.md Sec. 6.9 for why the two differ, e.g.
//  PhysRevD.94.082006's own Fig. 9).
// ===========================================================================

#pragma once

#include <string>
#include <vector>

namespace ccdarksens {

// ----------------------------------------------------------------------------
// BackgroundEfficiencyTable
//   Tabulated detection efficiency of the background versus energy.
// ----------------------------------------------------------------------------
struct BackgroundEfficiencyTable {
  std::vector<double> E_keV;  // energies [keV_ee], ascending
  std::vector<double> efficiency;  // efficiency at each energy, 0..1
};

/// Loads a two-column CSV (E_keV_ee, efficiency), skipping '#' comment lines.
/// Points must be sorted ascending in E_keV (not re-sorted here).
BackgroundEfficiencyTable LoadBackgroundEfficiencyTable(const std::string& csv_path);

/// Linear interpolation within the table's range; clamped flat to the
/// nearest endpoint value outside it (the table only covers what was
/// actually digitized/measured -- e.g. PhysRevD.94.082006's Fig. 9 only
/// plots to 2 keV_ee, and the paper's own text describes the background
/// efficiency as roughly constant at high energy, "dominated by the
/// contribution from Compton events" -- so a flat extrapolation is the
/// documented, deliberate choice here, not a silent default).
double EvaluateBackgroundEfficiency(const BackgroundEfficiencyTable& table, double E_keV);

}  // namespace ccdarksens
