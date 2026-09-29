// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  BackgroundEfficiencyTable.cc --
// ===========================================================================

#include "ccdarksens/response/BackgroundEfficiencyTable.hh"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ccdarksens {

// ----------------------------------------------------------------------------
// LoadBackgroundEfficiencyTable
//   I read a two-column CSV (E_keV_ee, efficiency). Blank and '#' lines and
//   rows that do not parse as two numbers are skipped. I do not sort the rows;
//   the file has to be ascending in energy.
//   Throws std::runtime_error if the file cannot be opened or has no valid rows.
// ----------------------------------------------------------------------------
BackgroundEfficiencyTable LoadBackgroundEfficiencyTable(const std::string& csv_path) {
  std::ifstream in(csv_path);
  if (!in.is_open()) {
    throw std::runtime_error("BackgroundEfficiencyTable: failed to open " + csv_path);
  }
  BackgroundEfficiencyTable table;
  std::string line;
  while (std::getline(in, line)) {
    while (!line.empty() && std::isspace(static_cast<unsigned char>(line.front()))) line.erase(line.begin());
    if (line.empty() || line[0] == '#') continue;

    std::stringstream ss(line);
    std::string col;
    if (!std::getline(ss, col, ',')) continue;
    double E_keV = 0.0;
    try {
      E_keV = std::stod(col);
    } catch (...) {
      continue;
    }
    if (!std::getline(ss, col, ',')) continue;
    double eff = 0.0;
    try {
      eff = std::stod(col);
    } catch (...) {
      continue;
    }
    table.E_keV.push_back(E_keV);
    table.efficiency.push_back(eff);
  }
  if (table.E_keV.empty()) {
    throw std::runtime_error("BackgroundEfficiencyTable: no rows parsed from " + csv_path);
  }
  return table;
}

// ----------------------------------------------------------------------------
// EvaluateBackgroundEfficiency
//   Efficiency at E_keV: linear interpolation inside the table, and the
//   first/last value held constant below/above its range.
// ----------------------------------------------------------------------------
double EvaluateBackgroundEfficiency(const BackgroundEfficiencyTable& table, double E_keV) {
  const std::size_t n = table.E_keV.size();
  if (E_keV <= table.E_keV.front()) return table.efficiency.front();
  if (E_keV >= table.E_keV.back()) return table.efficiency.back();

  const auto it = std::lower_bound(table.E_keV.begin(), table.E_keV.end(), E_keV);
  const std::size_t hi = static_cast<std::size_t>(it - table.E_keV.begin());
  const std::size_t lo = hi - 1;
  const double x0 = table.E_keV[lo], x1 = table.E_keV[hi];
  const double y0 = table.efficiency[lo], y1 = table.efficiency[hi];
  const double t = (x1 > x0) ? (E_keV - x0) / (x1 - x0) : 0.0;
  return y0 + t * (y1 - y0);
}

}  // namespace ccdarksens
