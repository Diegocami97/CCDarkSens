// ============================================================================
//  CCDarkSens — PatternRates
//  Header for folding n_e spectra to pattern rates and migrating true-pattern backgrounds to identified bins.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#pragma once

#include <map>
#include <utility>
#include <vector>

class TH1D;

namespace ccdarksens {

/**
 * Fold S(n_e) or B(n_e) histogram into rates per pattern bin.
 *
 * For each pattern p in pattern_roi:
 *   rate_p = sum_{ne = ne_min}^{ne_max} h_ne(ne) * epsilon(p, ne)
 *
 * Efficiency table: map (pattern_id, ne) -> efficiency (e.g. from CSV
 * columns pattern, ne, Efficiency). Missing (p, ne) entries are treated as 0.
 */
std::vector<double> FoldNeToPatternRates(
    const TH1D& h_ne,
    int ne_min,
    int ne_max,
    const std::vector<int>& pattern_roi,
    const std::map<std::pair<int, int>, double>& pattern_eff_map);

/**
 * Fold background from true-pattern space into observed pattern bins using
 * a migration matrix P(identified | true pattern).
 *
 * For each observed pattern iden_pat in pattern_roi:
 *   B_pat[iden_pat] = sum_{true_pat} B_true_pat[true_pat] * P(iden_pat | true_pat)
 *
 * Use when background is modeled as rate per true (ideal) pattern (e.g. DC
 * single-pixel (0)-(5)) and Background_efficiencies.csv gives the migration.
 */
std::vector<double> FoldBackgroundWithMigration(
    const std::map<int, double>& B_true_pat,
    const std::vector<int>& pattern_roi,
    const std::map<std::pair<int, int>, double>& migration_map);

}  // namespace ccdarksens
