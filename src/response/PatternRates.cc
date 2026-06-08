// ============================================================================
//  CCDarkSens — PatternRates
//  Folds n_e histograms into per-pattern expected rates using efficiency maps and applies background migration between true and identified patterns.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#include "ccdarksens/response/PatternRates.hh"
#include <TH1D.h>

namespace ccdarksens {

std::vector<double> FoldNeToPatternRates(
    const TH1D& h_ne,
    int ne_min,
    int ne_max,
    const std::vector<int>& pattern_roi,
    const std::map<std::pair<int, int>, double>& pattern_eff_map) {
  std::vector<double> out;
  out.reserve(pattern_roi.size());

  for (int pattern_id : pattern_roi) {
    double rate = 0.0;
    for (int ne = ne_min; ne <= ne_max; ++ne) {
      int bin = const_cast<TH1D&>(h_ne).FindBin(static_cast<double>(ne));
      double s_ne = h_ne.GetBinContent(bin);
      auto it = pattern_eff_map.find({pattern_id, ne});
      double eff = (it != pattern_eff_map.end()) ? it->second : 0.0;
      rate += s_ne * eff;
    }
    out.push_back(rate);
  }
  return out;
}

std::vector<double> FoldBackgroundWithMigration(
    const std::map<int, double>& B_true_pat,
    const std::vector<int>& pattern_roi,
    const std::map<std::pair<int, int>, double>& migration_map) {
  std::vector<double> out;
  out.reserve(pattern_roi.size());
  for (int iden_pat : pattern_roi) {
    double rate = 0.0;
    for (const auto& kv : B_true_pat) {
      int true_pat = kv.first;
      double b_true = kv.second;
      if (b_true == 0.0) continue;
      auto it = migration_map.find({true_pat, iden_pat});
      double p_iden_given_true = (it != migration_map.end()) ? it->second : 0.0;
      rate += b_true * p_iden_given_true;
    }
    out.push_back(rate);
  }
  return out;
}

}  // namespace ccdarksens
