// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PatternResponseFold.hh -- Header for the pattern-space ResponseFold
//  implementation, wrapping ChargeIonization::FoldToNe +
//  FoldNeToPatternRates.
// ===========================================================================

#pragma once

#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "ccdarksens/response/ResponseFold.hh"

namespace ccdarksens {

class ChargeIonization;

// ----------------------------------------------------------------------------
// PatternResponseFoldConfig
//   Inputs of the pattern-space fold, built once by the caller.
// ----------------------------------------------------------------------------
struct PatternResponseFoldConfig {
  int ne_min = 0;  // lowest n_e of the ionization fold
  int ne_max = 0;  // highest n_e of the ionization fold
  std::vector<int> pattern_roi;  // pattern IDs the likelihood is defined over
  /// (pattern_id, ne) -> efficiency. Built ONCE by the caller (EfficiencyMC's
  /// PrecomputeEpsilon output, flattened) -- not owned or rebuilt here.
  std::map<std::pair<int, int>, double> pattern_eff_map;
};

/**
 * Wraps the reference app's own signal-path fold sequence
 * (ccdarksens_scan_dmelectron_pattern.cc, "Use S_true, not S_obs" comment):
 *   S_true = ion->FoldToNe(dRdE, exposure, ne_min, ne_max)
 *   S_pat  = FoldNeToPatternRates(*S_true, ne_min, ne_max, pattern_roi, pattern_eff_map)
 * No new physics -- the same two already-validated calls, given the shared
 * ResponseFold interface so a generic scan app can treat pattern-space and
 * cluster-energy-space interchangeably.
 */
class PatternResponseFold : public ResponseFold {
 public:
  // Constructor: ionization model and the fold configuration.
  PatternResponseFold(std::shared_ptr<ChargeIonization> ion, PatternResponseFoldConfig cfg);

  std::vector<double> Fold(const TH1D& dRdE_density, double exposure_kg_year) const override;
  std::size_t NumBins() const override { return cfg_.pattern_roi.size(); }
  std::vector<std::string> BinLabels() const override;
  std::vector<double> FoldNeBackground(const TH1D& B_ne_asimov) const override;
  bool SupportsDarkCurrentBackground() const override { return true; }

 private:
  std::shared_ptr<ChargeIonization> ion_;  // P(n_e | E) folding
  PatternResponseFoldConfig cfg_;  // configuration
};

}  // namespace ccdarksens
