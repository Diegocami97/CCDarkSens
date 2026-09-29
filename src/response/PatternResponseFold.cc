// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PatternResponseFold.cc -- Pattern-space ResponseFold:
//  ChargeIonization::FoldToNe + FoldNeToPatternRates.
// ===========================================================================

#include "ccdarksens/response/PatternResponseFold.hh"

#include <TH1D.h>

#include "ccdarksens/response/ChargeIonization.hh"
#include "ccdarksens/response/PatternRates.hh"

namespace ccdarksens {

// Constructor: keep the ionization model and the configuration.
PatternResponseFold::PatternResponseFold(std::shared_ptr<ChargeIonization> ion,
                                          PatternResponseFoldConfig cfg)
    : ion_(std::move(ion)), cfg_(std::move(cfg)) {}

// ----------------------------------------------------------------------------
// PatternResponseFold::Fold
//   Signal per pattern bin: dR/dE -> n_e with the ionization table, then FoldNeToPatternRates.
// ----------------------------------------------------------------------------
std::vector<double> PatternResponseFold::Fold(const TH1D& dRdE_density,
                                               double exposure_kg_year) const {
  auto S_true = ion_->FoldToNe(dRdE_density, exposure_kg_year, cfg_.ne_min, cfg_.ne_max);
  return FoldNeToPatternRates(*S_true, cfg_.ne_min, cfg_.ne_max, cfg_.pattern_roi, cfg_.pattern_eff_map);
}

// ----------------------------------------------------------------------------
// PatternResponseFold::FoldNeBackground
//   Dark-current histogram (already in n_e) folded into the pattern bins.
// ----------------------------------------------------------------------------
std::vector<double> PatternResponseFold::FoldNeBackground(const TH1D& B_ne_asimov) const {
  return FoldNeToPatternRates(B_ne_asimov, cfg_.ne_min, cfg_.ne_max, cfg_.pattern_roi, cfg_.pattern_eff_map);
}

// ----------------------------------------------------------------------------
// PatternResponseFold::BinLabels
//   The pattern IDs as strings, one per bin.
// ----------------------------------------------------------------------------
std::vector<std::string> PatternResponseFold::BinLabels() const {
  std::vector<std::string> labels;
  labels.reserve(cfg_.pattern_roi.size());
  for (int pid : cfg_.pattern_roi) labels.push_back(std::to_string(pid));
  return labels;
}

}  // namespace ccdarksens
