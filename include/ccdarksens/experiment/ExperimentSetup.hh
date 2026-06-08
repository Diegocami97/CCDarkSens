// ============================================================================
//  CCDarkSens — ExperimentSetup
//  Header defining experiment modes, binning, pattern ROI, and exposure summary computation.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace ccdarksens {

enum class ExperimentMode { Observed, Asimov, Toys };

struct BinningNE {
  int ne_min = 0;
  int ne_max = 0; // inclusive
};

struct ExperimentConfig {
  ExperimentMode mode = ExperimentMode::Asimov;
  double livetime_days = 0.0;
  double duty_cycle    = 1.0;
  BinningNE binning;
  std::vector<int> roi_bins;
  /// Pattern IDs for likelihood when observable_bins == "pattern" (e.g. 11, 21, 111, 31, 22, 211).
  std::vector<int> pattern_roi;
  /// Observable for likelihood: "n_e" (default) or "pattern".
  std::string observable_bins = "n_e";
};

struct ExperimentSummary {
  double exposure_kg_year = 0.0;
  BinningNE binning;
  std::vector<int> roi_bins;
  std::vector<int> pattern_roi;
  std::string observable_bins = "n_e";
  uint64_t rng_seed_used = 0;
  std::string mode_string;
};

class ExperimentSetup {
public:
  ExperimentSetup(ExperimentConfig cfg, double detector_mass_kg, uint64_t rng_seed);
  ExperimentSummary prepare_summary() const;

private:
  ExperimentConfig cfg_;
  double           detector_mass_kg_;
  uint64_t         rng_seed_;
};

} // namespace ccdarksens
