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
};

struct ExperimentSummary {
  double exposure_kg_day = 0.0;
  BinningNE binning;
  std::vector<int> roi_bins;
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
