#pragma once

#include <string>
#include <vector>
#include <cstdint>

namespace ccdarksens::stats {

/**
 * Configuration for statistical methods, derived from JSON "run" (and
 * optionally a future "statistics" block).
 */
struct StatisticsConfig {
  // From "run"
  std::string   test_stat;   ///< e.g. "PLR"
  double        cl        = 0.90;  ///< confidence level (0 < cl < 1)
  int           n_toys    = 0;     ///< reserved for future toy-MC methods
  std::uint64_t rng_seed  = 0;
  int           verbosity = 0;     ///< 0=silent, 1=summary, 2+=debug

  // Bin selection: which n_e bins enter the likelihood
  enum class BinMode {
    ROI,    ///< use experiment.roi_bins
    ALL,    ///< use all [ne_min, ne_max]
    CUSTOM  ///< use custom_bins below
  };

  BinMode            bin_mode    = BinMode::ROI;
  std::vector<int>   custom_bins; ///< used only if bin_mode == CUSTOM
};

} // namespace ccdarksens::stats
