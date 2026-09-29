// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  TestStatisticFactory.cc -- Factory that instantiates the configured
//  binned test statistic (currently PoissonAsimovPLR) from StatisticsConfig.
// ===========================================================================

#include "ccdarksens/stats/TestStatisticFactory.hh"

#include <algorithm>
#include <iostream>

#include "ccdarksens/stats/PoissonAsimovPLR.hh"

namespace ccdarksens::stats {

namespace {

// Lower-case copy of s (for case-insensitive test-statistic names).
std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

} // namespace

// ----------------------------------------------------------------------------
// MakeTestStatistic
//   Build the test statistic named in cfg.test_stat ("PLR" / "poisson_plr" /
//   empty, case-insensitive). An unknown name prints a warning and falls back to
//   PoissonAsimovPLR.
// ----------------------------------------------------------------------------
std::unique_ptr<ITestStatistic>
MakeTestStatistic(const StatisticsConfig& cfg) {
  const std::string name = to_lower(cfg.test_stat);

  if (name == "plr" || name == "poisson_plr" || name.empty()) {
    if (cfg.verbosity > 0) {
      std::cout << "[stats] Using PoissonAsimovPLR test statistic (\"" << cfg.test_stat << "\")\n";
    }
    return std::make_unique<PoissonAsimovPLR>();
  }

  // Unknown name -> warn and fall back to PLR
  std::cerr << "[stats] WARNING: unknown test_stat=\"" << cfg.test_stat
            << "\". Falling back to PoissonAsimovPLR.\n";
  return std::make_unique<PoissonAsimovPLR>();
}

} // namespace ccdarksens::stats
