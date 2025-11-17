#include "ccdarksens/stats/TestStatisticFactory.hh"

#include <algorithm>
#include <iostream>

#include "ccdarksens/stats/PoissonAsimovPLR.hh"

namespace ccdarksens::stats {

namespace {

std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

} // namespace

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
