// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  TestStatisticFactory.hh -- Header declaring MakeTestStatistic for
//  constructing likelihood-ratio implementations from config.
// ===========================================================================

#pragma once

#include <memory>

#include "ccdarksens/stats/ITestStatistic.hh"
#include "ccdarksens/stats/StatisticsConfig.hh"

namespace ccdarksens::stats {

/**
 * Create an appropriate test statistic implementation based on StatisticsConfig.
 *
 * For now:
 *   test_stat = "PLR" (case-insensitive) -> PoissonAsimovPLR
 * Future:
 *   add other names and corresponding classes.
 */
std::unique_ptr<ITestStatistic> MakeTestStatistic(const StatisticsConfig& cfg);

} // namespace ccdarksens::stats
