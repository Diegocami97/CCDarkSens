// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  StatsUtils.hh -- I provide safe_log, the shared logarithm used in every
//  likelihood: it floors its argument at 1e-300 so that a zero or tiny
//  expected count never produces -inf or NaN. Use this instead of writing
//  another local lambda.
// ===========================================================================

#pragma once
#include <cmath>

namespace ccdarksens::stats {

// ln(x) with x floored at 1e-300.
inline double safe_log(double x) {
    constexpr double kFloor = 1e-300;
    return std::log(x < kFloor ? kFloor : x);
}

} // namespace ccdarksens::stats
