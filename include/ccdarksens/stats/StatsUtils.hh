#pragma once
#include <cmath>

namespace ccdarksens::stats {

inline double safe_log(double x) {
    constexpr double kFloor = 1e-300;
    return std::log(x < kFloor ? kFloor : x);
}

} // namespace ccdarksens::stats
