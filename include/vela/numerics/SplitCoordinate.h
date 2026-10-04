#pragma once
#include <cmath>
#include <stdexcept>

namespace vela::numerics {
// A normalized binary64 pair. Consumers must use both parts: hi + lo in
// binary64 alone would discard the correction this type exists to retain.
struct SplitCoordinate {
    double hi = 0.0;
    double lo = 0.0;

    static SplitCoordinate sum(double a, double b) {
        const double s = a + b;
        if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(s))
            throw std::invalid_argument("SplitCoordinate requires a finite sum");
        const double z = s - a;
        return {s, (a - (s - z)) + (b - z)};
    }

    SplitCoordinate shifted(double delta) const {
        const auto correction = sum(lo, delta);
        const auto primary = sum(hi, correction.hi);
        return sum(primary.hi, primary.lo + correction.lo);
    }
};
}
