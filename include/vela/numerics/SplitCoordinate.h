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

    SplitCoordinate added(const SplitCoordinate& rhs) const {
        const auto primary = sum(hi, rhs.hi);
        return sum(primary.hi, primary.lo + (lo + rhs.lo));
    }

    SplitCoordinate scaled(double factor) const {
        const double product = hi * factor;
        return sum(product, std::fma(hi, factor, -product) + lo * factor);
    }

    SplitCoordinate multiplied(const SplitCoordinate& rhs) const {
        const double product = hi * rhs.hi;
        return sum(product, std::fma(hi, rhs.hi, -product) + (hi * rhs.lo + lo * rhs.hi));
    }

    SplitCoordinate divided(double divisor) const {
        const double quotient = hi / divisor;
        return sum(quotient, (std::fma(-quotient, divisor, hi) + lo) / divisor);
    }
};
}
