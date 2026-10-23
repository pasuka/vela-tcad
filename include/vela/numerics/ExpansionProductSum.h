#pragma once

#include <cmath>
#include <stdexcept>

namespace vela::experimental {

// Form both parts from the original wide coefficient, before discarding it.
// Splitting an already rounded matrix cannot recover its missing coefficient.
struct ExpansionCoefficient {
    double high, low;
};

template<class Wide>
ExpansionCoefficient expansionCoefficient(const Wide& value) {
    const double high = static_cast<double>(value);
    if (!std::isfinite(high)) throw std::overflow_error("Unrepresentable expansion coefficient");
    const double low = static_cast<double>(value - Wide(high));
    if (!std::isfinite(low)) throw std::overflow_error("Unrepresentable expansion remainder");
    return {high, low};
}

// Two-part accumulation with FMA product errors. Requires ordinary IEEE
// round-to-nearest arithmetic, without fast-math reassociation. This is an
// experimental operator utility, not a solver acceptance or refinement policy.
class ExpansionProductSum {
public:
    void add(double value) {
        if (!std::isfinite(value)) throw std::overflow_error("Nonfinite expansion summand");
        const auto first = twoSum(high_, value);
        const auto second = twoSum(first.high, low_ + first.low);
        high_ = second.high;
        low_ = second.low;
        if (!std::isfinite(high_) || !std::isfinite(low_))
            throw std::overflow_error("Expansion sum overflow");
    }

    void addProduct(ExpansionCoefficient coefficient, double direction) {
        product(coefficient.high, direction);
        product(coefficient.low, direction);
    }

    // Corrections below the ULP of a Newton direction must survive as well.
    void addProduct(ExpansionCoefficient coefficient, ExpansionCoefficient direction) {
        addProduct(coefficient, direction.high);
        addProduct(coefficient, direction.low);
    }

    double value() const { return high_ + low_; }
    ExpansionCoefficient parts() const { return {high_, low_}; }

private:
    static ExpansionCoefficient twoSum(double a, double b) {
        const double s = a + b;
        const double z = s - a;
        return {s, (a - (s - z)) + (b - z)};
    }

    void product(double a, double b) {
        const double p = a * b;
        add(p);
        add(std::fma(a, b, -p));
    }

    double high_ = 0.0;
    double low_ = 0.0;
};

} // namespace vela::experimental
