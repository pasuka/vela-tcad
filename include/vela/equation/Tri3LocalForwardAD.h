#pragma once

#include "vela/core/ArithmeticPrecision.h"
#include "vela/core/Types.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>

namespace vela::detail {

constexpr std::size_t Tri3LocalPotentialDofCount = 9;

struct Tri3LocalForwardDual {
    Real value = 0.0;
    std::array<Real, Tri3LocalPotentialDofCount> derivative{};

    Tri3LocalForwardDual() = default;
    Tri3LocalForwardDual(Real scalar) : value(scalar) {}

    static Tri3LocalForwardDual variable(Real scalar, std::size_t localDof)
    {
        Tri3LocalForwardDual result(scalar);
        result.derivative.at(localDof) = 1.0;
        return result;
    }
};

inline Tri3LocalForwardDual operator+(
    const Tri3LocalForwardDual& lhs,
    const Tri3LocalForwardDual& rhs)
{
    Tri3LocalForwardDual result(lhs.value + rhs.value);
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = lhs.derivative[i] + rhs.derivative[i];
    return result;
}

inline Tri3LocalForwardDual operator-(
    const Tri3LocalForwardDual& lhs,
    const Tri3LocalForwardDual& rhs)
{
    Tri3LocalForwardDual result(lhs.value - rhs.value);
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = lhs.derivative[i] - rhs.derivative[i];
    return result;
}

inline Tri3LocalForwardDual operator-(const Tri3LocalForwardDual& operand)
{
    Tri3LocalForwardDual result(-operand.value);
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = -operand.derivative[i];
    return result;
}

inline Tri3LocalForwardDual operator*(
    const Tri3LocalForwardDual& lhs,
    const Tri3LocalForwardDual& rhs)
{
    Tri3LocalForwardDual result(lhs.value * rhs.value);
    for (std::size_t i = 0; i < result.derivative.size(); ++i) {
        result.derivative[i] =
            lhs.derivative[i] * rhs.value + lhs.value * rhs.derivative[i];
    }
    return result;
}

inline Tri3LocalForwardDual operator/(
    const Tri3LocalForwardDual& lhs,
    const Tri3LocalForwardDual& rhs)
{
    Tri3LocalForwardDual result(lhs.value / rhs.value);
    for (std::size_t i = 0; i < result.derivative.size(); ++i) {
        // d(a/b) = (da - (a/b) db)/b. Squaring b can underflow or
        // overflow even when the quotient and all derivatives are finite
        // (notably for depleted carrier-density ratios).
        result.derivative[i] =
            (lhs.derivative[i] - result.value * rhs.derivative[i]) / rhs.value;
    }
    return result;
}

inline Tri3LocalForwardDual dualAbs(const Tri3LocalForwardDual& operand)
{
    if (operand.value > 0.0)
        return operand;
    if (operand.value < 0.0)
        return -operand;
    return Tri3LocalForwardDual{};
}

// A structural zero stays zero even when exp(x) is outside double range.
inline Real dualExponentialDirection(Real x, Real exponential, Real direction)
{
    if (direction == 0.)
        return 0.;
    const Real ordinary = exponential * direction;
    if (exponential >= std::numeric_limits<Real>::min() &&
        std::isfinite(ordinary) && ordinary != 0.)
        return ordinary;
    if (!std::isfinite(x) || !std::isfinite(direction))
        return ordinary;

    // Outside this bound no nonzero finite double direction can recover a
    // representable result. It also bounds the exponent passed to ldexp.
    constexpr Real exponentLimit = 2048.;
    if (x <= -exponentLimit)
        return std::copysign(0., direction);
    if (x >= exponentLimit)
        return std::copysign(std::numeric_limits<Real>::infinity(), direction);

    // Reduce x by an integer multiple of ln(2). The high part has enough
    // trailing zero bits for an exact product with expPower in double.
    constexpr auto ln2High = VELA_ARITHMETIC_LITERAL(0x1.62e42fee00000p-1);
    constexpr auto ln2Low = VELA_ARITHMETIC_LITERAL(1.908214929270587816144265680755e-10);
    const auto argument = static_cast<ArithmeticReal>(x);
    const int expPower = static_cast<int>(argument * std::numbers::log2e_v<ArithmeticReal>);
    const auto reduced = (argument - expPower * ln2High) - expPower * ln2Low;
    int directionPower;
    const auto directionFraction = std::frexp(static_cast<ArithmeticReal>(direction), &directionPower);
    return static_cast<Real>(std::ldexp(std::exp(reduced) * directionFraction,
                                       expPower + directionPower));
}

inline Tri3LocalForwardDual dualExp(const Tri3LocalForwardDual& operand)
{
    Tri3LocalForwardDual result(std::exp(operand.value));
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = dualExponentialDirection(operand.value,result.value,operand.derivative[i]);
    return result;
}

inline Tri3LocalForwardDual dualExpm1(const Tri3LocalForwardDual& operand)
{
    Tri3LocalForwardDual result(std::expm1(operand.value));
    const Real slope = std::exp(operand.value);
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = dualExponentialDirection(operand.value,slope,operand.derivative[i]);
    return result;
}

inline Tri3LocalForwardDual dualSqrt(const Tri3LocalForwardDual& operand)
{
    if (operand.value <= 0.0)
        return Tri3LocalForwardDual{};
    Tri3LocalForwardDual result(std::sqrt(operand.value));
    const Real slope = 0.5 / result.value;
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = slope * operand.derivative[i];
    return result;
}

inline Tri3LocalForwardDual dualPow(
    const Tri3LocalForwardDual& operand,
    Real exponent)
{
    if (operand.value <= 0.0)
        return Tri3LocalForwardDual(std::pow(operand.value, exponent));
    Tri3LocalForwardDual result(std::pow(operand.value, exponent));
    const Real slope =
        exponent * std::pow(operand.value, exponent - 1.0);
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = slope * operand.derivative[i];
    return result;
}

} // namespace vela::detail
