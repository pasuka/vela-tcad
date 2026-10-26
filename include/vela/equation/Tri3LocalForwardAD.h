#pragma once

#include "vela/core/Types.h"

#include <array>
#include <cmath>
#include <cstddef>

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

inline Tri3LocalForwardDual dualExp(const Tri3LocalForwardDual& operand)
{
    Tri3LocalForwardDual result(std::exp(operand.value));
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = result.value * operand.derivative[i];
    return result;
}

inline Tri3LocalForwardDual dualExpm1(const Tri3LocalForwardDual& operand)
{
    Tri3LocalForwardDual result(std::expm1(operand.value));
    const Real slope = std::exp(operand.value);
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = slope * operand.derivative[i];
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
