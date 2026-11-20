#pragma once

#include "vela/core/ArithmeticPrecision.h"
#include "vela/core/Types.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <type_traits>

namespace vela::detail {

constexpr std::size_t Tri3LocalPotentialDofCount = 9;

template<class Scalar>
struct BasicTri3LocalForwardDual {
    Scalar value = 0.0;
    std::array<Scalar, Tri3LocalPotentialDofCount> derivative{};

    BasicTri3LocalForwardDual() = default;
    BasicTri3LocalForwardDual(Scalar scalar) : value(scalar) {}

    static BasicTri3LocalForwardDual variable(Scalar scalar, std::size_t localDof)
    {
        BasicTri3LocalForwardDual result(scalar);
        result.derivative.at(localDof) = 1.0;
        return result;
    }


    friend BasicTri3LocalForwardDual operator+(
        const BasicTri3LocalForwardDual& lhs,
        const BasicTri3LocalForwardDual& rhs)
    {
        BasicTri3LocalForwardDual result(lhs.value + rhs.value);
        for (std::size_t i = 0; i < result.derivative.size(); ++i)
            result.derivative[i] = lhs.derivative[i] + rhs.derivative[i];
        return result;
    }

    friend BasicTri3LocalForwardDual operator-(
        const BasicTri3LocalForwardDual& lhs,
        const BasicTri3LocalForwardDual& rhs)
    {
        BasicTri3LocalForwardDual result(lhs.value - rhs.value);
        for (std::size_t i = 0; i < result.derivative.size(); ++i)
            result.derivative[i] = lhs.derivative[i] - rhs.derivative[i];
        return result;
    }

    friend BasicTri3LocalForwardDual operator-(const BasicTri3LocalForwardDual& operand)
    {
        BasicTri3LocalForwardDual result(-operand.value);
        for (std::size_t i = 0; i < result.derivative.size(); ++i)
            result.derivative[i] = -operand.derivative[i];
        return result;
    }

    friend BasicTri3LocalForwardDual operator*(
        const BasicTri3LocalForwardDual& lhs,
        const BasicTri3LocalForwardDual& rhs)
    {
        BasicTri3LocalForwardDual result(lhs.value * rhs.value);
        for (std::size_t i = 0; i < result.derivative.size(); ++i) {
            result.derivative[i] =
                lhs.derivative[i] * rhs.value + lhs.value * rhs.derivative[i];
        }
        return result;
    }

    friend BasicTri3LocalForwardDual operator/(
        const BasicTri3LocalForwardDual& lhs,
        const BasicTri3LocalForwardDual& rhs)
    {
        BasicTri3LocalForwardDual result(lhs.value / rhs.value);
        for (std::size_t i = 0; i < result.derivative.size(); ++i) {
            // d(a/b) = (da - (a/b) db)/b. Squaring b can underflow or
            // overflow even when the quotient and all derivatives are finite
            // (notably for depleted carrier-density ratios).
            result.derivative[i] =
                (lhs.derivative[i] - result.value * rhs.derivative[i]) / rhs.value;
        }
        return result;
    }

};

using Tri3LocalForwardDual = BasicTri3LocalForwardDual<Real>;

template<class Scalar>
BasicTri3LocalForwardDual<Scalar> dualAbs(const BasicTri3LocalForwardDual<Scalar>& operand)
{
    if (operand.value > 0) return operand;
    if (operand.value < 0) return -operand;
    return {};
}

// A structural zero stays zero even when exp(x) is outside the scalar range.
template<class Scalar>
Scalar dualExponentialDirection(Scalar x, Scalar exponential, Scalar direction)
{
    if (direction == 0.)
        return 0.;
    const Scalar ordinary = exponential * direction;
    if (exponential >= std::numeric_limits<Scalar>::min() &&
        std::isfinite(ordinary) && ordinary != 0.)
        return ordinary;
    if (!std::isfinite(x) || !std::isfinite(direction))
        return ordinary;

    // Outside this bound no nonzero finite direction can recover a value.
    // It also bounds the exponent passed to ldexp.
    constexpr Scalar exponentLimit = 2 * std::numeric_limits<Scalar>::max_exponent;
    using Work = std::conditional_t<(sizeof(Scalar) > sizeof(ArithmeticReal)), Scalar, ArithmeticReal>;
    if (x <= -exponentLimit)
        return std::copysign(Scalar{0}, direction);
    if (x >= exponentLimit)
        return std::copysign(std::numeric_limits<Scalar>::infinity(), direction);

    // Reduce x by an integer multiple of ln(2). The high part has enough
    // trailing zero bits for an exact product with expPower in double.
    constexpr Work ln2High = 0x1.62e42fee00000p-1L;
    constexpr Work ln2Low = 1.908214929270587816144265680755e-10L;
    const Work argument = x;
    const int expPower = static_cast<int>(argument * std::numbers::log2e_v<Work>);
    const auto reduced = (argument - expPower * static_cast<Work>(ln2High)) - expPower * static_cast<Work>(ln2Low);
    int directionPower;
    const auto directionFraction = std::frexp(static_cast<Work>(direction), &directionPower);
    return static_cast<Scalar>(std::ldexp(std::exp(reduced) * directionFraction,
                                       expPower + directionPower));
}

template<class Scalar>
BasicTri3LocalForwardDual<Scalar> dualExp(const BasicTri3LocalForwardDual<Scalar>& operand)
{
    BasicTri3LocalForwardDual<Scalar> result(std::exp(operand.value));
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = dualExponentialDirection(operand.value,result.value,operand.derivative[i]);
    return result;
}

template<class Scalar>
BasicTri3LocalForwardDual<Scalar> dualExpm1(const BasicTri3LocalForwardDual<Scalar>& operand)
{
    BasicTri3LocalForwardDual<Scalar> result(std::expm1(operand.value));
    const Scalar slope = std::exp(operand.value);
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = dualExponentialDirection(operand.value,slope,operand.derivative[i]);
    return result;
}

template<class Scalar>
BasicTri3LocalForwardDual<Scalar> dualSqrt(const BasicTri3LocalForwardDual<Scalar>& operand)
{
    if (operand.value <= 0.0)
        return BasicTri3LocalForwardDual<Scalar>{};
    BasicTri3LocalForwardDual<Scalar> result(std::sqrt(operand.value));
    const Scalar slope = 0.5 / result.value;
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = slope * operand.derivative[i];
    return result;
}

template<class Scalar>
BasicTri3LocalForwardDual<Scalar> dualPow(
    const BasicTri3LocalForwardDual<Scalar>& operand,
    std::type_identity_t<Scalar> exponent)
{
    if (operand.value <= 0.0)
        return BasicTri3LocalForwardDual<Scalar>(std::pow(operand.value, exponent));
    BasicTri3LocalForwardDual<Scalar> result(std::pow(operand.value, exponent));
    const Scalar slope =
        exponent * std::pow(operand.value, exponent - 1.0);
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = slope * operand.derivative[i];
    return result;
}

} // namespace vela::detail
