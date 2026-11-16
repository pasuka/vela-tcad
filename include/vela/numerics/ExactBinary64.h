#pragma once

#include <boost/multiprecision/cpp_int.hpp>

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace vela::exact_binary64 {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::multiprecision::cpp_rational;

// Boost's floating-point rational constructor may promote binary64 inputs to
// long double. Build the dyadic fraction directly from the IEEE encoding.
inline Rational fromDouble(double value)
{
    if (!std::isfinite(value))
        throw std::invalid_argument("Expected a finite binary64 value");
    const auto bits = std::bit_cast<std::uint64_t>(value);
    const int exponent = static_cast<int>((bits >> 52) & 0x7ff);
    const auto mantissa = (bits & 0xfffffffffffffULL) | (exponent ? 0x10000000000000ULL : 0);
    if (mantissa == 0) return 0;
    Integer numerator = mantissa;
    if (bits >> 63) numerator = -numerator;
    const int power = (exponent == 0 ? 1 : exponent) - 1075;
    if (power >= 0) return Rational(numerator << power);
    return Rational(numerator, Integer(1) << -power);
}

// Linear refinement discards results below the binary64 normal range.
// Round the remaining quotient to nearest, ties to even, using integer
// arithmetic up to the final conversion of an exactly representable mantissa.
inline double toNormalDouble(const Rational& value)
{
    Integer numerator = boost::multiprecision::numerator(value);
    if (numerator == 0) return 0.;
    const bool negative = numerator < 0;
    if (negative) numerator = -numerator;
    Integer denominator = boost::multiprecision::denominator(value);
    int exponent = static_cast<int>(msb(numerator)) - static_cast<int>(msb(denominator));
    if (exponent >= 0 ? numerator < (denominator << exponent) : (numerator << -exponent) < denominator)
        --exponent;
    if (exponent < -1022) return 0.;
    if (exponent > 1023) {
        const double infinity = std::numeric_limits<double>::infinity();
        return negative ? -infinity : infinity;
    }
    const int shift = 52 - exponent;
    if (shift >= 0) numerator <<= shift;
    else denominator <<= -shift;
    Integer quotient, remainder;
    divide_qr(numerator, denominator, quotient, remainder);
    remainder <<= 1;
    if (remainder > denominator || (remainder == denominator && bit_test(quotient, 0)))
        ++quotient;
    const double result = std::scalbn(static_cast<double>(quotient.convert_to<std::uint64_t>()), exponent - 52);
    return negative ? -result : result;
}

} // namespace vela::exact_binary64
