#pragma once

#include <limits>
#include <stdexcept>
#include <string_view>

namespace vela {

enum class ArithmeticPrecision { Double, Extended80 };

#ifdef VELA_ARITHMETIC_EXTENDED80
using ArithmeticReal = long double;
static_assert(std::numeric_limits<ArithmeticReal>::radix == 2 &&
              std::numeric_limits<ArithmeticReal>::digits == 64 &&
              std::numeric_limits<ArithmeticReal>::max_exponent == 16384,
              "Extended80 requires an 80-bit floating-point format");
inline constexpr auto compiledArithmeticPrecision = ArithmeticPrecision::Extended80;
#define VELA_ARITHMETIC_LITERAL(value) value##L
#else
using ArithmeticReal = double;
inline constexpr auto compiledArithmeticPrecision = ArithmeticPrecision::Double;
#define VELA_ARITHMETIC_LITERAL(value) value
#endif

inline ArithmeticPrecision arithmeticPrecisionFromString(std::string_view value)
{
    if (value == "double") return ArithmeticPrecision::Double;
    if (value == "extended80") return ArithmeticPrecision::Extended80;
    throw std::invalid_argument("solver.arithmetic_precision must be double or extended80");
}

inline std::string_view arithmeticPrecisionName(ArithmeticPrecision precision)
{
    return precision == ArithmeticPrecision::Extended80 ? "extended80" : "double";
}

} // namespace vela
