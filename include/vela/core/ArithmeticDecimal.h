#pragma once

#include "vela/core/ArithmeticPrecision.h"
#include <cmath>
#include <locale>
#include <sstream>

namespace vela {

template<class Scalar = ArithmeticReal>
Scalar parseArithmeticDecimal(std::string_view text)
{
    // MinGW's floating-point from_chars can round long double through double.
    std::istringstream input{std::string(text)};
    input.imbue(std::locale::classic());
    Scalar value{};
    if (!(input >> std::noskipws >> value) || !std::isfinite(value) ||
        input.peek() != std::char_traits<char>::eof())
        throw std::invalid_argument("Arithmetic decimal is invalid or out of range");
    return value;
}

} // namespace vela
