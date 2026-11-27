#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace vela {

using Extended80Words = std::array<std::uint64_t, 2>;
inline constexpr int extended80ExponentOffset = 32768;

inline void requireExtended80Format() {
    using Limits = std::numeric_limits<long double>;
    if (Limits::radix != 2 || Limits::digits != 64 || Limits::max_exponent != 16384)
        throw std::invalid_argument("extended80 population archive requires a compatible long double");
}

// Store normalized significand and biased exponent plus sign. Encode numbers,
// rather than object bytes, so x87 padding, host byte order, and ABI are irrelevant.
inline Extended80Words encodeExtended80(long double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("nonfinite extended80 population state");
    const std::uint64_t sign = std::signbit(value) ? 1 : 0;
    if (value == 0) return {0, sign};
    int exponent = 0;
    const auto fraction = std::frexp(std::abs(value), &exponent);
    return {static_cast<std::uint64_t>(std::ldexp(fraction, 64)),
            (std::uint64_t(exponent + extended80ExponentOffset) << 1) | sign};
}

inline long double decodeExtended80(Extended80Words words) {
    const auto [significand, signExponent] = words;
    if (significand == 0) {
        if (signExponent > 1) throw std::invalid_argument("invalid extended80 zero");
        return signExponent == 1 ? -0.0L : 0.0L;
    }
    const auto biased = signExponent >> 1;
    if (biased < extended80ExponentOffset - 16444 || biased > extended80ExponentOffset + 16384 ||
        significand < (std::uint64_t{1} << 63))
        throw std::invalid_argument("invalid extended80 exponent or significand");
    const auto magnitude = std::ldexp(static_cast<long double>(significand),
                                      static_cast<int>(biased) - extended80ExponentOffset - 64);
    const auto value = (signExponent & 1) ? -magnitude : magnitude;
    if (encodeExtended80(value) != words) throw std::invalid_argument("noncanonical extended80 value");
    return value;
}

} // namespace vela
