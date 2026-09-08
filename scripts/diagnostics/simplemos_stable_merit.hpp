#pragma once
#include "vela/core/Types.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <bit>
#include <cstdint>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace simplemos_stable_merit {
struct Comparison {
    bool accepted = false;
    bool exactFallback = false;
    int sign = 0;
    long double delta = 0;
    long double roundoffBound = 0;
};

// Every finite binary64 square is an integer multiple of 2^-2148.
// Use this exact dyadic representation only when the fast sum cannot decide.
inline boost::multiprecision::cpp_int squaredInteger(double x) {
    const auto bits = std::bit_cast<std::uint64_t>(x);
    const unsigned exponent = static_cast<unsigned>((bits >> 52) & 0x7ffU);
    std::uint64_t mantissa = bits & ((std::uint64_t{1} << 52) - 1);
    if (exponent) mantissa |= std::uint64_t{1} << 52;
    boost::multiprecision::cpp_int square = mantissa;
    square *= mantissa;
    if (exponent) square <<= 2 * exponent - 2;
    return square;
}

inline Comparison compare(const vela::VectorXd& current, const vela::VectorXd& trial) {
    if (current.size() != trial.size()) throw std::invalid_argument("merit vector size mismatch");
    Comparison out;
    if (!current.allFinite() || !trial.allFinite()) return out;
    long double sum = 0, correction = 0, absolute = 0;
    bool currentZero = true, trialZero = true;
    for (int i = 0; i < current.size(); ++i) {
        const long double x = current(i), y = trial(i);
        currentZero = currentZero && x == 0;
        trialZero = trialZero && y == 0;
        const long double term = (y - x) * (y + x);
        const long double next = sum + term;
        correction += std::abs(sum) >= std::abs(term)
            ? (sum - next) + term : (term - next) + sum;
        sum = next;
        absolute += std::abs(term);
    }
    out.delta = sum + correction;
    out.roundoffBound = 32 * std::numeric_limits<long double>::epsilon() * absolute;
    if (std::isfinite(out.delta) && std::isfinite(out.roundoffBound)
        && std::abs(out.delta) > out.roundoffBound) {
        out.sign = out.delta < 0 ? -1 : 1;
    } else {
        out.exactFallback = true;
        boost::multiprecision::cpp_int exact = 0;
        for (int i = 0; i < current.size(); ++i)
            exact += squaredInteger(trial(i)) - squaredInteger(current(i));
        out.sign = exact < 0 ? -1 : exact > 0 ? 1 : 0;
    }
    // In this frozen profile alpha>0 and sufficientDecrease>0. The existing
    // (norm < base || norm <= target) rule is precisely strict decrease,
    // with the existing zero-to-zero exception retained.
    out.accepted = out.sign < 0 || (currentZero && trialZero);
    return out;
}
} // namespace simplemos_stable_merit
