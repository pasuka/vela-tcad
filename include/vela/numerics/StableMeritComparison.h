#pragma once
#include "vela/core/Types.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <bit>
#include <cstdint>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <array>

namespace vela::stable_merit {
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
    // For alpha>0 and sufficientDecrease>0, the existing
    // (norm < base || norm <= target) rule is precisely strict decrease,
    // with the existing zero-to-zero exception retained.
    out.accepted = out.sign < 0 || (currentZero && trialZero);
    return out;
}
// Preserve sum_b w_b * ||r_b||^2 / s_b^2 without rounding scaled vectors.
inline Comparison compareBlocks(const vela::VectorXd& current, const vela::VectorXd& trial,
                                const std::array<double,3>& scales,
                                const std::array<double,3>& weights) {
    if(current.size()!=trial.size() || current.size()%3!=0)
        throw std::invalid_argument("merit block size mismatch");
    std::array<double,3> safe{};
    for(int b=0;b<3;++b) {
        if(!std::isfinite(scales[b]) || !std::isfinite(weights[b]) || weights[b]<0.)
            throw std::invalid_argument("invalid merit block coefficient");
        safe[b]=std::max(std::abs(scales[b]),1e-300);
    }
    if(safe==std::array<double,3>{1.,1.,1.} && weights==std::array<double,3>{1.,1.,1.})
        return compare(current,trial);
    Comparison out;
    if(!current.allFinite() || !trial.allFinite())return out;
    const int n=current.size()/3;
    long double sum=0,correction=0,absolute=0;
    bool currentZero=true,trialZero=true;
    for(int b=0;b<3;++b) {
        if(weights[b]==0.)continue;
        const long double coefficient=static_cast<long double>(weights[b])/safe[b]/safe[b];
        for(int i=b*n;i<(b+1)*n;++i) {
            const long double x=current(i),y=trial(i);
            currentZero=currentZero && x==0;trialZero=trialZero && y==0;
            const long double term=((y-x)*(y+x))*coefficient;
            const long double next=sum+term;
            correction+=std::abs(sum)>=std::abs(term)?(sum-next)+term:(term-next)+sum;
            sum=next;absolute+=std::abs(term);
        }
    }
    out.delta=sum+correction;
    out.roundoffBound=64*std::numeric_limits<long double>::epsilon()*absolute;
    if(std::isfinite(out.delta) && std::isfinite(out.roundoffBound) && std::abs(out.delta)>out.roundoffBound)
        out.sign=out.delta<0?-1:1;
    else {
        out.exactFallback=true;
        using boost::multiprecision::cpp_int;
        using boost::multiprecision::cpp_rational;
        cpp_rational exact=0;
        for(int b=0;b<3;++b) {
            if(weights[b]==0.)continue;
            cpp_int difference=0;
            for(int i=b*n;i<(b+1)*n;++i)
                difference+=squaredInteger(trial(i))-squaredInteger(current(i));
            const cpp_rational denominator=cpp_rational(safe[b])*cpp_rational(safe[b]);
            exact+=cpp_rational(difference)*cpp_rational(weights[b])/denominator;
        }
        out.sign=exact<0?-1:exact>0?1:0;
    }
    out.accepted=out.sign<0 || (currentZero && trialZero);
    return out;
}
} // namespace vela::stable_merit
