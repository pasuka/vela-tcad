#pragma once
#include "vela/core/ArithmeticPrecision.h"
#include <cmath>
namespace vela {
// Preserve ordinary binary64 arithmetic; rescue only an avoidable intermediate
// overflow/underflow. A truly unrepresentable final result remains nonfinite/zero.
inline double scaledProduct(double a,double b,double c) {
    const double ab=a*b, value=ab*c;
    if(std::isfinite(ab)&&std::isfinite(value)&&
       (value!=0. || a==0. || b==0. || c==0.))return value;
    if constexpr (compiledArithmeticPrecision == ArithmeticPrecision::Extended80)
        return static_cast<double>(static_cast<ArithmeticReal>(a)*b*c);
    if (!(std::isfinite(a) && std::isfinite(b) && std::isfinite(c)))
        return value;
    int ea, eb, ec;
    const double ma = std::frexp(a, &ea);
    const double mb = std::frexp(b, &eb);
    const double mc = std::frexp(c, &ec);
    return std::scalbn(ma * mb * mc, ea + eb + ec);
}

inline long double scaledProduct(long double a,long double b,long double c) {
    const long double ab=a*b, value=ab*c;
    if (std::isfinite(ab) && std::isfinite(value) &&
        (value!=0.L || a==0.L || b==0.L || c==0.L)) return value;
    if (!(std::isfinite(a) && std::isfinite(b) && std::isfinite(c))) return value;
    int ea, eb, ec;
    const auto ma = std::frexp(a,&ea);
    const auto mb = std::frexp(b,&eb);
    const auto mc = std::frexp(c,&ec);
    return std::scalbn(ma*mb*mc,ea+eb+ec);
}
}
