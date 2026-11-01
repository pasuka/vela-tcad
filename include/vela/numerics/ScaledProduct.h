#pragma once
#include <cmath>
namespace vela {
// Preserve ordinary binary64 arithmetic; rescue only an avoidable intermediate
// overflow/underflow. A truly unrepresentable final result remains nonfinite/zero.
inline double scaledProduct(double a,double b,double c) {
    const double ab=a*b, value=ab*c;
    if(std::isfinite(ab)&&std::isfinite(value)&&
       (value!=0. || a==0. || b==0. || c==0.))return value;
    return static_cast<double>(static_cast<long double>(a)*b*c);
}
}
