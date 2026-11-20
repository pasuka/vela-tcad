#include "vela/discretization/Bernoulli.h"
#include <cmath>

namespace vela {

double bernoulli(double x)
{
    if (std::abs(x) < 1.0e-10) {
        // Taylor series: B(x) = 1 - x/2 + x^2/12 - ...
        return 1.0 - x * 0.5 + x * x / 12.0;
    }
    if (x > 500.0) {
        // exp(x) >> 1 -> B(x) = x / (exp(x)-1) ~= x * exp(-x)
        return x * std::exp(-x);
    }
    if (x < -500.0) {
        // exp(x) ~= 0 -> B(x) = x / (exp(x)-1) ~= x / (-1) = -x
        return -x;
    }
    return x / std::expm1(x);
}

double bernoulliDerivative(double x)
{
    if (std::abs(x) < 1.0e-8)
        return -0.5 + x / 6.0 - x * x * x / 180.0;
    if (x > 500.0)
        return (1.0 - x) * std::exp(-x);
    if (x < -500.0)
        return -1.0;
    const double em1=std::expm1(x);
    const double ex=em1+1.0;
    return (em1-x*ex)/(em1*em1);
}

} // namespace vela
