#include "vela/physics/SRHRecombination.h"
#include <cmath>

namespace vela {

double srhRate(double n, double p, double ni, double taun, double taup)
{
    // Trap level at mid-gap: n1 = p1 = ni
    const double n1  = ni;
    const double p1  = ni;
    const double num = n * p - ni * ni;
    const double den = taup * (n + n1) + taun * (p + p1);
    if (std::abs(den) < 1.0e-100)
        return 0.0;
    return num / den;
}

} // namespace vela
