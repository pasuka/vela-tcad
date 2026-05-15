#include "vela/physics/CarrierStatistics.h"
#include <cmath>
#include <algorithm>

namespace vela {

double electronDensity(double ni, double psi, double phin, double Vt)
{
    const double arg = std::clamp((psi - phin) / Vt, -500.0, 500.0);
    return ni * std::exp(arg);
}

double holeDensity(double ni, double psi, double phip, double Vt)
{
    const double arg = std::clamp((phip - psi) / Vt, -500.0, 500.0);
    return ni * std::exp(arg);
}


double intrinsicDensity(const Material& material, double temperature_K)
{
    return material.atTemperature(temperature_K).ni;
}

} // namespace vela
