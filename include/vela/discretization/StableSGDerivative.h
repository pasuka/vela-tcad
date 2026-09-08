#pragma once
#include <array>
#include <cmath>
#include <stdexcept>
#include "vela/discretization/Bernoulli.h"

namespace vela::stable_sg {
// Equal-ni Boltzmann transport: differentiate the already stable edge flux
// instead of subtracting large drift/diffusion derivatives near equilibrium.
inline std::array<double,2> psiDerivative(double flux,double B,double dB,double vt,bool electron) {
    if (!(B>0.) || !std::isfinite(B)) throw std::runtime_error("Unresolved Bernoulli logarithmic derivative");
    const double logDerivative=dB/B;
    const double f=flux/vt;
    if(electron)return {f*(1.+logDerivative),-f*logDerivative};
    return {-f*(1.+logDerivative),f*logDerivative};
}
}
