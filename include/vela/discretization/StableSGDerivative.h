#pragma once
#include <array>
#include <cmath>
#include <stdexcept>
#include "vela/discretization/Bernoulli.h"

namespace vela::stable_sg {
// Boltzmann transport with state-independent ni (including doping-only BGN):
// differentiate the stable edge flux instead of subtracting nearly cancelling
// drift/diffusion derivatives. B and dB must use the same effective SG argument
// as the flux: -eta for electrons, +eta for holes. ni and mobility are held
// fixed in these partials; any state-dependent chain terms are separate.
inline std::array<double,2> psiDerivative(double flux,double B,double dB,double vt,bool electron) {
    if (!(B>0.) || !std::isfinite(B)) throw std::runtime_error("Unresolved Bernoulli logarithmic derivative");
    const double logDerivative=dB/B;
    const double f=flux/vt;
    if(electron)return {f*(1.+logDerivative),-f*logDerivative};
    return {-f*(1.+logDerivative),f*logDerivative};
}
}
