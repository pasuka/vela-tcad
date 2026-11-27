#pragma once

#include "vela/physics/MobilityModel.h"

namespace vela {

struct ExtendedPhuMobEvaluation {
    long double mobility;
    ExtendedPhuMobLogDensityDerivatives slopes;
};

ExtendedPhuMobEvaluation evaluateExtendedPhuMobWithLogDensityDerivatives(
    CarrierType carrier, const ExtendedPhuMobScalarState& state,
    const ExtendedPhuMobParameters& parameters);

} // namespace vela
