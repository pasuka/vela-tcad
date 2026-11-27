#pragma once

#include "vela/core/ArithmeticPrecision.h"
#include <nlohmann/json_fwd.hpp>

namespace vela {

struct SemiconductorArithmetic {
    ArithmeticReal bandgap_eV;
    ArithmeticReal affinity_eV;
    ArithmeticReal electronsDOS_cm3;
    ArithmeticReal holesDOS_cm3;
    ArithmeticReal intrinsicDensity_cm3;
    ArithmeticReal referencePotential_V;
};

/// Evaluate the active Formula 1 parameters from their original PAR literals.
SemiconductorArithmetic evaluateSemiconductorArithmetic(
    const nlohmann::json& input, ArithmeticReal temperature_K);

} // namespace vela
