#pragma once

#include "vela/core/Types.h"
#include <cmath>
#include <stdexcept>

namespace vela {

/// One physical load-line evaluation shared by all numerical coordinates.
/// Retain this value for physical audits; never reconstruct it by subtracting
/// separately scaled voltage/current terms near a root.
struct LoadLineResidual {
    Real voltage_V;

    static LoadLineResidual evaluate(Real coefficient, Real innerVoltage_V,
                                     Real resistance_ohm_um,
                                     Real current_A_per_um, Real outerVoltage_V)
    {
        return {coefficient * innerVoltage_V +
                resistance_ohm_um * current_A_per_um - outerVoltage_V};
    }

    Real physicalRhs() const { return -voltage_V; }

    Real inVoltageUnits(Real voltageScale_V) const
    {
        if (!(voltageScale_V > 0) || !std::isfinite(voltageScale_V))
            throw std::invalid_argument("load-line voltage scale must be positive and finite");
        return voltage_V / voltageScale_V;
    }
};

} // namespace vela
