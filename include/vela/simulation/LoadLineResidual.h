#pragma once

#include "vela/core/ArithmeticPrecision.h"
#include "vela/core/Types.h"
#include <cmath>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace vela {

/// One physical load-line evaluation shared by all numerical coordinates.
/// Retain this value for physical audits; never reconstruct it by subtracting
/// separately scaled voltage/current terms near a root.
template<class Scalar>
struct BasicLoadLineResidual {
    Scalar voltage_V;

    static BasicLoadLineResidual evaluate(Scalar coefficient, Scalar innerVoltage_V,
        Scalar resistance_ohm_um, Scalar current_A_per_um, Scalar outerVoltage_V,
        std::optional<ArithmeticReal> decimalResistance_ohm_um = std::nullopt)
    {
        if (decimalResistance_ohm_um) {
            const auto drop = std::fma(static_cast<ArithmeticReal>(coefficient),
                                       static_cast<ArithmeticReal>(innerVoltage_V),
                                       -static_cast<ArithmeticReal>(outerVoltage_V));
            const auto residual = std::fma(*decimalResistance_ohm_um,
                                           static_cast<ArithmeticReal>(current_A_per_um), drop);
            return {static_cast<Scalar>(residual)};
        }
        if constexpr (!std::is_same_v<Scalar,Real>)
            return {std::fma(resistance_ohm_um, current_A_per_um,
                             std::fma(coefficient, innerVoltage_V, -outerVoltage_V))};
        return {coefficient * innerVoltage_V +
                resistance_ohm_um * current_A_per_um - outerVoltage_V};
    }

    static Scalar currentDerivative(Scalar resistance, Scalar derivative,
                                  std::optional<ArithmeticReal> decimalResistance)
    {
        return decimalResistance
            ? static_cast<Scalar>(*decimalResistance * static_cast<ArithmeticReal>(derivative))
            : resistance * derivative;
    }

    static Scalar voltageDerivative(Scalar coefficient, Scalar resistance, Scalar derivative,
                                  std::optional<ArithmeticReal> decimalResistance)
    {
        if (decimalResistance)
            return static_cast<Scalar>(std::fma(*decimalResistance, static_cast<ArithmeticReal>(derivative),
                                                static_cast<ArithmeticReal>(coefficient)));
        if constexpr (!std::is_same_v<Scalar,Real>)
            return std::fma(resistance, derivative, coefficient);
        return coefficient + resistance * derivative;
    }

    Scalar physicalRhs() const { return -voltage_V; }

    Scalar inVoltageUnits(Scalar voltageScale_V) const
    {
        if (!(voltageScale_V > 0) || !std::isfinite(voltageScale_V))
            throw std::invalid_argument("load-line voltage scale must be positive and finite");
        return voltage_V / voltageScale_V;
    }
};

using LoadLineResidual = BasicLoadLineResidual<Real>;
using ArithmeticLoadLineResidual = BasicLoadLineResidual<ArithmeticReal>;

} // namespace vela
