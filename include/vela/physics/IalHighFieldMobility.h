#pragma once
#include "vela/core/Types.h"

namespace vela {
/// T-2022.03 Caughey-Thomas, alpha=0, Vsat_Formula=1.
/// Explicit coefficients; no changes to existing coupled-solver defaults.
struct IalHighFieldParameters {
    Real saturationVelocity300_m_per_s = 1.07e5;
    Real beta300 = 1.109;
    Real saturationVelocityTemperatureExponent = .87;
    Real betaTemperatureExponent = .66;
};
struct IalHighFieldResult {
    Real mobility_m2_per_Vs;
    Real lowFieldDerivative;
    Real drivingFieldDerivative_m3_per_V2s;
    /// Partial at fixed low-field mobility and driving field.
    Real temperatureDerivative_m2_per_Vs_K;
};
/// Local diagnostic primitive. Caller owns interface/contact drive selection and
/// must chain the low-field temperature derivative for a total derivative.
/// At zero drive the selected drive slope is zero (a cusp if beta <= 1).
IalHighFieldResult evaluateIalHighFieldMobility(Real lowField_m2_per_Vs,
    Real drivingField_V_per_m, Real temperature_K, const IalHighFieldParameters&);
/// Independent explicit-partial candidate; reference AD entry above is retained.
IalHighFieldResult evaluateIalHighFieldMobilityExplicit(Real lowField_m2_per_Vs,
    Real drivingField_V_per_m, Real temperature_K, const IalHighFieldParameters&);
/// Same candidate arithmetic without partial evaluation.
Real evaluateIalHighFieldMobilityValue(Real lowField_m2_per_Vs,
    Real drivingField_V_per_m, Real temperature_K, const IalHighFieldParameters&);
} // namespace vela
