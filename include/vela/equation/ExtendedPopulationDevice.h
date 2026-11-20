#pragma once

#include "vela/equation/IalPopulationContinuity.h"
#include "vela/solver/ExtendedBoltzmannColdStart.h"

namespace vela {

struct ExtendedPopulationDeviceSystem {
    ExtendedSparseMatrix jacobian;
    ExtendedVector residual;
    // Derivative with respect to the resistor-side device voltage in volts.
    ExtendedVector voltageDerivative;
    std::vector<int> prescribedPopulationRows;
    std::vector<int> retainedIndependentPotentialRows;
    std::vector<bool> electronFree, holeFree;
};

struct ExtendedPopulationNormalization {
    long double concentration_cm3;
    long double poissonRowFactor;
    long double continuityRowFactor;
};

// Merge physical Poisson and continuity rows using the same state, then apply
// contact and insulation identities. The Poisson stiffness must be assembled
// without Dirichlet elimination. When requested, express the resulting rows
// and columns in equation coordinates. Inputs remain in TCAD basic units.
ExtendedPopulationDeviceSystem assembleExtendedPopulationDevice(
    const BasicPoissonAssembler<long double>& poisson, const ExtendedPopulationState& state,
    const ExtendedBoltzmannMaterials& material, const ExtendedOhmicTargets& contacts,
    std::span<const Index> drivenContactNodes, const ExtendedIalPopulationContinuity& continuity,
    const ExtendedAveragePoissonGeometry* geometry = nullptr,
    const ExtendedPopulationNormalization* normalization = nullptr);

} // namespace vela
