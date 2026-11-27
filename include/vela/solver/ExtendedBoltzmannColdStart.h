#pragma once

#include "vela/boundary/BoundaryCondition.h"
#include "vela/equation/ExtendedBoltzmannMaterials.h"
#include "vela/equation/PoissonAssembler.h"
#include "vela/solver/ExtendedPopulationState.h"
#include <optional>

namespace vela {

// Zero-bias Ohmic targets in internal potential and concentration coordinates.
// An empty potential denotes a node without a prescribed potential.
struct ExtendedOhmicTargets {
    std::vector<std::optional<long double>> potential;
    ExtendedVector electrons;
    ExtendedVector holes;
};

ExtendedOhmicTargets buildExtendedZeroBiasOhmicTargets(
    const DeviceMesh& mesh, const DopingModel& doping,
    const ExtendedEquilibriumMaterials& materials,
    const std::vector<ContactBoundarySpec>& contacts,
    long double thermalVoltage_V, bool dominantSignedContactMean);

struct ExtendedEquilibriumColdStartResult {
    // Populations are recovered from the final Poisson linearization, in the
    // same reference and units as the potential passed on to the coupled solve.
    ExtendedPopulationState state;
    long double residual;
    int updates;
    bool converged;
};

ExtendedEquilibriumColdStartResult solveExtendedEquilibriumColdStart(
    BasicPoissonAssembler<long double>& poisson, const DopingModel& doping,
    const ExtendedEquilibriumMaterials& materials, const ExtendedOhmicTargets& targets,
    long double thermalVoltage_V, long double concentrationScale_m3,
    int maxUpdates, long double absoluteResidualTolerance,
    const ExtendedAveragePoissonGeometry* geometry = nullptr);

// Seed a voltage stage in which each disconnected conducting component
// remains at a single electrochemical potential.
ExtendedEquilibriumColdStartResult solveExtendedEquilibriumPoissonSeed(
    BasicPoissonAssembler<long double>& poisson, const DopingModel& doping,
    const ExtendedEquilibriumMaterials& materials, const ExtendedOhmicTargets& targets,
    long double thermalVoltage_V, long double concentrationScale_m3,
    int maxUpdates, long double absoluteResidualTolerance,
    const ExtendedAveragePoissonGeometry* geometry,
    const ExtendedVector& initialPotential, const ExtendedVector& quasiFermi);

} // namespace vela
