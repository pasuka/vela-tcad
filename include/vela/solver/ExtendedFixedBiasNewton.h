#pragma once

#include "vela/equation/ExtendedPopulationDevice.h"
#include "vela/solver/DensityUpdateConvergence.h"
#include "vela/solver/PopulationState.h"
#include "vela/simulation/ExtendedBorderedLoadLine.h"
#include <functional>
#include <optional>
#include <string>

namespace vela {

struct ExtendedFixedBiasUpdateReferences {
    int digits = 5;
    long double potential_V = 0.025852L;
    long double electron_cm3 = 1e10L;
    long double hole_cm3 = 1e10L;
};

struct ExtendedRhsGrowthControl {
    long double factor;
    bool enforce;
};

struct ExtendedRhsGrowthReport {
    long double initialNorm;
    long double peakNorm;
    std::optional<int> firstExceededUpdate;
};

struct ExtendedFixedBiasNewtonConfig {
    int maxUpdates;
    int maxBacktracks;
    int refinementCorrections;
    long double damping;
    long double maxLinearBackwardError;
    long double maxScaledLinearDefect;
    ExtendedFixedBiasUpdateReferences updateError;
    PopulationUpdateConfig protection;
    bool lineSearch;
    std::optional<BasicPopulationTrialLimits<long double>> exactProtection;
    std::optional<ExtendedRhsGrowthControl> rhsGrowth = std::nullopt;
    bool refineToLinearGoal = false;
    ExtendedLinearDefectMethod linearDefectMethod = ExtendedLinearDefectMethod::mp100;
};

struct ExtendedFixedBiasNewtonResult {
    ExtendedPopulationState state;
    BasicDensityUpdateError<long double> updateError;
    int updates = 0;
    Eigen::Index factoredDimension = 0;
    long double scaledResidual = 0;
    long double linearBackwardError = 0;
    long double scaledLinearDefect = 0;
    long double acceptedAlpha = 0;
    bool converged = false;
    std::string failure;
    std::optional<ExtendedRhsGrowthReport> rhsGrowth;
};

using ExtendedFixedBiasAssembly =
    std::function<ExtendedPopulationDeviceSystem(const ExtendedPopulationState&)>;

/// Solve the population equations at the seed's fixed contact voltage.
/// Stop on the physical native-style update metric after an admissible trial.
/// The returned state is never narrowed and is distinct from publication of
/// an accepted sweep point by the calling stage.
ExtendedFixedBiasNewtonResult solveExtendedFixedBiasNewton(
    const ExtendedPopulationState& seed, const ExtendedFixedBiasNewtonConfig& config,
    const ExtendedFixedBiasAssembly& assemble);

struct ExtendedSeriesResistorNewtonConfig {
    long double resistance_ohm_um;
    long double tolerance_V;
    long double contactReference_V;
    std::size_t contactAliases;
    int direction;
    bool enforceLoadLine;
};

struct ExtendedSeriesResistorFrame {
    ExtendedPopulationDeviceSystem device;
    ExtendedVector portDerivative;
    long double current_A_per_um;
};

ExtendedFixedBiasNewtonResult solveExtendedSeriesResistorNewton(
    const ExtendedPopulationState& seed, const ExtendedFixedBiasNewtonConfig& config,
    const ExtendedSeriesResistorNewtonConfig& circuit,
    const std::function<ExtendedSeriesResistorFrame(const ExtendedPopulationState&)>& assemble);

} // namespace vela
