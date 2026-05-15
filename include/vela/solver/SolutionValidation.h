#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/solver/GummelSolver.h"
#include <string>
#include <unordered_map>
#include <vector>

namespace vela {

struct DDSolutionValidationOptions {
    Real carrierFloor = 1.0e-100;
    bool enforceMinimumCarrierDensity = false;
    Real minimumCarrierDensity = 0.0;
    bool checkContactQuasiFermiBias = true;
    Real contactPotentialAbsTolerance = 1.0e-8;
    Real contactPotentialRelTolerance = 1.0e-10;
};

struct DDSolutionFieldStats {
    Real min = 0.0;
    Real max = 0.0;
};

struct DDSolutionValidationResult {
    bool valid = true;
    std::vector<std::string> diagnostics;
    DDSolutionFieldStats psi;
    DDSolutionFieldStats phin;
    DDSolutionFieldStats phip;
    DDSolutionFieldStats n;
    DDSolutionFieldStats p;

    std::string diagnosticsString() const;
};

DDSolutionValidationResult validateDDSolution(
    const DDSolution& sol,
    const DeviceMesh& mesh,
    const std::unordered_map<std::string, Real>& contactBiases,
    const DDSolutionValidationOptions& options = {});

} // namespace vela
