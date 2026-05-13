#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/post/ContactCurrent.h"
#include "vela/solver/GummelSolver.h"
#include <string>
#include <vector>

namespace vela {

struct DCSweepConfig {
    std::string contact;
    Real start = 0.0;
    Real stop = 0.0;
    Real step = 0.0;
    Real minStep = 0.0;
    Real maxStep = 0.0;
    Real growthFactor = 1.0;
    Real shrinkFactor = 0.5;
    int maxRetries = 5;
    bool stopOnFailure = true;
    std::string currentContact;
    bool writeVtk = false;
    std::string vtkPrefix;
    std::string csvFile = "dc_sweep.csv";
};

struct DCSweepPoint {
    Real voltage = 0.0;
    Real electronCurrent = 0.0;
    Real holeCurrent = 0.0;
    Real totalCurrent = 0.0;
    bool converged = false;
    int iterations = 0;
    Real attemptedStep = 0.0;
    Real acceptedStep = 0.0;
    int retryCount = 0;
};

struct DCSweepResult {
    DeviceMesh mesh;
    std::vector<DCSweepPoint> points;
};

class DCSweep {
public:
    std::vector<DCSweepPoint> run(const std::string& configFile) const;
    DCSweepResult runWithResult(const std::string& configFile) const;
};

} // namespace vela
