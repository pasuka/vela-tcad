#pragma once

#include "vela/core/Types.h"
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
};

class DCSweep {
public:
    std::vector<DCSweepPoint> run(const std::string& configFile) const;
};

} // namespace vela
