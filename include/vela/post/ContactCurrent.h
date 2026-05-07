#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/solver/GummelSolver.h"
#include <string>

namespace vela {

struct ContactCurrentResult {
    Real electronCurrent = 0.0;
    Real holeCurrent = 0.0;
    Real totalCurrent = 0.0;
};

class ContactCurrent {
public:
    static ContactCurrentResult compute(const DeviceMesh& mesh,
                                        const MaterialDatabase& matdb,
                                        const DopingModel& doping,
                                        const DDSolution& solution,
                                        const std::string& contactName,
                                        const MobilityModelConfig& mobilityConfig = {});
};

} // namespace vela
