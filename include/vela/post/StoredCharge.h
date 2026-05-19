#pragma once

#include "vela/core/PhysicalConstants.h"
#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/solver/GummelSolver.h"
#include <string>
#include <vector>

namespace vela {

struct StoredChargeConfig {
    std::vector<std::string> regions;
    bool perMeter = true;
    Real depth_m = 1.0;
};

struct StoredChargeResult {
    Real charge = 0.0;
    bool perMeter = true;
};

class StoredCharge {
public:
    explicit StoredCharge(const DeviceMesh& mesh);

    StoredChargeResult compute(const DDSolution& solution,
                               const StoredChargeConfig& config) const;

private:
    const DeviceMesh& mesh_;
};

} // namespace vela
