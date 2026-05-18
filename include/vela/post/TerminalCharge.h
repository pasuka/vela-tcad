#pragma once

#include "vela/core/PhysicalConstants.h"
#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/GummelSolver.h"
#include <string>
#include <vector>

namespace vela {

struct TerminalChargeConfig {
    std::string name;
    std::string contact;
    std::vector<std::string> regions;
    Real contactRadius = 0.0;
    bool includeMobileCharge = true;
    bool includeIonizedDopants = true;
    bool perMeter = true;
    Real depth_m = 1.0;
};

struct TerminalChargeResult {
    Real charge = 0.0;
    bool perMeter = true;
};

class TerminalCharge {
public:
    TerminalCharge(const DeviceMesh& mesh, const DopingModel& doping);

    TerminalChargeResult compute(const DDSolution& solution,
                                 const TerminalChargeConfig& config) const;

    static TerminalChargeResult compute(const DeviceMesh& mesh,
                                        const DopingModel& doping,
                                        const DDSolution& solution,
                                        const TerminalChargeConfig& config);

private:
    const DeviceMesh& mesh_;
    const DopingModel& doping_;
};

} // namespace vela
