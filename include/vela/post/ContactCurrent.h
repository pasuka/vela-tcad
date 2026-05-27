#pragma once

#include "vela/core/PhysicalConstants.h"
#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/solver/GummelSolver.h"
#include "vela/equation/DDAssembler.h" // for DDScalingSpec
#include <string>
#include <vector>

namespace vela {

struct ContactCurrentResult {
    Real electronCurrent = 0.0;
    Real electronDriftCurrent = 0.0;
    Real electronDiffusionCurrent = 0.0;
    Real holeCurrent = 0.0;
    Real holeDriftCurrent = 0.0;
    Real holeDiffusionCurrent = 0.0;
    Real totalCurrent = 0.0;
};

class ContactCurrent {
public:
    ContactCurrent(const DeviceMesh& mesh,
                   const MaterialDatabase& matdb,
                   const DopingModel& doping,
                   MobilityModelConfig mobilityConfig = {},
                   Real temperature_K = constants::T0,
                   DDScalingSpec scaling = {},
                   BandgapNarrowingConfig bandgapNarrowingConfig = {});

    ContactCurrentResult compute(const DDSolution& solution,
                                 const std::string& contactName) const;

    static ContactCurrentResult compute(const DeviceMesh& mesh,
                                        const MaterialDatabase& matdb,
                                        const DopingModel& doping,
                                        const DDSolution& solution,
                                        const std::string& contactName,
                                        const MobilityModelConfig& mobilityConfig = {},
                                        Real temperature_K = constants::T0,
                                        DDScalingSpec scaling = {},
                                        const BandgapNarrowingConfig& bandgapNarrowingConfig = {});

private:
    const DeviceMesh& mesh_;
    const MaterialDatabase& matdb_;
    const DopingModel& doping_;
    std::vector<std::vector<Index>> edgeCells_;
    MobilityModelConfig mobilityConfig_;
    std::unique_ptr<MobilityModel> mobility_;
    Real thermalVoltage_;
    DDScalingSpec scaling_;
    std::vector<Real> ni_;
};

} // namespace vela
