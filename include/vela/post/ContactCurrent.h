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
#include <unordered_map>
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

struct ContactCurrentEdgeOverrides {
    std::unordered_map<Index, Real> holeQuasiFermiDropByEdge;
};

struct ContactCurrentEdgeDiagnostic {
    Index edgeId = -1;
    Index node0 = -1;
    Index node1 = -1;
    Real edgeLength_m = 0.0;
    Real edgeCouple_m = 0.0;
    Real outwardSign = 0.0;
    Real bernoulliU = 0.0;
    Real bernoulliBplus = 0.0;
    Real bernoulliBminus = 0.0;
    bool electronUsedQuasiFermi = false;
    bool holeUsedQuasiFermi = false;
    Real psi0 = 0.0;
    Real psi1 = 0.0;
    Real phin0 = 0.0;
    Real phin1 = 0.0;
    Real phip0 = 0.0;
    Real phip1 = 0.0;
    bool holeQfDropOverrideApplied = false;
    Real n0 = 0.0;
    Real n1 = 0.0;
    Real p0 = 0.0;
    Real p1 = 0.0;
    Real ni0 = 0.0;
    Real ni1 = 0.0;
    Real mun = 0.0;
    Real mup = 0.0;
    Real electronContinuityFlux = 0.0;
    Real holeContinuityFlux = 0.0;
    Real electronCurrent = 0.0;
    Real electronDriftCurrent = 0.0;
    Real electronDiffusionCurrent = 0.0;
    Real holeCurrent = 0.0;
    Real holeDriftCurrent = 0.0;
    Real holeDiffusionCurrent = 0.0;
    Real totalCurrent = 0.0;
};

struct ContactCurrentDetailedResult {
    ContactCurrentResult totals;
    std::vector<ContactCurrentEdgeDiagnostic> edges;
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
    ContactCurrentResult compute(const DDSolution& solution,
                                 const std::string& contactName,
                                 const ContactCurrentEdgeOverrides& overrides) const;

    ContactCurrentDetailedResult computeDetailed(const DDSolution& solution,
                                                 const std::string& contactName) const;
    ContactCurrentDetailedResult computeDetailed(const DDSolution& solution,
                                                 const std::string& contactName,
                                                 const ContactCurrentEdgeOverrides& overrides) const;

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
