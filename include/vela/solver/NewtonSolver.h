#pragma once

#include "vela/core/PhysicalConstants.h"
#include "vela/core/Types.h"
#include "vela/core/UnitScaling.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/numerics/LineSearch.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/solver/GummelSolver.h"
#include <nlohmann/json_fwd.hpp>
#include <string>
#include <unordered_map>
#include <vector>

namespace vela {

struct NewtonConfig {
    int maxIter = 20;
    Real reltol = 1.0e-8;
    Real abstol = 1.0e-18;
    Real temperature_K = constants::T0; ///< Lattice temperature [K]
    Real dampingFactor = 1.0;
    bool lineSearch = true;
    bool verbose = true;
    bool warmStart = false; ///< Preserve supplied quasi-Fermi potentials instead of resetting interiors.
    bool diagnostics = false; ///< Store detailed line-search diagnostics in NewtonResult::history.
    Real finiteDifferenceStep = 1.0e-6;
    std::string jacobian = "analytic"; ///< "analytic" or "finite_difference"
    std::string residualNorm = "block"; ///< "block" or "l2" convergence/line-search norm
    Real residualWeightPsi = 1.0;
    Real residualWeightPhin = 1.0;
    Real residualWeightPhip = 1.0;
    Real residualScalePsi = 0.0;  ///< <= 0 selects max(initial psi-block residual norm, 1)
    Real residualScalePhin = 0.0; ///< <= 0 selects max(initial electron-continuity residual norm, 1)
    Real residualScalePhip = 0.0; ///< <= 0 selects max(initial hole-continuity residual norm, 1)
    Real taun = 1.0e-7;
    Real taup = 1.0e-7;
    MobilityModelConfig mobility{}; ///< Mobility model configuration
    std::vector<std::string> recombination = {"srh"}; ///< e.g. {"srh", "auger"}
    ImpactIonizationModelConfig impactIonization; ///< Avalanche generation model.
    BandgapNarrowingConfig bandgapNarrowing; ///< Effective ni model for high doping.
};

struct NewtonIterationInfo {
    int iter = 0;
    Real residualNorm = 0.0;
    Real stepNorm = 0.0;
    Real dampingFactor = 0.0;
    Real relativeResidualNorm = 0.0;
    Real rawStepNorm = 0.0;
    int lineSearchAttempts = 0;
    bool lineSearchAccepted = false;
    std::vector<LineSearchIterationInfo> lineSearchHistory;
};

struct NewtonResult {
    DDSolution solution;
    bool converged = false;
    int iters = 0;
    Real initialResidualNorm = 0.0;
    Real finalResidualNorm = 0.0;
    std::vector<NewtonIterationInfo> history;
};

class NewtonSolver {
public:
    NewtonSolver(const DeviceMesh& mesh,
                 const MaterialDatabase& matdb,
                 const DopingModel& doping,
                 const std::unordered_map<std::string, Real>& contactBiases,
                 NewtonConfig cfg = {},
                 std::vector<RegionFixedChargeSpec> fixedCharges = {},
                 std::vector<InterfaceSheetChargeSpec> sheetCharges = {});

    NewtonResult solve() const;
    NewtonResult solve(const DDSolution& initial) const;

private:
    CoupledDDBoundaryConditions buildBoundaryConditions(
        const CoupledDDAssembler& assembler) const;
    DDSolution buildInitialGuess(const CoupledDDAssembler& assembler,
                                 const CoupledDDBoundaryConditions& bcs) const;
    DDSolution makeSolution(const CoupledDDAssembler& assembler,
                            const VectorXd& x,
                            int iters) const;

    const DeviceMesh& mesh_;
    const MaterialDatabase& matdb_;
    const DopingModel& doping_;
    std::unordered_map<std::string, Real> contactBiases_;
    NewtonConfig cfg_;
    std::vector<RegionFixedChargeSpec> fixedCharges_;
    std::vector<InterfaceSheetChargeSpec> sheetCharges_;
};

NewtonConfig newtonConfigFromJson(
    const nlohmann::json& cfg,
    UnitScalingConfig scaling = {});

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const NewtonConfig& cfg = {});

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const DDSolution& initial,
                       const NewtonConfig& cfg = {});


NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const NewtonConfig& cfg,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges);

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const DDSolution& initial,
                       const NewtonConfig& cfg,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges);

} // namespace vela
