#pragma once

#include "vela/core/Types.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/numerics/LineSearch.h"
#include "vela/physics/DopingModel.h"
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
    Real finiteDifferenceStep = 1.0e-6;
    std::string jacobian = "analytic"; ///< "analytic" or "finite_difference"
    Real taun = 1.0e-7;
    Real taup = 1.0e-7;
    std::string mobility = "constant"; ///< "constant" or "caughey_thomas"
    std::vector<std::string> recombination = {"srh"}; ///< e.g. {"srh", "auger"}
};

struct NewtonIterationInfo {
    int iter = 0;
    Real residualNorm = 0.0;
    Real stepNorm = 0.0;
    Real dampingFactor = 0.0;
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
                 NewtonConfig cfg = {});

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
};

NewtonConfig newtonConfigFromJson(const nlohmann::json& cfg);

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

} // namespace vela
