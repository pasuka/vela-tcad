#include "vela/solver/NewtonSolver.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/numerics/ResidualNorm.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/solver/LinearSolver.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace vela {
namespace {

inline double thermalVoltage(double T = 300.0)
{
    return constants::kb * T / constants::q;
}

inline double nEq(double Ndop, double ni)
{
    const double half = 0.5 * Ndop;
    return half + std::sqrt(half * half + ni * ni);
}

} // namespace


NewtonConfig newtonConfigFromJson(const nlohmann::json& json)
{
    NewtonConfig cfg;
    cfg.maxIter = json.value("max_iter", cfg.maxIter);
    cfg.reltol = json.value("reltol", cfg.reltol);
    cfg.abstol = json.value("abstol", cfg.abstol);
    cfg.dampingFactor = json.value("damping_factor", cfg.dampingFactor);
    cfg.lineSearch = json.value("line_search", cfg.lineSearch);
    cfg.verbose = json.value("verbose", cfg.verbose);
    cfg.finiteDifferenceStep = json.value("finite_difference_step", cfg.finiteDifferenceStep);
    cfg.taun = json.value("taun", cfg.taun);
    cfg.taup = json.value("taup", cfg.taup);
    cfg.mobility = json.value("mobility", cfg.mobility);

    if (json.contains("recombination")) {
        const auto& value = json.at("recombination");
        if (value.is_array())
            cfg.recombination = value.get<std::vector<std::string>>();
        else if (value.is_string())
            cfg.recombination = {value.get<std::string>()};
        else
            throw std::invalid_argument(
                "newtonConfigFromJson: recombination must be a string or string array.");
    }

    return cfg;
}

NewtonSolver::NewtonSolver(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const std::unordered_map<std::string, Real>& contactBiases,
    NewtonConfig cfg)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , contactBiases_(contactBiases)
    , cfg_(cfg)
{}

CoupledDDBoundaryConditions NewtonSolver::buildBoundaryConditions(
    const CoupledDDAssembler& assembler) const
{
    CoupledDDBoundaryConditions bcs;
    const auto& ni = assembler.intrinsicDensity();
    const double Vt = thermalVoltage();

    for (Index c = 0; c < mesh_.numContacts(); ++c) {
        const Contact& contact = mesh_.getContact(c);
        auto it = contactBiases_.find(contact.name);
        if (it == contactBiases_.end()) continue;

        const double Vbias = it->second;
        for (Index nid : contact.node_ids) {
            const double niNode = ni[nid];
            const double neq = nEq(doping_.netDoping(nid), niNode);
            double psiBuiltIn = 0.0;
            if (niNode > 0.0 && neq > 0.0)
                psiBuiltIn = Vt * std::log(neq / niNode);

            bcs.psi[nid] = Vbias + psiBuiltIn;
            bcs.phin[nid] = Vbias;
            bcs.phip[nid] = Vbias;
        }
    }
    return bcs;
}

DDSolution NewtonSolver::buildInitialGuess(
    const CoupledDDAssembler&, const CoupledDDBoundaryConditions&) const
{
    GummelConfig gcfg;
    gcfg.maxIter = 1;
    gcfg.reltol = 0.0;
    gcfg.dampingPsi = 0.5;
    gcfg.taun = cfg_.taun;
    gcfg.taup = cfg_.taup;
    gcfg.mobility = cfg_.mobility;
    gcfg.recombination = cfg_.recombination;
    DDSolution sol = runGummel(mesh_, matdb_, doping_, contactBiases_, gcfg);

    // The Gummel solver leaves tiny numerical noise (~1e-18 V) in the
    // quasi-Fermi potentials at interior nodes.  The balanced SG flux formula
    // multiplies exp(-phi/Vt) differences by coefficients of order exp(psi/Vt)
    // (~1e5 at a PN junction), so even 1-ULP rounding in exp(-tiny/Vt) creates
    // a spurious initial residual of O(0.05) that prevents Newton convergence.
    // Zeroing phi_n/phi_p is physically correct for equilibrium (flat
    // quasi-Fermi levels) and produces a small, well-conditioned initial
    // residual for the coupled Newton iteration.
    const int N = static_cast<int>(mesh_.numNodes());
    sol.phin = VectorXd::Zero(N);
    sol.phip = VectorXd::Zero(N);
    return sol;
}

DDSolution NewtonSolver::makeSolution(const CoupledDDAssembler& assembler,
                                      const VectorXd& x,
                                      int iters) const
{
    CoupledDDState state = assembler.unpack(x);
    DDSolution sol;
    sol.psi = state.psi;
    sol.phin = state.phin;
    sol.phip = state.phip;
    sol.n = assembler.electronDensity(x);
    sol.p = assembler.holeDensity(x);
    sol.iters = iters;
    return sol;
}

NewtonResult NewtonSolver::solve() const
{
    const double Vt = thermalVoltage();
    MobilityModelConfig mobilityConfig = mobilityModelConfig(cfg_.mobility);
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    CoupledDDAssembler assembler(mesh_, matdb_, doping_, Vt, mobilityConfig, recombinationConfig);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);
    return solve(buildInitialGuess(assembler, bcs));
}

NewtonResult NewtonSolver::solve(const DDSolution& initial) const
{
    const double Vt = thermalVoltage();
    MobilityModelConfig mobilityConfig = mobilityModelConfig(cfg_.mobility);
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    CoupledDDAssembler assembler(mesh_, matdb_, doping_, Vt, mobilityConfig, recombinationConfig);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);

    // The balanced Scharfetter-Gummel formula multiplies the quasi-Fermi
    // difference (expNegPhin[i]-expNegPhin[j]) by exp(+ψ[j]/Vt) for electrons
    // and exp(-ψ[i]/Vt) for holes.  At a PN junction these factors reach ~1e5,
    // so even sub-ULP noise in phin/phip from an external initial guess (e.g.
    // Gummel) produces O(1) residuals.  Zeroing the interior quasi-Fermi
    // potentials removes this amplification: at equilibrium the exact solution
    // has phin = phip = 0 everywhere, and for any bias Newton will converge to
    // the correct non-zero values from this well-conditioned start.
    VectorXd phinInit = initial.phin;
    VectorXd phipInit = initial.phip;
    const int N = static_cast<int>(mesh_.numNodes());
    for (int i = 0; i < N; ++i) {
        const Index nid = static_cast<Index>(i);
        if (bcs.phin.find(nid) == bcs.phin.end()) {
            phinInit(i) = 0.0;
            phipInit(i) = 0.0;
        }
    }

    VectorXd x = assembler.pack({initial.psi, phinInit, phipInit});
    VectorXd r = assembler.residual(x, bcs);
    const Real initialNorm = r.norm();

    NewtonResult result;
    result.solution = initial;
    result.initialResidualNorm = initialNorm;
    result.finalResidualNorm = initialNorm;

    if (cfg_.verbose) {
        std::cout << "Newton iter 0 residual=" << initialNorm
                  << " step=0 damping=0\n";
    }

    if (initialNorm <= cfg_.abstol) {
        result.converged = true;
        result.solution = makeSolution(assembler, x, 0);
        return result;
    }

    LinearSolver linearSolver;
    LineSearchConfig lscfg;
    lscfg.enabled = cfg_.lineSearch;
    lscfg.initialDamping = cfg_.dampingFactor;
    BacktrackingLineSearch lineSearch(lscfg);

    VectorXd acceptedX = x;
    VectorXd acceptedR = r;
    int acceptedIters = 0;

    for (int iter = 1; iter <= cfg_.maxIter; ++iter) {
        SparseMatrixd J = assembler.finiteDifferenceJacobian(
            x, bcs, cfg_.finiteDifferenceStep);
        VectorXd step;
        try {
            step = linearSolver.solve(J, -r);
        } catch (const std::runtime_error&) {
            result.finalResidualNorm = acceptedR.norm();
            result.iters = acceptedIters;
            result.solution = makeSolution(assembler, acceptedX, acceptedIters);
            return result;
        }
        const Real stepNorm = step.norm();

        const auto ls = lineSearch.search(
            x, step, r,
            [&](const VectorXd& candidate) { return assembler.residual(candidate, bcs); },
            [&](const VectorXd& candidate, const VectorXd&) {
                return assembler.hasPositiveFiniteCarriers(candidate);
            });

        if (!ls.accepted) {
            result.finalResidualNorm = acceptedR.norm();
            result.iters = acceptedIters;
            result.solution = makeSolution(assembler, acceptedX, acceptedIters);
            return result;
        }

        x = ls.x;
        r = ls.residual;
        acceptedX = x;
        acceptedR = r;
        acceptedIters = iter;

        // Record the norm of the actually applied update (damped step) so that
        // per-iteration metrics are consistent with the accepted solution.
        const Real appliedStepNorm = ls.damping * stepNorm;
        const Real residualNorm = ls.residualNorm;
        result.history.push_back({iter, residualNorm, appliedStepNorm, ls.damping});
        if (cfg_.verbose) {
            std::cout << "Newton iter " << iter
                      << " residual=" << residualNorm
                      << " step=" << appliedStepNorm
                      << " damping=" << ls.damping << '\n';
        }

        const Real rel = ResidualNorm::relative(residualNorm, initialNorm);
        if (residualNorm <= cfg_.abstol || rel <= cfg_.reltol) {
            result.converged = true;
            result.iters = iter;
            result.finalResidualNorm = residualNorm;
            result.solution = makeSolution(assembler, x, iter);
            return result;
        }
    }

    result.converged = false;
    result.iters = acceptedIters;
    result.finalResidualNorm = acceptedR.norm();
    result.solution = makeSolution(assembler, acceptedX, acceptedIters);
    return result;
}

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const NewtonConfig& cfg)
{
    return NewtonSolver(mesh, matdb, doping, contactBiases, cfg).solve();
}

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const DDSolution& initial,
                       const NewtonConfig& cfg)
{
    return NewtonSolver(mesh, matdb, doping, contactBiases, cfg).solve(initial);
}

} // namespace vela
