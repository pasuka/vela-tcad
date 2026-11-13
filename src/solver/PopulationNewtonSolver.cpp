#include "vela/solver/NewtonSolver.h"
#include "vela/solver/LinearSolver.h"
#include "vela/numerics/ExpansionProductSum.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/simulation/PopulationCoordinateNewton.h"
#include "vela/numerics/ResidualNorm.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace vela {
NewtonResult NewtonSolver::solvePopulation(const DDSolution& initial) const {
    if (cfg_.compensatedPopulationState)
        return solveCompensatedPopulation(initial);
    if(mesh_.contacts().empty())throw std::invalid_argument("population Newton requires an explicit contact");
    const std::string contact=mesh_.contacts().front().name;
    const auto bias=contactBiases_.find(contact);
    if(bias==contactBiases_.end())throw std::invalid_argument("population Newton requires the reference contact bias");
    auto system=makeArclengthSystem(contact);
    const auto y=packArclengthState(initial);
    const auto initialF=system.residual(y,bias->second);
    const auto blocks=ResidualNorm::computeBlocks(initialF,mesh_.numNodes());
    const ResidualBlockNormValue scales{
        cfg_.residualScalePsi>0.?cfg_.residualScalePsi:std::max(blocks.psi,1.),
        cfg_.residualScalePhin>0.?cfg_.residualScalePhin:std::max(blocks.phin,1.),
        cfg_.residualScalePhip>0.?cfg_.residualScalePhip:std::max(blocks.phip,1.),0.};
    const auto weights=ResidualBlockWeights{cfg_.residualWeightPsi,cfg_.residualWeightPhin,cfg_.residualWeightPhip};
    system.equationNorm=[this,scales,weights](const VectorXd& f) {
        return cfg_.residualNorm=="l2"?f.norm():ResidualNorm::normalizedBlockL2(
            ResidualNorm::computeBlocks(f,mesh_.numNodes()),scales,weights);
    };
    CoupledLoadLineConfig c;
    c.maxIterations=cfg_.maxIter;
    // Enforced block ceilings replace scalar relative/absolute stopping, as in
    // ordinary Newton. Do not add an unrelated initial-state-dependent gate.
    if(cfg_.blockAbsoluteConvergence.mode=="enforce") {
        const auto block=cfg_.blockAbsoluteConvergence;
        const int n=mesh_.numNodes();
        system.equationNorm=[block,n](const VectorXd& f) {
            if(!f.allFinite())return std::numeric_limits<Real>::infinity();
            return std::max({f.head(n).cwiseAbs().maxCoeff()/block.psiResidualCeiling,
                f.segment(n,n).cwiseAbs().maxCoeff()/block.electronResidualCeiling,
                f.tail(n).cwiseAbs().maxCoeff()/block.holeResidualCeiling});
        };
    }
    const double first=system.equationNorm(initialF);
    const auto acceptsLegacyFloor = [&](const ArclengthState& state, Real norm) {
        if (cfg_.coupledOperator != "legacy" || cfg_.populationStoppingRule != "legacy" ||
            cfg_.blockAbsoluteConvergence.mode == "enforce" ||
            !std::isfinite(norm) || norm > cfg_.stallResidualFloor ||
            std::abs(state.lambda - bias->second) > 1e-12)
            return false;
        if (system.convergenceSatisfied && !system.convergenceSatisfied(state.x, state.lambda))
            return false;
        const Real limit = cfg_.poissonLineSearchStallContactMajorityQfDropLimit_V;
        if (limit <= 0.0)
            return true;
        const Real drop = maxContactMajorityQuasiFermiDrop(unpackArclengthState(state.x));
        return std::isfinite(drop) && drop <= limit;
    };
    c.equationTolerance=cfg_.blockAbsoluteConvergence.mode=="enforce"?1.:
        std::max(cfg_.abstol,cfg_.reltol*first);
    c.loadLineTolerance_V=1e-12;
    c.dampingFactor=cfg_.dampingFactor;
    c.lineSearchMode=cfg_.lineSearch?"merit":"full_step";
    c.linearSolver="direct_bordered";
    c.directRefinementSteps=cfg_.linearRefinementIterations;
    c.updateConvergence = cfg_.populationUpdateConvergence;
    c.updateConvergence.mode="report";
    c.minimumOriginalCorrections=cfg_.primaryVariableCoordinates=="unit"?1:0;
    // The scalar fixes the selected contact voltage exactly. The zero
    // functional makes this the same device solve, without a physical resistor.
    ArclengthScalarFunctional zero;
    zero.value=[](const VectorXd&,Real){return 0.;};
    zero.linearize=[](const VectorXd& x,Real){
        ArclengthScalarLinearization out;out.stateDerivative=VectorXd::Zero(x.size());return out;
    };
    zero.linearizeInColumns=[](const VectorXd& x,Real,const VectorXd&){
        ArclengthScalarLinearization out;out.stateDerivative=VectorXd::Zero(x.size());return out;
    };
    CoupledLoadLineResult solved;
    const ArclengthState seed{y, bias->second};
    if (acceptsLegacyFloor(seed, first)) {
        solved.state = seed;
        solved.equationResidualNorm = first;
        solved.converged = true;
        solved.convergenceReason = "initial_stall_residual_floor";
    } else {
        solved = PopulationCoordinateNewton(system,zero.value,c,zero,
            makeArclengthUpdateErrorEvaluator(contact,c.updateConvergence,
                cfg_.coupledOperator != "legacy" || cfg_.populationStoppingRule == "legacy"),
            cfg_.primaryVariableCoordinates=="unit"?primaryVariableScales():PrimaryVariableScales{},
            cfg_.primaryVariableCoordinates=="unit",cfg_.populationStoppingRule=="native_approx",
            approximateNativeRhsScales()).solve(seed,bias->second,1.);
        if (!solved.converged &&
            (solved.failureReason.starts_with("coupled line search rejected update") ||
             solved.failureReason == "coupled Newton iteration budget exhausted") &&
            acceptsLegacyFloor(solved.state, solved.equationResidualNorm)) {
            solved.converged = true;
            solved.convergenceReason = "stall_residual_floor";
            solved.failureReason.clear();
        }
    }
    NewtonResult out;
    out.solution=unpackArclengthState(solved.state.x);
    out.converged=solved.converged;out.iters=solved.iterations;
    out.solution.converged=out.converged;out.solution.iters=out.iters;
    out.initialResidualNorm=first;out.finalResidualNorm=solved.equationResidualNorm;
    out.convergenceReason=solved.convergenceReason;
    out.failureDiagnostics.failureReason=solved.failureReason;
    const auto gates=makeArclengthConvergenceEvaluator(contact)(solved.state.x,solved.state.lambda);
    out.finalBlockNorms=gates.blockNorms;
    out.finalCarrierRowConvergence=gates.carrierRows;
    out.finalGlobalContinuityClosure=gates.globalContinuity;
    for(const auto& r:solved.trace) {
        NewtonIterationInfo info;info.iter=r.iteration;
        info.residualNorm=r.equationResidualNorm;info.rawStepNorm=r.rawStateUpdateNorm;
        info.dampingFactor=r.acceptedAlpha;
        info.lineSearchAttempts=solved.lineSearchAttempts;
        info.event="potential_density";out.trace.push_back(info);
        out.trace.back().populationUpdateError=r.updateError;
        out.trace.back().populationEffectiveAlpha=r.populationEffectiveAlpha;
        out.trace.back().protectedElectrons=r.protectedElectrons;
        out.trace.back().protectedHoles=r.protectedHoles;
        if(cfg_.diagnostics)out.history.push_back(out.trace.back());
    }
    return out;
}
NewtonResult NewtonSolver::solveCompensatedPopulation(const DDSolution& initial) const
{
    auto assembler = makeArclengthAssembler();
    const auto boundary = buildBoundaryConditions(*assembler);
    const auto scaling = buildScalingSpec();
    const int count = mesh_.numNodes();
    if (initial.psi.size() != count || initial.n.size() != count || initial.p.size() != count ||
        !initial.psi.allFinite() || !initial.n.allFinite() || !initial.p.allFinite() ||
        initial.n.minCoeff() <= 0. || initial.p.minCoeff() <= 0. ||
        (initial.populationLow.size() && initial.populationLow.size() != 3*count))
        throw std::invalid_argument("Invalid compensated population initial state");

    VectorXd state(3*count);
    state << initial.psi, initial.n, initial.p;
    VectorXd low = initial.populationLow.size() ? initial.populationLow : VectorXd::Zero(3*count);
    std::vector<Index> contacts;
    std::vector<int> prescribed;
    for (const auto& [node,target] : boundary.ohmicPopulation) {
        experimental::ExpansionProductSum potential;
        potential.add(boundary.ohmicBias.at(node));
        potential.add(target[0]);
        const auto value = potential.parts();
        state[node] = value.high;
        low[node] = value.low;
        state[count+node] = target[1];
        state[2*count+node] = target[2];
        low[count+node] = low[2*count+node] = 0.;
        contacts.push_back(node);
        for (int block = 0; block < 3; ++block)
            prescribed.push_back(block*count+node);
    }
    if (contacts.empty())
        throw std::invalid_argument("Compensated population Newton requires Ohmic contacts");

    // The PN2D precision validation uses absolute physical residual ceilings.
    // A small relative update alone cannot resolve the reverse leakage current.
    constexpr Real poissonCeiling_C_per_m = 8e-28;
    constexpr Real continuityCeiling_A_per_um = 5e-22;
    const Real poissonWeight = scaling.permittivityReference_F_per_m * scaling.V0 /
        poissonCeiling_C_per_m;
    const Real continuityScale = constants::q * assembler->continuityResidualScale() *
        scaling.currentDensityLineIntegralFactor * 1e-6;
    // Edge currents and source terms are binary64. Allow for their rounding
    // when a large forward flux enters and leaves the same control volume.
    // The absolute leakage floor still governs rows with small currents.
    constexpr Real continuityRoundoff = 16. * std::numeric_limits<Real>::epsilon();
    const auto merit = [&](const CoupledDDPhysicalDiagnostics& physical) {
        const VectorXd& residual = physical.residual;
        if (!residual.allFinite()) return std::numeric_limits<Real>::infinity();
        Real error = poissonWeight*residual.head(count).lpNorm<Eigen::Infinity>();
        for (int node = 0; node < count; ++node) {
            const auto& row = physical.constrainedTerms[node];
            const Real electronScale = row.electronFluxAbsSum +
                std::abs(row.electronRecombination) + std::abs(row.electronImpact);
            const Real holeScale = row.holeFluxAbsSum +
                std::abs(row.holeRecombination) + std::abs(row.holeImpact);
            error = std::max({error,
                continuityScale*std::abs(residual[count+node]) /
                    std::max(continuityCeiling_A_per_um,continuityRoundoff*continuityScale*electronScale),
                continuityScale*std::abs(residual[2*count+node]) /
                    std::max(continuityCeiling_A_per_um,continuityRoundoff*continuityScale*holeScale)});
        }
        return error;
    };
    const auto evaluate = [&](const VectorXd& x, const VectorXd& tails, SparseMatrixd* jacobian = nullptr) {
        return assembler->compensatedPopulationDiagnostics(
            x.head(count),x.segment(count,count),x.tail(count),tails,boundary,true,jacobian);
    };
    const auto blockInfo = [&](const VectorXd& residual) {
        const auto blocks = ResidualNorm::computeBlocks(residual,count);
        return NewtonBlockResidualInfo{blocks.psi,blocks.phin,blocks.phip,blocks.combined};
    };
    const auto gates = [&](const CoupledDDPhysicalDiagnostics& physical, NewtonResult& result) {
        result.finalBlockNorms = blockInfo(physical.residual);
        result.finalCarrierRowConvergence = evaluateCarrierRowConvergence(
            physical.constrainedTerms,cfg_.carrierRowConvergence);
        result.finalGlobalContinuityClosure = evaluateGlobalContinuityClosure(
            physical.unconstrainedTerms,contacts,contacts,cfg_.globalContinuityClosure);
        const auto& block = cfg_.blockAbsoluteConvergence;
        return (!result.finalCarrierRowConvergence.enforced || result.finalCarrierRowConvergence.satisfied) &&
               (!result.finalGlobalContinuityClosure.enforced || result.finalGlobalContinuityClosure.satisfied) &&
               (block.mode != "enforce" ||
                (physical.residual.head(count).lpNorm<Eigen::Infinity>() <= block.psiResidualCeiling &&
                 physical.residual.segment(count,count).lpNorm<Eigen::Infinity>() <= block.electronResidualCeiling &&
                 physical.residual.tail(count).lpNorm<Eigen::Infinity>() <= block.holeResidualCeiling));
    };

    NewtonResult result;
    LinearSolver linear;
    auto physical = evaluate(state,low);
    Real currentMerit = merit(physical);
    result.initialResidualNorm = currentMerit;
    result.failureDiagnostics.failureReason = "compensated population Newton iteration budget exhausted";
    for (int iteration = 0; iteration <= cfg_.maxIter; ++iteration) {
        if (currentMerit <= 1. && gates(physical,result)) {
            result.converged = true;
            result.convergenceReason = "compensated_population_residual";
            result.failureDiagnostics.failureReason.clear();
            break;
        }
        if (iteration == cfg_.maxIter)
            break;

        // Keep mobility, SRH and avalanche live in the residual. Freeze their
        // derivatives in this correction, which uses the scaled SG/Poisson stencil.
        SparseMatrixd jacobian;
        physical = evaluate(state,low,&jacobian);
        VectorXd rowMax = VectorXd::Zero(state.size());
        for (int column = 0; column < jacobian.outerSize(); ++column)
            for (SparseMatrixd::InnerIterator it(jacobian,column); it; ++it)
                rowMax[it.row()] = std::max(rowMax[it.row()],std::abs(it.value()));
        if ((rowMax.array() <= 0.).any() || !rowMax.allFinite()) {
            result.failureDiagnostics.failureReason = "invalid compensated population Jacobian row";
            break;
        }
        for (int column = 0; column < jacobian.outerSize(); ++column)
            for (SparseMatrixd::InnerIterator it(jacobian,column); it; ++it)
                it.valueRef() /= rowMax[it.row()];
        VectorXd rhs = -physical.residual.cwiseQuotient(rowMax);
        for (int row : prescribed) rhs[row] = 0.;
        VectorXd delta;
        try {
            delta = linear.solve(jacobian,rhs,prescribed);
        } catch (const std::runtime_error& error) {
            result.failureDiagnostics.failureReason = error.what();
            break;
        }
        delta.head(count) *= scaling.V0;
        delta.tail(2*count).array() *= state.tail(2*count).array();
        if (!delta.allFinite()) {
            result.failureDiagnostics.failureReason = "nonfinite compensated population correction";
            break;
        }

        bool accepted = false;
        Real alpha = cfg_.dampingFactor;
        for (int attempt = 1; attempt <= 20; ++attempt, alpha *= 0.5) {
            VectorXd trial(state.size()), trialLow(state.size());
            for (Eigen::Index i = 0; i < state.size(); ++i) {
                experimental::ExpansionProductSum sum;
                sum.add(state[i]);
                sum.add(low[i]);
                sum.add(alpha*delta[i]);
                const auto value = sum.parts();
                trial[i] = value.high;
                trialLow[i] = value.low;
            }
            if (trial.tail(2*count).minCoeff() <= 0.)
                continue;
            auto candidate = evaluate(trial,trialLow);
            const Real candidateMerit = merit(candidate);
            if (candidateMerit >= currentMerit)
                continue;
            state.swap(trial);
            low.swap(trialLow);
            physical = std::move(candidate);
            currentMerit = candidateMerit;
            result.iters = iteration+1;
            NewtonIterationInfo info;
            info.iter = result.iters;
            info.event = "compensated_population";
            info.residualNorm = currentMerit;
            info.blockResiduals = blockInfo(physical.residual);
            info.dampingFactor = alpha;
            info.lineSearchAttempts = attempt;
            info.lineSearchAccepted = true;
            result.trace.push_back(info);
            if (cfg_.diagnostics)
                result.history.push_back(info);
            accepted = true;
            break;
        }
        if (!accepted) {
            result.failureDiagnostics.failureReason = "compensated population line search rejected correction";
            break;
        }
    }
    gates(physical,result);
    result.finalResidualNorm = currentMerit;
    result.solution = makeSolution(*assembler,assembler->packPopulation(
        state.head(count),state.segment(count,count),state.tail(count)),result.iters);
    result.solution.psi = state.head(count);
    result.solution.populationPsi.resize(0);
    result.solution.populationPotentialScale_V = 0.;
    result.solution.populationLow = std::move(low);
    result.solution.converged = result.converged;
    return result;
}
} // namespace vela
