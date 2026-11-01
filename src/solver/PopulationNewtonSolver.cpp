#include "vela/solver/NewtonSolver.h"
#include "vela/simulation/CoupledLoadLine.h"
#include "vela/numerics/ResidualNorm.h"
#include <algorithm>
#include <stdexcept>

namespace vela {
NewtonResult NewtonSolver::solvePopulation(const DDSolution& initial) const {
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
    c.equationTolerance=cfg_.blockAbsoluteConvergence.mode=="enforce"?1.:
        std::max(cfg_.abstol,cfg_.reltol*first);
    c.loadLineTolerance_V=1e-12;
    c.dampingFactor=cfg_.dampingFactor;
    c.lineSearchMode=cfg_.lineSearch?"merit":"full_step";
    c.linearSolver="direct_bordered";
    c.directRefinementSteps=cfg_.linearRefinementIterations;
    c.updateConvergence.mode="report";
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
    const auto solved=CoupledLoadLineNewton(system,zero.value,c,zero,
        makeArclengthUpdateErrorEvaluator(contact,c.updateConvergence)).solve({y,bias->second},bias->second,1.);
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
} // namespace vela
