#include "vela/solver/NewtonSolver.h"
#include "vela/solver/LinearSolver.h"
#include "vela/numerics/ExpansionProductSum.h"
#include "vela/numerics/SplitCoordinate.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/core/RuntimeLog.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/simulation/PopulationCoordinateNewton.h"
#include "vela/numerics/ResidualNorm.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <set>

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
    c.reuseBorderedAnalysis = cfg_.reuseFixedVoltageBorderedAnalysis;
    c.directRefinementSteps=cfg_.linearRefinementIterations;
    c.bvMethodsEftDefect = cfg_.coupledOperator == "bvmethods_tri3_qf" &&
        linear_refinement::useBvMethodsEftDefect();
    c.eftLinearDefectWithoutMp100=cfg_.populationLinearDefect=="eft";
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
    for (const auto& [node, target] : boundary.psi) {
        if (boundary.ohmicPopulation.contains(node))
            continue;
        state[node] = target * scaling.V0;
        low[node] = 0.;
        prescribed.push_back(node);
    }
    const auto& intrinsic = assembler->intrinsicDensity();
    const auto positiveCarriers = [&](const VectorXd& values, const VectorXd& tails) {
        for (int node = 0; node < count; ++node) {
            if (intrinsic[node] > 0. &&
                (!(values[count + node] + tails[count + node] > 0.) ||
                 !(values[2 * count + node] + tails[2 * count + node] > 0.)))
                return false;
        }
        return true;
    };
    for (int node = 0; node < count; ++node) {
        if (intrinsic[node] > 0.)
            continue;
        for (int block = 1; block < 3; ++block) {
            state[block * count + node] = low[block * count + node] = 0.;
            prescribed.push_back(block * count + node);
        }
    }
    if (!positiveCarriers(state, low))
        throw std::invalid_argument("Invalid compensated semiconductor populations");
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
            const Real electronScale = row.electronFluxAbsSum + std::abs(row.electronBoundary) +
                std::abs(row.electronRecombination) + std::abs(row.electronImpact);
            const Real holeScale = row.holeFluxAbsSum + std::abs(row.holeBoundary) +
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

        // Evaluate the current transport coefficients before forming the correction.
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
            if (!positiveCarriers(trial, trialLow))
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
NewtonResult NewtonSolver::solveMonolithicDensityGradientBox(
    const DDSolution& initial,
    std::vector<DensityGradientCellMaterial> cells,
    const VectorXd& nodeMaterialDrive_V,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletLambda_V) const
{
    auto assembler = makeArclengthAssembler();
    const bool population = assembler->usesPopulationState();
    if (!population) {
        assembler->setFermiNodeCache(cfg_.diagnosticFermiNodeCache);
        if (cfg_.quasiFermiRecenterOnInitialState)
            assembler->setQuasiFermiReferenceFields(initial.phin, initial.phip);
    }
    const auto boundary = buildBoundaryConditions(*assembler);
    if (!boundary.thermionic.empty())
        throw std::invalid_argument("The four-field DG box currently requires non-thermionic contacts.");
    const int nodes = static_cast<int>(mesh_.numNodes());
    const int ddSize = 3 * nodes;
    const Real vt = constants::kb * cfg_.temperature_K / constants::q;
    const Real potentialScale =
        assembler->usesScaledState() ? assembler->potentialScale() : 1.;
    const auto unitSystem = cfg_.inputScaling.unitSystem();
    const auto& dg = cfg_.electronQuantumPotential;
    const auto statistics = carrierStatisticsModel(cfg_.carrierStatistics);
    const auto& Nc = assembler->electronDensityOfStates();
    const auto& Nv = assembler->holeDensityOfStates();
    const auto& ni = assembler->intrinsicDensity();
    std::vector<bool> electronFree(nodes, false);
    std::vector<bool> holeFree(nodes, false);
    for (int node = 0; node < nodes; ++node) {
        electronFree[node] = ni[node] > 0. && !boundary.phin.contains(node);
        holeFree[node] = ni[node] > 0. && !boundary.phip.contains(node);
    }
    const auto pinCarriers = [&](VectorXd& state) {
        for (int node = 0; node < nodes; ++node) {
            if (ni[node] <= 0.) {
                state[nodes + node] = 0.;
                state[2 * nodes + node] = 0.;
            }
        }
    };

    VectorXd maxDrive = nodeMaterialDrive_V;
    std::vector<bool> seen(nodes, false);
    for (const auto& cell : cells) {
        const auto& ids = mesh_.getCell(cell.cellId).node_ids;
        for (int i = 0; i < 3; ++i) {
            const Index node = ids[i];
            if (!seen[node] || cell.materialBandDrive_V[i] > maxDrive[node])
                maxDrive[node] = cell.materialBandDrive_V[i];
            seen[node] = true;
        }
    }
    // Edge SG, its cell-sized mobility stencils, and local SRH only touch
    // the original Tri3 one-ring. Columns sharing any affected row are
    // perturbed separately for the transport/quantum mixed block.
    const bool localDgTangent = assembler->supportsLocalDgTransportTangent();
    std::vector<std::vector<int>> support(nodes);
    for (int i = 0; i < nodes; ++i)
        support[i].push_back(i);
    for (Index i = 0; i < mesh_.numEdges(); ++i) {
        const auto& edge = mesh_.getEdge(i);
        support[edge.n0].push_back(edge.n1);
        support[edge.n1].push_back(edge.n0);
    }
    std::vector<std::set<int>> usedAtRow(nodes);
    std::vector<std::vector<int>> groups;
    for (int column = 0; !localDgTangent && column < nodes; ++column) {
        if (ni[column] <= 0.)
            continue;
        std::set<int> forbidden;
        for (int row : support[column])
            forbidden.insert(usedAtRow[row].begin(), usedAtRow[row].end());
        int color = 0;
        while (forbidden.contains(color))
            ++color;
        if (color == static_cast<int>(groups.size()))
            groups.emplace_back();
        groups[color].push_back(column);
        for (int row : support[column])
            usedAtRow[row].insert(color);
    }
    const auto materialBase = cells;
    const auto packLambda = [&](const VectorXd& state) -> VectorXd {
        return (state.head(nodes) + state.tail(nodes)) * potentialScale + maxDrive;
    };
    struct System {
        VectorXd residual;
        SparseMatrixd jacobian;
    };
    const auto evaluate = [&](const VectorXd& state, bool assembleJacobian,
                              bool includeTransportQuantumDerivative = true) {
        const VectorXd dd = state.head(ddSize);
        const VectorXd lambda = packLambda(state);
        assembler->setElectronQuantumPotential(lambda);
        const auto fields = assembler->unpack(dd);
        const VectorXd psi_V = fields.psi * potentialScale;
        const VectorXd phin_V = fields.phin * potentialScale;
        for (std::size_t c = 0; c < cells.size(); ++c) {
            auto& data = cells[c];
            const auto& ids = mesh_.getCell(data.cellId).node_ids;
            for (int i = 0; i < 3; ++i) {
                const int node = ids[i];
                data.materialBandDrive_V[i] = psi_V[node] + materialBase[c].materialBandDrive_V[i];
                data.dynamicDrivingPotential_V[i] = data.isTransport ? -phin_V[node] : -psi_V[node];
            }
        }
        const auto box = linearizeElectronDensityGradientBox(
            mesh_, cells, nodeMaterialDrive_V + psi_V, activeNodes, dirichletLambda_V,
            vt, unitSystem, dg, state.tail(nodes) * potentialScale, assembleJacobian);
        System result;
        result.residual.resize(4 * nodes);
        result.residual << assembler->residual(dd, boundary), box.residual;
        if (!assembleJacobian)
            return result;

        std::vector<Eigen::Triplet<Real>> entries;
        SparseMatrixd quantumTransport;
        const auto ddJacobian = localDgTangent && includeTransportQuantumDerivative
            ? assembler->assembleJacobianWithDgTransportTangent(dd, boundary, quantumTransport)
            : assembler->assembleJacobian(dd, boundary);
        const auto append = [&](const SparseMatrixd& matrix, int rowOffset, int columnOffset,
                                const VectorXd& scales) {
            for (int c = 0; c < matrix.outerSize(); ++c)
                for (SparseMatrixd::InnerIterator entry(matrix, c); entry; ++entry)
                    entries.emplace_back(rowOffset + entry.row(), columnOffset + entry.col(),
                                         entry.value() * scales[entry.col()]);
        };
        append(ddJacobian, 0, 0, VectorXd::Ones(ddSize));
        if (!population) {
            // Direct band-edge coupling gives dn/dLambda = dn/dphin.
            // Assemble this dependence in the Poisson rows without
            // subtracting two nearly identical charge residuals.
            for (int column = nodes; column < 2 * nodes; ++column) {
                const int node = column - nodes;
                for (SparseMatrixd::InnerIterator entry(ddJacobian, column); entry; ++entry) {
                    if (entry.row() >= nodes)
                        break;
                    entries.emplace_back(entry.row(), node, entry.value());
                    entries.emplace_back(entry.row(), ddSize + node, entry.value());
                }
            }
        }
        const VectorXd volts = VectorXd::Constant(nodes, potentialScale);
        if (quantumTransport.cols()) {
            const VectorXd normalized = VectorXd::Ones(nodes);
            append(quantumTransport, 0, 0, normalized);
            append(quantumTransport, 0, ddSize, normalized);
        }
        append(box.potentialLike, ddSize, ddSize, volts);
        if (population)
            append(box.transportDrive, ddSize, ddSize, volts);
        append(box.insulatorDrive, ddSize, 0, -volts);
        VectorXd electronDrive = VectorXd::Zero(nodes);
        for (int node = 0; node < nodes; ++node) {
            entries.emplace_back(ddSize + node, node, box.materialDriveDiagonal[node] * potentialScale);
            if (population && ni[node] > 0.) {
                Real densitySlope = fields.n[node];
                if (usesFermiDirac(statistics)) {
                    const Real eta = inverseFermiDiracHalf(fields.n[node] / Nc[node], statistics);
                    densitySlope = Nc[node] * fermiDiracHalfDerivative(eta, statistics);
                }
                electronDrive[node] = vt / densitySlope;
            }
        }
        if (population)
            append(box.transportDrive, ddSize, nodes, electronDrive);
        else
            append(box.transportDrive, ddSize, nodes, -volts);

        constexpr Real perturbation_V = 1e-6;
        for (const auto& group : groups) {
            if (!includeTransportQuantumDerivative)
                break;
            VectorXd offset = VectorXd::Zero(nodes);
            for (int column : group)
                offset[column] = perturbation_V;
            assembler->setElectronQuantumPotential(lambda + offset);
            const VectorXd plus = assembler->residual(dd, boundary);
            assembler->setElectronQuantumPotential(lambda - offset);
            const VectorXd diff = (plus - assembler->residual(dd, boundary))
                * (potentialScale / (2. * perturbation_V));
            for (int column : group) {
                for (int node : support[column]) {
                    for (int block = population ? 0 : 1; block < 3; ++block) {
                        const int row = block * nodes + node;
                        if (diff[row] == 0.)
                            continue;
                        // At fixed Phi/q, Lambda also shifts with psi.
                        entries.emplace_back(row, column, diff[row]);
                        entries.emplace_back(row, ddSize + column, diff[row]);
                    }
                }
            }
        }
        assembler->setElectronQuantumPotential(lambda);
        result.jacobian.resize(4 * nodes, 4 * nodes);
        result.jacobian.setFromTriplets(entries.begin(), entries.end());
        return result;
    };

    VectorXd state(4 * nodes);
    state.head(nodes) = initial.psi / potentialScale;
    const bool compensatePotential = cfg_.compensatedPopulationState;
    VectorXd potentialLow = VectorXd::Zero(nodes);
    if (compensatePotential && initial.dgElectrostaticPotentialLow_V.size() == nodes) {
        for (int node = 0; node < nodes; ++node) {
            const auto normalized = numerics::SplitCoordinate{
                initial.psi[node], initial.dgElectrostaticPotentialLow_V[node]}
                .divided(potentialScale);
            state[node] = normalized.hi;
            potentialLow[node] = normalized.lo;
        }
    }
    if (population) {
        state.segment(nodes, nodes) = initial.n;
        state.segment(2 * nodes, nodes) = initial.p;
        pinCarriers(state);
        for (const auto& [node, target] : boundary.ohmicPopulation) {
            state[nodes + node] = target[1];
            state[2 * nodes + node] = target[2];
        }
    } else {
        state.head(ddSize) = assembler->pack({
            initial.psi / potentialScale,
            cfg_.warmStart ? VectorXd(initial.phin / potentialScale) : VectorXd::Zero(nodes).eval(),
            cfg_.warmStart ? VectorXd(initial.phip / potentialScale) : VectorXd::Zero(nodes).eval()});
        for (int node = 0; node < nodes; ++node) {
            if (cfg_.warmStart && initial.hasReferencedElectronQuasiFermi()) {
                const auto increment = numerics::SplitCoordinate::sum(
                    initial.electronQuasiFermiReferenceAt(node), -assembler->electronQuasiFermiReferenceAt(node))
                    .added({initial.phinIncrement[node], 0.}).divided(potentialScale);
                state[nodes + node] = increment.hi;
            }
            if (cfg_.warmStart && initial.hasReferencedHoleQuasiFermi()) {
                const auto increment = numerics::SplitCoordinate::sum(
                    initial.holeQuasiFermiReferenceAt(node), -assembler->holeQuasiFermiReferenceAt(node))
                    .added({initial.phipIncrement[node], 0.}).divided(potentialScale);
                state[2 * nodes + node] = increment.hi;
            }
        }
        for (const auto& [node, target] : boundary.phin)
            state[nodes + node] = (target * potentialScale - assembler->electronQuasiFermiReferenceAt(node)) / potentialScale;
        for (const auto& [node, target] : boundary.phip)
            state[2 * nodes + node] = (target * potentialScale - assembler->holeQuasiFermiReferenceAt(node)) / potentialScale;
    }
    for (const auto& [node, target] : boundary.psi) {
        state[node] = target;
        potentialLow[node] = 0.;
    }
    if (initial.electronQuantumPotentialLike.size() == nodes) {
        state.tail(nodes) = initial.electronQuantumPotentialLike / potentialScale;
    } else {
        const VectorXd lambda = initial.electronQuantumPotential.size() == nodes
            ? initial.electronQuantumPotential : VectorXd::Zero(nodes);
        state.tail(nodes) = (lambda - maxDrive) / potentialScale - state.head(nodes);
    }
    for (const auto& [node, lambda] : dirichletLambda_V)
        state[ddSize + node] = (lambda - maxDrive[node]) / potentialScale - state[node];

    VectorXd columns = VectorXd::Ones(state.size());
    if (population)
        columns.head(ddSize) = primaryVariableScales().device;
    const auto solveScaled = [](LinearSolver& solver,
                                 SparseMatrixd& jacobian, const VectorXd& residual,
                                 const VectorXd& columnScale) -> VectorXd {
        VectorXd rowMax = VectorXd::Zero(residual.size());
        for (int c = 0; c < jacobian.outerSize(); ++c)
            for (SparseMatrixd::InnerIterator entry(jacobian, c); entry; ++entry) {
                entry.valueRef() *= columnScale[c];
                rowMax[entry.row()] = std::max(rowMax[entry.row()], std::abs(entry.value()));
            }
        if (!(rowMax.array() > 0.).all() || !rowMax.allFinite())
            throw std::runtime_error("Invalid DG Jacobian row");
        for (int c = 0; c < jacobian.outerSize(); ++c)
            for (SparseMatrixd::InnerIterator entry(jacobian, c); entry; ++entry)
                entry.valueRef() /= rowMax[entry.row()];
        return solver.solve(jacobian, -residual.cwiseQuotient(rowMax)).cwiseProduct(columnScale);
    };

    // The CMD starts with Poisson+DG. At fixed reservoir QFs, n and p
    // follow their statistics as functions of the two potential fields.
    // Use this smaller system to initialize the cold four-field state.
    const bool coldQuantum = initial.electronQuantumPotentialLike.size() != nodes &&
        (initial.electronQuantumPotential.size() != nodes || initial.electronQuantumPotential.isZero(0.));
    if (coldQuantum && initial.phin.size() == nodes && initial.phip.size() == nodes) {
        // The two-field initializer must not evict the sweep's four-field pattern.
        LinearSolver equilibriumSolver;
        const auto equilibriumPopulations = [&](VectorXd& x, VectorXd* slopes) {
            for (int node = 0; node < nodes; ++node) {
                if (ni[node] <= 0.)
                    continue;
                const Real psi = x[node] * potentialScale;
                const Real transportPotential = -x[ddSize + node] * potentialScale - maxDrive[node];
                x[nodes + node] = electronDensity(
                    ni[node], Nc[node], transportPotential, initial.phin[node], vt, statistics);
                x[2 * nodes + node] = holeDensity(
                    ni[node], Nv[node], psi, initial.phip[node], vt, statistics);
                if (slopes != nullptr) {
                    (*slopes)[node] = -potentialScale / vt * electronDensityDerivativeEta(
                        ni[node], Nc[node], transportPotential, initial.phin[node], vt, statistics);
                    (*slopes)[nodes + node] = -potentialScale / vt * holeDensityDerivativeEta(
                        ni[node], Nv[node], psi, initial.phip[node], vt, statistics);
                }
            }
        };
        const auto reducedResidual = [&](const VectorXd& full) -> VectorXd {
            VectorXd residual(2 * nodes);
            residual << full.head(nodes), full.tail(nodes);
            return residual;
        };
        if (population)
            equilibriumPopulations(state, nullptr);
        VectorXd residual = reducedResidual(evaluate(state, false).residual);
        const Real psiScale = std::max(1., residual.head(nodes).stableNorm());
        const Real quantumScale = std::max(1., residual.tail(nodes).stableNorm());
        const auto equilibriumMerit = [&](const VectorXd& f) {
            return std::hypot(f.head(nodes).stableNorm() / psiScale,
                f.tail(nodes).stableNorm() / quantumScale);
        };
        Real measure = equilibriumMerit(residual);
        for (int iter = 0; iter < dg.maxIterations && measure > dg.relativeTolerance; ++iter) {
            VectorXd slopes = VectorXd::Zero(2 * nodes);
            if (population)
                equilibriumPopulations(state, &slopes);
            const auto full = evaluate(state, true, false);
            std::vector<Eigen::Triplet<Real>> entries;
            for (int c = 0; c < full.jacobian.outerSize(); ++c) {
                for (SparseMatrixd::InnerIterator entry(full.jacobian, c); entry; ++entry) {
                    if ((entry.row() >= nodes && entry.row() < ddSize) ||
                        (!population && c >= nodes && c < ddSize))
                        continue;
                    const int row = entry.row() < nodes ? entry.row() : entry.row() - 2 * nodes;
                    const bool fromElectron = c >= nodes && c < 2 * nodes;
                    const bool fromHole = c >= 2 * nodes && c < ddSize;
                    const int column = c % nodes + ((c >= ddSize || fromElectron) ? nodes : 0);
                    const Real derivative = (fromElectron || fromHole)
                        ? entry.value() * slopes[c - nodes] : entry.value();
                    entries.emplace_back(row, column, derivative);
                }
            }
            SparseMatrixd jacobian(2 * nodes, 2 * nodes);
            jacobian.setFromTriplets(entries.begin(), entries.end());
            VectorXd delta;
            try {
                delta = solveScaled(equilibriumSolver, jacobian, residual, VectorXd::Ones(2 * nodes));
            } catch (const std::runtime_error&) {
                break;
            }
            if (!delta.allFinite())
                break;
            bool accepted = false;
            for (Real alpha = 1.; alpha >= 1.e-6; alpha *= 0.5) {
                VectorXd trial = state;
                trial.head(nodes) += alpha * delta.head(nodes);
                trial.tail(nodes) += alpha * delta.tail(nodes);
                if (population)
                    equilibriumPopulations(trial, nullptr);
                VectorXd next = reducedResidual(evaluate(trial, false).residual);
                const Real value = equilibriumMerit(next);
                if (!(value < measure))
                    continue;
                state.swap(trial);
                residual.swap(next);
                measure = value;
                accepted = true;
                runtimeLogInfo("Two-field DG equilibrium update " + std::to_string(iter + 1) +
                               " merit=" + std::to_string(measure) + " alpha=" + std::to_string(alpha));
                break;
            }
            if (!accepted)
                break;
        }
    }
    const auto compensate = [&](VectorXd& residual, const VectorXd& state,
                                const VectorXd& low, const SparseMatrixd& tangent) {
        if (!compensatePotential)
            return;
        residual.head(nodes) = assembler->compensatedPopulationPoissonResidual(
            state.head(ddSize), low, boundary);
        // Only apply the sub-ULP electrostatic remainder about the ordinary
        // transport/DG evaluation. The remainder is far below the truncation
        // scale of their Jacobians. Accumulate it before narrowing to double.
        std::vector<experimental::ExpansionProductSum> correction(3 * nodes);
        for (int c = 0; c < std::min<int>(nodes, tangent.outerSize()); ++c) {
            if (low[c] == 0.)
                continue;
            for (SparseMatrixd::InnerIterator entry(tangent, c); entry; ++entry)
                if (entry.row() >= nodes)
                    correction[entry.row() - nodes].addProduct(entry.value(), low[c]);
        }
        for (int row = nodes; row < residual.size(); ++row) {
            auto& value = correction[row - nodes];
            value.add(residual[row]);
            residual[row] = value.value();
        }
    };
    auto initialSystem = evaluate(state, compensatePotential && !potentialLow.isZero(0.));
    VectorXd residual = std::move(initialSystem.residual);
    compensate(residual, state, potentialLow, initialSystem.jacobian);
    std::array<Real, 4> tolerances{};
    for (int block = 0; block < 4; ++block) {
        const Real initialNorm = residual.segment(block * nodes, nodes).lpNorm<Eigen::Infinity>();
        tolerances[block] = block == 3
            ? dg.absoluteTolerance_V + dg.relativeTolerance * std::max(initialNorm, 1.)
            : std::max(cfg_.abstol, cfg_.reltol * std::max(initialNorm, 1.));
    }
    const auto merit = [&](const VectorXd& residual) {
        if (!residual.allFinite())
            return std::numeric_limits<Real>::infinity();
        Real value = 0.;
        for (int block = 0; block < 4; ++block)
            value = std::max(value, residual.segment(block * nodes, nodes).lpNorm<Eigen::Infinity>() / tolerances[block]);
        return value;
    };
    const auto globalizationMerit = [&](const VectorXd& residual) {
        if (!residual.allFinite())
            return std::numeric_limits<Real>::infinity();
        Real value = 0.;
        for (int block = 0; block < 4; ++block)
            value = std::hypot(value,
                residual.segment(block * nodes, nodes).stableNorm() / tolerances[block]);
        return value;
    };
    NewtonResult result;
    Real lineMerit = globalizationMerit(residual);
    Real norm = merit(residual);
    result.initialResidualNorm = norm;
    result.failureDiagnostics.failureReason = "Four-field DG Newton iteration budget exhausted";
    SparseMatrixd cachedTangent;
    bool reuseTangent = false;
    auto localSolver = cfg_.sequentialLinearSolver ? nullptr : std::make_unique<LinearSolver>();
    LinearSolver& coupledSolver = cfg_.sequentialLinearSolver ? *cfg_.sequentialLinearSolver : *localSolver;
    for (int iter = 0; iter < cfg_.maxIter && norm > 1.; ++iter) {
        const bool reusedTangent = reuseTangent;
        reuseTangent = false;
        auto system = reusedTangent ? System{residual, cachedTangent} : evaluate(state, true);
        if (!population && !reusedTangent)
            cachedTangent = system.jacobian;
        compensate(system.residual, state, potentialLow, system.jacobian);
        auto scaledJacobian = compensatePotential ? system.jacobian : std::move(system.jacobian);
        VectorXd direction;
        try {
            direction = solveScaled(coupledSolver, scaledJacobian, system.residual, columns);
        } catch (const std::runtime_error& error) {
            result.failureDiagnostics.failureReason = error.what();
            break;
        }
        bool accepted = false;
        Real alpha = 1.;
        if (!direction.allFinite()) {
            result.failureDiagnostics.failureReason = "Nonfinite four-field DG correction";
            break;
        }
        // Integrate the population tangent in the local statistics coordinate.
        // Its derivative at zero is the same density correction as in the
        // linear solve. This avoids changing depleted rows to an unrelated
        // clipped density when a trial crosses zero.
        VectorXd electronEta = VectorXd::Zero(nodes);
        VectorXd holeEta = VectorXd::Zero(nodes);
        VectorXd electronEtaDirection = VectorXd::Zero(nodes);
        VectorXd holeEtaDirection = VectorXd::Zero(nodes);
        if (!population && cfg_.quasiFermiUpdateLimit_V > 0.) {
            const Real maxUpdate = direction.segment(nodes, 2 * nodes).lpNorm<Eigen::Infinity>() * potentialScale;
            if (maxUpdate > cfg_.quasiFermiUpdateLimit_V)
                alpha = cfg_.quasiFermiUpdateLimit_V / maxUpdate;
        }
        for (int node = 0; population && node < nodes; ++node) {
            if (electronFree[node]) {
                const Real n = state[nodes + node];
                electronEta[node] = usesFermiDirac(statistics)
                    ? inverseFermiDiracHalf(n / Nc[node], statistics) : std::log(n);
                const Real slope = usesFermiDirac(statistics)
                    ? Nc[node] * fermiDiracHalfDerivative(electronEta[node], statistics) : n;
                electronEtaDirection[node] = direction[nodes + node] / slope;
            }
            if (holeFree[node]) {
                const Real p = state[2 * nodes + node];
                holeEta[node] = usesFermiDirac(statistics)
                    ? inverseFermiDiracHalf(p / Nv[node], statistics) : std::log(p);
                const Real slope = usesFermiDirac(statistics)
                    ? Nv[node] * fermiDiracHalfDerivative(holeEta[node], statistics) : p;
                holeEtaDirection[node] = direction[2 * nodes + node] / slope;
            }
        }
        const Real previousLineMerit = lineMerit;
        for (int attempt = 0; attempt < 20; ++attempt, alpha *= 0.5) {
            VectorXd trial = state + alpha * direction;
            VectorXd trialLow = potentialLow;
            if (compensatePotential) {
                for (int node = 0; node < nodes; ++node) {
                    experimental::ExpansionProductSum sum;
                    sum.add(state[node]);
                    sum.add(potentialLow[node]);
                    sum.addProduct(alpha, direction[node]);
                    const auto value = sum.parts();
                    trial[node] = value.high;
                    trialLow[node] = value.low;
                }
            }
            for (int node = 0; population && node < nodes; ++node) {
                if (electronFree[node]) {
                    const Real eta = electronEta[node] +
                        std::clamp(alpha * electronEtaDirection[node], -6., 6.);
                    trial[nodes + node] = usesFermiDirac(statistics)
                        ? Nc[node] * fermiDiracHalf(eta, statistics) : std::exp(eta);
                }
                if (holeFree[node]) {
                    const Real eta = holeEta[node] +
                        std::clamp(alpha * holeEtaDirection[node], -6., 6.);
                    trial[2 * nodes + node] = usesFermiDirac(statistics)
                        ? Nv[node] * fermiDiracHalf(eta, statistics) : std::exp(eta);
                }
            }
            if (population)
                pinCarriers(trial);
            bool positive = trial.allFinite();
            for (int node = 0; population && positive && node < nodes; ++node) {
                if (ni[node] > 0.)
                    positive = trial[nodes + node] > 0. && trial[2 * nodes + node] > 0.;
            }
            if (!positive)
                continue;
            auto candidate = evaluate(trial, false);
            compensate(candidate.residual, trial, trialLow, system.jacobian);
            const Real next = merit(candidate.residual);
            const Real trialMerit = globalizationMerit(candidate.residual);
            if (next > 1. && !(trialMerit < lineMerit))
                continue;
            lineMerit = trialMerit;
            state.swap(trial);
            potentialLow.swap(trialLow);
            residual.swap(candidate.residual);
            norm = next;
            result.iters = iter + 1;
            NewtonIterationInfo info;
            info.iter = result.iters;
            info.residualNorm = norm;
            info.dampingFactor = alpha;
            info.lineSearchAttempts = attempt + 1;
            info.lineSearchAccepted = true;

            info.event = "four_field_dg_newton";
            result.trace.push_back(info);
            if (cfg_.diagnostics)
                result.history.push_back(info);
            runtimeLogInfo("Four-field DG Newton update " + std::to_string(result.iters) +
                           " merit=" + std::to_string(norm) + " alpha=" + std::to_string(alpha) +
                           " normalized rows psi/n/p/DG=" +
                           std::to_string(residual.head(nodes).lpNorm<Eigen::Infinity>() / tolerances[0]) + "/" +
                           std::to_string(residual.segment(nodes, nodes).lpNorm<Eigen::Infinity>() / tolerances[1]) + "/" +
                           std::to_string(residual.segment(2 * nodes, nodes).lpNorm<Eigen::Infinity>() / tolerances[2]) + "/" +
                           std::to_string(residual.tail(nodes).lpNorm<Eigen::Infinity>() / tolerances[3]));
            // A well-predicted full step permits one chord step. A stale
            // tangent never replaces the convergence test on the live residual.
            reuseTangent = !population && !reusedTangent && alpha == 1.
                && lineMerit <= 0.2 * previousLineMerit;
            accepted = true;
            break;
        }
        if (!accepted) {
            if (reusedTangent) {
                --iter;
                continue;
            }
            result.failureDiagnostics.failureReason = "Four-field DG line search rejected the direction";
            break;
        }
    }
    result.converged = norm <= 1.;
    if (result.converged) {
        result.convergenceReason = "four_field_dg_newton+electron_density_gradient";
        result.failureDiagnostics.failureReason.clear();
    }
    result.finalResidualNorm = norm;
    const auto ddBlocks = ResidualNorm::computeBlocks(residual.head(ddSize), nodes);
    result.finalBlockNorms = {ddBlocks.psi, ddBlocks.phin, ddBlocks.phip, ddBlocks.combined};
    assembler->setElectronQuantumPotential(packLambda(state));
    result.solution = makeSolution(*assembler, state.head(ddSize), result.iters);
    result.solution.electronQuantumPotentialLike = state.tail(nodes) * potentialScale;
    if (compensatePotential) {
        result.solution.populationPsi.resize(0);
        result.solution.populationPotentialScale_V = 0.;
        auto& low = result.solution.dgElectrostaticPotentialLow_V;
        low.resize(nodes);
        for (int node = 0; node < nodes; ++node) {
            const auto physical = numerics::SplitCoordinate{state[node], potentialLow[node]}
                .scaled(potentialScale);
            result.solution.psi[node] = physical.hi;
            low[node] = physical.lo;
        }
    }
    result.solution.converged = result.converged;
    return result;
}

} // namespace vela
