#include "vela/simulation/CoupledLoadLine.h"
#include "vela/simulation/CarrierReducedLoadLine.h"
#include "vela/simulation/CommonRelativeLinearSolve.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace vela;

TEST_CASE("Common-relative solve preserves affine direction and tiny relative part", "[coupled_load_line][common_relative]") {
    SparseMatrixd a(3,3);std::vector<Eigen::Triplet<Real>> e={{0,0,2},{0,1,-1},{1,0,-1},{1,1,3},{1,2,1},{2,2,1}};a.setFromTriplets(e.begin(),e.end());
    VectorXd b(3);b<<1,4,2;VectorXd common=a*VectorXd((VectorXd(3)<<1,1,0).finished());
    auto z=solveCommonRelative(a,b,{{0,1}},{common},"sparselu");
    REQUIRE((a*(z.shared+z.relative)-b).norm()<1e-13);
    REQUIRE((z.matrix*z.direction-b).norm()<1e-13);
    REQUIRE_THROWS(solveCommonRelative(a,b,{{0,1},{1,2}},{common,common},"sparselu"));
    // The relative unknown remains available even when adding it to a common
    // binary64 value would round it away. Split-state consumers apply separately.
    SparseMatrixd identity(2,2);identity.setIdentity();VectorXd rhs(2);rhs<<1,1e-30;
    VectorXd col(2);col<<1,0;
    auto tiny=solveCommonRelative(identity,rhs,{{0,1}},{col},"sparselu");
    REQUIRE(tiny.shared[1]==1);REQUIRE(tiny.relative[1]==Catch::Approx(1e-30));
    REQUIRE(tiny.shared[1]+tiny.relative[1]==1);
}

TEST_CASE("Carrier-reduced real circuit border preserves affine constraints and psi rows",
          "[coupled_load_line][carrier_reduced_border]")
{
    // Two nodes: inactive carrier constraints have nonzero prescribed updates.
    // Retain both psi constraints, both active carrier equations and the circuit.
    SparseMatrixd j(6,6);
    std::vector<Eigen::Triplet<Real>> e={{0,0,1},{1,1,1},{2,2,3},{2,3,2},
        {3,3,2},{4,4,4},{4,5,1},{5,5,5}};
    j.setFromTriplets(e.begin(),e.end());
    VectorXd r(6); r<<-.2,.3,1,-.8,-2,.5;
    VectorXd fv=VectorXd::Zero(6); fv[2]=-.7; fv[4]=.4;
    VectorXd ig(6); ig<<.01,.02,.03,.04,.05,.06;
    const auto a=assembleCarrierReducedLoadLine(j,r,fv,ig,.1,.03,0,0,10,2);
    REQUIRE(a.reduced.rows()==5);
    REQUIRE(a.map[0]>=0); REQUIRE(a.map[1]>=0); REQUIRE(a.map[6]>=0);
    REQUIRE(a.map[3]==-1); REQUIRE(a.map[5]==-1);
    const auto d=a.solve("sparselu");
    REQUIRE((a.full*d-a.rhs).lpNorm<Eigen::Infinity>()<1e-14);
    REQUIRE(d[3]==Catch::Approx(.4)); REQUIRE(d[5]==Catch::Approx(-.1));
    // A voltage-dependent carrier constraint requires affine substitution into
    // both the device equations and the circuit row, then a voltage-aware lift.
    fv[3]=.5;
    const auto b=assembleCarrierReducedLoadLine(j,r,fv,ig,.1,.03,0,0,10,2);
    REQUIRE(b.reduced.rows()==5);
    REQUIRE(b.map[3]==-1);
    REQUIRE(b.voltageSlope[3]==Catch::Approx(-.25));
    REQUIRE((b.full*b.solve("sparselu")-b.rhs).lpNorm<Eigen::Infinity>()<1e-14);
    REQUIRE_THROWS_AS(a.solve("pardiso_mkl_eliminated_keep_psi"),std::invalid_argument);
}

TEST_CASE("Carrier-reduced border rejects invalid layout and nonfinite coefficients",
          "[coupled_load_line][carrier_reduced_border]")
{
    SparseMatrixd j(3,3); j.setIdentity();
    VectorXd z=VectorXd::Zero(3);
    REQUIRE_THROWS_AS(assembleCarrierReducedLoadLine(j,z,z,z,0,0,0,0,1,2),std::invalid_argument);
    j.coeffRef(0,0)=std::numeric_limits<Real>::infinity();
    REQUIRE_THROWS_AS(assembleCarrierReducedLoadLine(j,z,z,z,0,0,0,0,1,1),std::runtime_error);
}

namespace {

ArclengthSystem foldedDeviceSystem()
{
    // Device branch V = x^3 - x has folds at x = +/-1/sqrt(3).
    ArclengthSystem system;
    system.residual = [](const VectorXd& x, Real voltage) {
        VectorXd residual(1);
        residual(0) = x(0) * x(0) * x(0) - x(0) - voltage;
        return residual;
    };
    system.parameterDerivative = [](const VectorXd&, Real) {
        VectorXd derivative(1);
        derivative(0) = -1.0;
        return derivative;
    };
    system.jacobian = [](const VectorXd& x, Real) {
        SparseMatrixd jacobian(1, 1);
        jacobian.insert(0, 0) = 3.0 * x(0) * x(0) - 1.0;
        return jacobian;
    };
    system.solveJacobian = [](const VectorXd& x, Real, const VectorXd& rhs,
                              VectorXd& solution) {
        const Real jacobian = 3.0 * x(0) * x(0) - 1.0;
        if (std::abs(jacobian) < 1.0e-12)
            return false;
        solution.resize(1);
        solution(0) = rhs(0) / jacobian;
        return true;
    };
    return system;
}

ArclengthScalarFunctional foldedCurrentFunctional()
{
    ArclengthScalarFunctional current;
    current.value = [](const VectorXd& x, Real) { return x(0); };
    current.linearize = [](const VectorXd& x, Real) {
        ArclengthScalarLinearization result;
        result.value = x(0);
        result.stateDerivative = VectorXd::Ones(1);
        result.parameterDerivative = 0.0;
        return result;
    };
    return current;
}

} // namespace

TEST_CASE("Bounded diagnostic updates preserve gates and observe only accepted states",
          "[coupled_load_line][diagnostic]")
{
    auto system = foldedDeviceSystem();
    auto current = foldedCurrentFunctional();
    CoupledLoadLineConfig config;
    config.linearSolver = "direct_bordered_qr";
    config.equationTolerance = 1e-3;
    config.loadLineTolerance_V = 1e-3;
    ArclengthState initial{VectorXd::Constant(1, 0.1), -0.0989};
    const Real outer = initial.lambda + 3.0 * initial.x(0);
    CoupledLoadLineNewton early(system, current.value, config, current);
    REQUIRE(early.solve(initial, outer, 3.0).iterations == 0);
    config.diagnosticMinUpdates = 2;
    config.maxIterations = 2;
    CoupledLoadLineNewton deeper(system, current.value, config, current);
    int observed = 0;
    const auto result = deeper.solve(initial, outer, 3.0,
        [&](int iteration, const ArclengthState& state) {
            REQUIRE(iteration == ++observed);
            REQUIRE(state.x.allFinite());
        });
    REQUIRE(result.converged);
    REQUIRE(result.iterations == 2);
    REQUIRE(observed == 2);
    REQUIRE(result.equationResidualNorm <= config.equationTolerance);
    REQUIRE(std::abs(result.loadLineResidual_V) <= config.loadLineTolerance_V);
    config.diagnosticMinUpdates = 4;
    REQUIRE_THROWS_AS(CoupledLoadLineNewton(system, current.value, config, current), std::invalid_argument);
    config.diagnosticMinUpdates = 1;
    config.maxIterations = 1;
    config.equationTolerance = 1e-16;
    config.loadLineTolerance_V = 1e-16;
    CoupledLoadLineNewton limited(system, current.value, config, current);
    const auto failed = limited.solve(initial, outer, 3.0);
    REQUIRE_FALSE(failed.converged);
    system.trialState = [](const VectorXd&, const VectorXd&, Real) -> VectorXd {
        throw std::runtime_error("diagnostic test rejects every trial");
    };
    CoupledLoadLineNewton rejected(system, current.value, config, current);
    observed = 0;
    const auto noUpdate = rejected.solve(initial, outer, 3.0,
        [&](int, const ArclengthState&) { ++observed; });
    REQUIRE_FALSE(noUpdate.converged);
    REQUIRE(observed == 0);
}

TEST_CASE("Trial observation preserves rejection and propagates export errors",
          "[coupled_load_line][diagnostic]")
{
    auto system = foldedDeviceSystem();
    const auto current = foldedCurrentFunctional();
    CoupledLoadLineConfig config;
    config.linearSolver = "direct_bordered_qr";
    config.maxLineSearchSteps = 3;
    ArclengthState initial{VectorXd::Constant(1, 0.1), -0.09};
    const Real outer = initial.lambda + 3.0 * initial.x(0);
    // Every trial remains on the same equation defect and has a non-decreasing
    // circuit defect: the observer must not accidentally accept a stalled step.
    system.trialState = [](const VectorXd& x, const VectorXd&, Real) { return x; };
    system.limitUpdate = [](const VectorXd&, VectorXd&, Real& dv) { dv = 0; };
    CoupledLoadLineNewton solver(system, current.value, config, current);
    const auto baseline = solver.solve(initial, outer, 3.0);
    int observed = 0;
    const auto result = solver.solve(initial, outer, 3.0, {},
        [&](const CoupledLoadLineTrial& info, const ArclengthState& state) {
            REQUIRE(info.lineSearch == observed++);
            REQUIRE_FALSE(info.accepted);
            REQUIRE(info.equationResidual == Catch::Approx(0.009));
            REQUIRE(info.equationResidualIndex == 0);
            REQUIRE(state.lambda == initial.lambda);
        });
    REQUIRE(observed == 4);
    REQUIRE_FALSE(result.converged);
    REQUIRE(result.failureReason == baseline.failureReason);
    REQUIRE(result.state.x == baseline.state.x);
    REQUIRE_THROWS_AS(solver.solve(initial, outer, 3.0, {},
        [](const CoupledLoadLineTrial&, const ArclengthState&) {
            throw std::runtime_error("export failed");
        }), std::runtime_error);

    system = foldedDeviceSystem();
    CoupledLoadLineNewton normal(system, current.value, config, current);
    int accepted = 0;
    const auto solved = normal.solve(initial, outer, 3.0, {},
        [&](const CoupledLoadLineTrial& info, const ArclengthState&) {
            if (info.accepted) ++accepted;
        });
    REQUIRE(solved.converged);
    REQUIRE(accepted == solved.iterations);
}

TEST_CASE("Supply tangent predicts the full state through a device fold",
          "[coupled_load_line][tangent]")
{
    auto system = foldedDeviceSystem();
    auto current = foldedCurrentFunctional();
    CoupledLoadLineConfig config;
    config.linearSolver = "direct_bordered";
    CoupledLoadLineNewton solver(system, current.value, config, current);
    ArclengthState parent{VectorXd::Constant(1, 1.0 / std::sqrt(3.0)), 0.0};
    parent.lambda = std::pow(parent.x(0), 3) - parent.x(0);
    const Real outer = parent.lambda + 2.0 * parent.x(0);
    auto p = solver.tangentPredict(parent, outer + 1e-4, 2.0);
    REQUIRE(p.applied);
    REQUIRE(p.alpha == 1.0);
    REQUIRE(p.state.x(0) - parent.x(0) == Catch::Approx(5e-5).margin(1e-13));
    REQUIRE(std::abs(p.voltageIncrement) < 1e-12);
    auto half = solver.tangentPredict(parent, outer + 5e-5, 2.0);
    REQUIRE(half.applied);
    REQUIRE(p.deviceResidual / half.deviceResidual == Catch::Approx(4.0).epsilon(1e-3));
    auto solved = solver.solve(p.state, outer + 1e-4, 2.0);
    REQUIRE(solved.converged);
}

TEST_CASE("Supply tangent respects voltage trust region and rejects unqualified parents",
          "[coupled_load_line][tangent]")
{
    auto system = foldedDeviceSystem();
    auto current = foldedCurrentFunctional();
    CoupledLoadLineConfig config;
    config.linearSolver = "direct_bordered";
    config.maxVoltageUpdate_V = 1e-4;
    CoupledLoadLineNewton solver(system, current.value, config, current);
    ArclengthState parent{VectorXd::Zero(1), 0.0};
    const auto p = solver.tangentPredict(parent, 0.01, 2.0);
    REQUIRE(p.applied);
    REQUIRE(std::abs(p.voltageIncrement) <= config.maxVoltageUpdate_V);
    REQUIRE(p.state.x(0) == Catch::Approx(-p.voltageIncrement));
    parent.lambda = 0.1;
    const auto rejected = solver.tangentPredict(parent, 0.2, 2.0);
    REQUIRE_FALSE(rejected.applied);
    REQUIRE(rejected.reason == "parent_outside_device_tolerance");
    REQUIRE(rejected.state.lambda == parent.lambda);
    REQUIRE(rejected.state.x == parent.x);
    system.jacobian = {};
    CoupledLoadLineConfig fallbackConfig;
    CoupledLoadLineNewton missing(system, current.value, fallbackConfig, current);
    parent.lambda = 0.0;
    REQUIRE_FALSE(missing.tangentPredict(parent, 0.1, 2.0).applied);
}

TEST_CASE("Supply tangent includes explicit voltage current derivative and device limits",
          "[coupled_load_line][tangent]")
{
    auto system = foldedDeviceSystem();
    system.residual = [](const VectorXd& x, Real v) { return VectorXd(x.array() - 2*v); };
    system.parameterDerivative = [](const VectorXd&, Real) { return VectorXd::Constant(1,-2); };
    system.jacobian = [](const VectorXd&, Real) {
        SparseMatrixd j(1,1); j.insert(0,0)=1.; return j;
    };
    system.limitUpdate = [](const VectorXd&, VectorXd& dx, Real& dv) { dx*=.5; dv*=.5; };
    ArclengthScalarFunctional current;
    current.value = [](const VectorXd& x, Real v) { return 3*x(0)+4*v; };
    current.linearize = [current](const VectorXd& x, Real v) {
        return ArclengthScalarLinearization{current.value(x,v),VectorXd::Constant(1,3),4.};
    };
    CoupledLoadLineConfig config; config.linearSolver="direct_bordered";
    CoupledLoadLineNewton solver(system,current.value,config,current);
    ArclengthState parent{VectorXd::Zero(1),0.};
    auto p=solver.tangentPredict(parent,.21,2.);
    REQUIRE(p.applied);
    REQUIRE(p.voltageIncrement == Catch::Approx(.005));
    REQUIRE(p.state.x(0) == Catch::Approx(.01));
    REQUIRE(p.deviceResidual < 1e-12);
    auto q=solver.solve(p.state,.21,2.);
    REQUIRE(q.converged);
    REQUIRE(q.state.lambda == Catch::Approx(.01).margin(1e-7));
}

TEST_CASE("Coupled residual filter applies independent block envelopes",
          "[coupled_load_line][filter]")
{
    CoupledLoadLineConfig config;
    config.equationTolerance = 1.0e-6;
    config.loadLineTolerance_V = 1.0e-6;
    config.filterGamma = 1.0e-4;
    config.filterEnvelopeFactor = 2.0;

    // Regression for the SLOT-LDMOS trace: the load block improves slightly,
    // but the already-converged device block jumps by more than six orders of
    // magnitude.  A load-dominated shared envelope used to accept this step.
    REQUIRE_FALSE(coupledResidualFilterAccepts(
        config, 8.2476318868628675e-7, -1.6685485839818739,
        1.9599596827576904, -1.6620307947662241, 1.0 / 256.0));

    // Improving one block is sufficient only while the other remains inside
    // its own bounded envelope.
    REQUIRE(coupledResidualFilterAccepts(
        config, 8.0e-7, -1.0,
        1.5e-6, -0.9, 0.5));
    REQUIRE_FALSE(coupledResidualFilterAccepts(
        config, 1.0, -5.0e-7,
        0.5, -1.0, 0.5));

    // The solve-entry envelope is immutable. A sequence cannot ratchet an
    // initially converged PDE block from 1x to 2x, then 4x, merely because the
    // circuit row improves a little on every accepted step.
    REQUIRE_FALSE(coupledResidualFilterAccepts(
        config, 1.9e-6, -0.9,
        2.1e-6, -0.8, 0.5,
        1.0, 1.0e6));

    config.inexactDeviceForcingEnabled = true;
    config.inexactDeviceToleranceMax = 5.0e-6;
    config.inexactLoadActivationRatio = 100.0;
    REQUIRE(coupledEffectiveEquationTolerance(config, 0.0) ==
            Catch::Approx(5.0e-6));
    REQUIRE(coupledEffectiveEquationTolerance(config, 1.0e-4) ==
            Catch::Approx(1.0e-6));

    // Closing the scalar row activates the forcing allowance in the trial
    // state.  The same trial remains illegal with strict block scaling.
    REQUIRE(coupledResidualFilterAccepts(
        config, 5.0e-7, 2.0e-5,
        4.0e-6, 0.0, 1.0));
    config.inexactDeviceForcingEnabled = false;
    REQUIRE_FALSE(coupledResidualFilterAccepts(
        config, 5.0e-7, 2.0e-5,
        4.0e-6, 0.0, 1.0));
}

TEST_CASE("CoupledLoadLineNewton crosses a voltage-controlled fold",
          "[coupled_load_line]")
{
    CoupledLoadLineConfig config;
    config.maxIterations = 40;
    config.equationTolerance = 1.0e-11;
    config.loadLineTolerance_V = 1.0e-11;
    config.currentDirectionalStep = 1.0e-6;
    config.maxVoltageUpdate_V = 0.25;
    CoupledLoadLineNewton solver(
        foldedDeviceSystem(),
        [](const VectorXd& x, Real) { return x(0); },
        config);

    ArclengthState state;
    state.x = VectorXd::Zero(1);
    state.lambda = 0.0;
    for (Real outer = 0.25; outer <= 3.0; outer += 0.25) {
        const CoupledLoadLineResult result = solver.solve(state, outer, 3.0);
        REQUIRE(result.converged);
        REQUIRE(std::abs(result.loadLineResidual_V) <= 1.0e-11);
        state = result.state;
    }

    REQUIRE(state.x(0) == Catch::Approx(1.0).margin(1.0e-9));
    REQUIRE(state.lambda == Catch::Approx(0.0).margin(1.0e-9));
}

TEST_CASE("CoupledLoadLineNewton rejects invalid circuit input",
          "[coupled_load_line]")
{
    CoupledLoadLineConfig config;
    CoupledLoadLineNewton solver(
        foldedDeviceSystem(),
        [](const VectorXd& x, Real) { return x(0); },
        config);
    ArclengthState state;
    state.x = VectorXd::Zero(1);
    state.lambda = 0.0;
    REQUIRE_THROWS_AS(solver.solve(state, 1.0, 0.0), std::invalid_argument);
}

TEST_CASE("Coupled load-line backtracks nonlinear trial maps and closes both equations",
          "[coupled_load_line][positive_density]")
{
    auto system = foldedDeviceSystem();
    int calls = 0;
    int rejected = 0;
    system.trialState = [&](const VectorXd& x, const VectorXd& dx, Real alpha) {
        ++calls;
        // A domain-limited nonlinear path has the original Newton tangent.
        // Rejected trial evaluations must backtrack, not abort the solve.
        if (alpha > 0.5) {
            ++rejected;
            throw std::domain_error("Trial outside admissible domain");
        }
        VectorXd trial = x + alpha * dx;
        trial(0) += 0.1 * std::pow(alpha * dx(0), 2);
        return trial;
    };
    CoupledLoadLineConfig config;
    config.maxIterations = 100;
    config.equationTolerance = 1.0e-10;
    config.loadLineTolerance_V = 1.0e-10;
    auto current = foldedCurrentFunctional();
    const auto value = current.value;
    CoupledLoadLineNewton solver(system, value, config, std::move(current));
    ArclengthState initial;
    initial.x = VectorXd::Zero(1);
    initial.lambda = 0.0;
    const auto result = solver.solve(initial, 3.0, 3.0);
    REQUIRE(calls > 0);
    REQUIRE(rejected > 0);
    REQUIRE(result.converged);
    REQUIRE(system.residual(result.state.x, result.state.lambda).norm() <= 1.0e-10);
    REQUIRE(std::abs(result.state.lambda + 3.0 * result.state.x(0) - 3.0) <= 1.0e-10);
    REQUIRE(result.state.x(0) == Catch::Approx(1.0).margin(1.0e-9));
}

TEST_CASE("direct bordered load-line solve remains regular at a device fold",
          "[coupled_load_line][direct_bordered]")
{
    const Real xFold = 1.0 / std::sqrt(3.0);
    const Real innerAtFold = xFold * xFold * xFold - xFold;
    const Real resistance = 3.0;
    ArclengthState initial;
    initial.x = VectorXd::Constant(1, xFold);
    initial.lambda = innerAtFold;
    const Real targetOuter = innerAtFold + resistance * xFold + 0.05;

    CoupledLoadLineConfig schurConfig;
    schurConfig.equationTolerance = 1.0e-11;
    schurConfig.loadLineTolerance_V = 1.0e-11;
    CoupledLoadLineNewton schurSolver(
        foldedDeviceSystem(),
        [](const VectorXd& x, Real) { return x(0); },
        schurConfig);
    const CoupledLoadLineResult schur =
        schurSolver.solve(initial, targetOuter, resistance);
    REQUIRE_FALSE(schur.converged);
    REQUIRE(schur.failureReason == "device Jacobian solve failed");

    CoupledLoadLineConfig directConfig = schurConfig;
    directConfig.linearSolver = "direct_bordered";
    directConfig.maxVoltageUpdate_V = 0.25;
    ArclengthScalarFunctional current = foldedCurrentFunctional();
    const auto currentValue = current.value;
    CoupledLoadLineNewton directSolver(
        foldedDeviceSystem(), currentValue, directConfig, std::move(current));
    const CoupledLoadLineResult direct =
        directSolver.solve(initial, targetOuter, resistance);
    REQUIRE(direct.converged);
    REQUIRE(direct.equationResidualNorm <= directConfig.equationTolerance);
    REQUIRE(std::abs(direct.loadLineResidual_V) <=
            directConfig.loadLineTolerance_V);
}

TEST_CASE("direct bordered scalar row supports current control",
          "[coupled_load_line][direct_bordered][current_control]")
{
    CoupledLoadLineConfig config;
    config.linearSolver = "direct_bordered";
    config.voltageCoefficient = 0.0;
    config.equationTolerance = 1.0e-11;
    config.loadLineTolerance_V = 1.0e-11;
    config.maxVoltageUpdate_V = 0.5;
    ArclengthScalarFunctional current = foldedCurrentFunctional();
    const auto currentValue = current.value;
    CoupledLoadLineNewton solver(
        foldedDeviceSystem(), currentValue, config, std::move(current));

    ArclengthState initial;
    initial.x = VectorXd::Zero(1);
    initial.lambda = 0.0;
    const CoupledLoadLineResult result = solver.solve(initial, 0.5, 1.0);

    REQUIRE(result.converged);
    REQUIRE(result.state.x(0) == Catch::Approx(0.5).margin(1.0e-10));
    REQUIRE(result.state.lambda == Catch::Approx(-0.375).margin(1.0e-10));
    REQUIRE(std::abs(result.loadLineResidual_V) <=
            config.loadLineTolerance_V);
}

TEST_CASE("load-closed state can use configured inexact device forcing",
          "[coupled_load_line][inexact_forcing]")
{
    CoupledLoadLineConfig config;
    config.linearSolver = "direct_bordered";
    config.voltageCoefficient = 0.0;
    config.equationTolerance = 1.0e-6;
    config.loadLineTolerance_V = 1.0e-12;
    config.inexactDeviceForcingEnabled = true;
    config.inexactDeviceToleranceMax = 5.0e-6;
    config.inexactLoadActivationRatio = 10.0;
    ArclengthScalarFunctional current = foldedCurrentFunctional();
    const auto currentValue = current.value;
    CoupledLoadLineNewton solver(
        foldedDeviceSystem(), currentValue, config, std::move(current));

    ArclengthState initial;
    initial.x = VectorXd::Zero(1);
    initial.lambda = -2.0e-6;
    const CoupledLoadLineResult result = solver.solve(initial, 0.0, 1.0);

    REQUIRE(result.converged);
    REQUIRE(result.iterations == 0);
    REQUIRE(result.inexactConvergenceUsed);
    REQUIRE(result.equationResidualNorm == Catch::Approx(2.0e-6));
    REQUIRE(result.effectiveEquationTolerance == Catch::Approx(5.0e-6));
}

TEST_CASE("inexact forcing can follow the device tangent without polishing its residual",
          "[coupled_load_line][inexact_forcing][direct_bordered]")
{
    CoupledLoadLineConfig config;
    config.linearSolver = "direct_bordered";
    config.voltageCoefficient = 0.0;
    config.equationTolerance = 1.0e-6;
    config.loadLineTolerance_V = 1.0e-12;
    config.inexactDeviceForcingEnabled = true;
    config.inexactDeviceToleranceMax = 5.0e-6;
    config.inexactLoadActivationRatio = 100.0;
    config.zeroConvergedDeviceResidual = true;
    config.linearizationAuditEnabled = true;
    ArclengthScalarFunctional current = foldedCurrentFunctional();
    const auto currentValue = current.value;
    CoupledLoadLineNewton solver(
        foldedDeviceSystem(), currentValue, config, std::move(current));

    ArclengthState initial;
    initial.x = VectorXd::Zero(1);
    initial.lambda = -5.0e-7;
    const Real targetCurrent = 2.0e-11;
    REQUIRE(coupledEffectiveEquationTolerance(config, -targetCurrent) >
            5.0e-7);
    REQUIRE(coupledLinearizedDeviceTolerance(config, -targetCurrent) ==
            Catch::Approx(
                coupledEffectiveEquationTolerance(config, -targetCurrent)));

    const CoupledLoadLineResult result =
        solver.solve(initial, targetCurrent, 1.0);

    REQUIRE(result.converged);
    REQUIRE(result.iterations == 1);
    REQUIRE(result.trace.front().deviceResidualZeroed);
    REQUIRE(result.trace.front().deviceLinearizationResidualNorm <= 1.0e-14);
    REQUIRE(std::abs(result.trace.front().circuitLinearizationResidual) <=
            1.0e-14);
    REQUIRE(result.trace.front().currentDirectionalRelativeError <= 1.0e-10);
    REQUIRE(result.trace.front().deviceJvpErrorNorm <= 1.0e-12);
    REQUIRE(result.state.x(0) == Catch::Approx(targetCurrent));
    REQUIRE(result.equationResidualNorm == Catch::Approx(5.0e-7));
    REQUIRE(std::abs(result.loadLineResidual_V) <=
            config.loadLineTolerance_V);
}

TEST_CASE("secant predictor rejects an off-manifold extrapolation from a converged state",
          "[coupled_load_line][predictor][globalization]")
{
    CoupledLoadLineConfig config;
    config.voltageCoefficient = 0.0;
    config.equationTolerance = 1.0e-6;
    config.inexactDeviceForcingEnabled = true;
    config.inexactDeviceToleranceMax = 5.0e-6;
    config.filterEnvelopeFactor = 2.0;
    CoupledLoadLineNewton solver(
        foldedDeviceSystem(),
        [](const VectorXd& x, Real) { return x(0); },
        config);

    ArclengthState previous;
    previous.x = VectorXd::Zero(1);
    previous.lambda = 0.0;
    ArclengthState current;
    current.x = VectorXd::Constant(1, 0.5);
    current.lambda = -0.375;

    // Linear secant extrapolation to x=0.6 has an O(1e-2) nonlinear device
    // defect although the current point is exact, so it must be discarded.
    const ArclengthState predicted =
        solver.secantPredict(previous, current, 0.6, 1.0);
    REQUIRE(predicted.x(0) == Catch::Approx(current.x(0)));
    REQUIRE(predicted.lambda == Catch::Approx(current.lambda));
}

TEST_CASE("load-line secant predictor limits the complete bordered increment",
          "[coupled_load_line][predictor]")
{
    ArclengthSystem system = foldedDeviceSystem();
    system.limitUpdate = [](const VectorXd&, VectorXd& deltaState,
                            Real& deltaVoltage) {
        const Real scale = std::min(
            Real{1.0}, 0.5 / deltaState.lpNorm<Eigen::Infinity>());
        deltaState *= scale;
        deltaVoltage *= scale;
    };
    CoupledLoadLineNewton solver(
        std::move(system),
        [](const VectorXd& x, Real) { return x(0); },
        CoupledLoadLineConfig{});

    ArclengthState previous;
    previous.x = VectorXd::Zero(1);
    previous.lambda = 0.0;
    ArclengthState current;
    current.x = VectorXd::Constant(1, 2.0);
    current.lambda = 1.0;

    // With R=1, the previous and current outer voltages are 0 and 3 V.
    // A 6 V target gives a raw ratio of one and increments (dx,dV)=(2,1).
    // The trust region scales both components by 1/4.
    const ArclengthState predicted =
        solver.secantPredict(previous, current, 6.0, 1.0);
    REQUIRE(predicted.x(0) == Catch::Approx(2.5));
    REQUIRE(predicted.lambda == Catch::Approx(1.25));
}
