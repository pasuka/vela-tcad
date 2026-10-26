#include "vela/simulation/CoupledLoadLine.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace vela;

TEST_CASE("Report-only backtracking records damping but applies the full Newton step",
          "[coupled_load_line][report_only]")
{
    ArclengthSystem s;
    s.residual=[](const VectorXd& x,Real){return VectorXd::Constant(1,x(0)*x(0)-2.);};
    s.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Zero(1);};
    s.solveJacobian=[](const VectorXd& x,Real,const VectorXd& b,VectorXd& y){y=b/(2*x(0));return true;};
    auto current=[](const VectorXd&,Real){return 0.;};
    ArclengthState seed{VectorXd::Constant(1,.1),0.};
    CoupledLoadLineConfig c;c.maxIterations=1;c.equationTolerance=1e-10;
    const auto damped=CoupledLoadLineNewton(s,current,c).solve(seed,0.,1.);
    REQUIRE(damped.trace.front().acceptedAlpha<1.);
    c.lineSearchMode="report_only";c.dampingFactor=.01;
    const auto full=CoupledLoadLineNewton(s,current,c).solve(seed,0.,1.);
    REQUIRE_FALSE(full.converged);
    REQUIRE(full.state.x(0)==Catch::Approx(10.05));
    REQUIRE(full.trace.front().acceptedAlpha==1.);
    REQUIRE(full.trace.front().reportedAlpha==damped.trace.front().acceptedAlpha);
    REQUIRE(full.trace.front().fullStepMerit>full.trace.front().merit);
    REQUIRE(seed.x(0)==.1);
    c.lineSearchMode="full_step";
    const auto undamped=CoupledLoadLineNewton(s,current,c).solve(seed,0.,1.);
    REQUIRE(undamped.state.x(0)==full.state.x(0));
    REQUIRE(undamped.lineSearchAttempts==1);
    REQUIRE(undamped.trace.front().acceptedAlpha==1.);
    REQUIRE(undamped.trace.front().reportedAlpha==0.);
    c.lineSearchMode="report_only";
    c.maxLineSearchSteps=0;
    const auto noSuggestion=CoupledLoadLineNewton(s,current,c).solve(seed,0.,1.);
    REQUIRE(noSuggestion.trace.front().reportedAlpha==0.);
    REQUIRE(noSuggestion.state.x(0)==full.state.x(0));
    c.maxIterations=20;
    const auto root=CoupledLoadLineNewton(s,current,c).solve(seed,0.,1.);
    REQUIRE(root.converged);
    REQUIRE(root.state.x(0)==Catch::Approx(std::sqrt(2.)).margin(1e-10));
    REQUIRE(root.equationResidualNorm<c.equationTolerance);

    // Domain failures must not silently select a smaller trial in report mode.
    s.residual=[](const VectorXd& x,Real){
        if(x(0)>3.)throw std::domain_error("model domain");
        return VectorXd::Constant(1,x(0)*x(0)-2.);};
    const auto invalid=CoupledLoadLineNewton(s,current,c).solve(seed,0.,1.);
    REQUIRE_FALSE(invalid.converged);
    REQUIRE(invalid.state.x(0)==seed.x(0));
    REQUIRE(invalid.trace.front().acceptedAlpha==0.);
    REQUIRE(invalid.failureReason.find("full step invalid")!=std::string::npos);
    s.residual=[](const VectorXd&,Real){return VectorXd::Constant(1,std::numeric_limits<Real>::infinity());};
    REQUIRE_FALSE(CoupledLoadLineNewton(s,current,c).solve(seed,0.,1.).converged);
}

TEST_CASE("Density update stopping keeps raw direction and requires load closure",
          "[coupled_load_line][density_update_convergence]")
{
    // Constant row rescaling keeps the root but makes its raw residual large.
    ArclengthSystem system;
    system.residual=[](const VectorXd& x, Real){return VectorXd::Constant(1,1e12*(x(0)*x(0)-2.));};
    system.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Zero(1);};
    system.jacobian=[](const VectorXd& x,Real){SparseMatrixd J(1,1);J.insert(0,0)=2e12*x(0);return J;};
    system.solveJacobian=[](const VectorXd& x,Real,const VectorXd& b,VectorXd& y){y=b/(2e12*x(0));return true;};
    auto current=[](const VectorXd&,Real){return 0.;};
    CoupledLoadLineConfig cfg; cfg.equationTolerance=1e-20; cfg.maxIterations=8;
    DensityUpdateErrorEvaluator evaluator=[](const VectorXd&,const VectorXd& dx,const VectorXd& after,Real,Real){
        DensityUpdateError e;e.valid=true;e.rms=1e5*std::abs(dx(0))/(std::abs(after(0))+1.);return e;};
    ArclengthState seed{VectorXd::Constant(1,1.),0.};
    const auto strict = CoupledLoadLineNewton(system,current,cfg).solve(seed,.1,1.);
    REQUIRE_FALSE(strict.converged);
    cfg.updateConvergence.mode="report";
    const auto report = CoupledLoadLineNewton(system,current,cfg,{},evaluator).solve(seed,.1,1.);
    REQUIRE(report.converged == strict.converged);
    REQUIRE(report.state.x(0) == strict.state.x(0));
    REQUIRE(report.iterations == strict.iterations);
    cfg.updateConvergence.mode="enforce";
    const auto update = CoupledLoadLineNewton(system,current,cfg,{},evaluator).solve(seed,.1,1.);
    REQUIRE(update.converged); REQUIRE(update.convergenceReason=="coupled_density_update");
    REQUIRE(update.updateError.valid); REQUIRE(update.updateError.rms<1.);
    REQUIRE(update.equationResidualNorm>cfg.equationTolerance);
    REQUIRE(update.state.x(0)==Catch::Approx(std::sqrt(2.)).epsilon(1e-9));
    REQUIRE(std::abs(update.loadLineResidual_V)<=cfg.loadLineTolerance_V);
    // A tiny damping factor cannot manufacture a small raw update metric.
    cfg.dampingFactor=1e-8; cfg.maxIterations=1; cfg.equationTolerance=1e-20;
    Real capturedRaw=0.,capturedApplied=0.;
    auto observe=[&](const VectorXd& before,const VectorXd& dx,const VectorXd& after,Real,Real){
        capturedRaw=dx(0);capturedApplied=after(0)-before(0);
        DensityUpdateError e;e.valid=true;e.rms=1e5*std::abs(dx(0));return e;};
    const auto damped=CoupledLoadLineNewton(system,current,cfg,{},observe).solve(seed,.1,1.);
    REQUIRE_FALSE(damped.converged); REQUIRE(std::abs(capturedRaw)>.1);
    REQUIRE(std::abs(capturedApplied)<1e-7);
    // Even a small valid update never waives the external load equation.
    auto tiny=[](const VectorXd&,const VectorXd&,const VectorXd&,Real,Real){
        DensityUpdateError e;e.valid=true;e.rms=0.;return e;};
    const auto openLoad=CoupledLoadLineNewton(system,current,cfg,{},tiny).solve(seed,.1,1.);
    REQUIRE_FALSE(openLoad.converged);
    REQUIRE_THROWS(CoupledLoadLineNewton(system,current,cfg));
    cfg.updateConvergence.digits=16;
    REQUIRE_THROWS(CoupledLoadLineNewton(system,current,cfg,{},tiny));
}

TEST_CASE("resistor analytic diode audits raw and representable applied directions",
          "[coupled_load_line][resistor_contract]")
{
    // x is junction voltage; the boundary equality couples it to the port.
    constexpr Real thermal = .02585, saturation = 1e-9, resistance = 1e7;
    ArclengthSystem device;
    device.residual = [](const VectorXd& x, Real v) { return VectorXd::Constant(1, x(0)-v); };
    device.jacobian = [](const VectorXd&, Real) {
        SparseMatrixd j(1,1); j.insert(0,0)=1; return j;
    };
    device.parameterDerivative = [](const VectorXd&, Real) { return VectorXd::Constant(1,-1); };
    device.solveJacobian = [](const VectorXd&,Real,const VectorXd& b,VectorXd& y) { y=b; return true; };
    ArclengthScalarFunctional port;
    port.value = [](const VectorXd& x, Real) { return saturation*std::expm1(x(0)/thermal); };
    port.linearize = [port](const VectorXd& x, Real v) {
        ArclengthScalarLinearization l;
        l.value=port.value(x,v);
        l.stateDerivative=VectorXd::Constant(1,saturation/thermal*std::exp(x(0)/thermal));
        l.parameterDerivative=0; return l;
    };
    const Real targetInner=.06;
    const Real targetOuter=targetInner+resistance*saturation*std::expm1(targetInner/thermal);
    ArclengthState seed{VectorXd::Zero(1),0};
    CoupledLoadLineConfig cfg;
    cfg.linearSolver="direct_bordered";
    cfg.equationTolerance=1e-12; cfg.loadLineTolerance_V=1e-12;
    cfg.maxVoltageUpdate_V=.012; // Forces protection of the whole direction.
    const auto baseline=CoupledLoadLineNewton(device,port.value,cfg,port).solve(seed,targetOuter,resistance);
    cfg.linearizationAuditEnabled=true;
    const auto report=CoupledLoadLineNewton(device,port.value,cfg,port).solve(seed,targetOuter,resistance);
    REQUIRE(baseline.converged); REQUIRE(report.converged);
    REQUIRE(report.state.x(0)==baseline.state.x(0));
    REQUIRE(report.state.lambda==baseline.state.lambda);
    REQUIRE(report.iterations==baseline.iterations);
    REQUIRE(report.lineSearchAttempts==baseline.lineSearchAttempts);
    REQUIRE(report.state.lambda==Catch::Approx(targetInner).margin(1e-11));
    REQUIRE(seed.x(0)==0); REQUIRE(seed.lambda==0);
    bool sawLimited=false, sawNonlinear=false;
    for (std::size_t k=0;k<report.trace.size();++k) {
        const auto& t=report.trace[k];
        REQUIRE(t.acceptedAlpha==baseline.trace[k].acceptedAlpha);
        if(t.rawStateUpdate.size()==0) continue; // terminal convergence record
        REQUIRE(t.componentwiseBackwardError<=1e-10);
        REQUIRE(t.rawPredictedDeviceResidual.norm()<=1e-12);
        sawLimited |= std::abs(t.rawVoltageUpdate_V)>cfg.maxVoltageUpdate_V;
        for(const auto& trial:t.trials) {
            REQUIRE(trial.appliedStateUpdate(0)==Catch::Approx(trial.appliedVoltageUpdate_V).margin(1e-15));
            REQUIRE(trial.predictedDeviceResidual(0)==Catch::Approx(trial.actualDeviceResidual(0)).margin(1e-15));
            sawNonlinear |= std::abs(trial.predictedLoadResidual_V-trial.actualLoadResidual_V)>1e-8;
        }
    }
    REQUIRE(sawLimited); REQUIRE(sawNonlinear);
    cfg.directRefinementSteps=2;
    const auto refined=CoupledLoadLineNewton(device,port.value,cfg,port).solve(seed,targetOuter,resistance);
    REQUIRE(refined.converged);
    // Already accurate directions must not be perturbed by unnecessary corrections.
    REQUIRE(refined.state.x(0)==report.state.x(0));
    REQUIRE(refined.state.lambda==report.state.lambda);
    REQUIRE(refined.iterations==report.iterations);
    REQUIRE(refined.lineSearchAttempts==report.lineSearchAttempts);
    for(const auto& t:refined.trace)
        if(t.rawStateUpdate.size()) REQUIRE(t.componentwiseBackwardError<=1e-10);
    cfg.directRefinementSteps=-1;
    REQUIRE_THROWS_AS(CoupledLoadLineNewton(device,port.value,cfg,port),std::invalid_argument);
    cfg.directRefinementSteps=9;
    REQUIRE_THROWS_AS(CoupledLoadLineNewton(device,port.value,cfg,port),std::invalid_argument);
    cfg.linearSolver="schur";
    cfg.directRefinementSteps=1;
    REQUIRE_THROWS_AS(CoupledLoadLineNewton(device,port.value,cfg,port),std::invalid_argument);
    cfg.directRefinementSteps=0;
    const auto schur=CoupledLoadLineNewton(device,port.value,cfg,port).solve(seed,targetOuter,resistance);
    REQUIRE(schur.converged);
    REQUIRE(schur.state.lambda==Catch::Approx(report.state.lambda).margin(1e-11));
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
