#include "vela/numerics/WeightedSparseRows.h"
#include "vela/simulation/CoupledLoadLine.h"
#include "vela/solver/PopulationState.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace vela;

TEST_CASE("Automatic derivative recovery rejects an invalid full direction then restores Newton", "[coupled_load_line][population_core]") {
    ArclengthSystem s;
    s.residual=[](const VectorXd& x,Real){return VectorXd::Constant(1,x[0]*x[0]-2.);};
    s.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Zero(1);};
    s.solveJacobian=[](const VectorXd& x,Real,const VectorXd& b,VectorXd& y){y=b/(2*x[0]);return true;};
    int fullCalls=0;
    s.jacobian=[&](const VectorXd& x,Real){SparseMatrixd j(1,1);j.insert(0,0)=++fullCalls==1?0.:2*x[0];return j;};
    s.bootstrapJacobian=[](const VectorXd&,Real){SparseMatrixd j(1,1);j.insert(0,0)=2.;return j;};
    s.automaticDerivativeRecovery=true;
    ArclengthScalarFunctional port;
    port.value=[](const VectorXd&,Real){return 0.;};
    port.linearize=[](const VectorXd& x,Real){ArclengthScalarLinearization p;p.stateDerivative=VectorXd::Zero(x.size());return p;};
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";
    c.equationTolerance=1e-12;c.maxIterations=25;
    const ArclengthState initial{VectorXd::Constant(1,.5),0.};
    const auto result=CoupledLoadLineNewton(s,port.value,c,port).solve(initial,0.,1.);
    INFO(result.failureReason);REQUIRE(result.converged);
    REQUIRE(result.state.x[0]==Catch::Approx(std::sqrt(2.)).margin(1e-12));
    REQUIRE(initial.x[0]==.5);
    bool bootstrap=false,trial=false,rejected=false;
    for(const auto& t:result.trace){bootstrap|=t.derivativeMode=="bootstrap";trial|=t.derivativeMode=="trial";rejected|=t.iteration>0&&t.acceptedAlpha==0.;}
    REQUIRE(bootstrap);REQUIRE(trial);REQUIRE(rejected);REQUIRE(fullCalls>1);
    c.lineSearchMode="merit";REQUIRE_THROWS(CoupledLoadLineNewton(s,port.value,c,port));
}

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
    cfg.updateConvergence.mode="enforce_joint";
    const auto joint=CoupledLoadLineNewton(system,current,cfg,{},evaluator).solve(seed,.1,1.);
    REQUIRE_FALSE(joint.converged);
    cfg.updateConvergence.mode="enforce";
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
    INFO(baseline.failureReason);
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
        INFO(result.failureReason);
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
    INFO(direct.failureReason);
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

    INFO(result.failureReason);
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

    INFO(result.failureReason);
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

    INFO(result.failureReason);
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

TEST_CASE("Coupled continuation uses positive population predictor and rolls back invalid predictions", "[coupled_load_line][population_core]") {
    ArclengthSystem s;
    s.residual=[](const VectorXd& x,Real){return VectorXd::Zero(x.size()).eval();};
    s.parameterDerivative=[](const VectorXd& x,Real){return VectorXd::Zero(x.size()).eval();};
    s.solveJacobian=[](const VectorXd&,Real,const VectorXd& b,VectorXd& y){y=b;return true;};
    s.predictState=[](const VectorXd& a,const VectorXd& b,Real t){return predictPopulation(a,b,t);};
    CoupledLoadLineConfig cfg;
    CoupledLoadLineNewton solver(s,[](const VectorXd&,Real){return 0.;},cfg);
    VectorXd a(6),b(6);a<<0.,0.,100.,0.,1e-100,0.;b<<.01,0.,10.,0.,1e-100,0.;
    ArclengthState previous{a,0.},current{b,.01};
    const auto prediction=solver.secantPredict(previous,current,.02,1.);
    REQUIRE(prediction.x[2]==Catch::Approx(1.).epsilon(1e-13));REQUIRE(prediction.lambda==.02);
    REQUIRE(prediction.x[0]==.02);REQUIRE(prediction.x[3]==0.);
    const auto noHistory=solver.secantPredict(current,current,.02,1.);REQUIRE(noHistory.x==current.x);
    current.x[2]=1e-300;
    const auto invalid=solver.secantPredict(previous,current,.02,1.);REQUIRE(invalid.x==current.x);REQUIRE(invalid.lambda==current.lambda);
}

TEST_CASE("Terminal derivative uses only selected conservation rows", "[coupled_load_line][terminal_rows]") {
    SparseMatrixd j(3,2);j.insert(0,0)=2.;j.insert(0,1)=-3.;
    j.insert(1,0)=std::numeric_limits<Real>::quiet_NaN();
    j.insert(2,1)=std::numeric_limits<Real>::infinity();
    VectorXd w=VectorXd::Zero(3);w[0]=4.;Eigen::Index excluded=0;
    const auto d=weightedSparseRows(j,w,&excluded);
    REQUIRE(excluded==2);REQUIRE(d.allFinite());REQUIRE(d[0]==8.);REQUIRE(d[1]==-12.);
    w[1]=1.;REQUIRE_THROWS(weightedSparseRows(j,w));
}
TEST_CASE("Terminal row summation retains the finite sparse derivative", "[coupled_load_line][terminal_rows]") {
    SparseMatrixd j(3,2);j.insert(0,0)=2.;j.insert(1,0)=.5;j.insert(2,1)=-3.;
    VectorXd w(3);w<<1.,-2.,0.;
    REQUIRE((weightedSparseRows(j,w)-(j.transpose()*w)).isZero(0.));
}
TEST_CASE("Bordered equilibration preserves a finite product with overflowing intermediate", "[coupled_load_line][scaled_columns]") {
    ArclengthSystem s;const double goal=2e-320;
    s.residual=[goal](const VectorXd& x,Real){return VectorXd::Constant(1,1e303*(x[0]-goal));};
    s.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Zero(1);};
    s.jacobian=[](const VectorXd&,Real){SparseMatrixd j(1,1);j.insert(0,0)=1e303;return j;};
    s.solveJacobian=[](const VectorXd&,Real,const VectorXd& b,VectorXd& y){y=b/1e303;return true;};
    s.columnScales=[](const ArclengthState&){VectorXd s(2);s<<1e-320,1.;return s;};
    ArclengthScalarFunctional p;p.value=[](const VectorXd&,Real){return 0.;};
    p.linearize=[](const VectorXd&,Real){ArclengthScalarLinearization p;p.stateDerivative=VectorXd::Zero(1);return p;};
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";c.equationTolerance=1e-30;c.maxIterations=3;
    const auto r=CoupledLoadLineNewton(s,p.value,c,p).solve({VectorXd::Constant(1,1e-320),0.},0.,1.);
    INFO(r.failureReason);REQUIRE(r.converged);REQUIRE(std::abs(r.state.x[0]-goal)<=2*std::numeric_limits<double>::denorm_min());
    REQUIRE(std::isfinite(scaledProduct(1e17,1e303,1e-320)));
    REQUIRE(scaledProduct(1e17,1e303,1e-320)>0.);
}
TEST_CASE("Relative population callback publishes physical increments", "[coupled_load_line][scaled_columns]") {
    ArclengthSystem s;const double goal=2e-320;int relativeCalls=0;
    s.residual=[goal](const VectorXd& x,Real){return VectorXd::Constant(1,1e303*(x[0]-goal));};
    s.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Zero(1);};
    s.jacobian=[](const VectorXd&,Real)->SparseMatrixd{throw std::runtime_error("absolute assembly must not be used");};
    s.solveJacobian=[](const VectorXd&,Real,const VectorXd&,VectorXd&){return false;};
    s.jacobianInColumns=[&](const VectorXd&,Real,const VectorXd& t){++relativeCalls;SparseMatrixd j(1,1);j.insert(0,0)=1e303*t[0];return j;};
    s.bootstrapJacobianInColumns=s.jacobianInColumns;
    s.columnScales=[](const ArclengthState& st){VectorXd s(2);s<<st.x[0],1.;return s;};
    ArclengthScalarFunctional p;p.value=[](const VectorXd&,Real){return 0.;};
    p.linearize=[](const VectorXd&,Real)->ArclengthScalarLinearization{throw std::runtime_error("absolute terminal must not be used");};
    p.linearizeInColumns=[](const VectorXd&,Real,const VectorXd&){ArclengthScalarLinearization p;p.stateDerivative=VectorXd::Zero(1);return p;};
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";c.equationTolerance=1e-30;c.maxIterations=3;
    const auto r=CoupledLoadLineNewton(s,p.value,c,p).solve({VectorXd::Constant(1,1e-320),0.},0.,1.);
    INFO(r.failureReason);REQUIRE(r.converged);REQUIRE(relativeCalls>0);REQUIRE(std::abs(r.state.x[0]-goal)<=2*std::numeric_limits<double>::denorm_min());
}
TEST_CASE("Bordered diagnostics distinguish missing support from weighted overflow", "[coupled_load_line][scaled_columns]") {
    auto check=[](double diagonal,double physical){
        ArclengthSystem s;s.residual=[](const VectorXd&,Real){return VectorXd::Ones(1);};
        s.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Zero(1);};
        s.solveJacobian=[](const VectorXd&,Real,const VectorXd&,VectorXd&){return false;};
        s.jacobian=[=](const VectorXd&,Real){SparseMatrixd j(1,1);j.insert(0,0)=diagonal;return j;};
        s.columnScales=[=](const ArclengthState&){VectorXd a(2);a<<physical,1.;return a;};
        ArclengthScalarFunctional p;p.value=[](const VectorXd&,Real){return 0.;};
        p.linearize=[](const VectorXd&,Real){ArclengthScalarLinearization p;p.stateDerivative=VectorXd::Zero(1);return p;};
        CoupledLoadLineConfig c;c.linearSolver="direct_bordered";
        return CoupledLoadLineNewton(s,p.value,c,p).solve({VectorXd::Zero(1),0.},0.,1.).failureReason;
    };
    REQUIRE(check(0.,1.).find("structural zero row")!=std::string::npos);
    REQUIRE(check(1e300,1e300).find("non-finite weighted row magnitude")!=std::string::npos);
    REQUIRE(check(1e-300,1e-300).find("weighted row magnitude underflow")!=std::string::npos);
}

TEST_CASE("External resistor Quasistationary keeps the full inner voltage correction", "[coupled_load_line][quasistationary]") {
    // A linear device x=V with I=x/2 has the unique load-line solution
    // V + R I = 3, R=1, hence V=2. No nonlinear globalization is needed.
    ArclengthSystem system;
    system.residual=[](const VectorXd& x,Real v){return VectorXd::Constant(1,x[0]-v);};
    system.jacobian=[](const VectorXd&,Real){SparseMatrixd j(1,1);j.insert(0,0)=1.;return j;};
    system.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Constant(1,-1.);};
    system.solveJacobian=[](const VectorXd&,Real,const VectorXd& b,VectorXd& x){x=b;return true;};
    ArclengthScalarFunctional current;
    current.value=[](const VectorXd& x,Real){return .5*x[0];};
    current.linearize=[](const VectorXd&,Real){ArclengthScalarLinearization a;a.stateDerivative=VectorXd::Constant(1,.5);a.parameterDerivative=0.;return a;};
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";
    c.maxIterations=2;c.maxVoltageUpdate_V=.1;c.applyDeviceUpdateLimit=false;
    ArclengthState seed{VectorXd::Zero(1),0.};
    const auto limited=CoupledLoadLineNewton(system,current.value,c,current).solve(seed,3.,1.);
    REQUIRE_FALSE(limited.converged);
    c.applyVoltageUpdateLimit=false;
    const auto full=CoupledLoadLineNewton(system,current.value,c,current).solve(seed,3.,1.);
    REQUIRE(full.converged);
    REQUIRE(full.state.lambda==Catch::Approx(2.));
    REQUIRE(full.state.x[0]==Catch::Approx(2.));
    REQUIRE(std::abs(full.loadLineResidual_V)<1e-12);
}

TEST_CASE("Load-line report retains the complete residual branch", "[coupled_load_line][load_report]") {
    // Exact device root x=V, but the initial state does not satisfy V+I=3.
    ArclengthSystem system;
    system.residual=[](const VectorXd& x,Real v){return VectorXd::Constant(1,x[0]-v);};
    system.jacobian=[](const VectorXd&,Real){SparseMatrixd j(1,1);j.insert(0,0)=1.;return j;};
    system.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Constant(1,-1.);};
    system.solveJacobian=[](const VectorXd&,Real,const VectorXd& b,VectorXd& x){x=b;return true;};
    ArclengthScalarFunctional current;
    current.value=[](const VectorXd& x,Real){return .5*x[0];};
    current.linearize=[](const VectorXd&,Real){ArclengthScalarLinearization a;a.stateDerivative=VectorXd::Constant(1,.5);a.parameterDerivative=0.;return a;};
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";
    c.applyVoltageUpdateLimit=false;c.maxIterations=2;
    ArclengthState seed{VectorXd::Zero(1),0.};
    const auto closed=CoupledLoadLineNewton(system,current.value,c,current).solve(seed,3.,1.);
    REQUIRE(closed.converged);REQUIRE(closed.state.lambda==Catch::Approx(2.));
    c.enforceLoadLineConvergence=false;
    const auto reported=CoupledLoadLineNewton(system,current.value,c,current).solve(seed,3.,1.);
    REQUIRE(reported.converged);REQUIRE(reported.state.lambda==Catch::Approx(2.));
    REQUIRE(std::abs(reported.loadLineResidual_V)<1e-12);REQUIRE(reported.iterations==1);
    // An already closed complete system may still accept zero updates.
    const auto unchanged=CoupledLoadLineNewton(system,current.value,c,current).solve(reported.state,3.,1.);
    REQUIRE(unchanged.converged);REQUIRE(unchanged.iterations==0);
    // Every changed external target must re-close the circuit even from an exact PDE root.
    const auto next=CoupledLoadLineNewton(system,current.value,c,current).solve(reported.state,6.,1.);
    REQUIRE(next.converged);REQUIRE(next.iterations==1);
    REQUIRE(next.state.lambda==Catch::Approx(4.));
    // A cached-looking small update sample cannot authorize iter=0 in a new solve.
    c.updateConvergence.mode="enforce";
    DensityUpdateErrorEvaluator metric=[](const VectorXd&,const VectorXd&,const VectorXd&,Real,Real) {
        DensityUpdateError e;e.valid=true;e.rms=.01;return e;};
    const auto orNext=CoupledLoadLineNewton(system,current.value,c,current,metric).solve(next.state,9.,1.);
    REQUIRE(orNext.converged);REQUIRE(orNext.iterations==1);
    REQUIRE(orNext.state.lambda==Catch::Approx(6.));
    c.updateConvergence.mode="off";
    // Report does not bypass the device equation or another physical gate.
    seed.x[0]=1.;
    const auto corrected=CoupledLoadLineNewton(system,current.value,c,current).solve(seed,3.,1.);
    REQUIRE(corrected.converged);REQUIRE(corrected.iterations>0);
    REQUIRE(corrected.state.lambda==Catch::Approx(2.));
    system.convergenceSatisfied=[](const VectorXd&,Real){return false;};
    REQUIRE_FALSE(CoupledLoadLineNewton(system,current.value,c,current).solve(seed,3.,1.).converged);
    system.convergenceSatisfied={};
    current.value=[](const VectorXd&,Real){return std::numeric_limits<Real>::infinity();};
    REQUIRE_FALSE(CoupledLoadLineNewton(system,current.value,c,current).solve(seed,3.,1.).converged);
}

TEST_CASE("Complete residual and report-only load gate have independent update branches", "[coupled_load_line][load_report]") {
    // x=V, I=x^2. The first Newton step from V=1 toward V+I=2.001
    // leaves a quadratic circuit defect, while its relative voltage update is small.
    ArclengthSystem s;
    s.residual=[](const VectorXd& x,Real v){return VectorXd::Constant(1,x[0]-v);};
    s.jacobian=[](const VectorXd&,Real){SparseMatrixd j(1,1);j.insert(0,0)=1.;return j;};
    s.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Constant(1,-1.);};
    s.solveJacobian=[](const VectorXd&,Real,const VectorXd& b,VectorXd& x){x=b;return true;};
    ArclengthScalarFunctional current;
    current.value=[](const VectorXd& x,Real){return x[0]*x[0];};
    current.linearize=[](const VectorXd& x,Real){ArclengthScalarLinearization a;a.stateDerivative=VectorXd::Constant(1,2*x[0]);a.parameterDerivative=0.;return a;};
    DensityUpdateErrorEvaluator metric=[](const VectorXd&,const VectorXd& dx,const VectorXd& after,Real v,Real dv){
        DensityUpdateError e;e.valid=true;e.counts[0]=1;e.counts[3]=1;
        e.blockRms[0]=std::abs(dx[0])/(1+std::abs(after[0]));
        e.blockRms[3]=std::abs(dv)/(1+std::abs(v));
        e.rms=std::hypot(e.blockRms[0],e.blockRms[3])/std::sqrt(2.);return e;};
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";
    c.applyVoltageUpdateLimit=false;c.maxIterations=4;c.loadLineTolerance_V=1e-12;
    c.updateConvergence.mode="enforce";c.enforceLoadLineConvergence=false;
    ArclengthState seed{VectorXd::Constant(1,1.),1.};
    const auto report=CoupledLoadLineNewton(s,current.value,c,current,metric).solve(seed,2.001,1.);
    REQUIRE(report.converged);REQUIRE(report.iterations==1);
    REQUIRE(report.convergenceReason=="coupled_density_update");
    REQUIRE(std::abs(report.loadLineResidual_V)>c.loadLineTolerance_V);
    c.enforceLoadLineConvergence=true;
    const auto enforce=CoupledLoadLineNewton(s,current.value,c,current,metric).solve(seed,2.001,1.);
    REQUIRE(enforce.converged);REQUIRE(enforce.iterations>report.iterations);
    REQUIRE(std::abs(enforce.loadLineResidual_V)<=c.loadLineTolerance_V);
    c.enforceLoadLineConvergence=false;c.updateConvergence.mode="off";
    const auto residualOnly=CoupledLoadLineNewton(s,current.value,c,current).solve(seed,2.001,1.);
    REQUIRE(residualOnly.converged);REQUIRE(residualOnly.iterations>report.iterations);
}

TEST_CASE("Failed bootstrap solve restores full Jacobian at the unchanged target",
          "[coupled_load_line][population_core][linear_restore]") {
    ArclengthSystem s;
    s.residual=[](const VectorXd& x,Real){return VectorXd::Constant(1,x[0]-1.);};
    s.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Zero(1);};
    s.solveJacobian=[](const VectorXd&,Real,const VectorXd&,VectorXd&){return false;};
    int bootstrapCalls=0,fullCalls=0;bool fullFails=false,rejectTrial=false;
    const ArclengthState initial{VectorXd::Zero(1),.2};
    auto matrix=[&](bool full,const VectorXd& x,Real v,const VectorXd& columns){
        REQUIRE(x[0]==initial.x[0]);REQUIRE(v==initial.lambda);
        if(full)++fullCalls;else ++bootstrapCalls;
        SparseMatrixd j(1,1);j.insert(0,0)=full&&!fullFails?columns[0]:0.;return j;
    };
    s.jacobian=[&](const VectorXd& x,Real v){return matrix(true,x,v,VectorXd::Ones(1));};
    s.bootstrapJacobian=[&](const VectorXd& x,Real v){return matrix(false,x,v,VectorXd::Ones(1));};
    s.automaticDerivativeRecovery=true;
    s.makeTrial=[&](const ArclengthState& st,const VectorXd& dx,Real dv,Real alpha){
        if(rejectTrial)throw std::domain_error("diagnostic domain failure");
        ArclengthTrial t;t.state={st.x+alpha*dx,st.lambda+alpha*dv};return t;
    };
    ArclengthScalarFunctional p;p.value=[](const VectorXd&,Real){return 0.;};
    p.linearize=[](const VectorXd&,Real){ArclengthScalarLinearization a;a.stateDerivative=VectorXd::Zero(1);return a;};
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";
    c.equationTolerance=1e-12;c.maxIterations=1;c.applyVoltageUpdateLimit=false;
    c.updateConvergence.mode="enforce_joint";c.enforceLoadLineConvergence=true;
    double rms=.25;
    DensityUpdateErrorEvaluator metric=[&](const VectorXd&,const VectorXd&,const VectorXd&,Real,Real){
        DensityUpdateError e;e.valid=true;e.rms=rms;return e;
    };
    SECTION("Same iteration restores the ordinary bordered matrix") {}
    SECTION("Relative population columns also select the full callback") {
        s.columnScales=[](const ArclengthState&){VectorXd t(2);t<<1e-200,1.;return t;};
        s.jacobianInColumns=[&](const VectorXd& x,Real v,const VectorXd& t){return matrix(true,x,v,t);};
        s.bootstrapJacobianInColumns=[&](const VectorXd& x,Real v,const VectorXd& t){return matrix(false,x,v,t);};
        s.jacobian=[](const VectorXd&,Real)->SparseMatrixd{throw std::runtime_error("absolute full path forbidden");};
        s.bootstrapJacobian=[](const VectorXd&,Real)->SparseMatrixd{throw std::runtime_error("absolute bootstrap path forbidden");};
        p.linearizeInColumns=[](const VectorXd&,Real,const VectorXd&){ArclengthScalarLinearization a;a.stateDerivative=VectorXd::Zero(1);return a;};
    }
    SECTION("Two bad matrices return without mode cycling or state changes") {fullFails=true;}
    SECTION("An invalid restored candidate cannot renew the same-state probe") {rejectTrial=true;c.maxIterations=10;}
    SECTION("Restoration does not bypass the update stopping gate") {rms=2.;}
    SECTION("Restoration does not bypass an independent physical gate") {s.convergenceSatisfied=[](const VectorXd&,Real){return false;};}
    const auto result=CoupledLoadLineNewton(s,p.value,c,p,metric).solve(initial,.7,1.);
    INFO(result.failureReason);
    REQUIRE_FALSE(result.trace.front().bootstrapLinearFailureReason.empty());
    REQUIRE((result.trace.front().derivativeReason=="bootstrap_linear_failure_restore" || rejectTrial));
    if(fullFails||rejectTrial) {
        REQUIRE_FALSE(result.converged);REQUIRE(result.state.x[0]==initial.x[0]);
        REQUIRE(result.state.lambda==initial.lambda);REQUIRE(fullCalls==1);
        REQUIRE(bootstrapCalls==(rejectTrial?2:1));
        if(fullFails)REQUIRE(result.failureReason.find("full Jacobian retry failed")!=std::string::npos);
    } else {
        REQUIRE(bootstrapCalls==1);REQUIRE(fullCalls>=1);
        REQUIRE(result.state.x[0]==Catch::Approx(1.));REQUIRE(result.state.lambda==Catch::Approx(.7));
        REQUIRE(result.iterations==1);REQUIRE(result.trace.front().acceptedAlpha==1.);
        REQUIRE(result.equationResidualNorm<=c.equationTolerance);
        REQUIRE(std::abs(result.loadLineResidual_V)<=c.loadLineTolerance_V);
        REQUIRE(result.converged==(rms<1.&&!s.convergenceSatisfied));
    }
}

TEST_CASE("Protected bootstrap prediction is reviewed before changing the state",
          "[coupled_load_line][population_core][model_review]") {
    const ArclengthState seed{(VectorXd(2)<<0.,.1).finished(),.2};
    double bootSlope=-.001,rms=.25;bool rejectFull=false,fullFails=false;
    int fullCalls=0,bootCalls=0;
    ArclengthSystem s;
    s.residual=[](const VectorXd& x,Real){VectorXd f(2);f<<x[0]+x[1]-1.,x[1]-1.;return f;};
    s.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Zero(2);};
    s.solveJacobian=[](const VectorXd&,Real,const VectorXd&,VectorXd&){return false;};
    s.prescribedRows=[](const VectorXd&,Real){return std::vector<int>{};};
    auto matrix=[&](bool full,const VectorXd& x,Real v,const VectorXd& columns){
        REQUIRE((x-seed.x).norm()==0.);REQUIRE(v==seed.lambda);
        if(full)++fullCalls;else ++bootCalls;
        SparseMatrixd j(2,2);j.insert(0,0)=columns[0];j.insert(0,1)=columns[1];
        j.insert(1,1)=(full?(fullFails?0.:1.):bootSlope)*columns[1];return j;
    };
    s.jacobian=[&](const VectorXd& x,Real v){return matrix(true,x,v,VectorXd::Ones(2));};
    s.bootstrapJacobian=[&](const VectorXd& x,Real v){return matrix(false,x,v,VectorXd::Ones(2));};
    s.automaticDerivativeRecovery=true;
    s.makeTrial=[&](const ArclengthState& before,const VectorXd& dx,Real dv,Real alpha){
        if(rejectFull&&dx[1]>0)throw std::domain_error("full candidate outside domain");
        ArclengthTrial t;t.state={before.x+alpha*dx,before.lambda+alpha*dv};
        if(t.state.x[1]<=0){t.state.x[1]=.01*before.x[1];t.protectedHoles=1;}return t;
    };
    ArclengthScalarFunctional port;port.value=[](const VectorXd&,Real){return 0.;};
    port.linearize=[](const VectorXd&,Real){ArclengthScalarLinearization p;p.stateDerivative=VectorXd::Zero(2);return p;};
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";
    c.applyVoltageUpdateLimit=false;c.maxIterations=1;c.equationTolerance=1e-12;
    c.bootstrapModelGuardEnabled=true;c.updateConvergence.mode="enforce_joint";
    DensityUpdateErrorEvaluator metric=[&](const VectorXd&,const VectorXd&,const VectorXd&,Real,Real){DensityUpdateError e;e.valid=true;e.rms=rms;return e;};
    bool legacy=false,relative=false,feasible=false;
    SECTION("Review uses the original state, target and one applied-update budget"){}
    SECTION("Extreme relative columns retain their direct assembly contract"){
        relative=true;s.columnScales=[](const ArclengthState&){return (VectorXd(3)<<1e-200,1e-200,1.).finished();};
        s.jacobianInColumns=[&](const VectorXd& x,Real v,const VectorXd& t){return matrix(true,x,v,t);};
        s.bootstrapJacobianInColumns=[&](const VectorXd& x,Real v,const VectorXd& t){return matrix(false,x,v,t);};
        s.jacobian=[](const VectorXd&,Real)->SparseMatrixd{throw std::runtime_error("absolute full forbidden");};
        s.bootstrapJacobian=[](const VectorXd&,Real)->SparseMatrixd{throw std::runtime_error("absolute bootstrap forbidden");};
        port.linearizeInColumns=[](const VectorXd&,Real,const VectorXd&){ArclengthScalarLinearization a;a.stateDerivative=VectorXd::Zero(2);return a;};
    }
    SECTION("The explicit disabled control preserves the unsafe baseline") {legacy=true;c.bootstrapModelGuardEnabled=false;}
    SECTION("A rejected full candidate cannot reissue an unchanged-state probe") {rejectFull=true;c.maxIterations=10;}
    SECTION("An unusable full matrix exits without a mode cycle") {fullFails=true;c.maxIterations=10;}
    SECTION("Full review cannot bypass RMS") {rms=2.;}
    SECTION("Full review cannot bypass the independent physical gate") {s.convergenceSatisfied=[](const VectorXd&,Real){return false;};}
    SECTION("Feasible residuals keep the original stopping gate") {
        feasible=true;bootSlope=-100.;
        s.makeTrial=[](const ArclengthState& before,const VectorXd&,Real dv,Real alpha){ArclengthTrial t;t.state={(VectorXd(2)<<0.,1.).finished(),before.lambda+alpha*dv};return t;};
    }
    const auto a=CoupledLoadLineNewton(s,port.value,c,port,metric).solve(seed,.7,1.);
    INFO(a.failureReason);REQUIRE((seed.x-(VectorXd(2)<<0.,.1).finished()).norm()==0.);
    if(legacy){REQUIRE_FALSE(a.converged);REQUIRE(a.state.x[0]>100.);REQUIRE(fullCalls==0);return;}
    if(feasible){REQUIRE(a.converged);REQUIRE(fullCalls==0);REQUIRE_FALSE(a.trace.front().bootstrapModelRetry);return;}
    REQUIRE(a.trace.front().bootstrapModelRetry);REQUIRE(a.trace.front().candidateRejectedByRecovery);
    REQUIRE(a.trace.front().acceptedAlpha==0.);REQUIRE(a.trace.front().bootstrapModelAfter>a.trace.front().bootstrapModelBefore);
    REQUIRE(a.trace.front().iteration==0);REQUIRE(a.trace[1].iteration==0);
    if(rejectFull||fullFails){REQUIRE_FALSE(a.converged);REQUIRE((a.state.x-seed.x).norm()==0.);REQUIRE(a.state.lambda==seed.lambda);REQUIRE(fullCalls==1);REQUIRE(a.trace.size()<6);}
    else{REQUIRE(a.iterations==1);REQUIRE(a.state.x[0]==Catch::Approx(0.).margin(1e-12));REQUIRE(a.state.x[1]==Catch::Approx(1.));REQUIRE(a.state.lambda==Catch::Approx(.7));REQUIRE(a.converged==(rms<1.&&!s.convergenceSatisfied));}
    (void)relative;
}

TEST_CASE("Bootstrap review includes the circuit prediction even when device prediction improves",
          "[coupled_load_line][model_review]") {
    ArclengthSystem s;s.automaticDerivativeRecovery=true;
    s.residual=[](const VectorXd& x,Real){return VectorXd::Constant(1,x[0]-1.);};
    s.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Zero(1);};
    s.solveJacobian=[](const VectorXd&,Real,const VectorXd&,VectorXd&){return false;};
    s.prescribedRows=[](const VectorXd&,Real){return std::vector<int>{};};
    s.jacobian=[](const VectorXd&,Real){SparseMatrixd j(1,1);j.insert(0,0)=1.;return j;};
    s.bootstrapJacobian=[](const VectorXd&,Real){SparseMatrixd j(1,1);j.insert(0,0)=-.001;return j;};
    s.makeTrial=[](const ArclengthState& x,const VectorXd& dx,Real dv,Real alpha){ArclengthTrial t;t.state={x.x+alpha*dx,x.lambda+alpha*dv};if(t.state.x[0]<=0)t.state.x[0]=.01*x.x[0];return t;};
    ArclengthScalarFunctional p;p.value=[](const VectorXd& x,Real){return x[0];};p.linearize=[](const VectorXd&,Real){ArclengthScalarLinearization p;p.stateDerivative=VectorXd::Ones(1);return p;};
    CoupledLoadLineConfig c;c.bootstrapModelGuardEnabled=true;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";c.applyVoltageUpdateLimit=false;c.maxIterations=1;
    auto a=CoupledLoadLineNewton(s,p.value,c,p).solve({VectorXd::Constant(1,.1),.2},1.5,1.);
    REQUIRE(a.trace.front().bootstrapModelRetry);REQUIRE(a.converged);REQUIRE(a.iterations==1);
    REQUIRE(a.state.x[0]==Catch::Approx(1.));REQUIRE(a.state.lambda==Catch::Approx(.5));
}

TEST_CASE("Recovery history is explicit, copied and independent of the model guard",
          "[coupled_load_line][derivative_history]") {
    ArclengthSystem s;s.automaticDerivativeRecovery=true;
    s.residual=[](const VectorXd& x,Real v){return VectorXd::Constant(1,x[0]-v);};
    s.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Constant(1,-1.);};
    s.solveJacobian=[](const VectorXd&,Real,const VectorXd&,VectorXd&){return false;};
    int full=0,boot=0;s.jacobian=[&](const VectorXd&,Real){++full;SparseMatrixd j(1,1);j.insert(0,0)=1.;return j;};
    s.bootstrapJacobian=[&](const VectorXd&,Real){++boot;SparseMatrixd j(1,1);j.insert(0,0)=2.;return j;};
    ArclengthScalarFunctional p;p.value=[](const VectorXd&,Real){return 0.;};p.linearize=[](const VectorXd&,Real){ArclengthScalarLinearization p;p.stateDerivative=VectorXd::Zero(1);return p;};
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";c.applyVoltageUpdateLimit=false;c.maxIterations=2;
    derivative_recovery::DerivativeSwitch history({},true);auto reason=history.reason;
    CoupledLoadLineNewton solver(s,p.value,c,p);
    auto first=solver.solve({VectorXd::Zero(1),0.},.5,1.,&history);
    REQUIRE(first.converged);REQUIRE(boot==0);REQUIRE(first.derivativeHistory.has_value());
    REQUIRE(history.reason==reason);REQUIRE(history.mode()==derivative_recovery::DerivativeMode::Full);
    auto second=solver.solve(first.state,.6,1.,&*first.derivativeHistory);
    REQUIRE(second.converged);REQUIRE(boot==0);REQUIRE(full>0);
    s.automaticDerivativeRecovery=false;
    REQUIRE_THROWS_AS(CoupledLoadLineNewton(s,p.value,c,p).solve(first.state,.6,1.,&history),std::invalid_argument);
}


TEST_CASE("Full model quality reports do not select or alter candidates",
          "[coupled_load_line][model_report]") {
    ArclengthSystem system;
    system.residual=[](const VectorXd& x,Real){return VectorXd::Constant(1,x[0]*x[0]-2.);};
    system.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Zero(1);};
    system.jacobian=[](const VectorXd& x,Real){SparseMatrixd j(1,1);j.insert(0,0)=2.*x[0];return j;};
    system.solveJacobian=[](const VectorXd&,Real,const VectorXd&,VectorXd&){return false;};
    ArclengthScalarFunctional port;
    port.value=[](const VectorXd&,Real){return 0.;};
    port.linearize=[](const VectorXd&,Real){ArclengthScalarLinearization p;p.stateDerivative=VectorXd::Zero(1);return p;};
    CoupledLoadLineConfig config;config.maxIterations=1;config.linearSolver="direct_bordered";
    config.lineSearchMode="full_step";config.applyVoltageUpdateLimit=false;config.applyDeviceUpdateLimit=false;
    const auto initial=ArclengthState{VectorXd::Constant(1,.1),0.};
    const auto off=CoupledLoadLineNewton(system,port.value,config,port).solve(initial,0.,1.);
    config.reportModelQuality=true;
    const auto on=CoupledLoadLineNewton(system,port.value,config,port).solve(initial,0.,1.);
    REQUIRE(on.state.x[0]==off.state.x[0]);REQUIRE(on.state.lambda==off.state.lambda);
    REQUIRE(on.failureReason==off.failureReason);REQUIRE(on.iterations==off.iterations);
    REQUIRE(on.trace.front().modelReportValid);REQUIRE(on.trace.front().modelReportReliable);
    REQUIRE(on.trace.front().modelReportRatio<0.);
    REQUIRE_FALSE(on.trace.front().candidateRejectedByRecovery);
    REQUIRE(on.trace.front().acceptedAlpha==1.);
}

TEST_CASE("Full model report uses the coupled circuit prediction",
          "[coupled_load_line][model_report]") {
    ArclengthSystem system;
    system.residual=[](const VectorXd& x,Real v){return VectorXd::Constant(1,x[0]-v);};
    system.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Constant(1,-1.);};
    system.jacobian=[](const VectorXd&,Real){SparseMatrixd j(1,1);j.insert(0,0)=1.;return j;};
    system.solveJacobian=[](const VectorXd&,Real,const VectorXd&,VectorXd&){return false;};
    ArclengthScalarFunctional port;port.value=[](const VectorXd&,Real){return 0.;};
    port.linearize=[](const VectorXd&,Real){ArclengthScalarLinearization p;p.stateDerivative=VectorXd::Zero(1);return p;};
    system.makeTrial=[](const ArclengthState& x,const VectorXd& dx,Real dv,Real alpha){
        ArclengthTrial t;t.state={x.x+alpha*dx,x.lambda+alpha*dv};return t;};
    CoupledLoadLineConfig config;config.maxIterations=1;config.linearSolver="direct_bordered";
    config.lineSearchMode="full_step";config.applyVoltageUpdateLimit=false;config.reportModelQuality=true;
    const auto result=CoupledLoadLineNewton(system,port.value,config,port).solve({VectorXd::Zero(1),0.},.5,1.);
    REQUIRE(result.converged);REQUIRE(result.trace.front().modelReportReliable);
    REQUIRE(result.trace.front().modelReportRatio==Catch::Approx(1.));
    REQUIRE(result.trace.front().modelReportPredicted==Catch::Approx(0.));
    REQUIRE(result.trace.front().modelReportBefore>0.);
    REQUIRE(result.trace.front().linearGateAbsolute<=result.trace.front().linearGateAbsoluteTolerance);
    REQUIRE(result.trace.front().linearGateBackward<=result.trace.front().linearGateBackwardTolerance);
    REQUIRE(result.trace.front().linearGateBackwardTolerance==1e-10);
}
