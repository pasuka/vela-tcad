#include "vela/simulation/CarrierReducedLoadLine.h"
#include "vela/simulation/CommonRelativeLinearSolve.h"
#include "vela/simulation/CoupledLoadLine.h"
#include "vela/solver/PopulationState.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace vela;

TEST_CASE("Load-line observers see accepted trials before state promotion", "[coupled_load_line][observer]") {
    ArclengthSystem s;
    s.residual=[](const VectorXd& x,Real){return VectorXd::Constant(1,x[0]-2.);};
    s.parameterDerivative=[](const VectorXd&,Real){return VectorXd::Zero(1);};
    s.solveJacobian=[](const VectorXd&,Real,const VectorXd& b,VectorXd& y){y=b;return true;};
    const auto current=[](const VectorXd&,Real){return 0.;};
    for(const auto mode:{"merit","full_step","report_only"}) {
        CoupledLoadLineConfig c;c.lineSearchMode=mode;c.maxIterations=3;
        CoupledLoadLineNewton solver(s,current,c);
        const ArclengthState seed{VectorXd::Zero(1),0.};
        int calls=0;
        const auto observer=[&](const CoupledLoadLineTrial& t,const ArclengthState& candidate) {
            ++calls; REQUIRE(t.accepted); REQUIRE(candidate.x[0]==Catch::Approx(2.));
        };
        const auto result=solver.solve(seed,0.,1.,{},observer);
        REQUIRE(result.converged); REQUIRE(calls==1); REQUIRE(seed.x[0]==0.);
        REQUIRE_THROWS_AS(solver.solve(seed,0.,1.,{},
            [](const CoupledLoadLineTrial&,const ArclengthState&){throw std::runtime_error("export failed");}),
            std::runtime_error);
    }
}

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
