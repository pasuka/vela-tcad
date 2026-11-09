#include "vela/simulation/PopulationCoordinateNewton.h"
#include "vela/solver/PopulationState.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
using namespace vela;

namespace {
ArclengthSystem linearPopulationSystem() {
    ArclengthSystem s;
    s.residual=[](const VectorXd& x,Real v){VectorXd f(3);f<<x[0]-v,x[1]/1e20-2.,x[2]/.01-3.;return f;};
    s.parameterDerivative=[](const VectorXd&,Real){VectorXd f(3);f<<-1.,0.,0.;return f;};
    s.jacobian=[](const VectorXd&,Real){SparseMatrixd j(3,3);j.insert(0,0)=1.;j.insert(1,1)=1e-20;j.insert(2,2)=100.;return j;};
    s.solveJacobian=[](const VectorXd&,Real,const VectorXd& b,VectorXd& y){y=b;y[1]*=1e20;y[2]/=100.;return true;};
    s.columnScales=[](const ArclengthState& x){VectorXd d(4);d<<std::abs(x.x[0])+.025,std::abs(x.x[1]),std::abs(x.x[2]),std::abs(x.lambda)+.025;return d;};
    s.prescribedRows=[](const VectorXd&,Real){return std::vector<int>{};};
    s.makeTrial=[](const ArclengthState& before,const VectorXd& dx,Real dv,Real alpha){
        const auto p=makePopulationTrial(before.x,dx,alpha,{true},{true},PopulationUpdateConfig{});
        return ArclengthTrial{{p.state,before.lambda+p.alpha*dv},p.alpha,
            static_cast<int>(p.protectedElectrons.size()),static_cast<int>(p.protectedHoles.size())};
    };
    return s;
}
ArclengthScalarFunctional linearPort() {
    ArclengthScalarFunctional p;
    p.value=[](const VectorXd& x,Real v){return 1e-7*x[0]+1e-9*x[1]/1e20+1e-9*x[2]/.01+2e-8*v;};
    p.linearize=[](const VectorXd&,Real){ArclengthScalarLinearization a;a.stateDerivative.resize(3);a.stateDerivative<<1e-7,1e-29,1e-7;a.parameterDerivative=2e-8;return a;};
    return p;
}
PrimaryVariableScales refs(){VectorXd d(3);d<<.025,1e20,1e20;return {d,.025,3e-5};}
}

TEST_CASE("Unit primary coordinates preserve physical equations and circuit chain rule",
          "[coupled_load_line][primary_coordinates]") {
    const auto s=linearPopulationSystem();const auto p=linearPort();
    PrimaryVariableCoordinates b(refs());const auto mapped=b.system(s);const auto port=b.port(p);
    ArclengthState x{VectorXd(3),2.};x.x<<1.,1e20,.01;
    const auto u=b.scaled(x);const auto restored=b.physical(u);
    REQUIRE((mapped.residual(u.x,u.lambda)-s.residual(restored.x,restored.lambda)).norm()==0.);
    const auto J=mapped.jacobian(u.x,u.lambda);VectorXd du(3);du<<.5,.25,1e-22;
    const VectorXd expected=s.jacobian(restored.x,restored.lambda)*b.physical(du);
    REQUIRE((J*du-expected).norm()<1e-14);
    REQUIRE((mapped.parameterDerivative(u.x,u.lambda)-s.parameterDerivative(restored.x,restored.lambda)*.025).norm()==0.);
    const auto a=p.linearize(restored.x,restored.lambda),z=port.linearize(u.x,u.lambda);
    REQUIRE((z.stateDerivative.cwiseProduct(du)*3e-5-a.stateDerivative.cwiseProduct(b.physical(du))).norm()<1e-25);
    REQUIRE(z.parameterDerivative*3e-5==Catch::Approx(a.parameterDerivative*.025));
    VectorXd raw(3);raw<<0.,-2e20,-.02;
    const auto t=s.makeTrial(restored,raw,0.,1.);
    const auto mt=mapped.makeTrial(u,b.scaled(raw),0.,1.);
    REQUIRE(mt.protectedElectrons==t.protectedElectrons);
    REQUIRE(mt.protectedHoles==t.protectedHoles);
    const auto actual=b.physical(mt.state);
    REQUIRE(actual.x[1]==Catch::Approx(t.state.x[1]));
    REQUIRE(actual.x[2]==Catch::Approx(t.state.x[2]));
}

TEST_CASE("Unit population solve audits the original system and shares physical load RHS",
          "[coupled_load_line][primary_coordinates]") {
    const auto s=linearPopulationSystem();const auto p=linearPort();
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";
    c.applyDeviceUpdateLimit=false;c.applyVoltageUpdateLimit=false;c.voltageCoefficient=2.;
    c.maxIterations=4;c.equationTolerance=1e-10;c.minimumOriginalCorrections=1;
    ArclengthState seed{VectorXd(3),1.};seed.x<<1.,1e20,.01;
    const auto original=PopulationCoordinateNewton(s,p.value,c,p,{},refs(),false).solve(seed,10.,1e7);
    const auto unit=PopulationCoordinateNewton(s,p.value,c,p,{},refs(),true).solve(seed,10.,1e7);
    INFO(original.failureReason);INFO(unit.failureReason);
    REQUIRE(original.converged);REQUIRE(unit.converged);
    REQUIRE(unit.state.lambda==Catch::Approx((10.-.05)/3.2).margin(1e-12));
    REQUIRE(unit.state.lambda==Catch::Approx(original.state.lambda).margin(1e-12));
    REQUIRE(unit.state.x[1]/1e20==Catch::Approx(2.));
    REQUIRE(unit.state.x[2]==Catch::Approx(.03));
    const auto L=LoadLineResidual::evaluate(2.,unit.state.lambda,1e7,p.value(unit.state.x,unit.state.lambda),10.);
    REQUIRE(unit.loadLineResidual_V==L.voltage_V);
    int audits=0;
    for(const auto& t:unit.trace)if(std::isfinite(t.physicalLinearAbsolute)) {
        ++audits;REQUIRE(t.physicalLinearAbsolute<=1e-8);REQUIRE(t.physicalLinearBackward<=1e-10);
    }
    REQUIRE(audits>0);
    REQUIRE(seed.x[1]==1e20);
}

TEST_CASE("Unavailable native RHS mapping cannot authorize stopping without a valid E",
          "[coupled_load_line][primary_coordinates]") {
    auto s=linearPopulationSystem();auto p=linearPort();
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";c.maxIterations=1;
    c.equationTolerance=1e-30;
    ArclengthState seed{VectorXd(3),1.};seed.x<<1.,1e20,.01;
    DensityUpdateErrorEvaluator metric=[](const VectorXd&,const VectorXd&,const VectorXd&,Real,Real){DensityUpdateError e;e.valid=false;e.rms=0.;return e;};
    auto result=PopulationCoordinateNewton(s,p.value,c,p,metric,refs(),false,true,{}).solve(seed,10.,1e7);
    REQUIRE_FALSE(result.converged);
    metric=[](const VectorXd&,const VectorXd&,const VectorXd&,Real,Real){DensityUpdateError e;e.valid=true;e.rms=.5;return e;};
    result=PopulationCoordinateNewton(s,p.value,c,p,metric,refs(),false,true,{}).solve(seed,10.,1e7);
    REQUIRE(result.converged);REQUIRE(result.convergenceReason=="native_approx_update");
    REQUIRE(std::isnan(result.trace.back().mappedRhsNorm));
    REQUIRE_FALSE(ApproximateNativeRhsScales{}.available());
}

TEST_CASE("Unit stopping and physical audit share exact assemblies within one target only",
          "[coupled_load_line][primary_coordinates]") {
    auto s=linearPopulationSystem();const auto p=linearPort();
    int jacobians=0,residuals=0;
    const auto J=s.jacobian;const auto F=s.residual;
    s.jacobian=[&](const VectorXd& x,Real v){++jacobians;return J(x,v);};
    s.residual=[&](const VectorXd& x,Real v){++residuals;return F(x,v);};
    CoupledLoadLineConfig c;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";
    c.applyDeviceUpdateLimit=false;c.applyVoltageUpdateLimit=false;c.maxIterations=2;
    c.minimumOriginalCorrections=1;
    DensityUpdateErrorEvaluator metric=[](const VectorXd&,const VectorXd&,const VectorXd&,Real,Real){
        DensityUpdateError e;e.valid=true;e.rms=.5;return e;
    };
    ArclengthState seed{VectorXd(3),1.};seed.x<<1.,1e20,.01;
    PopulationCoordinateNewton solver(s,p.value,c,p,metric,refs(),true,true,{1.,1.,1.});
    for(double target:{10.,11.}) {
        const int beforeJ=jacobians,beforeF=residuals;
        const auto result=solver.solve(seed,target,1e7);
        INFO(result.failureReason);REQUIRE(result.converged);
        REQUIRE(result.iterations==1);
        // Initial and final physical borders plus the solver's separate direct
        // column assembly. The audit reuses the initial physical border/F.
        REQUIRE(jacobians-beforeJ==3);
        REQUIRE(residuals-beforeF==3);
        REQUIRE(result.trace.front().physicalLinearAbsolute<=1e-8);
        REQUIRE(result.trace.front().physicalLinearBackward<=1e-10);
        REQUIRE(std::abs(result.state.lambda+1e7*p.value(result.state.x,result.state.lambda)-target)<1e-12);
    }
}
