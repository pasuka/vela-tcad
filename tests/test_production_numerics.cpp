#include <catch2/catch_test_macros.hpp>
#include "vela/solver/LinearRefinement.h"
#include <limits>

using namespace vela::linear_refinement;

TEST_CASE("High precision dot residual retains a weak term lost in double summation") {
    vela::SparseMatrixd J(1,2);J.insert(0,0)=1e20;J.insert(0,1)=-1e20;
    vela::VectorXd F(1),x(2);F<<1.;x<<1.,1.;
    auto d=highPrecisionDefect(J,F,x);
    REQUIRE(d.values(0)==1.);
    REQUIRE(d.backwardMax>0.);
    REQUIRE(d.backwardMax<1e-19);
}

TEST_CASE("Refinement recovers a missing weak equation update with original row weights") {
    vela::SparseMatrixd J(3,3);J.insert(0,0)=2;J.insert(1,1)=1e-30;J.insert(2,2)=3;
    vela::VectorXd F(3),step(3),w(3);F<<-2.,-2.5e-47,-6.;step<<1.,0.,2.;w<<1.,1e12,1.;
    vela::LinearSolver solver;int observations=0;
    auto result=refine(J,F,w,step,solver,4,[&](int k,const auto&,const Defect& d){
        ++observations;if(k==0)REQUIRE(d.values(1)==F(1));
    });
    REQUIRE(observations==5);
    REQUIRE(std::abs(result(1)/2.5e-17-1)<1e-14);
    REQUIRE(result(0)==1.);REQUIRE(result(2)==2.);
    REQUIRE(std::abs(highPrecisionDefect(J,F,result).values(1)/F(1))<1e-14);
}

TEST_CASE("An already exact solve is invariant under fixed refinement") {
    vela::SparseMatrixd J(2,2);J.insert(0,0)=4.;J.insert(0,1)=-2.;J.insert(1,0)=1.;J.insert(1,1)=2.;
    vela::VectorXd F(2),step(2),w(2);step<<.5,.25;F<<-1.5,-1.;w<<8.,.125;
    vela::LinearSolver solver;
    auto result=refine(J,F,w,step,solver,4,[](int,const auto&,const Defect& d){REQUIRE(d.values.norm()==0.);});
    REQUIRE((result-step).norm()==0.);
}

TEST_CASE("Zero corrections preserves the supplied raw direction") {
    vela::SparseMatrixd J(2,2);J.setIdentity();
    vela::VectorXd F=vela::VectorXd::Ones(2),step=vela::VectorXd::Zero(2),w=F;
    vela::LinearSolver solver;int count=0;
    auto result=refine(J,F,w,step,solver,0,[&](int,const auto&,const auto&){++count;});
    REQUIRE(count==1);REQUIRE(result.norm()==0.);
}

#include <catch2/catch_test_macros.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "vela/discretization/StableSGDerivative.h"
#include "vela/discretization/ScharfetterGummel.h"
using HP=boost::multiprecision::cpp_dec_float_100;

namespace {
HP B(HP x) { if(x==0)return HP(1);return x/expm1(x); }
HP flux(HP a,HP b,HP q0,HP q1,bool electron) {
    const HP vt(".02585"),ni("1e10"),coef("2.4");const HP eta=(b-a)/vt;
    if(electron)return coef*ni*(exp((a-q0)/vt)*B(-eta)-exp((b-q1)/vt)*B(eta));
    return coef*ni*(exp((q0-a)/vt)*B(eta)-exp((q1-b)/vt)*B(-eta));
}
}

TEST_CASE("stable SG psi derivatives agree with independent 100-digit density-flux differences", "[simplemos][stable_sg]") {
    for(bool electron:{false,true})for(double potential:{-.43,.13})for(double delta:{0.,1e-18,1e-10,.03}) {
        const double a=potential,b=potential+.011,vt=.02585;const double eta=(b-a)/vt;
        const double f=electron?vela::sgElectronBoltzmannContinuityFlux(1e10,1e10,a,b,0.,delta,vt,2.4,{false,true})
                               :vela::sgHoleBoltzmannContinuityFlux(1e10,1e10,a,b,0.,delta,vt,2.4,{false,true});
        const HP x=electron?-HP(eta):HP(eta);const HP em1=expm1(x);const double db=static_cast<double>((em1-x*(em1+1))/(em1*em1));
        const auto derivative=vela::stable_sg::psiDerivative(f,vela::bernoulli(static_cast<double>(x)),db,vt,electron);
        const HP step("1e-20");
        for(int end=0;end<2;++end) {
            const HP plus=flux(HP(a)+(end==0?step:HP(0)),HP(b)+(end==1?step:HP(0)),HP(0),HP(delta),electron);
            const HP minus=flux(HP(a)-(end==0?step:HP(0)),HP(b)-(end==1?step:HP(0)),HP(0),HP(delta),electron);
            const HP ref=(plus-minus)/(2*step);
            if(delta==0.)REQUIRE(derivative[end]==0.);
            else REQUIRE(static_cast<double>(abs((HP(derivative[end])-ref)/ref))<2e-12);
        }
    }
}
