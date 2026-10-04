#include <catch2/catch_test_macros.hpp>
#include "vela/solver/LinearRefinement.h"
#include "vela/numerics/StableMeritComparison.h"
#include "vela/equation/ExtendedPoissonResidual.h"
#include "vela/numerics/SplitCoordinate.h"
#include "vela/solver/NewtonSolver.h"
#include <limits>
#include <nlohmann/json.hpp>

using namespace vela::linear_refinement;

TEST_CASE("Split coordinates accumulate reversible updates below one ULP") {
    using vela::numerics::SplitCoordinate;
    using Q=vela::poisson_precision::Q;
    const double base=-0x1.e53059760fd30p+3, step=0x1p-60;
    REQUIRE(base+step==base);
    SplitCoordinate x{base,0.};
    for(int i=0;i<1024;++i)x=x.shifted(step);
    REQUIRE(Q(x.hi)+Q(x.lo)==Q(base)+1024*Q(step));
    for(int i=0;i<1024;++i)x=x.shifted(-step);
    REQUIRE(x.hi==base); REQUIRE(x.lo==0.);
    REQUIRE_THROWS_AS(x.shifted(std::numeric_limits<double>::infinity()),std::invalid_argument);
}

TEST_CASE("Extended Poisson consumes split potential in flux and Boltzmann charge") {
    namespace pp=vela::poisson_precision;
    using vela::numerics::SplitCoordinate;
    pp::Node n{}; n.xpsi=1.;n.potentialScale=1.;n.vt=1.;n.ni=1.;
    n.q=1.;n.area=1.;n.scale=1.;n.voln=1.;n.volp=1.;
    auto other=n;other.xpsi=0.;std::vector<pp::Node> nodes{n,other};
    std::vector<pp::Edge> edges{{0,1,2.}};
    const double low=0x1p-55;
    const auto base=pp::evaluate(nodes,edges,true);
    const auto shifted=pp::evaluate(nodes,edges,true,{low,0.});
    const pp::Q t=low;
    const pp::Q expected=2*t+exp(pp::Q(1))*expm1(t)-exp(pp::Q(-1))*expm1(-t);
    REQUIRE(abs((shifted[0]-base[0])/expected-1)<pp::Q("1e-16"));
    REQUIRE(shifted[1]-base[1]==-2*t);
    // Same exact value with a different pair partition must give the same result.
    nodes[0].xpsi=std::nextafter(1.,2.);
    const auto repartitioned=pp::evaluate(nodes,edges,true,{low-(nodes[0].xpsi-1.),0.});
    REQUIRE(repartitioned==shifted);
    REQUIRE_THROWS_AS(pp::evaluate(nodes,edges,true,{low}),std::invalid_argument);
    REQUIRE_THROWS_AS(pp::evaluate(nodes,edges,false,{low,0.}),std::invalid_argument);
}

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

namespace {
// Independent, unfactorized density-flux reference. All intermediates,
// including eta and carrier populations, stay at 100 decimal digits.
HP variableNiDensityFlux(HP ni0, HP ni1, HP a, HP b, HP q0, HP q1,
                         HP vt, bool electron)
{
    const HP eta = (b-a)/vt + (electron ? log(ni1/ni0) : log(ni0/ni1));
    if (electron)
        return ni0*exp((a-q0)/vt)*B(-eta) - ni1*exp((b-q1)/vt)*B(eta);
    return ni0*exp((q0-a)/vt)*B(eta) - ni1*exp((q1-b)/vt)*B(-eta);
}

double bernoulliPrimeReference(double value)
{
    const HP x(value);
    if (x==0) return -.5;
    const HP e=expm1(x);
    return static_cast<double>((e-x*(e+1))/(e*e));
}
}

TEST_CASE("BGN SG potential derivatives preserve equilibrium and agree with 100-digit differences",
          "[sg][bgn][stable_sg]")
{
    const double vt = .025851999786435;
    for (bool electron : {false, true})
    for (double a : {-.43, .43})
    for (double difference : {-.6, 0., .011, .6})
    for (double ratio : {.001, 3., 1000.})
    for (double dqf : {0., -1e-18, 1e-18, -1e-10, .03}) {
        const double ni0 = 1e16, ni1 = ni0*ratio, b = a+difference;
        const double eta = (b-a)/vt + (electron ? std::log(ni1/ni0) : std::log(ni0/ni1));
        const double x = electron ? -eta : eta;
        const auto fluxFunction = electron ? vela::sgElectronContinuityFluxFromQuasiFermiVariableNi
                                           : vela::sgHoleContinuityFluxFromQuasiFermiVariableNi;
        const double f = fluxFunction(ni0,ni1,a,b,0.,dqf,vt,1.,true);
        const auto derivative = vela::stable_sg::psiDerivative(
            f,vela::bernoulli(x),bernoulliPrimeReference(x),vt,electron);
        const HP step("1e-25");
        for (int end=0; end<2; ++end) {
            const HP da = end==0 ? step : HP(0), db = end==1 ? step : HP(0);
            const HP plus = variableNiDensityFlux(HP(ni0),HP(ni1),HP(a)+da,HP(b)+db,HP(0),HP(dqf),HP(vt),electron);
            const HP minus = variableNiDensityFlux(HP(ni0),HP(ni1),HP(a)-da,HP(b)-db,HP(0),HP(dqf),HP(vt),electron);
            const HP reference = (plus-minus)/(2*step);
            CAPTURE(electron,a,difference,ratio,dqf,end);
            if (dqf==0.) {
                REQUIRE(derivative[end]==0.);
                REQUIRE(abs(reference)<HP("1e-35"));
            } else {
                REQUIRE(static_cast<double>(abs((HP(derivative[end])-reference)/reference))<2e-11);
            }
        }
        // Reversing the edge changes the sign and swaps endpoint derivatives.
        const double reverseFlux = fluxFunction(ni1,ni0,b,a,dqf,0.,vt,1.,true);
        const auto reverse = vela::stable_sg::psiDerivative(
            reverseFlux,vela::bernoulli(-x),bernoulliPrimeReference(-x),vt,electron);
        for (int end=0; end<2; ++end) {
            const double scale = std::max(std::abs(derivative[end]),std::abs(reverse[1-end]));
            REQUIRE(std::abs(derivative[end]+reverse[1-end]) <= 2e-11*scale);
        }
    }
}

TEST_CASE("Stable merit resolves physical decreases below the dominant residual norm") {
    vela::VectorXd base(2),trial(2);
    base << 1e-10,1e-30; trial << 1e-10,0.;
    REQUIRE(base.norm()==trial.norm());
    REQUIRE(vela::stable_merit::compare(base,trial).accepted);
    REQUIRE_FALSE(vela::stable_merit::compare(trial,base).accepted);
    REQUIRE_FALSE(vela::stable_merit::compare(base,base).accepted);
    base.setZero();trial.setZero();
    REQUIRE(vela::stable_merit::compare(base,trial).accepted);
}

TEST_CASE("Stable merit handles subnormal energy and rejects invalid residuals") {
    vela::VectorXd base(2),trial(2);
    base << std::numeric_limits<double>::denorm_min(),0.;trial.setZero();
    REQUIRE(base.squaredNorm()==0.);
    REQUIRE(vela::stable_merit::compare(base,trial).accepted);
    REQUIRE_FALSE(vela::stable_merit::compare(trial,base).accepted);
    trial(0)=std::numeric_limits<double>::infinity();
    REQUIRE_FALSE(vela::stable_merit::compare(base,trial).accepted);
    trial(0)=std::numeric_limits<double>::quiet_NaN();
    REQUIRE_FALSE(vela::stable_merit::compare(base,trial).accepted);
    REQUIRE_THROWS_AS(vela::stable_merit::compare(base,vela::VectorXd::Zero(3)),std::invalid_argument);
}

TEST_CASE("Extended Poisson preserves the charge response of a tiny referenced increment") {
    namespace pp=vela::poisson_precision;
    pp::Node node{};
    node.psi=1.;node.n=1.;node.p=1.;node.xpsi=10.;node.potentialScale=.1;
    node.eref=1.;node.href=1.;node.ni=1.;node.vt=1.;
    node.voln=1.;node.volp=1.;node.q=1.;node.area=1.;node.scale=1.;
    const auto before=pp::evaluate({node},{},true)[0];
    const auto legacy=pp::evaluate({node},{},false)[0];
    const double h=1e-18;
    node.xn=h;
    REQUIRE(node.eref+node.xn*node.potentialScale==node.eref);
    REQUIRE(pp::evaluate({node},{},false)[0]==legacy);
    const auto change=pp::evaluate({node},{},true)[0]-before;
    const HP exactPsi=HP(10.)*HP(.1);
    const HP expected=exp(exactPsi-HP(1.))*expm1(-HP(h)*HP(.1));
    REQUIRE(change<0);
    REQUIRE(std::abs(change.convert_to<double>()/expected.convert_to<double>()-1)<1e-12);
}

TEST_CASE("Extended Poisson retains signed edge geometry and pair conservation") {
    namespace pp=vela::poisson_precision;
    pp::Node left{},right{};
    left.psi=1.;right.psi=3.;left.scale=right.scale=1.;
    auto residual=pp::evaluate({left,right},{{0,1,-2.5}},false);
    REQUIRE(residual[0]==5.);
    REQUIRE(residual[1]==-5.);
    REQUIRE(residual[0]+residual[1]==0.);
}

TEST_CASE("Numerical precision options are explicit and reject unsupported combinations") {
    const auto defaults=vela::newtonConfigFromJson(nlohmann::json::object(),{});
    REQUIRE(defaults.poissonResidualPrecision=="double");
    REQUIRE_FALSE(defaults.stableMeritComparison);
    REQUIRE_FALSE(defaults.exactDirichletUpdates);
    const nlohmann::json valid={{"poisson_residual_precision","binary128"},
        {"stable_merit_comparison",true},{"exact_dirichlet_updates",true},
        {"residual_norm","l2"},{"line_search_mode","merit"}};
    auto cfg=vela::newtonConfigFromJson(valid,{});
    REQUIRE(cfg.poissonResidualPrecision=="binary128");
    REQUIRE(cfg.stableMeritComparison);REQUIRE(cfg.exactDirichletUpdates);
    auto invalid=valid;invalid["poisson_residual_precision"]="unknown";
    REQUIRE_THROWS_AS(vela::newtonConfigFromJson(invalid,{}),std::invalid_argument);
    invalid=valid;invalid["carrier_statistics"]="fermi_dirac";
    REQUIRE_THROWS_AS(vela::newtonConfigFromJson(invalid,{}),std::invalid_argument);
    invalid=valid;invalid["electron_quantum_potential"]=true;
    REQUIRE_THROWS_AS(vela::newtonConfigFromJson(invalid,{}),std::invalid_argument);
    auto block=valid;block["residual_norm"]="block";
    REQUIRE_NOTHROW(vela::newtonConfigFromJson(block,{}));
    block["residual_weights"]={{"psi",2.}};
    block["residual_scales"]={{"psi",2.}};
    REQUIRE_NOTHROW(vela::newtonConfigFromJson(block,{}));
    invalid=block;invalid["global_continuity_closure"]={{"mode","enforce"}};
    REQUIRE_THROWS_AS(vela::newtonConfigFromJson(invalid,{}),std::invalid_argument);
    invalid=valid;invalid["carrier_regularization_scale"]=1e-3;
    REQUIRE_THROWS_AS(vela::newtonConfigFromJson(invalid,{}),std::invalid_argument);
}

TEST_CASE("Stable block merit preserves weights and resolves exact cancellation") {
    vela::VectorXd x(3),y(3);x<<1.,0.,1e-150;y<<0.,2.,0.;
    std::array<double,3> scales{1.,4.,2.},weights{.25,1.,3.};
    auto c=vela::stable_merit::compareBlocks(x,y,scales,weights);
    REQUIRE(c.exactFallback);REQUIRE(c.sign==-1);REQUIRE(c.accepted);
    REQUIRE_FALSE(vela::stable_merit::compareBlocks(y,x,scales,weights).accepted);
    weights[2]=0.;REQUIRE_FALSE(vela::stable_merit::compareBlocks(x,y,scales,weights).accepted);
    x<<0.,0.,1.;y<<0.,0.,3.;weights={1.,1.,0.};
    REQUIRE(vela::stable_merit::compareBlocks(x,y,scales,weights).accepted);
    weights={.3,1.7,2.3};scales={13.1,.7,9.3};
    for(int k=1;k<=40;++k) {
        x<<k*.123,k*.77,k*.003;y<<k*.119,k*.779,k*.002;
        HP exact=0;
        for(int b=0;b<3;++b)exact+=(HP(y(b))*HP(y(b))-HP(x(b))*HP(x(b)))*HP(weights[b])/(HP(scales[b])*HP(scales[b]));
        const auto decision=vela::stable_merit::compareBlocks(x,y,scales,weights);
        REQUIRE(decision.sign==(exact<0?-1:exact>0?1:0));
        auto scaledX=x,scaledY=y;auto scaledS=scales;
        scaledX(1)*=8.;scaledY(1)*=8.;scaledS[1]*=8.;
        REQUIRE(vela::stable_merit::compareBlocks(scaledX,scaledY,scaledS,weights).sign==decision.sign);
    }
}
