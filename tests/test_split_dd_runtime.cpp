#include <catch2/catch_test_macros.hpp>
#include "vela/equation/SplitDDRuntime.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/post/ContactCurrent.h"
#include "vela/io/DDSolutionCsv.h"
#include "vela/numerics/LineSearch.h"
#include "vela/numerics/SplitCoordinate.h"
#include "vela/numerics/SplitDDState.h"
#include "vela/simulation/DCSweepPredictor.h"
#include "vela/core/PhysicalConstants.h"
#include <nlohmann/json.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
using namespace vela;
namespace {
nlohmann::json fixture() {
    nlohmann::json d={{"schema","vela.split-dd-operator.v1"},{"mesh_fingerprint","fixture"},
      {"temperature_K",300.},{"potential_scale",.02585},{"gauge_policy","encoded_reference"},
      {"phumob",{{"electron",{.1417,.00522,2.285,9.68e22,.68}},{"hole",{.04705,.00449,2.247,2.23e23,.719}},
       {"common",{4e26,7.2e26,.21,.5,1.,1.258,2.459,3.828,1.,1.,.89233,.41372,.005978,.28227,.72169,.19778,1.80618,1e-6,1e4}}}},
      {"poisson_edges",{{0,1,2e-10}}},{"current_factor",1.602176634e-19},
      {"boundary",nlohmann::json::array({nlohmann::json::array(),nlohmann::json::array(),nlohmann::json::array()})}};
    for(int i=0;i<2;++i){d["nodes"].push_back({0.,0.,0.,0.,0.,0.,.02585,0.,0.,1e16,.02585,1e20,1e-21,1e-21,1e-21,1.602176634e-19,1.,1.,0.});d["node_physics"].push_back({1e20,0.,1e-5,3e-6,1e-21});}
    d["transport_edges"].push_back({{"i",0},{"j",1},{"coefficient",.01},{"weights",{{0,.5},{1,.5}}},{"ports",{{"drain",1.},{"source",-1.}}}});return d;
}
struct Scratch {
    std::filesystem::path path;
    Scratch(){auto base=std::filesystem::temp_directory_path()/"vela_split_runtime";for(unsigned i=0;;++i){path=base.string()+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+"_"+std::to_string(i);if(std::filesystem::create_directory(path))break;}}
    ~Scratch(){std::error_code ec;std::filesystem::remove_all(path,ec);}
};
}
TEST_CASE("Split runtime checkpoint preserves sub ULP physics and refuses stale fields") {
    SplitDDRuntime runtime(fixture());VectorXd x(6);x<<1.,2.,.125,.125,.125,.125;
    auto low=runtime.low();VectorXd step=VectorXd::Zero(6);step(0)=0x1p-60;step(2)=-0x1p-60;step(4)=0x1p-60;
    auto r0=runtime.residual(x,{});auto next=runtime.candidate(x,step,1.,low);auto r1=runtime.residual(next,{});
    REQUIRE((next.array()==x.array()).all());REQUIRE(!runtime.low().isZero(0.));REQUIRE((r0.array()!=r1.array()).any());
    DDSolution state;runtime.save(next,state);REQUIRE(state.hasConsistentPackedState());
    Scratch scratch;auto csv=scratch.path/"state.csv";writeDDSolutionStateCsv(csv,state);auto loaded=readDDSolutionStateCsv(csv,2);
    REQUIRE((loaded.packedLow.array()==state.packedLow.array()).all());
    SplitDDRuntime fresh(fixture());auto restored=fresh.restore(loaded);REQUIRE((fresh.residual(restored,{}).array()==r1.array()).all());
    auto stale=loaded;stale.psi(0)=std::nextafter(stale.psi(0),INFINITY);REQUIRE_FALSE(stale.hasConsistentPackedState());REQUIRE_THROWS(writeDDSolutionStateCsv(scratch.path/"bad.csv",stale));
    auto wrong=loaded;wrong.packedMeshFingerprint="different";REQUIRE_THROWS_AS(fresh.restore(wrong),std::invalid_argument);
    auto model=fixture();model["nodes"][0][7]=1.;SplitDDRuntime wrongFrame(model);REQUIRE_THROWS_AS(wrongFrame.restore(loaded),std::invalid_argument);
}
TEST_CASE("Line search candidate builder retains a decrease hidden below coordinate ULP") {
    BacktrackingLineSearch search;VectorXd x(1),step(1),r(1);x<<1.;step<<5e-19;r<<-1.;double low=0;
    auto result=search.search(x,step,r,[&](const VectorXd&){VectorXd f(1);f<<-1.+1e18*low;return f;},{},{},{},
        [&](const VectorXd& base,const VectorXd& delta,Real alpha){auto p=numerics::SplitCoordinate::sum(base(0),alpha*delta(0));low=p.lo;VectorXd z(1);z<<p.hi;return z;});
    REQUIRE(result.accepted);REQUIRE(result.x(0)==1.);REQUIRE(result.residual(0)==-.5);REQUIRE(result.damping==1.);
    auto ordinary=search.search(x,step,r,[](const VectorXd&){VectorXd f(1);f<<-1.;return f;});REQUIRE_FALSE(ordinary.accepted);
}

TEST_CASE("Split runtime retains the original encoded inactive carrier constraints") {
    auto model=fixture();model["transport_edges"][0]["weights"]=nlohmann::json::array();
    for(int i=0;i<2;++i){model["nodes"][i][9]=0.;model["nodes"][i][7]=.05;model["nodes"][i][8]=1.;}
    SplitDDRuntime runtime(model);VectorXd x=VectorXd::Zero(6);x.segment(2,2).setConstant(-(.05/.02585));x.segment(4,2).setConstant(-(1./.02585));
    REQUIRE(runtime.residual(x,{}).tail(4).isZero(0.));
}

TEST_CASE("Split runtime ports and diagnostic rows share conservative fluxes") {
    auto model=fixture();model["continuity_scale"]=1.;model["field_factor"]=1.;
    auto& e=model["transport_edges"][0];e["id"]=0;e["length"]=1.;e["couple"]=1.;e["length_m"]=1.;e["couple_m"]=1.;
    SplitDDRuntime runtime(model);VectorXd x(6);x<<1.,2.,.125,.25,.2,.1;
    auto r=runtime.residual(x,{});auto terms=runtime.terms(x,{});auto edges=runtime.edges(x);
    REQUIRE(edges.size()==1);REQUIRE(terms[0].electronFlux==-terms[1].electronFlux);
    REQUIRE(terms[0].electronFlux==edges[0].electronFlux);
    REQUIRE(terms[0].holeFlux==edges[0].holeFlux);
    auto drain=runtime.contact(x,"drain"),source=runtime.contact(x,"source");
    REQUIRE(drain.totals.totalCurrent==-source.totals.totalCurrent);
    REQUIRE(std::abs(drain.totals.totalCurrent*1e-6/runtime.current(x,"drain")-1)<1e-14);
    CoupledDDBoundaryConditions bc;bc.phin[0]=x(2);
    auto constrained=runtime.terms(x,bc);REQUIRE_FALSE(constrained[0].electronContinuityActive);
    REQUIRE(constrained[0].electronFlux==0.);REQUIRE(constrained[0].electronRecombination==0.);
}

TEST_CASE("Split symmetric residual difference resolves increments below projected residual ULP") {
    SplitDDRuntime runtime(fixture());VectorXd x(6);x<<1.,2.,.125,.25,.2,.1;
    VectorXd dx=VectorXd::Zero(6);dx(0)=1e-30;VectorXd a,b;
    const auto fd=runtime.symmetricDifference(x,dx,{},a,b);
    REQUIRE(a(0)==b(0));REQUIRE(fd(0)!=0.);REQUIRE(fd.allFinite());
    REQUIRE(runtime.low().isZero(0.));
}

TEST_CASE("Split current decomposition has no diffusion for uniform populations") {
    auto model=fixture();model["continuity_scale"]=1.;model["field_factor"]=1.;
    auto& edge=model["transport_edges"][0];edge["id"]=0;edge["length"]=1.;edge["couple"]=1.;edge["length_m"]=1.;edge["couple_m"]=1.;
    SplitDDRuntime runtime(model);VectorXd x(6);x<<1.,2.,1.,2.,1.,2.;
    const auto current=runtime.contact(x,"drain").totals;
    REQUIRE(current.electronCurrent<0.);REQUIRE(current.holeCurrent>0.);
    REQUIRE(std::abs(current.electronDiffusionCurrent/current.electronCurrent)<1e-80);
    REQUIRE(std::abs(current.holeDiffusionCurrent/current.holeCurrent)<1e-80);
}

namespace {
split_dd::Wide exactPotential(const DDSolution& s,int b,int i) {
    using split_dd::Wide;const int k=b*s.psi.size()+i;
    return (Wide(s.packedState(k))+Wide(s.packedLow(k)))*Wide(s.packedPotentialScale_V)+
        (b==1?Wide(s.electronQfReference(i)):b==2?Wide(s.holeQfReference(i)):Wide(0));
}
}
TEST_CASE("Split continuation rebases references within pair rounding and keeps strict restore", "[split_continuation]") {
    using split_dd::Wide;
    auto m=fixture();SplitDDRuntime original(m);VectorXd x(6);x<<1.,2.,1e-24,2e-24,-1e-24,-2e-24;
    VectorXd lo=VectorXd::Constant(6,1e-42);lo(0)=0x1p-60;original.setLow(lo);
    DDSolution s;original.save(x,s);
    auto same=original.restoreForContinuation(s);
    REQUIRE((same.array()==x.array()).all());REQUIRE((original.low().array()==lo.array()).all());
    auto changed=m;changed["nodes"][0][7]=.1;changed["nodes"][1][8]=-.05;
    SplitDDRuntime next(changed);REQUIRE_THROWS(next.restore(s));
    auto y=next.restoreForContinuation(s);DDSolution t;next.save(y,t);
    const Wide eps=std::numeric_limits<double>::epsilon();
    for(int b=0;b<3;++b)for(int i=0;i<2;++i) {
        const Wide bound=8*eps*eps*(1+abs(Wide(t.packedState(b*2+i))))*Wide(s.packedPotentialScale_V);
        REQUIRE(abs(exactPotential(t,b,i)-exactPotential(s,b,i))<bound);
    }
    REQUIRE(t.packedLow(0)==s.packedLow(0));
    REQUIRE(t.packedState(3)==s.packedState(3));REQUIRE(t.packedLow(3)==s.packedLow(3));
    REQUIRE(t.packedState(4)==s.packedState(4));REQUIRE(t.packedLow(4)==s.packedLow(4));
    REQUIRE(std::abs(next.current(y,"drain")/original.current(x,"drain")-1)<1e-7);
    auto back=original.restoreForContinuation(t);DDSolution roundtrip;original.save(back,roundtrip);
    for(int b=0;b<3;++b)for(int i=0;i<2;++i)
        REQUIRE(abs(exactPotential(roundtrip,b,i)-exactPotential(s,b,i))<Wide("1e-32"));
    const auto before=next.low();auto wrong=s;wrong.packedMeshFingerprint="wrong";
    REQUIRE_THROWS(next.restoreForContinuation(wrong));REQUIRE((next.low().array()==before.array()).all());
    auto wrongScale=changed;wrongScale["potential_scale"]=1.;SplitDDRuntime scale(wrongScale);
    REQUIRE_THROWS(scale.restoreForContinuation(s));
    wrong=s;wrong.psi(0)+=.1;REQUIRE_THROWS(next.restoreForContinuation(wrong));
    CoupledDDBoundaryConditions bc;bc.psi[0]=.25;bc.phin[0]=.1/.02585;bc.phip[1]=-.05/.02585;
    const auto untouched=y(1);const auto untouchedLow=next.low()(1);
    next.applyBoundary(y,bc);const auto r=next.residual(y,bc);
    REQUIRE(r(0)==0.);REQUIRE(r(2)==0.);REQUIRE(r(5)==0.);
    REQUIRE(next.low()(0)==0.);REQUIRE(next.low()(2)==0.);REQUIRE(next.low()(5)==0.);
    REQUIRE(y(1)==untouched);REQUIRE(next.low()(1)==untouchedLow);
}

TEST_CASE("Split sweep predictor updates full state across references and retains retry state", "[split_continuation]") {
    using split_dd::Wide;auto m=fixture();for(int i=0;i<2;++i)m["nodes"][i][10]=constants::Vt_300;
    SplitDDRuntime first(m);VectorXd x(6);x<<1.,2.,.1,.2,.3,.4;DDSolution a;first.save(x,a);
    m["nodes"][0][7]=.1;m["nodes"][1][8]=-.05;SplitDDRuntime second(m);
    auto y=second.restoreForContinuation(a);VectorXd step=VectorXd::Constant(6,1e-4);step(0)=1e-22;
    y=second.candidate(y,step,1.,second.low());DDSolution b;second.save(y,b);
    SweepPredictorConfig cfg;cfg.mode="linear";cfg.fields={"psi","phin","phip"};cfg.maxExtrapolationRatio=3.;
    const auto predicted=detail::predictDCSweepInitialState(cfg,&a,b,0.,.1,.25);
    REQUIRE(predicted.hasConsistentSplitPackedState());
    const double ratio=(.25-.1)/.1;
    for(int block=0;block<3;++block)for(int i=0;i<2;++i) {
        const Wide target=exactPotential(b,block,i)+Wide(ratio)*(exactPotential(b,block,i)-exactPotential(a,block,i));
        REQUIRE(abs(exactPotential(predicted,block,i)-target)<Wide("1e-32"));
    }
    REQUIRE(predicted.packedState(0)==b.packedState(0));
    REQUIRE(predicted.packedLow(0)!=b.packedLow(0));
    auto z=second.restore(predicted);DDSolution recomputed;second.save(z,recomputed);
    REQUIRE((predicted.n-recomputed.n).norm()/recomputed.n.norm()<1e-14);
    REQUIRE((predicted.p-recomputed.p).norm()/recomputed.p.norm()<1e-14);
    auto retry=detail::predictDCSweepInitialState(cfg,&a,b,0.,.1,.25,1);
    REQUIRE((retry.packedState.array()==b.packedState.array()).all());
    REQUIRE((retry.packedLow.array()==b.packedLow.array()).all());
    cfg.fields={"psi"};auto partial=detail::predictDCSweepInitialState(cfg,&a,b,0.,.1,.25);
    REQUIRE((partial.packedState.tail(4).array()==b.packedState.tail(4).array()).all());
    REQUIRE((partial.packedLow.tail(4).array()==b.packedLow.tail(4).array()).all());
    auto invalid=a;invalid.packedLow.resize(0);
    REQUIRE_THROWS(detail::predictDCSweepInitialState(cfg,&invalid,b,0.,.1,.25));
}
