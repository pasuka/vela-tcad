#include <catch2/catch_test_macros.hpp>
#include "vela/equation/SplitDDRuntime.h"
#include "vela/equation/SplitDDOperator.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/post/ContactCurrent.h"
#include "vela/io/DDSolutionCsv.h"
#include "vela/numerics/LineSearch.h"
#include "vela/numerics/SplitCoordinate.h"
#include "vela/numerics/SplitDDState.h"
#include "vela/simulation/DCSweepPredictor.h"
#include "vela/core/PhysicalConstants.h"
#include "../src/equation/SplitMeshFingerprint.h"
#include <nlohmann/json.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
using namespace vela;

TEST_CASE("Split mesh fingerprint preserves SHA-1 block padding and checkpoint format") {
    using detail::splitMeshFingerprint;
    REQUIRE(splitMeshFingerprint("") == "sha1:da39a3ee5e6b4b0d3255bfef95601890afd80709");
    REQUIRE(splitMeshFingerprint(std::string(55, 'a')) == "sha1:c1c8bbdc22796e28c0e15163d20899b65621d65a");
    REQUIRE(splitMeshFingerprint(std::string(56, 'a')) == "sha1:c2db330f6083854c99d4b5bfb6e8f29f201be699");
    REQUIRE(splitMeshFingerprint(std::string(64, 'a')) == (BOOST_VERSION < 108600
        ? "sha1:98ba824b5c16427bd7a1122a5a442a25ec644d"
        : "sha1:0098ba824b5c16427bd7a1122a5a442a25ec644d"));
    REQUIRE(splitMeshFingerprint(std::string(65, 'a')) == "sha1:11655326c708d70319be2610e8a57d9a5b959d3b");
}

TEST_CASE("Split DD field norm retains the software precision over a wide exponent range") {
    using split_dd::Wide;
    REQUIRE(split_dd::fieldSquareRoot(Wide(0)) == 0);
    for (int exponent : {-3001, -45, 0, 73, 3001}) {
        const Wide squared = ldexp(Wide("0.79"), exponent);
        const auto root = split_dd::fieldSquareRoot(squared);
        REQUIRE(root > 0);
        REQUIRE(abs(root * root / squared - 1) < Wide("1e-97"));
    }
}

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

namespace {
nlohmann::json sgFixture() {
    auto model=fixture();model["continuity_scale"]=1.;model["field_factor"]=1.;
    auto& edge=model["transport_edges"][0];edge["id"]=0;edge["length"]=1.;edge["couple"]=1.;edge["length_m"]=1.;edge["couple_m"]=1.;
    // Deliberately asymmetric carrier prefactors and common ni: using common
    // ni in either SG weight would create a fictitious equilibrium current.
    model["population_prefactors"]={{2e16,7e15},{9e16,3e16}};
    return model;
}
}
TEST_CASE("Split SG diagnostics preserve equilibrium with carrier specific prefactors", "[split_sg_diagnostics]") {
    SplitDDRuntime runtime(sgFixture());VectorXd x(6);x<<1.,2.,0.,0.,0.,0.;
    const auto before=runtime.residual(x,{});const auto current=runtime.current(x,"drain");const auto low=runtime.low();
    const auto edge=runtime.edges(x).at(0);
    REQUIRE(edge.electronFlux==0.);REQUIRE(edge.holeFlux==0.);
    REQUIRE(std::abs(edge.electronSgReconstructedFluxScaled)/edge.electronSgHighPrecisionReferenceTermScaleScaled<1e-90);
    REQUIRE(std::abs(edge.holeSgReconstructedFluxScaled)/edge.holeSgHighPrecisionReferenceTermScaleScaled<1e-90);
    REQUIRE(std::abs(edge.electronSgHighPrecisionReferenceFluxScaled)/edge.electronSgHighPrecisionReferenceTermScaleScaled<1e-180);
    REQUIRE(std::abs(edge.holeSgHighPrecisionReferenceFluxScaled)/edge.holeSgHighPrecisionReferenceTermScaleScaled<1e-180);
    REQUIRE(edge.electronSgPopulationPrefactor0Internal==2e16);
    REQUIRE(edge.holeSgPopulationPrefactor0Internal==7e15);
    REQUIRE(edge.splitSgReferenceDecimalDigits==200);
    REQUIRE((runtime.residual(x,{}).array()==before.array()).all());
    REQUIRE(runtime.current(x,"drain")==current);REQUIRE((runtime.low().array()==low.array()).all());
}
TEST_CASE("Split SG independent references retain sub ULP drive and reverse orientation", "[split_sg_diagnostics]") {
    auto model=sgFixture();SplitDDRuntime runtime(model);VectorXd x(6);x<<1.,2.,.125,.125,.125,.125;
    auto low=runtime.low();low(3)=1e-28;low(5)=-2e-28;runtime.setLow(low);
    const auto a=runtime.edges(x).at(0);
    REQUIRE(a.phin0_V==a.phin1_V);REQUIRE(a.phip0_V==a.phip1_V);
    REQUIRE(a.electronFlux!=0.);REQUIRE(a.holeFlux!=0.);
    REQUIRE(std::abs(a.electronSgHighPrecisionReferenceFluxScaled/a.electronFlux-1)<1e-14);
    REQUIRE(std::abs(a.holeSgHighPrecisionReferenceFluxScaled/a.holeFlux-1)<1e-14);
    const auto packet=nlohmann::json::parse(a.splitSgOperandsJson);
    REQUIRE(packet["primitives"]["packed_low"][3].get<double>()==low(3));
    REQUIRE(packet["primitives"]["electron_prefactors_internal"][0]==2e16);
    // Direct cache serialization must retain a drive hidden by double output.
    using split_dd::Wide;
    const auto& cache=packet.at("production_cache");
    const Wide q0(cache["phin_V"][0].get<std::string>()),q1(cache["phin_V"][1].get<std::string>());
    REQUIRE(q0!=q1);
    REQUIRE(static_cast<double>(q0)==static_cast<double>(q1));
    const Wide vt(model["nodes"][0][10].get<double>());
    const Wide psi0(cache["psi_V"][0].get<std::string>());
    const Wide n0(cache["n_internal"][0].get<std::string>());
    REQUIRE(abs(log(n0/Wide(2e16))-(psi0-q0)/vt)<Wide("1e-90"));
    REQUIRE(static_cast<double>(Wide(cache["flux_scaled"][0].get<std::string>()))==a.electronFlux);
    // Reversing the edge must reverse both particle fluxes without changing
    // any state or constitutive parameter.
    model["transport_edges"][0]["i"]=1;model["transport_edges"][0]["j"]=0;
    SplitDDRuntime reversed(model);reversed.setLow(low);const auto b=reversed.edges(x).at(0);
    REQUIRE(std::abs(b.electronSgHighPrecisionReferenceFluxScaled/a.electronFlux+1)<1e-14);
    REQUIRE(std::abs(b.holeSgHighPrecisionReferenceFluxScaled/a.holeFlux+1)<1e-14);
}
TEST_CASE("Split SG references cover finite drive and common intrinsic fallback", "[split_sg_diagnostics]") {
    auto model=sgFixture();model.erase("population_prefactors");SplitDDRuntime runtime(model);
    VectorXd x(6);x<<1.,2.,.125,.25,.2,.1;const auto a=runtime.edges(x).at(0);
    REQUIRE(a.electronSgPopulationPrefactor0Internal==1e16);
    REQUIRE(a.holeSgPopulationPrefactor1Internal==1e16);
    REQUIRE(std::abs(a.electronSgHighPrecisionReferenceFluxScaled/a.electronFlux-1)<1e-14);
    REQUIRE(std::abs(a.holeSgHighPrecisionReferenceFluxScaled/a.holeFlux-1)<1e-14);
}

TEST_CASE("Split SG absent transport support is explicitly unavailable", "[split_sg_diagnostics]") {
    auto model=sgFixture();model["transport_edges"][0]["weights"]=nlohmann::json::array();
    SplitDDRuntime runtime(model);VectorXd x(6);x<<1.,2.,.125,.25,.2,.1;
    const auto edge=runtime.edges(x).at(0);const auto packet=nlohmann::json::parse(edge.splitSgOperandsJson);
    REQUIRE_FALSE(edge.electronSgBoltzmannDecompositionAvailable);
    REQUIRE(edge.splitSgReferenceDecimalDigits==0);
    REQUIRE_FALSE(packet["production_transport_active"].get<bool>());
    REQUIRE_FALSE(packet.contains("reference_eta"));
    REQUIRE(edge.electronFlux==0.);REQUIRE(edge.holeFlux==0.);
}

TEST_CASE("Split boundary retains tiny bias and is idempotent", "[split_boundary]") {
    using split_dd::Wide;
    constexpr double scale=.02585,bias=1e-22;
    auto model=fixture();model["nodes"][0][7]=.8;model["nodes"][0][8]=-.7;
    SplitDDRuntime runtime(model);VectorXd x=VectorXd::Zero(6);
    CoupledDDBoundaryConditions zero,tiny;
    for(int b=0;b<3;++b) {
        setSplitBoundaryTarget(zero,b,0,b==0?.8:0.,0.,scale);
        setSplitBoundaryTarget(tiny,b,0,b==0?.8:0.,bias,scale);
    }
    runtime.applyBoundary(x,zero);auto base=runtime.residual(x,zero);auto before=runtime.residual(x,tiny);
    for(int b=0;b<3;++b)REQUIRE(std::abs(before[b*2]-base[b*2]+bias/scale)<1e-28);
    runtime.applyBoundary(x,tiny);auto r=runtime.residual(x,tiny);
    for(int b=0;b<3;++b)REQUIRE(std::abs(r[b*2])<1e-30);
    const auto hi=x,lo=runtime.low();runtime.applyBoundary(x,tiny);
    REQUIRE((hi.array()==x.array()).all());REQUIRE((lo.array()==runtime.low().array()).all());
    REQUIRE(tiny.psiLow.at(0)!=0.);
    DDSolution saved;runtime.save(x,saved);SplitDDRuntime again(model);
    const auto loaded=again.restore(saved);REQUIRE((again.low().array()==lo.array()).all());
    REQUIRE((again.residual(loaded,tiny).array()==r.array()).all());
    REQUIRE_THROWS(setSplitBoundaryTarget(tiny,0,0,0.,0.,0.));
}

#include "vela/numerics/StableBoltzmannSg.h"
#include "vela/numerics/SplitContactVoltage.h"
TEST_CASE("Contact low voltage is retained by boundary targets circuit and archive", "[split_boundary]") {
    using split_dd::Wide;
    constexpr double voltage=.075,scale=.02585,dv=-4.4e-25;
    const numerics::SplitContactVoltage base{{voltage,0}};
    const auto changed=base.shifted(1.,dv);
    REQUIRE(changed.value.hi==voltage);REQUIRE(changed.value.lo==dv);
    REQUIRE(changed.loadLine(1e12,0.,voltage)==dv);
    REQUIRE(base.shifted(0.,dv).value.lo==0.);
    SplitDDRuntime runtime(fixture());VectorXd x=VectorXd::Zero(6);
    CoupledDDBoundaryConditions before,after;
    for(int b=0;b<3;++b) {
        setSplitBoundaryTarget(before,b,0,b==0?.8:0.,voltage,scale);
        setSplitBoundaryTarget(after,b,0,b==0?.8:0.,changed.value.hi,scale,changed.value.lo);
    }
    runtime.applyBoundary(x,before);const auto low=runtime.low();
    const auto beforeR=runtime.residual(x,before),afterR=runtime.residual(x,after);
    for(int b=0;b<3;++b)REQUIRE(std::abs(afterR[2*b]-beforeR[2*b]+dv/scale)<1e-30);
    auto shifted=before;shiftSplitBoundaryVoltage(shifted,0,changed.value.lo,scale);
    for(int b=0;b<3;++b) {
        const auto& hi=b==0?after.psi:b==1?after.phin:after.phip;
        const auto& lo=b==0?after.psiLow:b==1?after.phinLow:after.phipLow;
        const auto& shi=b==0?shifted.psi:b==1?shifted.phin:shifted.phip;
        const auto& slo=b==0?shifted.psiLow:b==1?shifted.phinLow:shifted.phipLow;
        REQUIRE(abs(Wide(hi.at(0))+Wide(lo.at(0))-Wide(shi.at(0))-Wide(slo.at(0)))<Wide("1e-30"));
    }
    runtime.applyBoundary(x,after);const auto r=runtime.residual(x,after);
    DDSolution saved;runtime.save(x,saved);SplitDDRuntime again(fixture());auto restored=again.restore(saved);
    REQUIRE((again.residual(restored,after).array()==r.array()).all());
    for(int b=0;b<3;++b)REQUIRE(std::abs(r[2*b])<1e-30);
    REQUIRE_THROWS(changed.loadLine(-1,0,voltage));
}

TEST_CASE("Constant Boltzmann SG derivative preserves equilibrium and relative directions", "[split_sg_diagnostics]") {
    using split_dd::Wide;
    const auto flux=[](Wide cp,Wide eta,Wide q)->Wide {
        return cp*(eta==0?Wide(1):Wide(eta/expm1(eta)))*expm1(q);
    };
    for(double eta:{-8.,-1e-10,0.,1e-10,2.}) {
        const Wide cp=1e12,q=Wide("1e-30"),common=Wide("1e-18");
        const Wide d=numerics::stableBoltzmannSgDirection<Wide>(cp,Wide(eta),q,-common,0,0);
        REQUIRE(abs(d+flux(cp,Wide(eta),q)*common)<abs(d)*Wide("1e-30"));
        REQUIRE(numerics::stableBoltzmannSgDirection<Wide>(cp,Wide(eta),0,1,2,0)==0);
        const Wide finiteQ=.13,dl=.2,de=-.3,dq=.7,h=Wide("1e-8");
        const Wide analytic=numerics::stableBoltzmannSgDirection<Wide>(cp,Wide(eta),finiteQ,dl,de,dq);
        const Wide plus=flux(cp*exp(h*dl),Wide(eta)+h*de,finiteQ+h*dq);
        const Wide minus=flux(cp*exp(-h*dl),Wide(eta)-h*de,finiteQ-h*dq);
        REQUIRE(abs((plus-minus)/(2*h)-analytic)<abs(analytic)*Wide("1e-14"));
    }
}

TEST_CASE("Rounded nonzero voltage projection removes boundary tail but preserves free updates", "[split_boundary]") {
    constexpr double scale=.02585,voltage=.075,dv=-4.4e-25;
    REQUIRE(voltage+dv==voltage);
    SplitDDRuntime runtime(fixture());VectorXd x=VectorXd::Zero(6);
    CoupledDDBoundaryConditions bcs;
    for(int block=0;block<3;++block)setSplitBoundaryTarget(bcs,block,0,0.,voltage,scale);
    runtime.applyBoundary(x,bcs);const auto base=runtime.low();
    VectorXd step=VectorXd::Constant(6,dv/scale);
    auto trial=runtime.candidate(x,step,1.,base);
    const auto fullDelta=runtime.coordinateDelta(trial,x,base);
    REQUIRE_THROWS(runtime.coordinateDelta(trial,x,VectorXd::Zero(5)));
    const double tailRoundoff=4*std::numeric_limits<double>::epsilon()*
        (base.lpNorm<Eigen::Infinity>()+step.lpNorm<Eigen::Infinity>());
    REQUIRE((fullDelta-step).lpNorm<Eigen::Infinity>()<=tailRoundoff);
    REQUIRE((fullDelta-step).lpNorm<Eigen::Infinity>()<step.lpNorm<Eigen::Infinity>()*1e-6);
    const auto before=runtime.residual(trial,bcs);const auto hi=trial,lo=runtime.low();
    REQUIRE(std::abs(before[0])>1e-24);
    runtime.applyBoundary(trial,bcs);const auto after=runtime.residual(trial,bcs);
    const auto projectedDelta=runtime.coordinateDelta(trial,x,base);
    for(int block=0;block<3;++block) {
        REQUIRE(std::abs(after[2*block])<1e-30);
        REQUIRE(trial[2*block+1]==hi[2*block+1]);
        REQUIRE(runtime.low()[2*block+1]==lo[2*block+1]);
        REQUIRE(projectedDelta[2*block]==0.);
        REQUIRE(projectedDelta[2*block+1]==fullDelta[2*block+1]);
    }
    DDSolution saved;runtime.save(trial,saved);SplitDDRuntime reread(fixture());
    const auto restored=reread.restore(saved);
    REQUIRE((reread.residual(restored,bcs).array()==after.array()).all());
}
