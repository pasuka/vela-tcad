#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/physics/SiliconThermalPhysics.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/RecombinationModel.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/equation/LatticeBandEdgeWork.h"
#include "vela/discretization/ThermalSgCurrent.h"
#include "vela/discretization/ScharfetterGummel.h"
#include <cmath>
#include <limits>
#include <stdexcept>
using namespace vela;
TEST_CASE("Thermal preparations preserve values and partials and reject stale identities", "[thermal][silicon][preparation]") {
    SiliconThermalPhysics model,other;
    const std::array<ThermalQuantity SiliconThermalResult::*,18> quantities{
        &SiliconThermalResult::bandgap_eV,&SiliconThermalResult::affinity_eV,&SiliconThermalResult::Nc_m3,&SiliconThermalResult::Nv_m3,
        &SiliconThermalResult::ni_m3,&SiliconThermalResult::effectiveNi_m3,&SiliconThermalResult::conductionBand_eV,&SiliconThermalResult::valenceBand_eV,
        &SiliconThermalResult::electrons_m3,&SiliconThermalResult::holes_m3,&SiliconThermalResult::electronEta,&SiliconThermalResult::holeEta,
        &SiliconThermalResult::srhRate_m3_per_s,&SiliconThermalResult::augerRate_m3_per_s,&SiliconThermalResult::electronLifetime_s,&SiliconThermalResult::holeLifetime_s,
        &SiliconThermalResult::augerElectron_m6_per_s,&SiliconThermalResult::augerHole_m6_per_s};
    const auto doping=model.prepareDoping(1e23,2e22);
    CHECK_THROWS_AS(other.prepareTemperature(300.,doping),std::invalid_argument);
    for(Real temperature:{100.,300.,401.,514.,600.})for(Real shift:{0.,40.}) {
        const auto prepared=model.prepareTemperature(temperature,doping);
        SiliconThermalState state{shift+.5,.01,.03,temperature,1e23,2e22,shift,shift};
        const auto plain=model.evaluate(state),cached=model.evaluate(state,prepared);
        for(auto member:quantities){CHECK((plain.*member).value==(cached.*member).value);CHECK((plain.*member).derivative==(cached.*member).derivative);}
        CHECK(plain.bandgapNarrowing_eV==cached.bandgapNarrowing_eV);
        const auto densities=model.carrierDensities(state,prepared);
        CHECK(densities[0].value==plain.electrons_m3.value);CHECK(densities[0].derivative==plain.electrons_m3.derivative);
        CHECK(densities[1].value==plain.holes_m3.value);CHECK(densities[1].derivative==plain.holes_m3.derivative);
        auto changed=state;changed.temperature_K+=.01;
        CHECK_THROWS_AS(model.evaluate(changed,prepared),std::invalid_argument);
        changed=state;changed.donors_m3*=1.01;
        CHECK_THROWS_AS(model.evaluate(changed,prepared),std::invalid_argument);
        CHECK_THROWS_AS(other.evaluate(state,prepared),std::invalid_argument);
        auto copy=model;CHECK(copy.evaluate(state,prepared).electrons_m3.value==plain.electrons_m3.value);
    }
}
TEST_CASE("Thermal silicon recovers audited 300 K material and Fermi statistics", "[thermal][silicon]") {
    SiliconThermalPhysics model;
    auto r=model.evaluate({.15,.03,.07,300.,2e23,1e22});
    CHECK(r.bandgap_eV.value==Catch::Approx(1.12416).epsilon(1e-14));
    CHECK(r.Nc_m3.value==Catch::Approx(2.856679069182081e25).epsilon(1e-13));
    CHECK(r.Nv_m3.value==Catch::Approx(3.1046570334558958e25).epsilon(1e-13));
    CHECK(r.ni_m3.value==Catch::Approx(1.0750016577953705e16).epsilon(1e-13));
    CHECK(r.electrons_m3.value==Catch::Approx(electronDensity(r.effectiveNi_m3.value,r.Nc_m3.value,.15,.03,constants::Vt_300,CarrierStatisticsModel::FermiDirac)).epsilon(1e-12));
    CHECK(r.holes_m3.value==Catch::Approx(holeDensity(r.effectiveNi_m3.value,r.Nv_m3.value,.15,.07,constants::Vt_300,CarrierStatisticsModel::FermiDirac)).epsilon(1e-12));
    CHECK(r.electronLifetime_s.value==Catch::Approx(1e-5/22.));
}
TEST_CASE("Thermal silicon preserves equilibrium source and potential gauge", "[thermal][silicon]") {
    SiliconThermalPhysics model;
    for(Real t:{100.,300.,401.,514.,600.}){
        auto a=model.evaluate({.1,.02,.02,t,1e23,3e22});
        CHECK(a.srhRate_m3_per_s.value==0.);CHECK(a.augerRate_m3_per_s.value==0.);
        CHECK(a.srhRate_m3_per_s.derivative[3]==0.);
        auto b=model.evaluate({2.1,2.02,2.02,t,1e23,3e22});
        CHECK(a.electrons_m3.value==Catch::Approx(b.electrons_m3.value).epsilon(1e-12));
        CHECK(a.holes_m3.value==Catch::Approx(b.holes_m3.value).epsilon(1e-12));
        CHECK(a.bandgapNarrowing_eV==b.bandgapNarrowing_eV);
    }
}
TEST_CASE("Thermal silicon all four nodal partials match independent perturbations", "[thermal][silicon][jacobian]") {
    SiliconThermalPhysics model;
    const std::array<Real SiliconThermalState::*,4> fields{&SiliconThermalState::potential_V,&SiliconThermalState::electronQf_V,&SiliconThermalState::holeQf_V,&SiliconThermalState::temperature_K};
    const std::array<ThermalQuantity SiliconThermalResult::*,10> quantities{&SiliconThermalResult::bandgap_eV,&SiliconThermalResult::affinity_eV,&SiliconThermalResult::Nc_m3,&SiliconThermalResult::Nv_m3,&SiliconThermalResult::ni_m3,&SiliconThermalResult::effectiveNi_m3,&SiliconThermalResult::electrons_m3,&SiliconThermalResult::holes_m3,&SiliconThermalResult::srhRate_m3_per_s,&SiliconThermalResult::augerRate_m3_per_s};
    for(Real t:{200.,300.,401.,514.,600.})for(Real psi:{-.5,.1,.7})for(Real split:{-.03,.06}) {
        SiliconThermalState s{psi,0.,split,t,1e23,2e22};auto actual=model.evaluate(s);
        for(int k=0;k<4;++k)for(Real fraction:{1.,.25}) {
            const Real step=(k==3?.002:1e-6)*fraction;auto a=s,b=s;a.*fields[k]+=step;b.*fields[k]-=step;
            auto plus=model.evaluate(a),minus=model.evaluate(b);
            for(auto quantity:quantities){auto q=actual.*quantity;Real scale=std::max(1e-100,std::abs(q.value));Real unit=k==3?t:1.;
                INFO("T="<<t<<" psi="<<psi<<" split="<<split<<" column="<<k);
                CHECK(q.derivative[k]*unit/scale==Catch::Approx(((plus.*quantity).value-(minus.*quantity).value)*unit/(2.*step*scale)).epsilon(3e-6).margin(2e-7));
            }
        }
    }
}
TEST_CASE("Thermal silicon rejects invalid temperatures and coefficients", "[thermal][silicon]") {
    SiliconThermalPhysics model;
    CHECK_THROWS_AS(model.evaluate({0.,0.,0.,49.,0.,0.}),std::invalid_argument);
    CHECK_THROWS_AS(model.evaluate({0.,0.,0.,300.,-1.,0.}),std::invalid_argument);
    SiliconThermalParameters p;p.bandgapBeta_K=std::numeric_limits<Real>::quiet_NaN();
    CHECK_THROWS_AS(SiliconThermalPhysics(p),std::invalid_argument);
}

TEST_CASE("Audited Auger defaults to recombination only and preserves explicit generation", "[thermal][silicon][auger]") {
    SiliconThermalPhysics audited;
    SiliconThermalParameters p;p.augerWithGeneration=true;SiliconThermalPhysics signedModel(p);
    for(Real t:{300.,515.})for(Real split:{-.03,0.,.03}){
        SiliconThermalState state{.1,0.,split,t,1e23,2e22};
        const auto a=audited.evaluate(state),b=signedModel.evaluate(state);
        CHECK(a.srhRate_m3_per_s.value==b.srhRate_m3_per_s.value);
        if(split<=0.){
            CHECK(a.augerRate_m3_per_s.value==0.);
            for(Real derivative:a.augerRate_m3_per_s.derivative)CHECK(derivative==0.);
        }else{
            CHECK(a.augerRate_m3_per_s.value>0.);
            CHECK(a.augerRate_m3_per_s.value==b.augerRate_m3_per_s.value);
            CHECK(a.augerRate_m3_per_s.derivative==b.augerRate_m3_per_s.derivative);
        }
        if(split<0.){
            CHECK(b.augerRate_m3_per_s.value<0.);
            auto hi=state,lo=state;hi.temperature_K+=.002;lo.temperature_K-=.002;
            CHECK(b.augerRate_m3_per_s.derivative[3]==Catch::Approx((signedModel.evaluate(hi).augerRate_m3_per_s.value-signedModel.evaluate(lo).augerRate_m3_per_s.value)/.004).epsilon(3e-6));
        }
    }
}

TEST_CASE("Default lattice edge work preserves signs and discrete power balance", "[thermal][silicon][heat]") {
    auto r=latticeBandEdgeWork(-2.,-.1,1.,.8,0.,-.2);
    CHECK(r.power_W_per_m==Catch::Approx(.42));
    CHECK(latticeBandEdgeWork(2.,.1,.8,1.,-.2,0.).power_W_per_m==r.power_W_per_m);
    CHECK(latticeBandEdgeWork(-2.,0.,0.,1.,0.,1.).power_W_per_m==-2.);
    // Constant through-current: work telescopes to terminal electrical power.
    Real sum=0.;const std::array<Real,4> ec{1.,.8,.3,-.2},ev{0.,-.2,-.7,-1.2};
    for(int i=0;i<3;++i)sum+=latticeBandEdgeWork(-2.,-.1,ec[i],ec[i+1],ev[i],ev[i+1]).power_W_per_m;
    CHECK(sum==Catch::Approx(2.1*1.2));
    std::array<Real,6> input{-2.,-.1,1.,.8,0.,-.2};
    for(int k=0;k<6;++k){auto a=input,b=input;a[k]+=1e-5;b[k]-=1e-5;
        const auto eval=[](const auto&v){return latticeBandEdgeWork(v[0],v[1],v[2],v[3],v[4],v[5]).power_W_per_m;};
        CHECK(r.derivative[k]==Catch::Approx((eval(a)-eval(b))/2e-5).margin(1e-10));
    }
}

TEST_CASE("Thermal SG recovers uniform-temperature generalized SG and flat QF", "[thermal][silicon][sg]"){
    SiliconThermalPhysics model;
    for(bool electron:{false,true})for(Real psi:{-.5,.1,.7}){
        SiliconThermalState a{psi,.01,.04,300.,1e23,1e22},b{psi+.05,.015,.047,300.,2e23,1e22};
        auto pa=model.evaluate(a),pb=model.evaluate(b);Real mu=.03,w=.7,vt=constants::Vt_300;
        Real n0=electron?pa.electrons_m3.value:pa.holes_m3.value,n1=electron?pb.electrons_m3.value:pb.holes_m3.value;
        Real eta0=electron?(-a.electronQf_V-pa.conductionBand_eV.value)/vt:(pa.valenceBand_eV.value+a.holeQf_V)/vt;
        Real eta1=electron?(-b.electronQf_V-pb.conductionBand_eV.value)/vt:(pb.valenceBand_eV.value+b.holeQf_V)/vt;
        Real q0=electron?a.electronQf_V:a.holeQf_V,q1=electron?b.electronQf_V:b.holeQf_V;
        Real drift=q1-q0+(electron?1.:-1.)*vt*(eta1-eta0);
        Real expected=electron?-constants::q*sgElectronFermiDiracContinuityFlux(n0,n1,eta0,eta1,drift,q0,q1,vt,mu*vt*w):
            constants::q*sgHoleFermiDiracContinuityFlux(n0,n1,eta0,eta1,drift,q0,q1,vt,mu*vt*w);
        CHECK(thermalSgCurrent(a,pa,b,pb,mu,w,electron).current_A_per_m==Catch::Approx(expected).epsilon(2e-12));
        auto forward=thermalSgCurrent(a,pa,b,pb,mu,w,electron);auto reverse=thermalSgCurrent(b,pb,a,pa,mu,w,electron);
        CHECK(forward.current_A_per_m==Catch::Approx(-reverse.current_A_per_m).epsilon(2e-12));
        b.electronQf_V=a.electronQf_V;b.holeQf_V=a.holeQf_V;b.temperature_K=450.;pb=model.evaluate(b);
        CHECK(thermalSgCurrent(a,pa,b,pb,mu,w,electron).current_A_per_m==0.);
    }
}
TEST_CASE("Thermal SG eight state columns include temperature and equilibrium conductance", "[thermal][silicon][sg][jacobian]"){
    SiliconThermalPhysics model;std::array<Real SiliconThermalState::*,4> fields{&SiliconThermalState::potential_V,&SiliconThermalState::electronQf_V,&SiliconThermalState::holeQf_V,&SiliconThermalState::temperature_K};
    for(bool electron:{false,true})for(bool flat:{false,true})for(bool same:{false,true}){
        std::array<SiliconThermalState,2> s{{{.12,.01,.04,350.,1e23,1e22},{.17,flat?.01:.015,flat?.04:.047,450.,1e23,1e22}}};
        if(same)s[1]=s[0];
        const auto eval=[&](const auto& a){return thermalSgCurrent(a[0],model.evaluate(a[0]),a[1],model.evaluate(a[1]),.03,.7,electron);};
        auto r=eval(s);
        for(int k=0;k<8;++k)for(Real fraction:{1.,.25}){
            const Real step=(k%4==3?.001:2e-7)*fraction;auto a=s,b=s;a[k/4].*fields[k%4]+=step;b[k/4].*fields[k%4]-=step;
            Real fd=(eval(a).current_A_per_m-eval(b).current_A_per_m)/(2.*step);
            const Real scale=std::max({std::abs(r.derivative[k]),std::abs(fd),1e-100});
            CHECK(r.derivative[k]/scale==Catch::Approx(fd/scale).epsilon(4e-5).margin(4e-7));
        }
    }
}


TEST_CASE("Thermal SG preserves QF increments below the absolute potential spacing", "[thermal][silicon][sg][reference]") {
    SiliconThermalPhysics physics;
    for(bool electron:{false,true}) {
        SiliconThermalState a{40.6,0.,0.,400.,1e23,1e22,40.,40.},b=a;
        if(electron)b.electronQf_V=1e-18;else b.holeQf_V=1e-18;
        auto pa=physics.evaluate(a),pb=physics.evaluate(b);
        auto zero=thermalSgCurrent(a,pa,a,pa,.03,.7,electron);
        auto small=thermalSgCurrent(a,pa,b,pb,.03,.7,electron);
        REQUIRE(small.current_A_per_m!=0.);
        CHECK(small.current_A_per_m==Catch::Approx(zero.derivative[electron?5:6]*1e-18).epsilon(1e-11));
        CHECK(small.current_A_per_m==Catch::Approx(-thermalSgCurrent(b,pb,a,pa,.03,.7,electron).current_A_per_m).epsilon(1e-11));
        CHECK(pb.srhRate_m3_per_s.value!=0.);
        auto shifted=a;shifted.potential_V-=40.;shifted.electronQfReference_V=0.;shifted.holeQfReference_V=0.;
        auto ps=physics.evaluate(shifted);
        CHECK(pa.electrons_m3.value==Catch::Approx(ps.electrons_m3.value).epsilon(1e-12));
        CHECK(pa.holes_m3.value==Catch::Approx(ps.holes_m3.value).epsilon(1e-12));
    }
}
