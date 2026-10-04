#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/physics/IalMobility.h"
#include "vela/physics/IalMobilityJson.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "vela/core/IalKernelProfiling.h"
#include <chrono>
#include <iostream>
#include <vector>

using namespace vela;

TEST_CASE("IALMob omitted screening method uses safeguarded Halley and keeps explicit legacy", "[mobility][ialmob][screening_default]") {
    IalKernelProfilingScope profile(false);
    const auto root=ialScreeningMinimum(1.258,401.);
    CHECK(ialKernelProfile.screeningCandidateCalls==1);
    CHECK(root==ialScreeningMinimum(1.258,401.,IalScreeningMethod::Halley));
    IalScreeningCache cache;
    CHECK(cache.minimum(1.258,401.)==root);
    CHECK(cache.minimum(1.258,401.,IalScreeningMethod::Halley)==root);
    CHECK(cache.size()==1);
    const auto before=ialKernelProfile.screeningCandidateCalls;
    cache.minimum(1.258,401.,IalScreeningMethod::Legacy);
    CHECK(ialKernelProfile.screeningCandidateCalls==before);
    CHECK(cache.size()==2);
    IalMobility model(IalMobility::siliconDefaults(true),true);
    CHECK(model.screeningMethod()==IalScreeningMethod::Halley);
    CHECK(model.withScreeningMethod(IalScreeningMethod::Legacy).screeningMethod()==IalScreeningMethod::Legacy);
    using nlohmann::json;
    json cfg={{"mobility_SI",{{"ialmob",json::object()}}}};
    CHECK(ial_json::electrothermalScreeningMethod(cfg)==IalScreeningMethod::Halley);
    cfg["mobility_SI"]["ialmob"]["screening_method"]="legacy";
    CHECK(ial_json::electrothermalScreeningMethod(cfg)==IalScreeningMethod::Legacy);
    cfg["diagnostic_ialmob_screening_method"]="halley";
    CHECK_THROWS_AS(ial_json::electrothermalScreeningMethod(cfg),std::invalid_argument);
    cfg["mobility_SI"]["ialmob"].erase("screening_method");
    cfg["diagnostic_ialmob_screening_method"]="legacy";
    CHECK(ial_json::electrothermalScreeningMethod(cfg)==IalScreeningMethod::Legacy);
}

TEST_CASE("Screening candidate scalar timing", "[.][screening_benchmark]") {
    const std::array methods{IalScreeningMethod::Legacy,IalScreeningMethod::Newton,IalScreeningMethod::Halley,IalScreeningMethod::Toms748};
    const std::array names{"legacy","newton","halley","toms748"};
    std::vector<std::pair<Real,Real>> inputs;
    for(Real mass:{1.,1.258})for(int i=0;i<257;++i)inputs.emplace_back(mass,50.+950.*i/256.);
    for(int round=0;round<3;++round)for(int order=0;order<4;++order){
        const int index=round%2?3-order:order;IalKernelProfilingScope scope(false);
        Real checksum=0.;const auto start=std::chrono::steady_clock::now();
        for(int repeat=0;repeat<256;++repeat)for(const auto& [mass,temp]:inputs)
            checksum+=ialScreeningMinimum(mass,temp,methods[index]);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        REQUIRE(std::isfinite(checksum));
        std::cout<<"SCREENING_BENCH {\"round\":"<<round<<",\"method\":\""<<names[index]
            <<"\",\"calls\":"<<inputs.size()*256<<",\"seconds\":"<<seconds
            <<",\"function_evaluations\":"<<ialKernelProfile.screeningFunctionEvaluations
            <<",\"fallbacks\":"<<ialKernelProfile.screeningFallbacks<<",\"checksum\":"<<checksum<<"}\n";
    }
}

TEST_CASE("Screening candidates preserve positive stationary roots and method cache identity", "[mobility][ialmob][screening_candidate]") {
    using MP=boost::multiprecision::cpp_dec_float_50;
    for(Real mass:{.01,.258,1.,1.258,100.})for(Real temperature:{50.,299.,300.,301.,600.,1000.,100000.}){
        CAPTURE(mass,temperature);
        const MP a=pow(MP(temperature)/(MP(300)*MP(mass)),MP(".28227"));
        const MP b=pow(MP(mass)*MP(300)/MP(temperature),MP(".72169"));
        MP lo("1e-12"),hi("1e12");
        for(int i=0;i<210;++i){
            const MP p=sqrt(lo*hi);
            const MP derivative=MP(".89233")*MP(".19778")*a/pow(MP(".41372")+a*p,MP("1.19778"))
                -MP(".005978")*MP("1.80618")*b/pow(b*p,MP("2.80618"));
            if(derivative<0)lo=p;else hi=p;
        }
        const Real reference=static_cast<Real>(sqrt(lo*hi));
        IalScreeningCache cache;
        for(auto method:{IalScreeningMethod::Legacy,IalScreeningMethod::Newton,IalScreeningMethod::Halley,IalScreeningMethod::Toms748}){
            CAPTURE(static_cast<int>(method));
            const Real root=ialScreeningMinimum(mass,temperature,method);
            REQUIRE(std::isfinite(root));REQUIRE(root>0.);
            CHECK(std::abs(root/reference-1.)<1e-12);
            const MP p(root);
            const MP left=MP(".89233")*MP(".19778")*a/pow(MP(".41372")+a*p,MP("1.19778"));
            const MP right=MP(".005978")*MP("1.80618")*b/pow(b*p,MP("2.80618"));
            CHECK(static_cast<Real>(abs((left-right)/left))<3e-12);
            CHECK(cache.minimum(mass,temperature,method)==root);
            CHECK(cache.minimum(mass,temperature,method)==root);
        }
        CHECK(cache.size()==4);
    }
    for(auto method:{IalScreeningMethod::Newton,IalScreeningMethod::Halley,IalScreeningMethod::Toms748}){
        IalKernelProfilingScope profile(false);
        CHECK(ialScreeningMinimum(.001,400.,method)==ialScreeningMinimum(.001,400.,IalScreeningMethod::Legacy));
        CHECK(ialKernelProfile.screeningFallbacks==1);
        CHECK_THROWS_AS(ialScreeningMinimum(1.,49.,method),std::invalid_argument);
        CHECK_THROWS_AS(ialScreeningMinimum(0.,300.,method),std::invalid_argument);
    }
    CHECK_THROWS_AS(ialScreeningMethod("invalid"),std::invalid_argument);
}

TEST_CASE("Screening candidates preserve mobility and coupled temperature response", "[mobility][ialmob][screening_candidate]") {
    for(bool electron:{true,false})for(Real t:{50.,299.,300.,301.,400.,1000.})for(Real density:{1e18,1e23,1e28}){
        const IalMobility base(IalMobility::siliconDefaults(electron),electron,IalScreeningMethod::Legacy);
        const IalMobilityState state{density,.8*density,.3*density,.7*density,1e6,2e-8,t};
        const auto old=base.evaluateWithDerivatives(state);
        for(auto method:{IalScreeningMethod::Newton,IalScreeningMethod::Halley,IalScreeningMethod::Toms748}){
            CAPTURE(electron,t,density,static_cast<int>(method));
            const auto model=base.withScreeningMethod(method);
            const auto now=model.evaluateWithDerivatives(state);
            const Real scale=old.result.mobility_m2_per_Vs;
            CHECK(std::abs(now.result.mobility_m2_per_Vs/scale-1.)<2e-12);
            const std::array<Real,6> units{density,density,density,density,1e6,2e-8};
            for(int i=0;i<6;++i)CHECK(std::abs(now.derivative_SI[i]-old.derivative_SI[i])*units[i]/scale<2e-10);
            CHECK(std::abs(now.temperatureDerivative_m2_per_Vs_K-old.temperatureDerivative_m2_per_Vs_K)*t/scale<2e-10);
            if(t>50.){
                auto plus=state,minus=state;const Real step=t*1e-5;plus.temperature_K+=step;minus.temperature_K-=step;
                const Real fd=(model.evaluate(plus).mobility_m2_per_Vs-model.evaluate(minus).mobility_m2_per_Vs)/(2.*step);
                CHECK(std::abs(fd-now.temperatureDerivative_m2_per_Vs_K)*t/scale<2e-5);
            }
        }
    }
}

TEST_CASE("IALMob differentiated evaluation preserves every scalar component", "[mobility][ialmob][temperature]")
{
    const std::array<Real IalMobilityResult::*,6> fields{&IalMobilityResult::mobility_m2_per_Vs,
        &IalMobilityResult::coulomb3d_m2_per_Vs,&IalMobilityResult::coulomb2d_m2_per_Vs,
        &IalMobilityResult::coulomb_m2_per_Vs,&IalMobilityResult::phonon_m2_per_Vs,
        &IalMobilityResult::roughness_m2_per_Vs};
    for(bool electron:{false,true}){
        const IalMobility model(IalMobility::siliconDefaults(electron),electron);
        IalScreeningCache cache;
        for(Real temperature:{50.,299.,300.,401.,600.,1000.})
        for(Real density:{0.,1e18,1e23,1e28})for(Real field:{0.,1e6,1e8}){
            CAPTURE(electron,temperature,density,field);
            const IalMobilityState state{density,.8*density,.3*density,.7*density,field,2e-8,temperature};
            IalMobilityResult scalar;
            try {scalar=model.evaluate(state);}
            catch(const std::runtime_error&){
                CHECK_THROWS_AS(model.evaluateWithDerivatives(state,&cache),std::runtime_error);
                continue;
            }
            IalMobilityResult dual;
            REQUIRE_NOTHROW(dual=model.evaluateWithDerivatives(state,&cache).result);
            for(auto member:fields)CHECK(scalar.*member==dual.*member);
        }
    }
}

TEST_CASE("IALMob screening optimization retains the frozen bisection root", "[mobility][ialmob][temperature]")
{
    IalScreeningCache cache;
    // Compare the previous fixed-budget oracle across both carrier masses and
    // a broad temperature range spanning 300 K.
    for (Real mass : {.1, .258, .5, 1., 1.5, 3.}) {
        for (int index=0;index<=100;++index) {
            const Real temperature=50.+9.5*index;
            Real lower=1e-12,upper=1e12;
            for(int iteration=0;iteration<100;++iteration) {
                const Real p=std::sqrt(lower*upper);
                const Real a=std::pow(temperature/(300.*mass),.28227);
                const Real b=std::pow(mass*300./temperature,.72169);
                const Real derivative=.89233*.19778*a/std::pow(.41372+a*p,1.19778)
                    -.005978*1.80618*b/std::pow(b*p,2.80618);
                if(derivative<0.)lower=p;else upper=p;
            }
            CHECK(cache.minimum(mass,temperature,IalScreeningMethod::Legacy)==std::sqrt(lower*upper));
        }
    }
}

TEST_CASE("IALMob predicted screening bracket preserves roots and stationarity", "[mobility][ialmob][temperature]")
{
    IalScreeningCache cache;
    // Includes the optimization boundaries and the unaccelerated fallback.
    for (Real mass : {.001,.01,.1,.258,.5,1.,1.258,3.,10.,100.,1000.}) {
        for (int i=0;i<=100;++i) {
            const Real temperature=50.*std::pow(4000.,i/100.);
            const Real a=std::pow(temperature/(300.*mass),.28227);
            const Real b=std::pow(mass*300./temperature,.72169);
            const auto terms=[&](Real p) {
                return std::array<Real,2>{.89233*.19778*a/std::pow(.41372+a*p,1.19778),
                    .005978*1.80618*b/std::pow(b*p,2.80618)};
            };
            Real lower=1e-12,upper=1e12;
            for(int k=0;k<100;++k) {
                const Real p=std::sqrt(lower*upper);const auto t=terms(p);
                if(t[0]-t[1]<0.)lower=p;else upper=p;
            }
            const Real result=cache.minimum(mass,temperature,IalScreeningMethod::Legacy);
            CHECK(result==std::sqrt(lower*upper));
            const auto t=terms(result);
            CHECK(std::abs(t[0]-t[1])/std::max(t[0],t[1])<1e-12);
        }
    }
}

TEST_CASE("IALMob local preparation tracks all node inputs and model identity", "[mobility][ialmob][temperature]")
{
    IalMobilityPreparationCache prepared;IalScreeningCache screening;
    const std::array<IalMobility,2> models{IalMobility(IalMobility::siliconDefaults(true),true),
        IalMobility(IalMobility::siliconDefaults(false),false)};
    const std::array fields{&IalMobilityResult::mobility_m2_per_Vs,&IalMobilityResult::coulomb3d_m2_per_Vs,
        &IalMobilityResult::coulomb2d_m2_per_Vs,&IalMobilityResult::coulomb_m2_per_Vs,
        &IalMobilityResult::phonon_m2_per_Vs,&IalMobilityResult::roughness_m2_per_Vs};
    const std::array inputs{&IalMobilityState::donors_m3,&IalMobilityState::acceptors_m3,
        &IalMobilityState::electrons_m3,&IalMobilityState::holes_m3,
        &IalMobilityState::interfaceDistance_m,&IalMobilityState::temperature_K};
    const IalMobilityState base{2e22,1e23,1e24,1e16,1e7,2e-9,367.};
    for(const auto& model:models)for(std::size_t varied=0;varied<=inputs.size();++varied){
        auto state=base;if(varied<inputs.size())state.*inputs[varied]*=1.01;
        for(Real field:{0.,1e5,1e7,1e8,1e5}){
            state.normalField_V_per_m=field;
            const auto expected=model.evaluateWithDerivatives(state);
            const auto actual=model.evaluateWithDerivatives(state,&screening,&prepared);
            const auto scalar=model.evaluate(state,&screening,&prepared);
            for(auto member:fields){CHECK(actual.result.*member==expected.result.*member);CHECK(scalar.*member==expected.result.*member);}
            CHECK(actual.derivative_SI==expected.derivative_SI);
            CHECK(actual.temperatureDerivative_m2_per_Vs_K==expected.temperatureDerivative_m2_per_Vs_K);
        }
    }
    CHECK(prepared.hits()>0);
}

TEST_CASE("IALMob intrinsic zero-field mobility is the lattice limit in SI", "[mobility][ialmob]")
{
    for (bool electron : {false,true}) {
        const auto parameters=IalMobility::siliconDefaults(electron);
        const IalMobility model(parameters,electron);
        for (double density : {0.,1e16,1e22}) {
            IalMobilityState state;
            if (electron) state.electrons_m3=density;
            else state.holes_m3=density;
            const auto result=model.evaluate(state);
            CHECK(result.mobility_m2_per_Vs == Catch::Approx(parameters.muMax*1e-4));
            CHECK(std::isinf(result.coulomb_m2_per_Vs));
            CHECK(std::isinf(result.roughness_m2_per_Vs));
        }
    }
}

TEST_CASE("IALMob recovers the bulk majority-carrier doping curve far from interfaces", "[mobility][ialmob]")
{
    for (bool electron : {false,true}) {
        auto parameters=IalMobility::siliconDefaults(electron);
        // Suppress clustering to isolate the analytical majority-carrier limit.
        parameters.nRefD=parameters.nRefA=1e100;
        parameters.lCrit=parameters.lCritC=1e-6;
        const IalMobility model(parameters,electron);
        double previous=parameters.muMax*1e-4;
        for (double doping_cm3 : {1e14,1e16,1e18,1e20}) {
            IalMobilityState state;
            state.interfaceDistance_m=1.;
            state.normalField_V_per_m=1e8;
            if (electron) state.donors_m3=state.electrons_m3=doping_cm3*1e6;
            else state.acceptors_m3=state.holes_m3=doping_cm3*1e6;
            const double value=model.evaluate(state).mobility_m2_per_Vs;
            const double expected=(parameters.muMin+(parameters.muMax-parameters.muMin)/
                (1.+std::pow(doping_cm3/parameters.nRef,parameters.alpha)))*1e-4;
            CHECK(value == Catch::Approx(expected).epsilon(1e-12));
            CHECK(value < previous);
            previous=value;
        }
    }
}

TEST_CASE("IALMob surface scattering strengthens with field and decays with distance", "[mobility][ialmob]")
{
    for (bool electron : {false,true}) {
        auto parameters=IalMobility::siliconDefaults(electron);
        parameters.lCrit=parameters.lCritC=1e-6;
        const IalMobility model(parameters,electron);
        IalMobilityState state;
        double previous=parameters.muMax*1e-4;
        for (double field : {1e3,1e5,1e7,1e8}) {
            state.normalField_V_per_m=field;
            const auto result=model.evaluate(state);
            CHECK(result.mobility_m2_per_Vs < previous);
            CHECK(result.mobility_m2_per_Vs > 0.);
            previous=result.mobility_m2_per_Vs;
        }
        for (double distance : {1e-9,1e-8,1e-7,1e-5}) {
            state.interfaceDistance_m=distance;
            const double value=model.evaluate(state).mobility_m2_per_Vs;
            CHECK(value > previous);
            previous=value;
        }
        CHECK(previous == Catch::Approx(parameters.muMax*1e-4).epsilon(1e-12));
    }
}

TEST_CASE("IALMob retains compensated impurities and opposite-carrier scattering", "[mobility][ialmob]")
{
    for (bool electron : {false,true}) {
        auto parameters=IalMobility::siliconDefaults(electron);
        parameters.lCrit=parameters.lCritC=1e-6;
        const IalMobility model(parameters,electron);
        IalMobilityState state;
        // Isolate bulk scattering: the 2D impurity term diverges at zero doping.
        state.interfaceDistance_m=1.;
        state.electrons_m3=state.holes_m3=1e16;
        const double intrinsic=model.evaluate(state).mobility_m2_per_Vs;
        state.donors_m3=state.acceptors_m3=1e24;
        const double compensated=model.evaluate(state).mobility_m2_per_Vs;
        CHECK(compensated < intrinsic);
        state.donors_m3=state.acceptors_m3=0.;
        if (electron) state.holes_m3=1e24;
        else state.electrons_m3=1e24;
        CHECK(model.evaluate(state).mobility_m2_per_Vs < intrinsic);
    }
}

TEST_CASE("IALMob remains finite across the isothermal LDMOS state domain", "[mobility][ialmob]")
{
    for (bool electron : {false,true}) {
        auto parameters=IalMobility::siliconDefaults(electron);
        parameters.lCrit=parameters.lCritC=1e-6;
        const IalMobility model(parameters,electron);
        for (double nd : {0.,1e20,1e24,1e27})
        for (double na : {0.,1e20,1e24,1e27})
        for (double n : {0.,1e16,1e24,1e27})
        for (double p : {0.,1e16,1e24,1e27})
        for (double field : {0.,1e4,1e7,1e8}) {
            const auto result=model.evaluate({nd,na,n,p,field,1e-9});
            REQUIRE(std::isfinite(result.mobility_m2_per_Vs));
            REQUIRE(result.mobility_m2_per_Vs > 0.);
            REQUIRE(result.mobility_m2_per_Vs <= parameters.muMax*1e-4);
        }
    }
}

TEST_CASE("IALMob cubic family selection respects signs permutations and scale", "[mobility][ialmob]")
{
    CHECK(IalMobility::orientationFamily({1.,0.,0.}) == 100);
    CHECK(IalMobility::orientationFamily({0.,-8.,0.}) == 100);
    CHECK(IalMobility::orientationFamily({1.,1.,0.}) == 110);
    CHECK(IalMobility::orientationFamily({0.,-2.,2.}) == 110);
    CHECK(IalMobility::orientationFamily({1.,1.,1.}) == 111);
    CHECK(IalMobility::orientationFamily({-1e300,1e300,-1e300}) == 111);
    CHECK(IalMobility::orientationFamily({1.,.1,0.}) == 100);
    CHECK(IalMobility::orientationFamily({1.,.8,.1}) == 110);
    CHECK(IalMobility::orientationFamily({1.,.8,.7}) == 111);
    CHECK_THROWS_AS(IalMobility::orientationFamily({0.,0.,0.}), std::invalid_argument);
    CHECK_THROWS_AS(IalMobility::orientationFamily({0.,INFINITY,0.}), std::invalid_argument);
}

TEST_CASE("IALMob rejects invalid input before numerical evaluation", "[mobility][ialmob]")
{
    const IalMobility model(IalMobility::siliconDefaults(true),true);
    IalMobilityState state;
    state.donors_m3=-1.;
    CHECK_THROWS_AS(model.evaluate(state), std::invalid_argument);
    state.donors_m3=std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(model.evaluate(state), std::invalid_argument);
    auto parameters=IalMobility::siliconDefaults(true);
    parameters.lCrit=0.;
    CHECK_THROWS_AS(IalMobility(parameters,true), std::invalid_argument);
}

TEST_CASE("IALMob temperature lattice limit and its analytic slope", "[mobility][ialmob][temperature]")
{
    for (bool electron : {false,true}) {
        const auto p=IalMobility::siliconDefaults(electron);
        const IalMobility model(p,electron);
        for (Real t : {50.,200.,300.,400.,514.,600.}) {
            IalMobilityState state; state.temperature_K=t;
            const auto r=model.evaluateWithDerivatives(state);
            const Real expected=p.muMax*1e-4*std::pow(t/300.,-p.theta);
            CHECK(r.result.mobility_m2_per_Vs==Catch::Approx(expected).epsilon(1e-13));
            CHECK(r.temperatureDerivative_m2_per_Vs_K==Catch::Approx(-p.theta*expected/t).epsilon(1e-13));
        }
    }
}

TEST_CASE("IALMob hot-state derivatives include moving screening clamp", "[mobility][ialmob][temperature][jacobian]")
{
    for (bool electron : {false,true}) {
        auto p=IalMobility::siliconDefaults(electron);
        p.alpha1Inv=.7; p.alpha2Inv=-.4; p.alpha1Acc=.3; p.alpha2Acc=.8;
        p.lCrit=p.lCritC=1e-6;
        const IalMobility model(p,electron);
        // Dilute and strongly screened states cover either side of Pmin.
        for (Real density : {1e18,1e23,1e28})
        for (Real field : {0.,1e6,1e8})
        for (Real t : {100.,299.,300.,401.,514.,600.}) {
            IalMobilityState state{density,.8*density,.3*density,.7*density,field,2e-8,t};
            const auto r=model.evaluateWithDerivatives(state);
            for (Real dt : {1e-2,1e-3}) {
                auto a=state,b=state; a.temperature_K+=dt; b.temperature_K-=dt;
                const Real fd=(model.evaluate(a).mobility_m2_per_Vs-model.evaluate(b).mobility_m2_per_Vs)/(2.*dt);
                INFO("electron="<<electron<<" T="<<t<<" density="<<density<<" field="<<field);
                CHECK(r.temperatureDerivative_m2_per_Vs_K*t/r.result.mobility_m2_per_Vs==
                    Catch::Approx(fd*t/r.result.mobility_m2_per_Vs).epsilon(2e-6).margin(2e-8));
            }
        }
    }
}

TEST_CASE("IALMob rejects unsupported temperature and nonfinite thermal exponents", "[mobility][ialmob][temperature]")
{
    auto p=IalMobility::siliconDefaults(true);
    const IalMobility model(p,true);
    IalMobilityState state;
    for (Real t : {0.,49.,-300.,std::numeric_limits<Real>::infinity(),std::numeric_limits<Real>::quiet_NaN()}) {
        state.temperature_K=t;
        CHECK_THROWS_AS(model.evaluate(state),std::invalid_argument);
    }
    for (auto member : {&IalMobilityParameters::theta,&IalMobilityParameters::k,&IalMobilityParameters::alpha1Inv}) {
        auto bad=p; bad.*member=std::numeric_limits<Real>::quiet_NaN();
        CHECK_THROWS_AS(IalMobility(bad,true),std::invalid_argument);
    }
}
