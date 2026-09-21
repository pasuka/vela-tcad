#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/physics/IalHighFieldMobility.h"
#include "vela/physics/IalMobility.h"
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace vela;
TEST_CASE("Generated zero surface offset keeps selected power slopes finite", "[ialmob][generated][jacobian]") {
    for(bool electron:{false,true})for(Real field:{0.,1e6}) {
        auto p=IalMobility::siliconDefaults(electron);
        p.N2=0.;p.lambda=p.lambdaSr=0.;
        IalMobility model(p,electron);
        const IalMobilityState state{0.,0.,1e20,1e20,field,5e-9,400.};
        const auto a=model.evaluateWithDerivatives(state),b=model.evaluateWithDerivatives(state,nullptr,nullptr,true);
        CHECK(b.result.mobility_m2_per_Vs==Catch::Approx(a.result.mobility_m2_per_Vs).epsilon(2e-12));
        for(int k=0;k<6;++k)CHECK(b.derivative_SI[k]==Catch::Approx(a.derivative_SI[k]).epsilon(2e-10).margin(1e-30));
        CHECK(b.temperatureDerivative_m2_per_Vs_K==Catch::Approx(a.temperatureDerivative_m2_per_Vs_K).epsilon(2e-10));
    }
}
TEST_CASE("Generated screening envelope agrees across its moving clamp", "[ialmob][generated][jacobian]") {
    for(bool electron:{false,true})for(Real temperature:{200.,300.,514.,1000.}) {
        const auto p=IalMobility::siliconDefaults(electron);
        IalMobility model(p,electron);IalScreeningCache roots;
        const Real minimum=roots.minimum(p.mass,temperature),t=temperature/300.;
        IalMobilityState state{1e20,2e19,1e14,1e14,1e6,5e-9,temperature};
        const auto star=[](Real n,Real ref,Real c){const Real x=n/ref;return n*(1+x*x/(1+c*x*x));};
        const Real other=1e8,nsc=star(1e14,p.nRefD,p.cRefD)+star(2e13,p.nRefA,p.cRefA)+other;
        const Real boundary=(t*t/minimum-(2.459/3.97e13)*std::pow(nsc,2./3.))*(p.mass*1.36e20/3.828)-other;
        REQUIRE(boundary>0.);
        for(Real factor:{1.-1e-6,1.,1.+1e-6})for(Real distance:{0.,5e-9,1.}) {
            (electron?state.electrons_m3:state.holes_m3)=boundary*1e6*factor;
            state.interfaceDistance_m=distance;
            const auto a=model.evaluateWithDerivatives(state),b=model.evaluateWithDerivatives(state,nullptr,nullptr,true);
            CHECK(b.result.mobility_m2_per_Vs==Catch::Approx(a.result.mobility_m2_per_Vs).epsilon(2e-12));
            CHECK(model.evaluate(state,nullptr,nullptr,true).mobility_m2_per_Vs==b.result.mobility_m2_per_Vs);
            CHECK(b.temperatureDerivative_m2_per_Vs_K*temperature/b.result.mobility_m2_per_Vs==
                Catch::Approx(a.temperatureDerivative_m2_per_Vs_K*temperature/a.result.mobility_m2_per_Vs).epsilon(2e-10).margin(2e-12));
            auto plus=state,minus=state;const Real step=temperature*1e-5;
            plus.temperature_K+=step;minus.temperature_K-=step;
            const Real fd=(model.evaluate(plus,nullptr,nullptr,true).mobility_m2_per_Vs-
                model.evaluate(minus,nullptr,nullptr,true).mobility_m2_per_Vs)/(2*step);
            CHECK(b.temperatureDerivative_m2_per_Vs_K*temperature/b.result.mobility_m2_per_Vs==
                Catch::Approx(fd*temperature/b.result.mobility_m2_per_Vs).epsilon(2e-5).margin(2e-7));
        }
    }
}
TEST_CASE("Generated low-field graph matches AD values slopes and zero branches", "[ialmob][generated][jacobian]") {
    const std::array<Real IalMobilityState::*,7> inputs{&IalMobilityState::donors_m3,&IalMobilityState::acceptors_m3,
        &IalMobilityState::electrons_m3,&IalMobilityState::holes_m3,&IalMobilityState::normalField_V_per_m,
        &IalMobilityState::interfaceDistance_m,&IalMobilityState::temperature_K};
    const std::array<Real IalMobilityResult::*,6> components{&IalMobilityResult::mobility_m2_per_Vs,
        &IalMobilityResult::coulomb3d_m2_per_Vs,&IalMobilityResult::coulomb2d_m2_per_Vs,
        &IalMobilityResult::coulomb_m2_per_Vs,&IalMobilityResult::phonon_m2_per_Vs,&IalMobilityResult::roughness_m2_per_Vs};
    for(bool electron:{false,true})for(bool special:{false,true}) {
        auto params=IalMobility::siliconDefaults(electron);
        // Density in this exponent is in cm^-3: keep the nonzero perturbation
        // representable while exercising its carrier and temperature chains.
        if(special){params.alphaSr=1e-15;params.nu=.2;params.alpha1Inv=.5;params.alpha2Acc=-.3;params.lCrit=1e-6;params.lCritC=1e-6;}
        IalMobility model(params,electron);IalMobilityPreparationCache cache;
        for(Real t:{50.,200.,300.,514.,1000.})for(Real field:{0.,1e-8,1e2,1e6,1e10})for(int mode=0;mode<5;++mode) {
            INFO("electron="<<electron<<" special="<<special<<" T="<<t<<" field="<<field<<" mode="<<mode);
            IalMobilityState state{1e23,2e22,3e23,1e19,field,5e-9,t};
            if(mode==1)state.donors_m3=0.;if(mode==2)state.acceptors_m3=0.;
            if(mode==3)state.donors_m3=state.acceptors_m3=state.electrons_m3=state.holes_m3=0.;
            if(mode==4)state.electrons_m3=state.holes_m3=0.;
            INFO("reference evaluation");
            const auto a=model.evaluateWithDerivatives(state);
            INFO("generated evaluation");
            const auto b=model.evaluateWithDerivatives(state,nullptr,&cache,true);
            const auto value=model.evaluate(state,nullptr,&cache,true);
            for(auto component:components) {
                if(std::isinf(a.result.*component))CHECK(std::isinf(b.result.*component));
                else CHECK(b.result.*component==Catch::Approx(a.result.*component).epsilon(2e-12));
                CHECK(value.*component==b.result.*component);
            }
            for(int k=0;k<7;++k){
                const Real da=k==6?a.temperatureDerivative_m2_per_Vs_K:a.derivative_SI[k];
                const Real db=k==6?b.temperatureDerivative_m2_per_Vs_K:b.derivative_SI[k];
                const Real scale=std::max(state.*inputs[k],std::array<Real,7>{1e20,1e20,1e20,1e20,1.,1e-10,50.}[k])/a.result.mobility_m2_per_Vs;
                CHECK(db*scale==Catch::Approx(da*scale).epsilon(2e-10).margin(2e-12));
                if(mode==0 && field>=1e2 && t>50.){
                    auto plus=state,minus=state;const Real step=state.*inputs[k]*1e-5;
                    plus.*inputs[k]+=step;minus.*inputs[k]-=step;
                    const Real fd=(model.evaluate(plus,nullptr,nullptr,true).mobility_m2_per_Vs-model.evaluate(minus,nullptr,nullptr,true).mobility_m2_per_Vs)/(2*step);
                    CHECK(db*scale==Catch::Approx(fd*scale).epsilon(2e-5).margin(2e-7));
                }
            }
            // Same cache can safely alternate candidate and reference kernels.
            const auto again=model.evaluateWithDerivatives(state,nullptr,&cache);
            CHECK(again.result.mobility_m2_per_Vs==a.result.mobility_m2_per_Vs);
            CHECK(again.derivative_SI==a.derivative_SI);
        }
    }
}
TEST_CASE("Explicit IALMob HFS partials preserve reference and independent perturbations", "[ialmob][temperature][hfs][explicit]") {
    for(auto p:{IalHighFieldParameters{},IalHighFieldParameters{8.37e4,1.213,.52,.17}})
    for(Real t:{50.,200.,300.,514.,1000.})for(Real m:{1e-5,.03,1.})
    for(Real e:{0.,1e-100,1e-10,1e2,1e6,1e9,1e14}) {
        const auto a=evaluateIalHighFieldMobility(m,e,t,p),b=evaluateIalHighFieldMobilityExplicit(m,e,t,p);
        CHECK(b.mobility_m2_per_Vs==a.mobility_m2_per_Vs);
        CHECK(evaluateIalHighFieldMobilityValue(m,e,t,p)==b.mobility_m2_per_Vs);
        CHECK(b.lowFieldDerivative*m/b.mobility_m2_per_Vs==Catch::Approx(a.lowFieldDerivative*m/a.mobility_m2_per_Vs).epsilon(1e-11).margin(1e-13));
        CHECK(b.drivingFieldDerivative_m3_per_V2s*e/b.mobility_m2_per_Vs==Catch::Approx(a.drivingFieldDerivative_m3_per_V2s*e/a.mobility_m2_per_Vs).epsilon(1e-11).margin(1e-13));
        CHECK(b.temperatureDerivative_m2_per_Vs_K*t/b.mobility_m2_per_Vs==Catch::Approx(a.temperatureDerivative_m2_per_Vs_K*t/a.mobility_m2_per_Vs).epsilon(1e-11).margin(1e-13));
        if(e>=1e2 && t>50.)for(int k=0;k<3;++k) {
            std::array<Real,3> u{m,e,t},v=u;const Real step=u[k]*1e-5;u[k]+=step;v[k]-=step;
            const Real fd=(evaluateIalHighFieldMobilityValue(u[0],u[1],u[2],p)-evaluateIalHighFieldMobilityValue(v[0],v[1],v[2],p))/(2.*step);
            const Real slope=std::array{b.lowFieldDerivative,b.drivingFieldDerivative_m3_per_V2s,b.temperatureDerivative_m2_per_Vs_K}[k];
            CHECK(slope*std::array{m,e,t}[k]/b.mobility_m2_per_Vs==Catch::Approx(fd*std::array{m,e,t}[k]/b.mobility_m2_per_Vs).epsilon(2e-6).margin(1e-8));
        }
    }
    CHECK_THROWS(evaluateIalHighFieldMobilityExplicit(.03,-1.,300.,{}));
    CHECK_THROWS(evaluateIalHighFieldMobilityValue(.03,1.,49.,{}));
}
TEST_CASE("IALMob hot HFS preserves zero drive and saturation velocity", "[ialmob][temperature][hfs]")
{
    for (auto p : {IalHighFieldParameters{},IalHighFieldParameters{8.37e4,1.213,.52,.17}})
    for (Real t : {200.,300.,401.,514.,600.}) {
        CHECK(evaluateIalHighFieldMobility(.03,0.,t,p).mobility_m2_per_Vs==.03);
        const auto high=evaluateIalHighFieldMobility(.03,1e12,t,p);
        CHECK(high.mobility_m2_per_Vs*1e12==Catch::Approx(p.saturationVelocity300_m_per_s*
            std::pow(t/300.,-p.saturationVelocityTemperatureExponent)).epsilon(1e-4));
        const auto actual=evaluateIalHighFieldMobility(.03,1e7,t,p);
        CHECK(actual.mobility_m2_per_Vs>0.);
        CHECK(actual.mobility_m2_per_Vs<.03);
        CHECK(actual.lowFieldDerivative>=0.);
        CHECK(actual.drivingFieldDerivative_m3_per_V2s<0.);
    }
}
TEST_CASE("IALMob hot HFS three partials match perturbations", "[ialmob][temperature][hfs][jacobian]")
{
    for (auto p : {IalHighFieldParameters{},IalHighFieldParameters{8.37e4,1.213,.52,.17}})
    for (Real t : {200.,300.,401.,514.,600.})
    for (Real field : {1e2,1e6,1e9}) {
        const auto r=evaluateIalHighFieldMobility(.03,field,t,p);
        for (Real fraction : {1e-4,1e-5}) {
            const Real dm=.03*fraction, df=field*fraction, dt=t*fraction;
            const auto value=[&](Real m,Real f,Real temp){return evaluateIalHighFieldMobility(m,f,temp,p).mobility_m2_per_Vs;};
            CHECK(r.lowFieldDerivative*.03/r.mobility_m2_per_Vs==Catch::Approx(
                (value(.03+dm,field,t)-value(.03-dm,field,t))/(2.*fraction*r.mobility_m2_per_Vs)).epsilon(2e-6).margin(1e-8));
            CHECK(r.drivingFieldDerivative_m3_per_V2s*field/r.mobility_m2_per_Vs==Catch::Approx(
                (value(.03,field+df,t)-value(.03,field-df,t))/(2.*fraction*r.mobility_m2_per_Vs)).epsilon(2e-6).margin(1e-8));
            CHECK(r.temperatureDerivative_m2_per_Vs_K*t/r.mobility_m2_per_Vs==Catch::Approx(
                (value(.03,field,t+dt)-value(.03,field,t-dt))/(2.*fraction*r.mobility_m2_per_Vs)).epsilon(2e-6).margin(1e-8));
        }
    }
}

TEST_CASE("IALMob hot HFS rejects invalid inputs", "[ialmob][temperature][hfs]")
{
    IalHighFieldParameters p;
    CHECK_THROWS_AS(evaluateIalHighFieldMobility(.03,-1.,300.,p),std::invalid_argument);
    CHECK_THROWS_AS(evaluateIalHighFieldMobility(.03,1.,49.,p),std::invalid_argument);
    CHECK_THROWS_AS(evaluateIalHighFieldMobility(0.,1.,300.,p),std::invalid_argument);
    p.betaTemperatureExponent=std::numeric_limits<Real>::quiet_NaN();
    CHECK_THROWS_AS(evaluateIalHighFieldMobility(.03,1.,300.,p),std::invalid_argument);
}
TEST_CASE("IALMob screening reuse preserves values and all derivatives with exact mass and temperature keys", "[ialmob][temperature][preparation]") {
    IalScreeningCache cache;
    const std::array<Real IalMobilityResult::*,6> fields{
        &IalMobilityResult::mobility_m2_per_Vs,&IalMobilityResult::coulomb3d_m2_per_Vs,
        &IalMobilityResult::coulomb2d_m2_per_Vs,&IalMobilityResult::coulomb_m2_per_Vs,
        &IalMobilityResult::phonon_m2_per_Vs,&IalMobilityResult::roughness_m2_per_Vs};
    for(bool electron:{true,false}) {
        IalMobility model(IalMobility::siliconDefaults(electron),electron);
        for(Real temperature:{300.,350.,350.,401.,514.,350.,std::nextafter(350.,351.)}) {
            IalMobilityState s{1e23,2e22,3e23,1e19,2e6,5e-9,temperature};
            const auto plain=model.evaluateWithDerivatives(s),cached=model.evaluateWithDerivatives(s,&cache);
            for(auto field:fields)CHECK(plain.result.*field==cached.result.*field);
            CHECK(plain.derivative_SI==cached.derivative_SI);
            CHECK(plain.temperatureDerivative_m2_per_Vs_K==cached.temperatureDerivative_m2_per_Vs_K);
        }
    }
    CHECK(cache.size()==8); // Four distinct non-300 K values for each mass.
    CHECK_THROWS_AS(cache.minimum(0.,400.),std::invalid_argument);
    CHECK_THROWS_AS(cache.minimum(1.,49.),std::invalid_argument);
    CHECK_THROWS_AS(cache.minimum(1.,std::numeric_limits<Real>::quiet_NaN()),std::invalid_argument);
}
