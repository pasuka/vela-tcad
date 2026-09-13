#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/physics/IalHighFieldMobility.h"
#include "vela/physics/IalMobility.h"
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace vela;
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
