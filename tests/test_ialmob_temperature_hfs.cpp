#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/physics/IalHighFieldMobility.h"
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
