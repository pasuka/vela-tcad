#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/physics/IalMobility.h"
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace vela;

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
