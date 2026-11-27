#include <catch2/catch_test_macros.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "vela/solver/SlotPopulationNumerics.h"
#include <bit>

using namespace vela;
using namespace vela::slot_population;
using WidePopulation = boost::multiprecision::cpp_dec_float_100;

TEST_CASE("Population coordinates survive an unrepresentable Fermi ratio", "[slot_population_coordinates]") {
    using W=WidePopulation;
    for(auto [n,dos]: {std::pair{1.368160756855898e-305,3.104628932474838e19},
                       std::pair{3.5192399979968014e-305,3.104628932474838e19},
                       std::pair{1e-300,1e8},std::pair{1e-280,1e28},
                       std::pair{std::numeric_limits<double>::denorm_min(),1e300}}) {
        CAPTURE(n,dos);
        const double eta=populationEta(n,dos);
        REQUIRE(std::isfinite(eta));
        CHECK(abs(exp(W(eta))*W(dos)/W(n)-1)<W("1e-12"));
        CHECK(populationLogSlope(eta)==1.);
        for(double rel: {0.,-.9,-.999,std::nextafter(-1.,0.),1e-12,1e5}) {
            const double d=relativeInversePopulation(n,dos,rel,true);
            CHECK(std::isfinite(d));
            // Compare in extended precision: the target need not fit in double.
            CHECK(abs(exp(W(d))/(W(1)+W(rel))-1)<W("1e-13"));
            if(rel==0.)CHECK(d==0.);
        }
    }
    CHECK_THROWS(inverseFermiDiracHalf(1.368160756855898e-305/3.104628932474838e19));
    CHECK(std::isnan(fermiDiracHalfDerivative(-747.)/fermiDiracHalf(-747.)));
}

TEST_CASE("Population coordinate repair preserves ordinary arithmetic and full Fermi endpoint", "[slot_population_coordinates]") {
    for(double z: {1e-280,1e-13,1e-5,1.,20.}) {
        const double dos=1e19,n=z*dos,ratio=n/dos;
        const double eta=populationEta(n,dos);
        CHECK(std::bit_cast<std::uint64_t>(eta)==std::bit_cast<std::uint64_t>(inverseFermiDiracHalf(ratio)));
        CHECK(populationLogSlope(eta)==fermiDiracHalfDerivative(eta)/fermiDiracHalf(eta));
        for(bool repair:{false,true})for(double rel:{-.5,.2}) {
            const double old=repair?relativeInverseDiluteLog(ratio,rel):relativeInverse(ratio,rel);
            CHECK(relativeInversePopulation(n,dos,rel,repair)==old);
            CHECK(usesDiluteLogInversePopulation(n,dos,rel)==usesDiluteLogInverse(ratio,rel));
        }
    }
    // Positive recovery can leave the dilute domain; exp-only inversion is wrong.
    const double n=1e-300,dos=1e8,rel=1e308;
    const double eta=populationEta(n,dos),d=relativeInversePopulation(n,dos,rel,true);
    const WidePopulation target=WidePopulation(n)/WidePopulation(dos)*(1+WidePopulation(rel));
    CHECK(abs(WidePopulation(fermiDiracHalf(eta+d))/target-1)<WidePopulation("2e-12"));
    CHECK_FALSE(usesDiluteLogInversePopulation(n,dos,rel));
}

TEST_CASE("Stable population coordinates retain invalid input protections", "[slot_population_coordinates]") {
    const double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
    for(double bad:{0.,-1.,nan,inf}) {
        CHECK_THROWS(populationEta(bad,1.));
        CHECK_THROWS(populationEta(1.,bad));
        CHECK_THROWS(relativeInversePopulation(bad,1.,-.5,true));
        CHECK_THROWS(relativeInversePopulation(1.,bad,-.5,true));
    }
    for(double bad:{-1.,-2.,nan,inf})CHECK_THROWS(relativeInversePopulation(1e-305,1e19,bad,true));
    CHECK_THROWS(populationLogSlope(nan));
    CHECK_THROWS(populationLogSlope(inf));
}
