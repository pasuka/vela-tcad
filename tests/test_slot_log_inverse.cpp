#include <catch2/catch_test_macros.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "vela/solver/SlotPopulationNumerics.h"
#include <bit>

namespace {
using MP = boost::multiprecision::cpp_dec_float_100;
MP fermiDilute(MP eta) {
    const MP t=exp(eta); MP sum=0, power=t;
    for(int k=1;k<=16;++k) {
        sum += (k%2 ? MP(1) : MP(-1))*power/(MP(k)*sqrt(MP(k)));
        power *= t;
    }
    return sum;
}
double candidate(double z,double rel) {
#ifdef VELA_LOG_INVERSE_RED
    return vela::slot_population::relativeInverse(z,rel);
#else
    return vela::slot_population::relativeInverseDiluteLog(z,rel);
#endif
}
}

TEST_CASE("Dilute inverse preserves final population near complete depletion", "[slot_log_inverse]") {
    for(double z: {1e-300,1e-30,5.948816250887021e-13,1e-12}) {
        MP eta=log(MP(z));
        for(int i=0;i<5;++i) eta += log(MP(z)/fermiDilute(eta));
        for(double rel: {-.9,-.999,-.9999993512311891,-(1.-1e-9),std::nextafter(-1.,0.),1e-12}) {
            CAPTURE(z,rel);
            const double d=candidate(z,rel);
            const MP error=fermiDilute(eta+MP(d))/(MP(z)*(MP(1)+MP(rel)))-1;
            REQUIRE(abs(error)<MP("2e-13"));
        }
    }
}

#ifndef VELA_LOG_INVERSE_RED
TEST_CASE("Dilute inverse also handles population increases within its domain", "[slot_log_inverse]") {
    for(auto [z,rel]: {std::pair{1e-300,1e287},std::pair{1e-30,5e17},
                       std::pair{5e-13,.5},std::pair{1e-13,1e-20}}) {
        CAPTURE(z,rel);
        REQUIRE(vela::slot_population::usesDiluteLogInverse(z,rel));
        MP eta=log(MP(z));
        for(int i=0;i<5;++i)eta+=log(MP(z)/fermiDilute(eta));
        const double d=candidate(z,rel);
        REQUIRE(abs(fermiDilute(eta+MP(d))/(MP(z)*(1+MP(rel)))-1)<MP("2e-13"));
    }
}

TEST_CASE("Dilute inverse has identity validation and unchanged outside-domain fallback", "[slot_log_inverse]") {
    using namespace vela::slot_population;
    REQUIRE(relativeInverseDiluteLog(1e-13,0)==0.);
    for(double z: {1e-10,1e-5,1.,20.}) for(double rel: {-.5,.2})
        REQUIRE(std::bit_cast<std::uint64_t>(relativeInverseDiluteLog(z,rel)) ==
                std::bit_cast<std::uint64_t>(relativeInverse(z,rel)));
    REQUIRE_FALSE(usesDiluteLogInverse(1e-13,100.));
    REQUIRE_FALSE(usesDiluteLogInverse(1e-11,-.999));
    REQUIRE(usesDiluteLogInverse(1e-13,-.999));
    REQUIRE_THROWS(relativeInverseDiluteLog(0.,-.5));
    REQUIRE_THROWS(relativeInverseDiluteLog(1e-13,-1.));
    REQUIRE_THROWS(relativeInverseDiluteLog(1e-13,std::numeric_limits<double>::infinity()));
}
#endif
