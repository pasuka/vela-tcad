#include <catch2/catch_test_macros.hpp>
#include "simplemos_compensated_step.hpp"
#include <boost/multiprecision/cpp_dec_float.hpp>
using simplemos_compensated_step::propose;
TEST_CASE("sub-ulp potential steps accumulate without changing rejected state", "[simplemos][step]") {
    const double step=std::ldexp(1.,-58);double x=1.,carry=0.;
    REQUIRE(x+step==x);
    for(int i=0;i<1000;++i) {const auto p=propose(x,step,1.,carry);x=p.value;carry=p.remainder;}
    using MP=boost::multiprecision::cpp_dec_float_100;
    REQUIRE(abs((MP(x)+MP(carry))-(MP(1)+MP(step)*1000))<MP("1e-90"));
    const double oldX=x,oldCarry=carry;const auto rejected=propose(x,1.,1.,carry);(void)rejected;
    REQUIRE(x==oldX);REQUIRE(carry==oldCarry);
}
TEST_CASE("candidate construction retains full and backtracked arithmetic", "[simplemos][step]") {
    using MP=boost::multiprecision::cpp_dec_float_100;
    for(double x:{-31.5,0.,0.125,30.25})for(double step:{-1e-14,1e-17,0.125})for(double alpha:{1.,0.5,0.000244140625}) {
        const auto p=propose(x,step,alpha,0.);
        const MP expected=MP(x)+MP(alpha)*MP(step);
        REQUIRE(abs(MP(p.value)+MP(p.remainder)-expected)<MP("1e-70"));
    }
}
