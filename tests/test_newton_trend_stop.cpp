#include <catch2/catch_test_macros.hpp>
#include "vela/solver/detail/NewtonTrendStop.h"
#include <limits>
using vela::detail::NewtonTrendSample;
using vela::detail::NewtonTrendStop;
static NewtonTrendSample sample(int k) {
    NewtonTrendSample s; s.iteration=k; s.rows={0,0,0,100,300,1.5};
    s.gates={.5,0,0,0,0,0}; s.rawDirection=100; return s;
}
TEST_CASE("Trend policy requires warmup and repeated windows", "[trend_stop]") {
    NewtonTrendStop t;
    for(int k=0;k<14;++k) REQUIRE(t.observe(sample(k)).empty());
    REQUIRE(t.observe(sample(14))=="stagnation_detected");
}
TEST_CASE("Trend policy preserves genuine row or failing gate progress", "[trend_stop]") {
    for(int channel=0;channel<4;++channel) {
        NewtonTrendStop t;
        for(int k=0;k<40;++k) {
            auto s=sample(k);
            if(channel==0) s.rows[3]=100*std::pow(.98,k);
            if(channel==1) s.rows[4]=300*std::pow(.98,k);
            if(channel==2) s.gates[1]=100*std::pow(.98,k);
            if(channel==3) s.rawDirection=100*std::pow(.90,k);
            REQUIRE(t.observe(s).empty());
        }
    }
}
TEST_CASE("Already satisfied global residual cannot mask frozen bad rows", "[trend_stop]") {
    NewtonTrendStop t;
    for(int k=0;k<14;++k) {auto s=sample(k);s.gates[0]=std::pow(.1,k);REQUIRE(t.observe(s).empty());}
    REQUIRE(t.observe(sample(14))=="stagnation_detected");
}
TEST_CASE("Tiny steps require three nonimproving transitions and success wins", "[trend_stop]") {
    NewtonTrendStop t;
    REQUIRE(t.observe(sample(0)).empty());
    for(int k=1;k<3;++k) {auto s=sample(k);s.alpha=1./1024;REQUIRE(t.observe(s).empty());}
    auto s=sample(3);s.alpha=1./1024;
    REQUIRE(t.observe(s)=="tiny_step_stagnation");
    REQUIRE(t.observe(sample(4),true).empty());
    NewtonTrendStop good;
    for(int k=0;k<10;++k) {auto q=sample(k);q.alpha=1./1024;q.rows[4]*=std::pow(.98,k);REQUIRE(good.observe(q).empty());}
}
TEST_CASE("One improving window resets confirmation and nonfinite metrics fail", "[trend_stop]") {
    NewtonTrendStop t;
    for(int k=0;k<13;++k) REQUIRE(t.observe(sample(k)).empty());
    auto a=sample(13);a.rows[4]=200;REQUIRE(t.observe(a).empty());
    a.iteration=14;REQUIRE(t.observe(a).empty());REQUIRE(t.consecutiveWindows()==0);
    a.rawDirection=std::numeric_limits<double>::infinity();
    REQUIRE(t.observe(a)=="nonfinite_trend_metric");
}
