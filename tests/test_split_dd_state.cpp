#include <catch2/catch_test_macros.hpp>
#include "vela/equation/SplitDDOperator.h"
#include "vela/physics/MobilityModel.h"

using namespace vela::split_dd;
static nlohmann::json parameters() {
    return {{"electron",{.1417,.00522,2.285,9.68e22,.68}},
        {"hole",{.04705,.00449,2.247,2.23e23,.719}},
        {"common",{4e26,7.2e26,.21,.5,1.,1.258,2.459,3.828,1.,1.,.89233,.41372,.005978,.28227,.72169,.19778,1.80618,1e-6,1e4}}};
}
TEST_CASE("Split DD checkpoint keeps every block and refuses incompatible metadata") {
    State s({{1,0},{2,0},{-1,0}}, {0.05},{0.05},.02585,"mesh-a");
    auto t=s.shifted({0x1p-60,-0x1p-61,0x1p-62},1.);
    for(int b=0;b<3;++b)REQUIRE(t.coordinate(b)!=s.coordinate(b));
    auto loaded=State::restore(nlohmann::json::parse(t.checkpoint().dump()),"mesh-a");
    for(int b=0;b<3;++b)REQUIRE(loaded.potential(b,0)==t.potential(b,0));
    REQUIRE_THROWS_AS(State::restore(t.checkpoint(),"mesh-b"),std::invalid_argument);
    auto malformed=t.checkpoint();malformed["coordinates"][1]=nlohmann::json::array({2.});
    REQUIRE_THROWS_AS(State::restore(malformed,"mesh-a"),std::invalid_argument);
    auto inverse=t.shifted({-0x1p-60,0x1p-61,-0x1p-62},1.);
    for(int b=0;b<3;++b)REQUIRE(inverse.coordinate(b)==s.coordinate(b));
    State tiny({{1.,0},{1e-25,0},{1e-25,0}},{1.},{1.},.02585,"mesh-a");
    auto changed=tiny.shifted({0.,1e-44,-1e-44},1.);
    REQUIRE(changed.potential(1,0)!=tiny.potential(1,0));
    REQUIRE(changed.potential(2,0)!=tiny.potential(2,0));
}
TEST_CASE("Wide PhuMob retains the calibrated double constitutive model") {
    PhuMob m(parameters());vela::PhuMobParameters p;
    for(double doping:{1e19,1e23,1e26})for(double n:{1e10,1e20,1e25})for(int b=0;b<2;++b) {
        auto expected=vela::evaluatePhuMobScalar(b?vela::CarrierType::Hole:vela::CarrierType::Electron,{doping,doping*.1,n,1e14,300.},p).mobility;
        Wide result=m.mobility(b,Wide(doping),Wide(doping*.1),Wide(n),Wide(1e14));
        REQUIRE(abs(result/Wide(expected)-1)<Wide("2e-13"));
    }
}
TEST_CASE("Split DD operator and ports are invariant under exact repartition") {
    nlohmann::json d={{"schema","vela.split-dd-operator.v1"},{"mesh_fingerprint","mesh-a"},
        {"temperature_K",300.},{"potential_scale",.02585},{"phumob",parameters()},
        {"poisson_edges",{{0,1,2e-10}}},{"current_factor",1.602176634e-19},
        {"boundary",nlohmann::json::array({nlohmann::json::array(),nlohmann::json::array(),nlohmann::json::array()})}};
    for(int i=0;i<2;++i) {
        d["nodes"].push_back({0.,0.,0.,0.,0.,0.,.02585,0.,0.,1e16,.02585,1e20,1e-21,1e-21,1e-21,1.602176634e-19,1.,1.,0.});
        d["node_physics"].push_back({1e20,0.,1e-5,3e-6,1e-21});
    }
    d["transport_edges"].push_back({{"i",0},{"j",1},{"coefficient",.01},
        {"weights",{{0,.5},{1,.5}}},{"ports",{{"drain",1.},{"source",-1.}}}});
    Operator op(d);State s({{1.,0},{2.,0},{.1,0},{.2,0},{.11,0},{.23,0}},{0.,0.},{0.,0.},.02585,"mesh-a");
    auto wrong=s.checkpoint();wrong["potential_scale_V"]=.026;
    REQUIRE_THROWS_AS(op.evaluate(State::restore(wrong,"mesh-a")),std::invalid_argument);
    wrong=s.checkpoint();wrong["electron_reference_V"][0]=.05;
    REQUIRE_THROWS_AS(op.evaluate(State::restore(wrong,"mesh-a")),std::invalid_argument);
    wrong=s.checkpoint();wrong["hole_reference_V"][1]=.05;
    REQUIRE_THROWS_AS(op.evaluate(State::restore(wrong,"mesh-a")),std::invalid_argument);
    auto a=op.evaluate(s);auto j=s.checkpoint();
    for(auto& pair:j["coordinates"]) {double old=pair[0];double hi=std::nextafter(old,INFINITY);pair={hi,old-hi};}
    auto b=op.evaluate(State::restore(j,"mesh-a"));
    REQUIRE(a.residual==b.residual);REQUIRE(a.mun==b.mun);REQUIRE(a.mup==b.mup);
    REQUIRE(a.n==b.n);REQUIRE(a.p==b.p);REQUIRE(a.srh==b.srh);REQUIRE(a.currents==b.currents);
    REQUIRE(a.currents.at("drain")==-a.currents.at("source"));
    REQUIRE(a.currents.at("drain")<0); // node 0 contact, both positive conventional-current contributions leave node 1
    auto t=op.evaluate(s.shifted({0x1p-60,0.,0x1p-60,0.,-0x1p-60,0.},1.));
    REQUIRE(t.n[0]==a.n[0]); // common psi/phin shift leaves electron population fixed
    REQUIRE(t.p[0]!=a.p[0]);REQUIRE(t.mun[0]!=a.mun[0]);REQUIRE(t.srh[0]!=a.srh[0]);
    REQUIRE(t.currents.at("drain")!=a.currents.at("drain"));
}
