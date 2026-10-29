#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "bvmethods_live_geometry.h"
#include "bvmethods_precision_sg.h"
#include <limits>
using Catch::Approx;
static nlohmann::json fixture() {
    return {{"volume_um2",{.1,.2,.3}}, {"stencil",{
        {{"row",0},{"column",0},{"Ex",1e5},{"Ey",1e5}},
        {{"row",0},{"column",1},{"Ex",-1e5},{"Ey",0}},
        {{"row",0},{"column",2},{"Ex",0},{"Ey",-1e5}}}}};
}
TEST_CASE("Private live geometry preserves affine field and gauge") {
    bvmethods_live::Geometry g(fixture());const auto a=g.field({3.,5.,7.});
    REQUIRE(a[0][0]==Approx(-2e5));REQUIRE(a[0][1]==Approx(-4e5));
    const auto b=g.field({103.,105.,107.});REQUIRE(a==b);
    REQUIRE(g.field({9.,9.,9.})[0][0]==0);
}
TEST_CASE("Private live source responds to actual potential") {
    bvmethods_live::Geometry g(fixture());auto rate=[](double F){return F*F;};
    const auto a=g.integrals({.1,.2,.4},rate),b=g.integrals({.1,.3,.4},rate);
    REQUIRE(a[0]!=b[0]);REQUIRE(a[1]==0);
}
TEST_CASE("Private live E2 chain matches same-direction finite differences") {
    bvmethods_live::Geometry g(fixture());const std::vector<double> psi{.1,2.,4.},direction{.02,-.03,.04};
    auto rate=[](double F){return 3.5e23*F*F*std::exp(-2.25e9/F);};
    auto partial=[&](double F){return rate(F)*(2/F+2.25e9/(F*F));};
    double analytic=0;for(const auto& [i,j,v]:g.chain(psi,partial))if(i==0)analytic+=v*direction[j];
    REQUIRE(std::abs(analytic)>1e-20);
    auto p=psi,m=psi;const double h=1e-5;for(int i=0;i<3;++i){p[i]+=h*direction[i];m[i]-=h*direction[i];}
    const double fd=(g.integrals(p,rate)[0]-g.integrals(m,rate)[0])/(2*h);
    REQUIRE(fd==Approx(analytic).epsilon(2e-7));
}
TEST_CASE("Private live chain has a constant-potential nullspace") {
    bvmethods_live::Geometry g(fixture());double sum=0,scale=0;
    for(const auto& [i,j,v]:g.chain({1.,2.,3.},[](double F){return 2*F;})){sum+=v;scale+=std::abs(v);}
    // Three floating-point contributions: bound summation roundoff relative
    // to the derivative scale, rather than an arbitrary dimensional floor.
    REQUIRE(std::abs(sum)<=4*std::numeric_limits<double>::epsilon()*scale);
}
TEST_CASE("Private live zero-field source and chain are finite") {
    bvmethods_live::Geometry g(fixture());for(const auto& [i,j,v]:g.chain({1.,1.,1.},[](double){return 0.;}))REQUIRE(v==0);
}
TEST_CASE("Private live geometry rejects invalid support") {
    auto j=fixture();j["volume_um2"][0]=-1;REQUIRE_THROWS_AS(bvmethods_live::Geometry(j),std::invalid_argument);
    j=fixture();j["stencil"].push_back(j["stencil"][0]);REQUIRE_THROWS_AS(bvmethods_live::Geometry(j),std::invalid_argument);
    j=fixture();j["stencil"][0]["column"]=9;REQUIRE_THROWS_AS(bvmethods_live::Geometry(j),std::invalid_argument);
}
TEST_CASE("Private live geometry rejects malformed states") {
    bvmethods_live::Geometry g(fixture());REQUIRE_THROWS_AS(g.field({1,2}),std::invalid_argument);
    REQUIRE_THROWS_AS(g.field({1,2,std::numeric_limits<double>::infinity()}),std::invalid_argument);
}
TEST_CASE("Private decimal SG retains a tiny resolved current and orientation") {
    const long double drop=1e-20L;
    const auto f=bvmethods_live::decimalSg(1e20,1.01e20,1.01e20,drop,.02585,1.,true);
    const auto reverse=bvmethods_live::decimalSg(1.01e20,1e20,1e20,-drop,.02585,1.,true);
    REQUIRE(std::isfinite(f));REQUIRE(f!=0);REQUIRE(double(reverse)==Approx(-double(f)).epsilon(1e-14));
}
TEST_CASE("Private decimal SG preserves zero drive and carrier signs") {
    REQUIRE(bvmethods_live::decimalSg(1e18,2e18,2e18,0,.02585,1,true)==0);
    REQUIRE(bvmethods_live::decimalSg(1e18,1e18,1e18,1e-20L,.02585,1,true)>0);
    REQUIRE(bvmethods_live::decimalSg(1e18,1e18,1e18,1e-20L,.02585,1,false)<0);
}
