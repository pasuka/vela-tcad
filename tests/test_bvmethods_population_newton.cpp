#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "../scripts/diagnostics/BvMethodsPopulationNewton.h"
using namespace vela;
using Catch::Approx;

TEST_CASE("Population Newton coordinates reproduce density tangents across scales", "[population_newton]") {
    const double vt=.025852,v0=.037;
    VectorXd psi(3),n(3),p(3);psi<<.2,-.1,0.;n<<1e18,1e-200,0.;p<<1e-180,2e16,0.;
    const VectorXd ns=n/vt,ps=p/vt;
    const auto c=diagnostics::populationCoordinates(psi,n,p,ns,ps,{true,true,false},v0,vt,.2);
    VectorXd z(10);z<<.1,-.3,.2,.7,-.2,.8,-.4,.6,-.5,.1;
    const VectorXd x=c.qfFromScaledPopulation*z;
    for(int i=0;i<2;++i){
        REQUIRE(ns[i]*(x[i]-x[3+i])*v0==Approx(c.electronScale_cm3[i]*z[3+i]).epsilon(1e-13));
        REQUIRE(ps[i]*(x[6+i]-x[i])*v0==Approx(c.holeScale_cm3[i]*z[6+i]).epsilon(1e-13));
    }
    REQUIRE(x[5]==z[5]);REQUIRE(x[8]==z[8]);
    REQUIRE(x[9]==Approx(c.voltageScale_V*z[9]));
}
TEST_CASE("Population fraction step preserves positivity without density projection", "[population_newton]") {
    VectorXd n(2),p(2),dn(2),dp(2);n<<100.,1e-200;p<<4.,2.;dn<<-120.,-2e-200;dp<<.3,-.5;
    const auto f=diagnostics::populationFractionToBoundary(n,p,dn,dp,{true,true},.99);
    REQUIRE(f.alpha==Approx(.495));REQUIRE(f.limitingNode==1);REQUIRE(f.limitingCarrier==0);
    for(int i=0;i<2;++i){REQUIRE(n[i]+f.alpha*dn[i]>=.009999*n[i]);REQUIRE(p[i]+f.alpha*dp[i]>=.009999*p[i]);}
    dn.setConstant(1.);dp.setConstant(1.);
    REQUIRE(diagnostics::populationFractionToBoundary(n,p,dn,dp,{true,true},.99).alpha==1.);
    REQUIRE_THROWS_AS(diagnostics::populationFractionToBoundary(n,p,dn,dp,{true,true},1.),std::invalid_argument);
}
TEST_CASE("Population equilibration preserves the original linear system", "[population_newton]") {
    SparseMatrixd a(3,3);std::vector<Eigen::Triplet<double>> t{{0,0,2e-200},{0,1,1e-100},{1,1,3e100},{1,2,-1e100},{2,0,1e-100},{2,2,4.}};a.setFromTriplets(t.begin(),t.end());
    const auto e=diagnostics::equilibratePopulationSystem(a);
    VectorXd z(3);z<<.2,-.7,.3;const VectorXd x=e.columnWeights.cwiseProduct(z);
    const VectorXd lhs=e.matrix*z,rhs=e.rowWeights.cwiseProduct(a*x);
    REQUIRE((lhs-rhs).norm()<1e-14*lhs.norm());
    REQUIRE((e.rowWeights.array()>0.).all());REQUIRE((e.columnWeights.array()>0.).all());
}
TEST_CASE("Population derivative restoration requires both residual and raw update improvement", "[population_newton]") {
    REQUIRE(diagnostics::restorationNoWorse(.5,.8,1.,1.));
    REQUIRE_FALSE(diagnostics::restorationNoWorse(.5,2.,1.,1.));
    REQUIRE_FALSE(diagnostics::restorationNoWorse(2.,.5,1.,1.));
    REQUIRE_FALSE(diagnostics::restorationNoWorse(std::numeric_limits<double>::quiet_NaN(),.5,1.,1.));
}
TEST_CASE("Population target-aware restoration uses the existing stopping targets", "[population_newton]") {
    REQUIRE(diagnostics::restorationWithinTargets(1e-10,1e-8,.15,1e-8,1e-6));
    REQUIRE_FALSE(diagnostics::restorationWithinTargets(2e-8,1e-8,.15,1e-8,1e-6));
    REQUIRE_FALSE(diagnostics::restorationWithinTargets(1e-10,2e-6,.15,1e-8,1e-6));
    REQUIRE_FALSE(diagnostics::restorationWithinTargets(1e-10,1e-8,1.,1e-8,1e-6));
}
TEST_CASE("Population autonomous restoration keeps converged candidates despite merit noise", "[population_newton]") {
    REQUIRE_FALSE(diagnostics::derivativeBootstrapRequired(true,1.,2e-11,1e-11,true));
    REQUIRE(diagnostics::derivativeBootstrapRequired(true,1.,2e-11,1e-11,false));
    REQUIRE(diagnostics::derivativeBootstrapRequired(false,1.,0.,1.,true));
    REQUIRE(diagnostics::derivativeBootstrapRequired(true,.01,.5,1.,false));
}
