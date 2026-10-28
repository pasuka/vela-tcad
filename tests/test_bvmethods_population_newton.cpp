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
TEST_CASE("Population relative columns preserve minority tangents below the old overflow threshold", "[population_newton]") {
    const double vt=.025852;
    VectorXd psi=VectorXd::Constant(1,-.1713951393747497);
    VectorXd n=VectorXd::Constant(1,6.981265876883068e-301),p=VectorXd::Constant(1,2.65e13);
    for(int step=0;step<5;++step) {
        const auto c=diagnostics::populationCoordinates(psi,n,p,n/vt,p/vt,{true},vt,vt,1.05);
        REQUIRE(c.electronScale_cm3[0]==n[0]);
        REQUIRE(c.qfFromScaledPopulation.coeff(1,1)==Approx(-1.));
        REQUIRE(c.qfFromScaledPopulation.coeff(2,2)==Approx(1.));
        VectorXd z=VectorXd::Zero(4);z[1]=-1210.4320286817224;
        const VectorXd dx=c.qfFromScaledPopulation*z;
        REQUIRE(dx.allFinite());
        const double dc=(n[0]/vt)*(dx[0]-dx[1])*vt;
        REQUIRE(dc/n[0]==Approx(z[1]).epsilon(1e-13));
        n[0]=diagnostics::localNonpositivePopulation(n[0],dc,1.).value;
    }
    REQUIRE(n[0]>0.);
}
TEST_CASE("Population coordinate validation rejects unusable restart slopes", "[population_newton]") {
    VectorXd psi=VectorXd::Zero(1),n=VectorXd::Ones(1),p=n,slope=n;
    slope[0]=std::numeric_limits<double>::infinity();
    REQUIRE_THROWS_AS(diagnostics::populationCoordinates(psi,n,p,slope,p,{true},1.,.025,0.),std::domain_error);
    slope[0]=0.;
    REQUIRE_THROWS_AS(diagnostics::populationCoordinates(psi,n,p,slope,p,{true},1.,.025,0.),std::domain_error);
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
TEST_CASE("Population local protection preserves unrelated finite updates", "[population_newton]") {
    VectorXd c(4),dc(4);c<<100.,4.,100.,1e-200;dc<<-120.,3.,-99.5,-1e-200;
    std::vector<diagnostics::LocalPopulationCandidate> candidates;
    for(int i=0;i<c.size();++i)candidates.push_back(diagnostics::localNonpositivePopulation(c[i],dc[i],1.));
    REQUIRE(candidates[0].protectedNonpositive);
    REQUIRE(candidates[0].value==Approx(.01*c[0]));
    REQUIRE_FALSE(candidates[1].protectedNonpositive);
    REQUIRE(candidates[1].value==c[1]+dc[1]);
    REQUIRE_FALSE(candidates[2].protectedNonpositive);
    REQUIRE(candidates[2].value==Approx(.005*c[2]));
    REQUIRE(candidates[3].protectedNonpositive);
    REQUIRE(candidates[3].value/c[3]==Approx(.01));
    for(const auto& a:candidates)REQUIRE(a.value>0.);
    const auto feasible=diagnostics::localNonpositivePopulation(100.,-120.,.5);
    REQUIRE_FALSE(feasible.protectedNonpositive);
    REQUIRE(feasible.value==40.);
    REQUIRE_THROWS_AS(diagnostics::localNonpositivePopulation(0.,1.,1.),std::domain_error);
    REQUIRE_THROWS_AS(diagnostics::localNonpositivePopulation(1.,std::numeric_limits<double>::infinity(),1.),std::domain_error);
    REQUIRE_THROWS_AS(diagnostics::localNonpositivePopulation(std::numeric_limits<double>::denorm_min(),-1.,1.),std::domain_error);
}
TEST_CASE("Population candidate selection uses physical quality and original joint gates", "[population_newton]") {
    using diagnostics::PopulationChoice;
    using diagnostics::selectPopulationCandidate;
    REQUIRE(selectPopulationCandidate({true,false,10.},{true,false,1.})==PopulationChoice::Local);
    REQUIRE(selectPopulationCandidate({true,false,1.},{true,false,10.})==PopulationChoice::Fraction);
    REQUIRE(selectPopulationCandidate({true,true,10.},{true,false,1.})==PopulationChoice::Fraction);
    REQUIRE(selectPopulationCandidate({true,false,1.},{true,true,10.})==PopulationChoice::Local);
    REQUIRE(selectPopulationCandidate({false,true,0.},{true,false,2.})==PopulationChoice::Local);
    REQUIRE(selectPopulationCandidate({true,false,2.},{true,false,2.})==PopulationChoice::Fraction);
    REQUIRE(selectPopulationCandidate({true,false,std::numeric_limits<double>::infinity()},{})==PopulationChoice::None);
    REQUIRE(selectPopulationCandidate({},{true,false,std::numeric_limits<double>::quiet_NaN()})==PopulationChoice::None);
    REQUIRE(selectPopulationCandidate({true,true,std::numeric_limits<double>::infinity()},{true,false,1.})==PopulationChoice::Fraction);
}
TEST_CASE("Population target merit preserves tolerance and unit scaling", "[population_newton]") {
    const auto merit=diagnostics::populationTargetMerit;
    REQUIRE(merit(2e-8,5e-7,1e-8,1e-6)==Approx(2.));
    REQUIRE(merit(2e-5,5e-4,1e-5,1e-3)==Approx(2.));
    REQUIRE(merit(1e-8,-3e-6,1e-8,1e-6)==Approx(3.));
    REQUIRE(merit(1e-8,1e-6,1e-8,1e-6)==Approx(1.));
    REQUIRE(std::isinf(merit(std::numeric_limits<double>::quiet_NaN(),0.,1e-8,1e-6)));
    REQUIRE_THROWS_AS(merit(1.,1.,0.,1.),std::invalid_argument);
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
