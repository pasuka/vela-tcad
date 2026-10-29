#include <catch2/catch_test_macros.hpp>
#include "vela/solver/SlotPopulationHold.h"
#include "vela/solver/SlotPopulationNumerics.h"

TEST_CASE("SLOT low density hold requires repeated nonpositive proposals", "[slot][population]") {
 SlotPopulationHoldState s;
 s.observe(1e-9,1e-6,-2.,2);REQUIRE_FALSE(s.held);
 s.observe(1e-9,1e-6,-.5,2);REQUIRE(s.consecutiveNonpositive==0);
 s.observe(1e-9,1e-6,-1.,2);REQUIRE_FALSE(s.held);
 s.observe(1e-9,1e-6,-2.,2);REQUIRE(s.held);REQUIRE(s.activations==1);
 s.observe(1e-9,1e-6,0.,2);REQUIRE(s.held);
 s.observe(1e-9,1e-6,.2,2);REQUIRE_FALSE(s.held);REQUIRE(s.releases==1);
 s.observe(1e-9,1e-6,-2.,2);REQUIRE_FALSE(s.held);
}
TEST_CASE("SLOT hold threshold is not a population floor", "[slot][population]") {
 SlotPopulationHoldState high;
 for(int i=0;i<4;++i)high.observe(1.,1e-6,-2.,2);
 REQUIRE_FALSE(high.held);
 SlotPopulationHoldState recovering;
 for(int i=0;i<4;++i)recovering.observe(1e-100,1e-6,.5,2);
 REQUIRE_FALSE(recovering.held);
}
TEST_CASE("SLOT concentration constraints preserve psi coupling and neighboring equations", "[slot][population]") {
 using namespace vela;
 SparseMatrixd a(6,6);
 std::vector<Eigen::Triplet<Real>> entries{{0,0,2},{1,1,2},{2,2,3},{3,3,3},{4,4,4},{4,0,-2},{5,5,4},{5,4,-1}};
 a.setFromTriplets(entries.begin(),entries.end());
 const auto constrained=slot_population::populationConstraintMatrix(a,{4},2);
 REQUIRE(constrained.coeff(4,4)==1.);
 REQUIRE(constrained.coeff(4,0)==-1.);
 REQUIRE(constrained.coeff(5,4)==a.coeff(5,4));
 VectorXd d=VectorXd::Zero(6);d[0]=.7;d[4]=.7;
 REQUIRE((constrained*d)[4]==0.);
 REQUIRE((a*d)[4]!=0.); // Real equation remains independently measurable.
}
