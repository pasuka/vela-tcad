#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "vela/solver/NewtonSolver.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "../src/solver/PopulationStepControl.h"
#include <nlohmann/json.hpp>
using namespace vela;
using Catch::Approx;

TEST_CASE("Global population policy is explicit and validates its fraction", "[population_step]") {
    CHECK(newtonConfigFromJson(nlohmann::json::object()).populationStepFraction==0.);
    CHECK(newtonConfigFromJson({{"population_step_control",{{"mode","legacy"}}}}).populationStepFraction==0.);
    CHECK(newtonConfigFromJson({{"population_step_control",{{"mode","fraction_to_boundary"},{"fraction",.99}}}}).populationStepFraction==.99);
    for(double t:{0.,1.,-1.,2.})CHECK_THROWS(newtonConfigFromJson({{"population_step_control",{{"mode","fraction_to_boundary"},{"fraction",t}}}}));
    CHECK_THROWS(newtonConfigFromJson({{"population_step_control",{{"mode","clip"}}}}));
}
TEST_CASE("Relative global boundary is unit invariant and preserves the whole direction", "[population_step]") {
    VectorXd rn(3),rp(3);rn<<-2.,1.,-1e5;rp<<1.,-4.,-1e6;
    const std::vector<bool> active{true,true,false};
    auto f=population_step::relativeFractionToBoundary(rn,rp,active,.99);
    REQUIRE(f.alpha==Approx(.2475));CHECK(f.limitingNode==1);CHECK(f.limitingCarrier==1);
    VectorXd direction=VectorXd::LinSpaced(10,-.3,.7),scaled=f.alpha*direction;
    // The final component represents an independent external contact voltage.
    for(int i=0;i<10;++i)CHECK(scaled[i]==f.alpha*direction[i]);
    for(double c:{1e-300,1.,1e200}) {
        const VectorXd n=VectorXd::Constant(3,c),p=n,dn=c*rn,dp=c*rp;
        CHECK(population_step::populationFractionToBoundary(n,p,dn,dp,active,.99).alpha==Approx(f.alpha));
    }
    rn[0]=std::numeric_limits<double>::infinity();
    CHECK_THROWS(population_step::relativeFractionToBoundary(rn,rp,active,.99));
}
TEST_CASE("Finite global population mapping follows the physical statistics and boundary mask", "[population_step]") {
    DeviceMesh mesh;
    for(Index i=0;i<3;++i){Node n;n.id=i;n.x=i==1?1e-6:0.;n.y=i==2?1e-6:0.;mesh.addNode(n);}
    Cell cell;cell.id=0;cell.type=CellType::Tri3;cell.region_id=0;cell.node_ids={0,1,2};mesh.addCell(cell);
    Region region;region.id=0;region.name="silicon";region.material="Si";region.cell_ids={0};mesh.addRegion(region);mesh.buildEdges();
    MaterialDatabase materials;DopingModel doping(3);CarrierStatisticsConfig stats;
    SECTION("Boltzmann") {stats.model="boltzmann";}
    SECTION("Fermi") {stats.model="fermi_dirac";}
    const double vt=.025851995266484913;const bool scaled=GENERATE(false,true);
    DDScalingSpec scaling;scaling.enabled=scaled;scaling.V0=vt;scaling.C0=1e21;
    CoupledDDAssembler a(mesh,materials,doping,vt,mobilityModelConfig("constant"),recombinationModelConfig({"none"}),{},{},{},{},scaling,{},stats);
    CoupledDDState state;state.psi=VectorXd::Constant(3,.02);state.phin=VectorXd::Constant(3,.01);state.phip=VectorXd::Constant(3,-.03);
    if(stats.model=="fermi_dirac")for(int i=0;i<3;++i) {
        state.phin[i]=state.psi[i]-vt*(3.-std::log(a.electronReferencePrefactor()[i]/a.electronDensityOfStates()[i]));
        state.phip[i]=state.psi[i]+vt*(2.5-std::log(a.holeReferencePrefactor()[i]/a.holeDensityOfStates()[i]));
    }
    if(scaled){state.psi/=vt;state.phin/=vt;state.phip/=vt;}
    VectorXd x=a.pack(state),d=VectorXd::Constant(9,.04);
    d.segment(3,3).array()+=5*vt;d.tail(3).array()+=.3*vt;
    if(scaled)d/=vt;
    CoupledDDBoundaryConditions bc;bc.phin[0]=0.;bc.phip[0]=0.;
    const auto f=population_step::frame(a,x,d,bc);
    REQUIRE(f.relativeN[0]==0.);REQUIRE(f.relativeP[0]==0.);
    const auto bound=population_step::relativeFractionToBoundary(f.relativeN,f.relativeP,f.active,.99);
    for(double alpha:{bound.alpha,bound.alpha*.5}) {
        const VectorXd mapped=population_step::finiteStep(a,f,d,alpha);
        CHECK(mapped[3]==alpha*d[3]);CHECK(mapped[6]==alpha*d[6]);
        const VectorXd xn=x+mapped,n=a.electronDensity(xn),p=a.holeDensity(xn);
        for(int i=1;i<3;++i){
            CHECK(n[i]>0.);CHECK(p[i]>0.);
            CHECK(n[i]/f.n[i]==Approx(1+alpha*f.relativeN[i]).epsilon(2e-10));
            CHECK(p[i]/f.p[i]==Approx(1+alpha*f.relativeP[i]).epsilon(2e-10));
        }
    }
}
TEST_CASE("Population logarithmic remainder retains tiny represented corrections", "[population_step]") {
    const double t=1e-15;
    REQUIRE(std::log1p(t)-t==Approx(-.5*t*t).margin(1e-30));
    CHECK(population_step::log1pRemainder(t)!=0.);
    CHECK(population_step::log1pRemainder(t)/(-.5*t*t)==Approx(1.).epsilon(1e-14));
}
