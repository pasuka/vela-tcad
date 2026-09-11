#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/equation/IalTransport.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/post/ContactCurrent.h"
#include "vela/solver/GummelSolver.h"
#include "vela/core/PhysicalConstants.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <chrono>

using namespace vela;
namespace {
struct Fixture {
    DeviceMesh mesh;MaterialDatabase materials;DopingModel doping{6};
    MobilityModelConfig mobility;
    std::filesystem::path path;
    Fixture(bool tcad=false,bool interfaceCorner=false) {
        const UnitScalingConfig inputScaling{tcad?UnitScalingMode::UnitScaling:UnitScalingMode::LegacySI};
        materials=MaterialDatabase(inputScaling);
        const Real L=tcad?.1:1e-7;
        const std::array<std::array<Real,2>,6> xy{{{0,0},{L,0},{L,L},{0,L},{0,2*L},{L,2*L}}};
        for (Index i=0;i<6;++i) {Node n;n.id=i;n.x=xy[i][0];n.y=xy[i][1];mesh.addNode(n);}
        const std::array<std::array<Index,3>,4> ns{{{0,1,2},{0,2,3},{3,2,5},{3,5,4}}};
        for (Index i=0;i<4;++i) {Cell c;c.id=i;c.region_id=i<2?0:1;c.node_ids.assign(ns[i].begin(),ns[i].end());mesh.addCell(c);}
        Region si;si.id=0;si.name="silicon";si.material="Si";si.cell_ids={0,1};mesh.addRegion(si);
        Region ox;ox.id=1;ox.name="oxide";ox.material="SiO2";ox.cell_ids={2,3};mesh.addRegion(ox);
        Contact left;left.id=0;left.name="left";left.node_ids={0,3};left.edge_node_ids={{{0,3}}};
        if (interfaceCorner) {left.node_ids={0};left.edge_node_ids.clear();}
        mesh.addContact(left);
        Contact right;right.id=1;right.name="right";right.node_ids={1,2};right.edge_node_ids={{{1,2}}};mesh.addContact(right);
        Contact bottom;bottom.id=2;bottom.name="non_electrode";bottom.node_ids={0,1};bottom.edge_node_ids={{{0,1}}};mesh.addContact(bottom);
        mesh.buildEdges();
        doping=DopingModel::fromMeshAndRegions(mesh,{{"silicon",tcad?2e17:2e23,tcad?8e16:8e22}});
        using json=nlohmann::json;
        json cells=json::array();
        std::map<std::pair<Index,Index>,Real> sums;
        for (Index i=0;i<2;++i) {
            const std::array<Real,3> coeff{.4+.1*i,.3,.8};
            cells.push_back({{"cell_id",i},{"node_ids",ns[i]},{"vertex_measure_m2",{1e-15,2e-15,2e-15}},{"edge_coefficients",coeff}});
            for (int k=0;k<3;++k) sums[std::minmax(ns[i][k],ns[i][(k+1)%3])]+=coeff[k];
        }
        for (const auto& e:mesh.edges()) mesh.setTransportCouple(e.id,sums[std::minmax(e.n0,e.n1)]*e.length);
        path=std::filesystem::temp_directory_path()/
            ("vela_ialmob_transport_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".json");
        std::ofstream(path)<<json({{"schema","vela.ialmob.transport_geometry.v1"},{"node_count",6},{"cell_count",4},{"cells",cells}});
        mobility=mobilityModelConfigFromJson({{"model","ialmob"},{"high_field_driving_force","quasi_fermi_gradient"},
            {"ialmob",{{"geometry_file",path.string()},{"effective_electrodes",{"left","right"}},
                {"crystal_x",{1,0,0}},{"crystal_y",{0,1,0}},
                {"electron_parameters_cm",{{"100",{{"l_crit",1e-6},{"l_crit_c",1e-6}}}}},
                {"hole_parameters_cm",{{"100",{{"l_crit",1e-6},{"l_crit_c",1e-6}}}}}}}},inputScaling);
    }
    ~Fixture() {std::error_code ec;std::filesystem::remove(path,ec);}
};
}

TEST_CASE("IALMob full residual Jacobian and terminal currents share live mobility", "[ialmob][jacobian][transport]")
{
    const Real vt=constants::kb*300./constants::q;
    for (bool tcad:{false,true}) {
    Fixture f(tcad);DDScalingSpec scaling;
    if (tcad) {
        scaling.enabled=true;scaling.V0=vt;scaling.C0=1e17;scaling.D0=1000.*vt;
        scaling.mu0=1000.;scaling.L0=.1;
        scaling.permittivityReference_F_per_m=constants::eps0*11.7;
        scaling.unitSystem=PhysicalUnitSystem::tcadInternal();
        scaling.chargeAreaFactor=scaling.unitSystem.chargeAreaFactor();
        scaling.chargeLineFactor=scaling.unitSystem.chargeLineFactor();
        scaling.fieldFromCoordinateDeltaFactor=scaling.unitSystem.fieldFromCoordinateDeltaFactor();
        scaling.currentDensityLineIntegralFactor=scaling.unitSystem.currentDensityAM2PerInternal()*scaling.unitSystem.lengthMPerInternal();
    }
    const Real potentialScale=tcad?vt:1.;
    for (bool fermi:{false,true}) {
        CarrierStatisticsConfig stats;stats.model=fermi?"fermi_dirac":"boltzmann";
        RecombinationModelConfig recombination;recombination.mechanisms.clear();
        CoupledDDAssembler assembler(f.mesh,f.materials,f.doping,vt,f.mobility,recombination, {}, {}, {}, {}, scaling, {},stats);
        CoupledDDState s;s.psi=VectorXd::Zero(6);s.phin=VectorXd::Zero(6);s.phip=VectorXd::Zero(6);
        // Both carrier populations must be resolved by double-precision
        // residual differences, including the cross-carrier mobility columns.
        for (int i=0;i<6;++i) {s.psi[i]=.48+.012*i;s.phin[i]=.003*i;s.phip[i]=.96-.002*i;}
        s.psi/=potentialScale;s.phin/=potentialScale;s.phip/=potentialScale;
        const auto x=assembler.pack(s);CoupledDDBoundaryConditions bc;
        const auto jac=assembler.assembleJacobian(x,bc);
        for (int k=0;k<x.size();++k) {
            auto plus=x,minus=x;const Real step=2e-7/potentialScale;plus[k]+=step;minus[k]-=step;
            const VectorXd fd=(assembler.residual(plus,bc)-assembler.residual(minus,bc))/(2.*step);
            const VectorXd column=VectorXd(jac.col(k));
            for (int block=0;block<3;++block) {
                INFO("TCAD="<<tcad<<" Fermi="<<fermi<<" column="<<k<<" rowblock="<<block);
                const auto a=column.segment(block*6,6).eval(),b=fd.segment(block*6,6).eval();
                CHECK((a-b).norm()/std::max({a.norm(),b.norm(),1e-30})<3e-5);
            }
        }
        // A cross-carrier column must be present even without recombination
        // at a node outside the SG edge: it enters the live IALMob mobility.
        DDSolution solution;solution.psi=s.psi*potentialScale;solution.phin=s.phin*potentialScale;solution.phip=s.phip*potentialScale;
        solution.n=assembler.electronDensity(x);solution.p=assembler.holeDensity(x);
        ContactCurrent current(f.mesh,f.materials,f.doping,f.mobility,300.,scaling, {},stats);
        for (const std::string name:{"left","right"}) {
            const auto edge=current.compute(solution,name);
            const auto residual=current.computeFromResidual(assembler,x,name);
            CHECK(edge.electronCurrent==Catch::Approx(residual.electronCurrent).epsilon(1e-10).margin(1e-22));
            CHECK(edge.holeCurrent==Catch::Approx(residual.holeCurrent).epsilon(1e-10).margin(1e-22));
        }
        // Boundary rows and the enlarged fixed pattern must remain consistent.
        bc.psi[0]=s.psi[0];bc.phin[0]=s.phin[0];bc.phip[0]=s.phip[0];
        const auto constrained=assembler.assembleJacobian(x,bc);
        for (int row:{0,6,12}) {
            CHECK(constrained.coeff(row,row)==1.);
            for (int col=0;col<18;++col) if (col!=row) CHECK(constrained.coeff(row,col)==0.);
        }
    }
}

}

TEST_CASE("IALMob transport rejects geometry coupling mismatch", "[ialmob][transport]")
{
    Fixture f;f.mesh.setTransportCouple(0,123.);
    VectorXd psi=VectorXd::Zero(6),n=VectorXd::Constant(6,1e20),p=n;
    CHECK_THROWS_AS(updateIalTransportState(f.mobility,f.mesh,f.doping,psi,n,p,psi,psi),std::invalid_argument);
}

TEST_CASE("IALMob interface termination uses material interface before exterior face", "[ialmob][transport]")
{
    Fixture f(false,true);
    VectorXd psi=VectorXd::Zero(6),n=VectorXd::Constant(6,1e20),p=n;
    updateIalTransportState(f.mobility,f.mesh,f.doping,psi,n,p,psi,psi);
    const auto& g=f.mobility.ialmobGeometry->cells.at(1).geometry;
    CHECK(g.partialBoundaryLayer);
    CHECK(std::abs(g.boundaryTangent[0])==Catch::Approx(1.));
    CHECK(g.boundaryTangent[1]==0.);
}
