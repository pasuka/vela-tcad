#include "vela/equation/ElectrothermalAssembler.h"
#include "vela/simulation/ElectrothermalSimulation.h"
#include "vela/core/IalKernelProfiling.h"
#include <sstream>
#include "vela/solver/ElectrothermalTangent.h"
#include <Eigen/SparseLU>
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
TEST_CASE("Sweep preparation observes exact sources and never caches carrier state", "[electrothermal][preparation_context]") {
    using J=nlohmann::json;
    const auto path=std::filesystem::temp_directory_path()/
        ("vela_preparation_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".json");
    struct Cleanup {std::filesystem::path p;~Cleanup(){std::error_code e;std::filesystem::remove(p,e);}} cleanup{path};
    J mesh={{"regions",J::array({{{"id",0},{"name","silicon"},{"material","Silicon"},{"cell_ids",{0,1}}}})},
        {"nodes",J::array()},{"triangles",J::array({{{"id",0},{"region_id",0},{"node_ids",{0,1,2}}},
            {{"id",1},{"region_id",0},{"node_ids",{0,2,3}}}})},{"contacts",J::array()}};
    for(int i=0;i<4;++i){mesh["nodes"].push_back({{"id",i},{"x",i==1||i==2?1.:0.},{"y",i>=2?1.:0.}});
        mesh["contacts"].push_back({{"id",i},{"region_id",0},{"name",std::array{"source","drain","substrate","gate"}[i]},{"node_ids",{i}}});}
    const auto write=[&](){std::ofstream f(path);f<<mesh.dump();};write();
    J cfg={{"mesh_file",path.string()},{"coordinate_to_metres",1.},
        {"region_conductivity",J::array({{{"region_id",0},{"model","constant"},{"value_W_per_m_K",1.}}})},
        {"thermodes",J::array()},{"silicon_area_m2",{1./3,1./6,1./3,1./6}},
        {"fixed_charge_C_per_m",{0.,0.,0.,0.}},{"donors_m3",{0.,0.,0.,0.}},{"acceptors_m3",{0.,0.,0.,0.}},
        {"mobility_SI",{{"model","constant"}}},{"boundaries",J::array()},
        {"state_interleaved",J::array()},{"edge_geometry",J::array()},
        {"performance_profiling",true},{"diagnostic_newton_max_iterations",0}};
    for(int i=0;i<4;++i)for(int k=0;k<4;++k){const double v=k==3?300.:0.;cfg["state_interleaved"].push_back(v);
        cfg["boundaries"].push_back({{"node",i},{"kind",std::array{"psi","fn","fp","temperature"}[k]},{"value",v}});}
    for(auto e:std::array<std::array<int,2>,5>{{{0,1},{1,2},{0,2},{2,3},{0,3}}})
        cfg["edge_geometry"].push_back({{"nodes",e},{"poisson_F_per_m",1e-10},{"transport_weight",0.}});
    ElectrothermalPreparationContext context;std::ostringstream log;
    const auto compare=[&](bool hit){
        auto actual=solveElectrothermalPoint(cfg,log,&context),expected=solveElectrothermalPoint(cfg,log);
        CHECK(actual["performance"]["static_preparation_reused"]==hit);
        actual.erase("performance");expected.erase("performance");CHECK(actual==expected);
    };
    compare(false);compare(true);
    cfg["state_interleaved"][3]=350.;cfg["boundaries"][0]["value"]=.1;compare(true);
    cfg["donors_m3"][0]=1e20;compare(false);compare(true);
    cfg["coordinate_to_metres"]=2.;compare(false);
    cfg["edge_geometry"][0]["poisson_F_per_m"]=2e-10;compare(false);
    // Same path, byte count and timestamp: a metadata-only cache would be stale.
    const auto stamp=std::filesystem::last_write_time(path);const auto size=std::filesystem::file_size(path);
    mesh["nodes"][1]["x"]=2.;write();std::filesystem::last_write_time(path,stamp);
    REQUIRE(std::filesystem::file_size(path)==size);compare(false);compare(true);
    {std::ofstream f(path);f<<"bad mesh";}
    CHECK_THROWS(solveElectrothermalPoint(cfg,log,&context));write();compare(true);
}
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

TEST_CASE("IALMob Fermi Auger enhancement has a consistent coupled source Jacobian",
          "[ialmob][auger][density_dependence][jacobian]")
{
    Fixture f;
    const Real vt=constants::kb*300./constants::q;
    RecombinationModelConfig active, off;
    active.mechanisms={"auger"};off.mechanisms.clear();
    active.augerDensityDependence={true,3.46667,8.25688,1e24,1e24};
    // Resolve the source derivatives above transport subtraction roundoff.
    active.augerCn*=1e6;active.augerCp*=1e6;
    const CarrierStatisticsConfig stats{"fermi_dirac"};
    CoupledDDAssembler on(f.mesh,f.materials,f.doping,vt,f.mobility,active,{}, {}, {}, {}, {}, {},stats);
    CoupledDDAssembler zero(f.mesh,f.materials,f.doping,vt,f.mobility,off,{}, {}, {}, {}, {}, {},stats);
    CoupledDDState s;s.psi=VectorXd::Zero(6);s.phin=VectorXd::Zero(6);s.phip=VectorXd::Zero(6);
    for(int i=0;i<6;++i){s.psi[i]=.48+.003*i;s.phin[i]=.002*i;s.phip[i]=.96-.002*i;}
    const auto x=on.pack(s);CoupledDDBoundaryConditions bc;
    const auto jac=(on.assembleJacobian(x,bc)-zero.assembleJacobian(x,bc)).eval();
    for(int k=0;k<x.size();++k){
        auto plus=x,minus=x;constexpr Real h=2e-6;plus[k]+=h;minus[k]-=h;
        const VectorXd sourcePlus=on.residual(plus,bc)-zero.residual(plus,bc);
        const VectorXd sourceMinus=on.residual(minus,bc)-zero.residual(minus,bc);
        const VectorXd fd=(sourcePlus-sourceMinus)/(2*h);
        const VectorXd col=VectorXd(jac.col(k));
        INFO("column="<<k);
        CHECK((col-fd).norm()/std::max({col.norm(),fd.norm(),1e-20})<2e-5);
    }
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
        auto preparedMobility=current.prepareMobility(solution);
        ContactCurrent preparedCurrent(f.mesh,f.materials,f.doping,preparedMobility,300.,scaling, {},stats);
        for (const std::string name:{"left","right"}) {
            const auto edge=current.compute(solution,name);
            const auto cached=preparedCurrent.compute(solution,name);
            CHECK(cached.electronCurrent==edge.electronCurrent);
            CHECK(cached.holeCurrent==edge.holeCurrent);
            auto changed=solution;
            changed.phin[2]+=.01;
            changed.n*=2.;
            changed.electronQfReference_V=16.;
            changed.phinIncrement=changed.phin.array()-16.;
            const auto changedCold=current.compute(changed,name);
            const auto changedPrepared=preparedCurrent.compute(changed,name);
            CHECK(changedPrepared.electronCurrent==changedCold.electronCurrent);
            CHECK(changedPrepared.holeCurrent==changedCold.holeCurrent);
            ContactCurrent referencedCurrent(f.mesh,f.materials,f.doping,
                current.prepareMobility(changed),300.,scaling, {},stats);
            CHECK(referencedCurrent.compute(changed,name).electronCurrent==changedCold.electronCurrent);
            CHECK(referencedCurrent.compute(changed,name).holeCurrent==changedCold.holeCurrent);
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

TEST_CASE("IALMob prepared geometry preserves live state changes across consumers", "[ialmob][transport]")
{
    Fixture f;
    const auto cold = f.mobility;
    prepareIalTransportGeometry(f.mobility,f.mesh);
    REQUIRE(f.mobility.ialmobGeometry);
    CHECK_FALSE(f.mobility.ialmobState);
    auto first=f.mobility,second=f.mobility;
    VectorXd psi=VectorXd::LinSpaced(6,0.,.12);
    VectorXd qn=VectorXd::LinSpaced(6,0.,.03),qp=-qn;
    VectorXd n=VectorXd::Constant(6,2e23),p=VectorXd::Constant(6,8e22);
    updateIalTransportState(first,f.mesh,f.doping,psi,n,p,qn,qp);
    const auto before=first.ialmobState;
    n*=3.;qn*=2.;
    updateIalTransportState(second,f.mesh,f.doping,psi,n,p,qn,qp);
    auto independent=cold;
    updateIalTransportState(independent,f.mesh,f.doping,psi,n,p,qn,qp);
    CHECK(first.ialmobGeometry==second.ialmobGeometry);
    CHECK(first.ialmobState==before);
    CHECK(first.ialmobState!=second.ialmobState);
    bool changed=false;
    for (Index edge=0;edge<f.mesh.numEdges();++edge) {
        for (auto carrier:{CarrierType::Electron,CarrierType::Hole}) {
            const Real fresh=ialEdgeMobility(independent,edge,carrier);
            CHECK(ialEdgeMobility(second,edge,carrier)==fresh);
            changed |= ialEdgeMobility(first,edge,carrier)!=fresh;
        }
    }
    CHECK(changed);
}

TEST_CASE("IALMob temperature state participates in cache and edge mobility", "[ialmob][thermal][transport]") {
    Fixture f;VectorXd psi=VectorXd::LinSpaced(6,.1,.12),fn=VectorXd::LinSpaced(6,0.,.005),fp=VectorXd::Constant(6,.05);
    VectorXd n=VectorXd::Constant(6,1e23),p=VectorXd::Constant(6,1e22),dn=n/.04,dp=p/.04;
    VectorXd t=VectorXd::Constant(6,400.),nt=n*.003,pt=p*.002;
    auto cfg=f.mobility;
    updateIalTransportState(cfg,f.mesh,f.doping,psi,n,p,fn,fp,dn,dp,t,nt,pt);
    auto old=cfg.ialmobState;
    updateIalTransportState(cfg,f.mesh,f.doping,psi,n,p,fn,fp,dn,dp,t,nt,pt);
    CHECK(cfg.ialmobState==old);
    t[0]+=.01; n[0]+=.01*nt[0];p[0]+=.01*pt[0];
    updateIalTransportState(cfg,f.mesh,f.doping,psi,n,p,fn,fp,dn,dp,t,nt,pt);
    CHECK(cfg.ialmobState!=old);
    bool changed=false;
    for(int i=0;i<f.mesh.numEdges();++i)changed=changed||cfg.ialmobState->electronEdges[i]!=old->electronEdges[i];
    CHECK(changed);
    CHECK_THROWS_AS(updateIalTransportState(cfg,f.mesh,f.doping,psi,n,p,fn,fp,dn,dp,{},nt,pt),std::invalid_argument);
}


TEST_CASE("Four-equation operator couples live IALMob current and conservative heat", "[thermal][ialmob][electrothermal]") {
    Fixture f;
    ElectrothermalGeometry g;
    g.siliconArea_m2=VectorXd::Zero(6);g.fixedCharge_C_per_m=VectorXd::Zero(6);
    g.poissonEdge_F_per_m=VectorXd::Zero(f.mesh.numEdges());g.transportWeight=g.poissonEdge_F_per_m;
    for(const auto& cell:f.mesh.cells())if(cell.region_id==0)
        for(Index i:cell.node_ids)g.siliconArea_m2[i]+=1e-14/6.;
    for(const auto& edge:f.mesh.edges()){
        g.poissonEdge_F_per_m[edge.id]=constants::eps0*8.*edge.couple/edge.length;
        g.transportWeight[edge.id]=edge.transport_couple/edge.length;
    }
    // Charge geometry can differ from transport/recombination volumes. Exercise
    // nonuniform explicit native-style measures through every coupled partial.
    g.recombinationArea_m2=g.siliconArea_m2;
    SECTION("Existing charge geometry"){}
    SECTION("Independent nonuniform Poisson edge and charge-volume pair"){
        for(int i=0;i<6;++i)g.siliconArea_m2[i]*=1.+.2*i;
        for(const auto& edge:f.mesh.edges())g.poissonEdge_F_per_m[edge.id]*=1.+.1*edge.id;
    }
    SECTION("Reused thermal high-field partials preserve the full coupled Jacobian"){
        auto options=std::make_shared<IalTransportOptions>(*f.mobility.ialmob);
        options->element.reuseThermalHighField=true;f.mobility.ialmob=options;
    }
    SECTION("Explicit high-field partials preserve the full coupled Jacobian"){
        auto options=std::make_shared<IalTransportOptions>(*f.mobility.ialmob);
        options->element.reuseThermalHighField=true;options->element.explicitHighFieldPartials=true;
        f.mobility.ialmob=options;
    }
    SECTION("Generated low-field and explicit high-field preserve the coupled Jacobian"){
        auto options=std::make_shared<IalTransportOptions>(*f.mobility.ialmob);
        options->element.reuseThermalHighField=true;options->element.explicitHighFieldPartials=true;
        options->element.generatedLowFieldPartials=true;f.mobility.ialmob=options;
    }
    LatticeConductivity law;law.model=LatticeConductivity::Model::InverseQuadratic;
    law.numerator=100.;law.denominator={-.0393,.00155,1.82e-6};
    LatticeHeatAssembler heat(f.mesh,1.,{{0,law},{1,law}},{{{0,1},300.,2e6}});
    ElectrothermalAssembler coupled(f.mesh,f.doping,g,heat,f.mobility);
    VectorXd x(24);for(int i=0;i<6;++i){x[4*i]=.48+.012*i;x[4*i+1]=.003*i;x[4*i+2]=.96-.002*i;x[4*i+3]=350.+12.*i;}
    for(int mode:{0,1,2}){
        const bool constrained=mode!=0;
        ElectrothermalBoundary bc;
        if(constrained){bc.neutralContactBias_V={{0,0.}};bc.potential_V={{4,.2}};bc.temperature_K={{5,400.}};}
        if(mode==2)bc.holeRecombination={{0,{1.93e4,1e-7}}};

        const auto base=coupled.assemble(x,bc);
        std::vector<ElectrothermalHoleRowAudit> audit(1);audit[0].node=1;
        const auto observed=coupled.assemble(x,bc,{},{},true,false,false,false,&audit);
        CHECK((observed.residual.array()==base.residual.array()).all());
        CHECK((Eigen::MatrixXd(observed.jacobian).array()==Eigen::MatrixXd(base.jacobian).array()).all());
        Real sum=audit[0].source,absolute=0.;
        for(const auto& e:audit[0].edges){sum+=(e.a==1?1.:-1.)*e.current;absolute+=std::abs(e.current);}
        CHECK(sum==base.residual[6]);CHECK(absolute==base.holeFluxAbs_A_per_m[1]);
        CHECK(audit[0].residual==sum);REQUIRE(!audit[0].edges.empty());
        audit[0].node=6;CHECK_THROWS_AS(coupled.assemble(x,bc,{},{},false,false,false,false,&audit),std::invalid_argument);
        if(mode==2){
            audit[0].node=0;coupled.assemble(x,bc,{},{},false,false,false,false,&audit);
            REQUIRE(audit[0].finiteContact);Real sum=audit[0].source,absolute=0.;
            for(const auto& e:audit[0].edges){sum+=(e.a==0?1.:-1.)*e.current;absolute+=std::abs(e.current);}
            sum+=audit[0].contactOutward;absolute+=std::abs(audit[0].contactOutward);
            CHECK(sum==base.residual[2]);CHECK(absolute==base.holeFluxAbs_A_per_m[0]);
        }
        const auto frozenMu=coupled.assemble(x,bc,{},{},true,false,true,false);
        const auto frozenRec=coupled.assemble(x,bc,{},{},true,false,false,true);
        const auto frozenBoth=coupled.assemble(x,bc,{},{},true,false,true,true);
        for(const auto* frozen:{&frozenMu,&frozenRec,&frozenBoth}) {
            CHECK((frozen->residual.array()==base.residual.array()).all());
            CHECK(frozen->latticeSource_W_per_m==base.latticeSource_W_per_m);
        }
        CHECK((base.jacobian-frozenMu.jacobian).norm()>0.);
        CHECK((base.jacobian-frozenRec.jacobian).norm()>0.);
        CHECK((base.jacobian-frozenMu.jacobian-frozenRec.jacobian+frozenBoth.jacobian).norm()<1e-12*base.jacobian.norm());
        // Omitting a recombination derivative must not change Poisson/heat or
        // any row replaced by a boundary condition. Full Jacobian FD checks
        // below continue to use the complete operator.
        for(int i=0;i<24;++i)if(i%4==0 || i%4==3 || (constrained && i<3 && !(mode==2 && i==2)))
            for(int j=0;j<24;++j)CHECK(base.jacobian.coeff(i,j)==frozenRec.jacobian.coeff(i,j));
        VectorXd eRef(6),hRef(6),referenced=x;
        for(int i=0;i<6;++i){eRef[i]=x[4*i+1];hRef[i]=x[4*i+2];referenced[4*i+1]=0.;referenced[4*i+2]=0.;}
        const auto recentered=coupled.assemble(referenced,bc,eRef,hRef);
        for(int block=0;block<4;++block){Real error=0.,scale=0.;for(int i=0;i<6;++i){const int row=4*i+block;
            error+=std::pow(base.residual[row]-recentered.residual[row],2);scale+=base.residual[row]*base.residual[row];}
            CHECK(std::sqrt(error/std::max(scale,1e-100))<1e-10);}
        for(int k=0;k<24;++k){auto a=referenced,b=referenced;const Real h=k%4==3?.002:1e-6;a[k]+=h;b[k]-=h;
            const VectorXd fd=(coupled.assemble(a,bc,eRef,hRef).residual-coupled.assemble(b,bc,eRef,hRef).residual)/(2.*h);
            const VectorXd exact=recentered.jacobian.col(k);
            for(int block=0;block<4;++block){Real err=0.,scale=0.;for(int i=0;i<6;++i){const int row=4*i+block;
                err+=std::pow(fd[row]-exact[row],2);scale+=fd[row]*fd[row]+exact[row]*exact[row];}
                CHECK(std::sqrt(err/std::max(scale,1e-100))<3e-5);}
        }

        if(!constrained){Real sum=0.;for(int i=0;i<6;++i)sum+=base.residual[4*i+3];
            CHECK(sum==Catch::Approx(base.boundaryHeat_W_per_m-base.latticeSource_W_per_m).epsilon(1e-12));}
        CHECK(std::abs(base.electronOutflow_A_per_m.sum())<1e-10);
        if(mode!=2)CHECK(std::abs(base.holeOutflow_A_per_m.sum())<1e-10);
        for(int k=0;k<24;++k)for(Real fraction:{1.,.25}){
            auto a=x,b=x;Real h=(k%4==3?.002:1e-6)*fraction;a[k]+=h;b[k]-=h;
            const VectorXd fd=(coupled.assemble(a,bc).residual-coupled.assemble(b,bc).residual)/(2.*h);
            const VectorXd exact=base.jacobian.col(k);
            for(int block=0;block<4;++block){Real err=0.,scale=0.;
                for(int i=0;i<6;++i){const int row=4*i+block;err+=std::pow(fd[row]-exact[row],2);scale+=fd[row]*fd[row]+exact[row]*exact[row];}
                INFO("column="<<k<<" block="<<block<<" boundary="<<constrained);
                CHECK(std::sqrt(err/std::max(scale,1e-100))<3e-5);
            }
        }
    }
}

TEST_CASE("Electrothermal preparation reuse refreshes temperature and doping without changing the operator", "[thermal][electrothermal][preparation]") {
    Fixture f;ElectrothermalGeometry g;
    auto fastMobility=f.mobility;
    SECTION("Existing full local evaluation"){}
    SECTION("Node preparation and values-only residual candidates"){
        auto options=std::make_shared<IalTransportOptions>(*f.mobility.ialmob);
        options->element.reuseLocalPreparation=true;options->element.residualValuesOnly=true;
        fastMobility.ialmob=options;
    }
    SECTION("High-field work reused only within one element state"){
        auto options=std::make_shared<IalTransportOptions>(*f.mobility.ialmob);
        options->element.reuseLocalPreparation=true;options->element.residualValuesOnly=true;
        options->element.reuseThermalHighField=true;fastMobility.ialmob=options;
    }
    g.siliconArea_m2=VectorXd::Zero(6);g.fixedCharge_C_per_m=VectorXd::Zero(6);
    g.poissonEdge_F_per_m=VectorXd::Zero(f.mesh.numEdges());g.transportWeight=g.poissonEdge_F_per_m;
    for(const auto& cell:f.mesh.cells())if(cell.region_id==0)for(Index i:cell.node_ids)g.siliconArea_m2[i]+=1e-14/6.;
    for(const auto& edge:f.mesh.edges()){
        g.poissonEdge_F_per_m[edge.id]=constants::eps0*8.*edge.couple/edge.length;
        g.transportWeight[edge.id]=edge.transport_couple/edge.length;
    }
    LatticeConductivity law;law.constant_W_per_m_K=100.;
    LatticeHeatAssembler heat(f.mesh,1.,{{0,law},{1,law}},{{{0,1},300.,2e6}});
    ElectrothermalAssembler plain(f.mesh,f.doping,g,heat,f.mobility),cached(f.mesh,f.doping,g,heat,fastMobility,SiliconThermalPhysics{},.1,.04,true,true,true);
    ElectrothermalBoundary bc;bc.neutralContactBias_V={{0,0.}};bc.holeRecombination={{0,{1.93e4,1e-7}}};
    VectorXd x(24);
    for(Real temperature:{300.,401.,401.,299.,514.}){
        for(int i=0;i<6;++i){x[4*i]=.48+.012*i;x[4*i+1]=.003*i;x[4*i+2]=.96-.002*i;x[4*i+3]=temperature+3.*i;}
        IalKernelProfilingScope profile(true);
        const auto a=plain.assemble(x,bc);const auto beforeWork=ialKernelProfile;
        const auto b=cached.assemble(x,bc);const auto afterWork=ialKernelProfile;
        if(fastMobility.ialmob->element.reuseThermalHighField && afterWork.passes[1]>beforeWork.passes[1]) {
            CHECK(afterWork.highFieldReuses>beforeWork.highFieldReuses);
            CHECK(afterWork.highFieldEvaluations-beforeWork.highFieldEvaluations<beforeWork.highFieldEvaluations);
            CHECK(afterWork.seconds[2]>beforeWork.seconds[2]);
        }
        const auto residualOnly=cached.assemble(x,bc,{},{},false);
        CHECK(residualOnly.jacobian.rows()==0);
        CHECK((b.residual-residualOnly.residual).norm()==0.);
        CHECK((b.electronOutflow_A_per_m-residualOnly.electronOutflow_A_per_m).norm()==0.);
        CHECK((b.holeOutflow_A_per_m-residualOnly.holeOutflow_A_per_m).norm()==0.);
        CHECK((b.electronFluxAbs_A_per_m-residualOnly.electronFluxAbs_A_per_m).norm()==0.);
        CHECK((b.holeFluxAbs_A_per_m-residualOnly.holeFluxAbs_A_per_m).norm()==0.);
        CHECK((b.recombination_A_per_m-residualOnly.recombination_A_per_m).norm()==0.);
        CHECK(b.latticeSource_W_per_m==residualOnly.latticeSource_W_per_m);
        CHECK(b.boundaryHeat_W_per_m==residualOnly.boundaryHeat_W_per_m);
        CHECK((b.jacobian-cached.assemble(x,bc).jacobian).norm()==0.);
        ElectrothermalAssembler partial(f.mesh,f.doping,g,heat,fastMobility,SiliconThermalPhysics{},.1,.04,true,true);
        const auto firstResidual=partial.assemble(x,bc,{},{},false);
        CHECK((b.residual-firstResidual.residual).norm()==0.);
        // A residual cache may lack all derivative columns; the same-state full
        // request must rebuild them, not reuse incomplete derivatives.
        CHECK((b.jacobian-partial.assemble(x,bc).jacobian).norm()==0.);
        CHECK((a.residual-b.residual).norm()==0.);CHECK((a.jacobian-b.jacobian).norm()==0.);
        const auto before=cached.preparationCounts();const auto repeat=cached.assemble(x,bc);
        CHECK((b.residual-repeat.residual).norm()==0.);
        CHECK(cached.preparationCounts()[0]==before[0]);CHECK(cached.preparationCounts()[1]==before[1]);
        for(Real bias:{0.,40.}){
            const auto expected=plain.neutralPotential(0,bias,temperature),actual=cached.neutralPotential(0,bias,temperature);
            CHECK(expected==actual);
            const auto counts=cached.neutralRootCounts();
            CHECK(cached.neutralPotential(0,bias,temperature)==expected);
            CHECK(cached.neutralRootCounts()[0]==counts[0]);
            CHECK(cached.neutralRootCounts()[1]==counts[1]+1);
            const Real adjacent=std::nextafter(temperature,1000.);
            CHECK(cached.neutralPotential(0,bias,adjacent)==plain.neutralPotential(0,bias,adjacent));
            CHECK(cached.neutralRootCounts()[0]==counts[0]+1);
        }
    }
    auto poissonState=x;ElectrothermalBoundary poissonBoundary;
    for(int i=0;i<6;++i){
        poissonState[4*i+1]=poissonState[4*i+2]=0.;poissonState[4*i+3]=300.;
        poissonBoundary.electronQf_V[i]=poissonBoundary.holeQf_V[i]=0.;
        poissonBoundary.temperature_K[i]=300.;
    }
    for(int perturbation=0;perturbation<3;++perturbation){
        auto state=poissonState;
        if(perturbation==1)state[1]=.001;
        if(perturbation==2)state[3]=301.;
        const auto full=plain.assemble(state,poissonBoundary);
        const auto fast=cached.assemble(state,poissonBoundary,{},{},true,true);
        CHECK((full.residual-fast.residual).norm()==0.);
        CHECK((full.jacobian-fast.jacobian).norm()==0.);
        CHECK((full.electronOutflow_A_per_m-fast.electronOutflow_A_per_m).norm()==0.);
        CHECK((full.holeOutflow_A_per_m-fast.holeOutflow_A_per_m).norm()==0.);
        CHECK(full.latticeSource_W_per_m==fast.latticeSource_W_per_m);
        CHECK(full.boundaryHeat_W_per_m==fast.boundaryHeat_W_per_m);
    }
    const auto before=cached.preparationCounts();
    auto tinyQf=poissonState;const VectorXd largeReference=VectorXd::Constant(6,40.);
    for(int i=0;i<6;++i)tinyQf[4*i]+=40.;
    tinyQf[1]=1e-18;
    const auto tinyFull=plain.assemble(tinyQf,poissonBoundary,largeReference,largeReference);
    const auto tinyFast=cached.assemble(tinyQf,poissonBoundary,largeReference,largeReference,true,true);
    CHECK(tinyFull.electronOutflow_A_per_m.norm()>0.);
    CHECK((tinyFull.electronOutflow_A_per_m-tinyFast.electronOutflow_A_per_m).norm()==0.);
    CHECK((tinyFull.residual-tinyFast.residual).norm()==0.);
    f.doping.setNodeDoping(0,f.doping.donors(0)+1e22,f.doping.acceptors(0));
    CHECK(cached.neutralPotential(0,40.,514.)==plain.neutralPotential(0,40.,514.));
    const auto a=plain.assemble(x,bc),b=cached.assemble(x,bc);
    CHECK((a.residual-b.residual).norm()==0.);CHECK((a.jacobian-b.jacobian).norm()==0.);
    CHECK(cached.preparationCounts()[0]==before[0]+1);
    for(int node:{0,2}){
        auto hi=x,lo=x;hi[4*node+3]+=.002;lo[4*node+3]-=.002;
        const VectorXd fd=(cached.assemble(hi,bc).residual-cached.assemble(lo,bc).residual)/.004;
        const VectorXd exact=b.jacobian.col(4*node+3);
        for(int k=0;k<4;++k){Real error=0.,scale=0.;for(int i=0;i<6;++i){int row=4*i+k;
            error+=std::pow(fd[row]-exact[row],2);scale+=fd[row]*fd[row]+exact[row]*exact[row];}
            CHECK(std::sqrt(error/std::max(scale,1e-100))<3e-5);
        }
    }
}

TEST_CASE("Safeguarded neutral Newton preserves charge balance and temperature response", "[thermal][electrothermal][neutral_root]") {
    Fixture f;ElectrothermalGeometry g;
    g.siliconArea_m2=VectorXd::Zero(6);g.siliconArea_m2[0]=1e-14;
    g.fixedCharge_C_per_m=VectorXd::Zero(6);
    g.poissonEdge_F_per_m=VectorXd::Zero(f.mesh.numEdges());g.transportWeight=g.poissonEdge_F_per_m;
    LatticeConductivity law;law.constant_W_per_m_K=100.;
    LatticeHeatAssembler heat(f.mesh,1.,{{0,law},{1,law}},{});
    ElectrothermalAssembler plain(f.mesh,f.doping,g,heat,{},SiliconThermalPhysics{},.1,.04,true);
    ElectrothermalAssembler fast(f.mesh,f.doping,g,heat,{},SiliconThermalPhysics{},.1,.04,true,false,false,true);
    SiliconThermalPhysics physics;
    for(const auto& doping:std::array<std::array<Real,2>,7>{{{0.,0.},{1e19,0.},{1e23,0.},
            {1e26,0.},{0.,1e23},{0.,1e26},{1e24,1e24}}}) {
        f.doping.setNodeDoping(0,doping[0],doping[1]);
        for(Real t:{200.,300.,std::nextafter(300.,400.),450.,600.})for(Real bias:{-40.,-1.,-0.,0.,.375,40.}) {
            INFO("Nd="<<doping[0]<<" Na="<<doping[1]<<" T="<<t<<" bias="<<bias);
            const auto a=plain.neutralPotential(0,bias,t),b=fast.neutralPotential(0,bias,t);
            CHECK(std::abs(a.first-b.first)<=8.*std::numeric_limits<Real>::epsilon()*std::max(1.,std::abs(bias)));
            CHECK(std::abs(a.second-b.second)<=1e-12*std::max(1e-6,std::abs(a.second)));
            const auto state=physics.evaluate({b.first,bias,bias,t,doping[0],doping[1]});
            const Real scale=std::max({state.electrons_m3.value,state.holes_m3.value,doping[0],doping[1],1.});
            CHECK(std::abs(state.electrons_m3.value-state.holes_m3.value-doping[0]+doping[1])/scale<1e-11);
            const Real h=.002;
            const Real fd=(fast.neutralPotential(0,bias,t+h).first-fast.neutralPotential(0,bias,t-h).first)/(2.*h);
            CHECK(std::abs(fd-b.second)<3e-8*std::max(.001,std::abs(b.second)));
        }
    }
    const auto iterations=fast.neutralRootIterationCounts();
    CHECK(iterations[1]>0);CHECK(iterations[2]>0);CHECK(iterations[3]>0);
    // Counts include every derivative check and any original-method fallback.
    CHECK(iterations[0]<73*fast.neutralRootCounts()[0]);
    REQUIRE_THROWS(fast.neutralPotential(0,0.,-1.));
    REQUIRE_THROWS(fast.neutralPotential(4,0.,300.));
}

TEST_CASE("Four-equation silicon resistor closes self-heating and neutral contact temperature response", "[thermal][electrothermal][newton]") {
    bool finiteHoleContact=false;
    SECTION("Ideal contact"){}
    SECTION("Finite hole exchange"){finiteHoleContact=true;}
    DeviceMesh mesh;mesh.addRegion({0,"silicon","Si",{}});
    constexpr int nx=4,ny=2;constexpr Real L=1e-6;
    for(int j=0;j<=ny;++j)for(int i=0;i<=nx;++i)mesh.addNode({mesh.numNodes(),L*i/nx,L*j/ny,0.});
    for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){Index a=j*(nx+1)+i,b=a+1,c=a+nx+1,d=c+1;
        mesh.addCell({mesh.numCells(),CellType::Tri3,0,{a,b,d}});mesh.addCell({mesh.numCells(),CellType::Tri3,0,{a,d,c}});}
    mesh.buildEdges();const Index n=mesh.numNodes();DopingModel doping(n);
    ElectrothermalGeometry g;g.siliconArea_m2=VectorXd::Zero(n);g.fixedCharge_C_per_m=VectorXd::Zero(n);
    g.poissonEdge_F_per_m=VectorXd::Zero(mesh.numEdges());g.transportWeight=g.poissonEdge_F_per_m;
    for(Index i=0;i<n;++i){doping.setNodeDoping(i,1e23,0.);g.siliconArea_m2[i]=mesh.getNode(i).volume;}
    for(const auto& e:mesh.edges()){g.poissonEdge_F_per_m[e.id]=constants::eps0*11.7*e.couple/e.length;g.transportWeight[e.id]=e.couple/e.length;}
    LatticeConductivity law;law.model=LatticeConductivity::Model::InverseQuadratic;law.numerator=100.;law.denominator={-.0393,.00155,1.82e-6};
    std::vector<LatticeThermodeEdge> thermodes;for(int i=0;i<nx;++i)thermodes.push_back({{Index(i),Index(i+1)},300.,2e6});
    ElectrothermalAssembler coupled(mesh,doping,g,LatticeHeatAssembler(mesh,1.,{{0,law}},thermodes),{});
    ElectrothermalBoundary bc;
    for(int j=0;j<=ny;++j){bc.neutralContactBias_V[j*(nx+1)]=0.;bc.neutralContactBias_V[j*(nx+1)+nx]=0.;}
    if(finiteHoleContact)for(int j=0;j<=ny;++j){
        const Real length=L/ny*((j==0 || j==ny)?.5:1.);
        for(Index i:{Index(j*(nx+1)),Index(j*(nx+1)+nx)})bc.holeRecombination[i]={1.93e4,length};
    }
    VectorXd x(4*n);Real psi=coupled.neutralPotential(0,0.,300.).first;
    for(Index i=0;i<n;++i){x[4*i]=psi;x[4*i+1]=0.;x[4*i+2]=0.;x[4*i+3]=300.;}
    // Scale each physical equation separately. This test gate is for a
    // manufactured resistor, not a replacement for LDMOS electrical gates.
    VectorXd scales(4*n);for(Index i=0;i<n;++i){scales[4*i]=1e11;scales[4*i+1]=.01;scales[4*i+2]=.01;scales[4*i+3]=.01;}
    for(const auto& [i,v]:bc.neutralContactBias_V)for(int k=0;k<3;++k)
        if(!(finiteHoleContact && k==2))scales[4*i+k]=1.;
    auto norm=[&](const ElectrothermalAssembly& a){return (a.residual.array()*scales.array()).matrix().norm();};
    REQUIRE(norm(coupled.assemble(x,bc))<1e-8);
    for(Real bias:{.01,.05,.1}){
        for(int j=0;j<=ny;++j)bc.neutralContactBias_V[j*(nx+1)+nx]=bias;
        for(int iteration=0;iteration<35;++iteration){
            auto a=coupled.assemble(x,bc);if(norm(a)<1e-9)break;
            SparseMatrixd jac=a.jacobian;
            for(int k=0;k<jac.outerSize();++k)for(SparseMatrixd::InnerIterator it(jac,k);it;++it)it.valueRef()*=scales[it.row()];
            Eigen::SparseLU<SparseMatrixd> lu;lu.compute(jac);REQUIRE(lu.info()==Eigen::Success);
            VectorXd rhs=-(a.residual.array()*scales.array()).matrix();VectorXd delta=lu.solve(rhs);REQUIRE(lu.info()==Eigen::Success);
            Real alpha=1.;bool accepted=false;
            for(int trial=0;trial<20;++trial){VectorXd candidate=x+alpha*delta;
                try{if(norm(coupled.assemble(candidate,bc))<norm(a)){x=candidate;accepted=true;break;}}catch(const std::exception&){}
                alpha*=.5;
            }
            REQUIRE(accepted);
        }
        const auto a=coupled.assemble(x,bc);REQUIRE(norm(a)<1e-9);
        REQUIRE(a.latticeSource_W_per_m>0.);
        REQUIRE(std::abs(a.boundaryHeat_W_per_m-a.latticeSource_W_per_m)/a.latticeSource_W_per_m<1e-8);
        Real peak=300.;for(Index i=0;i<n;++i)peak=std::max(peak,x[4*i+3]);REQUIRE(peak>300.);
        Real left=0.,right=0.;for(int j=0;j<=ny;++j){Index l=j*(nx+1),r=l+nx;
            left+=a.electronOutflow_A_per_m[l]+a.holeOutflow_A_per_m[l];right+=a.electronOutflow_A_per_m[r]+a.holeOutflow_A_per_m[r];}
        REQUIRE(std::abs(left+right)/std::abs(left)<1e-9);
    }
}

TEST_CASE("Finite Ohmic hole contact preserves equilibrium and sub-ULP exchange", "[thermal][electrothermal][hrec]") {
    Fixture f;
    ElectrothermalGeometry g;g.siliconArea_m2=VectorXd::Zero(6);g.fixedCharge_C_per_m=VectorXd::Zero(6);
    for(const auto& cell:f.mesh.cells())if(cell.region_id==0)for(Index i:cell.node_ids)g.siliconArea_m2[i]+=1e-14/6.;
    g.poissonEdge_F_per_m=VectorXd::Zero(f.mesh.numEdges());g.transportWeight=g.poissonEdge_F_per_m;
    SiliconThermalParameters params;params.srhTau300_s={1e100,1e100};params.augerElectron_m6_per_s={0.,0.,0.};params.augerHole_m6_per_s={0.,0.,0.};
    SiliconThermalPhysics physics(params);LatticeConductivity law;
    ElectrothermalAssembler coupled(f.mesh,f.doping,g,LatticeHeatAssembler(f.mesh,1.,{{0,law},{1,law}},{}),{},physics);
    ElectrothermalBoundary bc;bc.neutralContactBias_V={{0,40.}};bc.holeRecombination={{0,{1.93e4,1e-7}}};
    for(Real t:{300.,515.}){
        VectorXd x=VectorXd::Zero(24),ref=VectorXd::Constant(6,40.);
        const Real psi=coupled.neutralPotential(0,40.,t).first;
        for(int i=0;i<6;++i){x[4*i]=psi;x[4*i+3]=t;}
        const auto eq=coupled.assemble(x,bc,ref,ref);
        CHECK(eq.residual[2]==0.);CHECK(eq.holeOutflow_A_per_m[0]==0.);
        // A finite contact must use its algebraic potential at the current
        // temperature, not a root cached from the previous temperature. At
        // equilibrium the constrained temperature direction also cancels the
        // free-hole residual derivative. Repairing psi must not clamp fp.
        auto heated=x;heated[3]=t+.01;
        const auto target=coupled.neutralPotential(0,40.,heated[3]);
        REQUIRE(target.first!=heated[0]);
        const auto stale=coupled.assemble(heated,bc,ref,ref);
        REQUIRE(stale.residual[0]!=0.);REQUIRE(stale.residual[2]!=0.);
        heated[0]=target.first;
        const auto consistent=coupled.assemble(heated,bc,ref,ref);
        CHECK(consistent.residual[0]==0.);CHECK(consistent.residual[2]==0.);
        CHECK(std::abs(consistent.jacobian.coeff(2,3)+consistent.jacobian.coeff(2,0)*target.second)
            <1e-12*std::max(std::abs(consistent.jacobian.coeff(2,3)),1e-100));
        heated[2]=1e-18;
        CHECK(coupled.assemble(heated,bc,ref,ref).residual[2]!=0.);
        for(Real increment:{-1e-18,1e-18}){
            auto trial=x;trial[2]=increment;const auto a=coupled.assemble(trial,bc,ref,ref);
            REQUIRE(a.holeOutflow_A_per_m[0]!=0.);
            CHECK(a.holeOutflow_A_per_m[0]*increment<0.);
            CHECK(a.residual[2]==Catch::Approx(-a.holeOutflow_A_per_m[0]).epsilon(1e-12));
            CHECK(a.residual[2]==Catch::Approx(eq.jacobian.coeff(2,2)*increment).epsilon(1e-11));
            auto blocked=bc;blocked.holeRecombination[0].velocity_m_per_s=0.;
            CHECK(coupled.assemble(trial,blocked,ref,ref).holeOutflow_A_per_m[0]==0.);
        }
    }
    VectorXd x=VectorXd::Zero(24);for(int i=0;i<6;++i)x[4*i+3]=300.;
    auto bad=bc;bad.holeRecombination[0].velocity_m_per_s=-1.;CHECK_THROWS_AS(coupled.assemble(x,bad),std::invalid_argument);
    bad=bc;bad.holeRecombination[0].boundaryLength_m=0.;CHECK_THROWS_AS(coupled.assemble(x,bad),std::invalid_argument);
    bad=bc;bad.holeQf_V[0]=40.;CHECK_THROWS_AS(coupled.assemble(x,bad),std::invalid_argument);
    bad=bc;bad.neutralContactBias_V.clear();CHECK_THROWS_AS(coupled.assemble(x,bad),std::invalid_argument);
}

TEST_CASE("Electrothermal tangent contact derivative matches fixed-state bias differences", "[thermal][electrothermal][predictor]") {
    Fixture f;
    ElectrothermalGeometry g;g.siliconArea_m2=VectorXd::Zero(6);g.fixedCharge_C_per_m=VectorXd::Zero(6);
    for(const auto& cell:f.mesh.cells())if(cell.region_id==0)for(Index i:cell.node_ids)g.siliconArea_m2[i]+=1e-14/6.;
    g.poissonEdge_F_per_m=VectorXd::Zero(f.mesh.numEdges());g.transportWeight=g.poissonEdge_F_per_m;
    LatticeConductivity law;
    ElectrothermalAssembler coupled(f.mesh,f.doping,g,LatticeHeatAssembler(f.mesh,1.,{{0,law},{1,law}},{}),{});
    for(bool finite:{false,true})for(Real t:{300.,515.})for(Real bias:{0.,40.}) {
        ElectrothermalBoundary bc;bc.neutralContactBias_V={{0,bias}};
        if(finite)bc.holeRecombination={{0,{1.93e4,1e-7}}};
        VectorXd x=VectorXd::Zero(24),ref=VectorXd::Constant(6,bias);
        const Real psi=coupled.neutralPotential(0,bias,t).first;
        for(int i=0;i<6;++i){x[4*i]=psi+.001;x[4*i+2]=.0001;x[4*i+3]=t;}
        const auto derivative=experimental::electrothermalContactBiasDerivative(6,bc,{0});
        constexpr Real h=1e-4;
        auto plus=bc,minus=bc;plus.neutralContactBias_V[0]+=h;minus.neutralContactBias_V[0]-=h;
        const VectorXd difference=(coupled.assemble(x,plus,ref,ref,false).residual-
            coupled.assemble(x,minus,ref,ref,false).residual)/(2*h);
        CHECK((difference-derivative).lpNorm<Eigen::Infinity>()<1e-7);
        CHECK(derivative[2]==(finite?0.:-1.));
        auto bad=bc;bad.potential_V[0]=psi;
        CHECK_THROWS_AS(experimental::electrothermalContactBiasDerivative(6,bad,{0}),std::invalid_argument);
    }
}
