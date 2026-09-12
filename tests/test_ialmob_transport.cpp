#include "vela/equation/ElectrothermalAssembler.h"
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
    LatticeConductivity law;law.model=LatticeConductivity::Model::InverseQuadratic;
    law.numerator=100.;law.denominator={-.0393,.00155,1.82e-6};
    LatticeHeatAssembler heat(f.mesh,1.,{{0,law},{1,law}},{{{0,1},300.,2e6}});
    ElectrothermalAssembler coupled(f.mesh,f.doping,g,heat,f.mobility);
    VectorXd x(24);for(int i=0;i<6;++i){x[4*i]=.48+.012*i;x[4*i+1]=.003*i;x[4*i+2]=.96-.002*i;x[4*i+3]=350.+12.*i;}
    for(bool constrained:{false,true}){
        ElectrothermalBoundary bc;
        if(constrained){bc.neutralContactBias_V={{0,0.}};bc.potential_V={{4,.2}};bc.temperature_K={{5,400.}};}

        const auto base=coupled.assemble(x,bc);
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
        CHECK(std::abs(base.holeOutflow_A_per_m.sum())<1e-10);
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

TEST_CASE("Four-equation silicon resistor closes self-heating and neutral contact temperature response", "[thermal][electrothermal][newton]") {
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
    VectorXd x(4*n);Real psi=coupled.neutralPotential(0,0.,300.).first;
    for(Index i=0;i<n;++i){x[4*i]=psi;x[4*i+1]=0.;x[4*i+2]=0.;x[4*i+3]=300.;}
    // Scale each physical equation separately. This test gate is for a
    // manufactured resistor, not a replacement for LDMOS electrical gates.
    VectorXd scales(4*n);for(Index i=0;i<n;++i){scales[4*i]=1e11;scales[4*i+1]=.01;scales[4*i+2]=.01;scales[4*i+3]=.01;}
    for(const auto& [i,v]:bc.neutralContactBias_V)for(int k=0;k<3;++k)scales[4*i+k]=1.;
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
