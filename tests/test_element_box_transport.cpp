#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>
#include "vela/mesh/DelaunayBox.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/simulation/ConfigParsing.h"
#include "vela/solver/NewtonSolver.h"
#include "vela/post/ContactCurrent.h"
#include "vela/equation/SplitDDRuntime.h"

using namespace vela;
namespace {
DeviceMesh pairMesh(double lower=.9,bool boundary=false,bool regionBoundary=false)
{
    DeviceMesh mesh;
    const std::array<std::array<double,2>,4> points={{{0,0},{1,0},{.5,.3},{.5,-lower}}};
    for(Index i=0;i<4;++i){Node n;n.id=i;n.x=points[i][0]*1e-6;n.y=points[i][1]*1e-6;mesh.addNode(n);}
    for(Index i=0;i<(boundary?1:2);++i){Cell c;c.id=i;c.type=CellType::Tri3;c.region_id=regionBoundary?i:0;c.node_ids=i==0?std::vector<Index>{0,1,2}:std::vector<Index>{1,0,3};mesh.addCell(c);}
    Region r;r.id=0;r.name="Si";r.material="Si";r.cell_ids={0};if(!boundary&&!regionBoundary)r.cell_ids.push_back(1);mesh.addRegion(r);
    if(regionBoundary){r.id=1;r.name="oxide";r.material="SiO2";r.cell_ids={1};mesh.addRegion(r);}
    mesh.buildEdges();return mesh;
}
BoxGeometryBuilder::Options options()
{
    BoxGeometryBuilder::Options o;o.poissonPermittivityPolicy=BoxGeometryBuilder::PoissonPermittivityPolicy::CellMaterial;
    o.cellBoxPolicy=BoxGeometryBuilder::CellBoxPolicy::DelaunayTransfer;return o;
}
MobilityModelConfig mobilityConfig()
{
    MobilityModelConfig c;c.model="masetti";c.dopingConcentrationBasis="total_impurity";c.edgeAveraging="element_box";return c;
}
DDScalingSpec scaling()
{
    DDScalingSpec s;s.regionResolvedInterfaceAssembly.transportEdgeGeometry="element_box";return s;
}
DeviceMesh surfaceMesh() {
    DeviceMesh mesh;const std::array<std::array<double,2>,5> xy{{{0,0},{1,0},{.5,.8},{.5,-.8},{-.5,.8}}};
    for(Index i=0;i<5;++i){Node n;n.id=i;n.x=xy[i][0]*1e-8;n.y=xy[i][1]*1e-8;mesh.addNode(n);}
    const std::array<std::array<Index,3>,3> triangles{{{0,1,2},{0,2,4},{1,0,3}}};
    for(Index i=0;i<3;++i){Cell cell;cell.id=i;cell.type=CellType::Tri3;cell.region_id=i==2?1:0;cell.node_ids.assign(triangles[i].begin(),triangles[i].end());mesh.addCell(cell);}
    Region si;si.id=0;si.name="Si";si.material="Si";si.cell_ids={0,1};mesh.addRegion(si);
    Region ox;ox.id=1;ox.name="oxide";ox.material="SiO2";ox.cell_ids={2};mesh.addRegion(ox);
    mesh.buildEdges();mesh.buildBoxGeometry(options());Contact contact;contact.name="test";contact.node_ids={1};mesh.addContact(contact);return mesh;
}
MobilityModelConfig elementLombardiConfig() {
    return mobilityModelConfigFromJson(nlohmann::json{{"model","phumob_lombardi"},{"edge_averaging","element_box_phumob"},{"doping_concentration_basis","total_impurity"},{"surface",{{"discretization","element_distance_gradient"}}}});
}
}

TEST_CASE("Element Lombardi uses the nonunit gradient of true interface distances", "[element_lombardi]") {
    auto mesh=surfaceMesh();MaterialDatabase db;auto cfg=elementLombardiConfig();auto mats=detail::buildCellMaterials(mesh,db,300.);
    VectorXd psi(5);psi<<.01,.04,-.1,.02,.2;
    detail::updateSurfaceMobilityCellGeometry(cfg,mesh,detail::buildEdgeCellMap(mesh),psi,1.,&mats);
    REQUIRE(cfg.surface.nodeDistances[0]==0.);REQUIRE(cfg.surface.nodeDistances[1]==0.);
    REQUIRE(cfg.surface.nodeDistances[4]==Catch::Approx(std::sqrt(.89)*1e-8).epsilon(1e-14));
    REQUIRE(std::abs(std::hypot(cfg.surface.cellNormalX[1],cfg.surface.cellNormalY[1])-1.)>1e-3);
    auto shifted=psi.array()+4.;VectorXd translated=shifted;
    REQUIRE(detail::liveSurfaceNormalFieldForCell(cfg,mesh,translated,1)==Catch::Approx(detail::liveSurfaceNormalFieldForCell(cfg,mesh,psi,1)).epsilon(1e-13));
    REQUIRE_THROWS(mobilityModelConfigFromJson(nlohmann::json{{"model","phumob"},{"surface",{{"discretization","element_distance_gradient"}}}}));
}

TEST_CASE("Element Lombardi complete Jacobian agrees with wide state differences", "[element_lombardi]") {
    auto mesh=surfaceMesh();MaterialDatabase db;DopingModel doping(5);
    for(Index i=0;i<5;++i)doping.setNodeDoping(i,1e23*(i+1),2e22*(5-i));
    auto cfg=elementLombardiConfig();CoupledDDAssembler assembler(mesh,db,doping,constants::Vt_300,cfg,recombinationModelConfig({"srh"}),{}, {},{}, {},scaling());
    assembler.enableSplitDDState(true);const Real ni=db.getMaterial("Si",300.).ni;
    for(bool minority:{false,true}) {
        CoupledDDState state;state.psi=VectorXd(5);state.psi<<.01,.04,-.1,.02,.2;state.phin=VectorXd(5);state.phip=VectorXd(5);
        for(int i=0;i<5;++i){state.phin(i)=state.psi(i)-constants::Vt_300*std::log(3e22*(i+1)/ni);state.phip(i)=state.psi(i)+constants::Vt_300*std::log((minority?1e8:2e22)*(5-i)/ni);}
        const auto x=assembler.pack(state);const auto J=assembler.assembleJacobian(x,{});auto runtime=assembler.splitRuntime();
        for(int col=0;col<15;++col)for(double step:{1e-7,1e-12,1e-25}) {
            VectorXd delta=VectorXd::Zero(15);delta(col)=step;VectorXd plus,minus;
            const VectorXd fd=runtime->symmetricDifference(x,delta,{},plus,minus)/step;
            for(int b=0;b<3;++b) {
                const VectorXd analytic=J.col(col).toDense().segment(5*b,5),numeric=fd.segment(5*b,5);
                if(analytic.norm()==0.){REQUIRE(numeric.norm()<1e-80);continue;}
                CAPTURE(minority,col,b,step);
                REQUIRE((analytic-numeric).norm()/analytic.norm()<2e-6);
                for(int row=0;row<5;++row) {
                    CAPTURE(row);
                    const Real magnitude=std::max(std::abs(analytic(row)),std::abs(numeric(row)));
                    if(magnitude==0.)continue;
                    REQUIRE(std::abs(analytic(row)-numeric(row))/magnitude<2e-6);
                }
            }
        }
        if(minority)REQUIRE(J.coeff(5,14)!=0.); // nonlocal hole column: no local SRH term can hide it
        auto r=runtime->residual(x,{});auto hi=x;VectorXd low=VectorXd::Zero(15);
        for(int i=0;i<15;++i){hi(i)=std::nextafter(x(i),INFINITY);low(i)=x(i)-hi(i);}
        runtime->setLow(low);REQUIRE((runtime->residual(hi,{}).array()==r.array()).all());runtime->setLow(VectorXd::Zero(15));
        DDSolution solution;runtime->save(x,solution);ContactCurrent port(mesh,db,doping,cfg,300.,scaling());
        auto actual=port.compute(solution,"test"),reference=port.computeFromResidual(assembler,x,"test");
        REQUIRE(actual.totalCurrent==Catch::Approx(reference.totalCurrent).epsilon(2e-12));
    }
}

TEST_CASE("Element HFS applies the cutoff and differentiates the active Canali branch", "[element_hfs]") {
    FieldMobilityParameters p{1.07e5,1.109};
    for(Real f:{0.,99.9999,100.0001,1e5,1e8}) {
        const Real m=.08;auto v=evaluateElementCanali(m,f,p,100.);
        if(f<100.) {REQUIRE(v.mobility==m);REQUIRE(v.lowDerivative==1.);REQUIRE(v.fieldDerivative==0.);}
        else {
            REQUIRE(v.mobility<m);
            const auto dm=(elementCanali((long double)m+1e-8L,(long double)f,p,100.)-elementCanali((long double)m-1e-8L,(long double)f,p,100.))/2e-8L;
            const auto df=(elementCanali((long double)m,(long double)f+1e-6L,p,100.)-elementCanali((long double)m,(long double)f-1e-6L,p,100.))/2e-6L;
            REQUIRE(std::abs(v.lowDerivative/static_cast<Real>(dm)-1)<1e-8);
            REQUIRE(std::abs(v.fieldDerivative/static_cast<Real>(df)-1)<1e-4);
        }
    }
    REQUIRE(elementCanali(.08,99.9999,p,100.)>elementCanali(.08,100.0001,p,100.));
}

TEST_CASE("Element HFS matrix includes contact, projected and interior cell chains", "[element_hfs]") {
    DeviceMesh mesh;
    for(int y=-1;y<4;++y)for(int x=0;x<4;++x){Node n;n.id=mesh.numNodes();n.x=x*1e-8;n.y=y*1e-8;mesh.addNode(n);}
    Region si;si.id=0;si.name="Si";si.material="Si";Region ox;ox.id=1;ox.name="oxide";ox.material="SiO2";
    for(int y=0;y<4;++y)for(int x=0;x<3;++x)for(int t=0;t<2;++t) {
        Index a=y*4+x,b=a+1,c=a+4,d=c+1;Cell cell;cell.id=mesh.numCells();cell.type=CellType::Tri3;cell.region_id=y==0?1:0;
        cell.node_ids=t?std::vector<Index>{a,d,c}:std::vector<Index>{a,b,d};
        (y==0?ox:si).cell_ids.push_back(cell.id);mesh.addCell(cell);
    }
    mesh.addRegion(si);mesh.addRegion(ox);mesh.buildEdges();mesh.buildBoxGeometry(options());Contact contact;contact.name="test";contact.node_ids={19};mesh.addContact(contact);
    auto cfg=elementLombardiConfig();cfg.model="phumob_field_lombardi";cfg.highFieldDrivingForce="quasi_fermi_gradient";cfg.highFieldGradientDiscretization="element_vertex_partial_layer";
    MaterialDatabase db;DopingModel doping(20);for(Index i=0;i<20;++i)doping.setNodeDoping(i,1e23*(i+1),2e22*(20-i));
    auto mats=detail::buildCellMaterials(mesh,db,300.);detail::updateElementHighFieldGeometry(cfg,mesh,detail::buildEdgeCellMap(mesh),mats);
    int contacts=0,zeros=0,full=0;
    for(Index cid:si.cell_ids){auto g=cfg.highFieldCells[cid];contacts+=g.contact;zeros+=g.stencil==std::array<Real,4>{};full+=!g.contact&&std::abs(g.stencil[0]*g.stencil[3]-g.stencil[1]*g.stencil[2])>1.;}
    REQUIRE(contacts>0);REQUIRE(zeros>0);REQUIRE(full>0);
    CoupledDDAssembler assembler(mesh,db,doping,constants::Vt_300,cfg,recombinationModelConfig({"srh"}),{}, {},{}, {},scaling());assembler.enableSplitDDState(true);
    const Real ni=db.getMaterial("Si",300.).ni;
    for(bool minority:{false,true}) {
        CoupledDDState state;state.psi=VectorXd(20);state.phin=VectorXd(20);state.phip=VectorXd(20);
        for(int i=0;i<20;++i){state.psi(i)=.001*(i*i+i+1);state.phin(i)=state.psi(i)-constants::Vt_300*std::log(3e22*(i+1)/ni);state.phip(i)=state.psi(i)+constants::Vt_300*std::log((minority?1e8:2e22)*(20-i)/ni);}
        auto x=assembler.pack(state);auto J=assembler.assembleJacobian(x,{});auto runtime=assembler.splitRuntime();
        for(int col=0;col<60;++col)for(double step:{1e-7,1e-12,1e-25}) {
            VectorXd delta=VectorXd::Zero(60);delta(col)=step;VectorXd plus,minus;const VectorXd fd=runtime->symmetricDifference(x,delta,{},plus,minus)/step;
            for(int row=0;row<60;++row) {
                const Real an=J.coeff(row,col),scale=std::max(std::abs(an),std::abs(fd(row)));if(scale==0.)continue;
                CAPTURE(minority,col,row,step,an,fd(row));REQUIRE(std::abs(an-fd(row))/scale<2e-6);
            }
        }
        auto r=runtime->residual(x,{});VectorXd hi=x,lo=VectorXd::Zero(60);for(int i=0;i<60;++i){hi(i)=std::nextafter(x(i),INFINITY);lo(i)=x(i)-hi(i);}
        runtime->setLow(lo);REQUIRE((runtime->residual(hi,{}).array()==r.array()).all());runtime->setLow(VectorXd::Zero(60));
        DDSolution solution;runtime->save(x,solution);ContactCurrent port(mesh,db,doping,cfg,300.,scaling());
        REQUIRE(port.compute(solution,"test").totalCurrent==Catch::Approx(port.computeFromResidual(assembler,x,"test").totalCurrent).epsilon(2e-12));
    }
}

TEST_CASE("Zero element surface scattering recovers the unchanged bulk profile", "[element_lombardi]") {
    auto mesh=surfaceMesh();MaterialDatabase db;DopingModel doping(5);for(Index i=0;i<5;++i)doping.setNodeDoping(i,1e23*(i+1),1e22);
    auto cfg=elementLombardiConfig();cfg.electronLombardi.acousticFactor=cfg.holeLombardi.acousticFactor=0.;cfg.electronLombardi.roughnessFactor=cfg.holeLombardi.roughnessFactor=0.;
    auto plain=cfg;plain.model="phumob";plain.surface.discretization="legacy_cell_centroid";
    CoupledDDAssembler surface(mesh,db,doping,constants::Vt_300,cfg,recombinationModelConfig({"srh"}),{}, {},{}, {},scaling());
    CoupledDDAssembler bulk(mesh,db,doping,constants::Vt_300,plain,recombinationModelConfig({"srh"}),{}, {},{}, {},scaling());
    surface.enableSplitDDState(true);bulk.enableSplitDDState(true);CoupledDDState s{VectorXd::Constant(5,.01),VectorXd::Constant(5,-.1),VectorXd::Constant(5,.1)};s.psi(2)=.025;s.phin(1)=-.08;
    auto x=surface.pack(s);auto y=bulk.pack(s);auto lhs=surface.residual(x,{}),rhs=bulk.residual(y,{});
    REQUIRE((lhs-rhs).norm()/rhs.norm()<1e-13);
}

TEST_CASE("Delaunay transfer preserves signed edge totals and nodal box measures", "[element_box]")
{
    auto mesh=pairMesh();const auto originalVolumes=detail::computeNodeVolumes(mesh);
    const auto raw=detail::computeTransportSignedAverageBoxNodeVolumes(mesh,detail::buildCellMaterials(mesh,MaterialDatabase{},300));
    mesh.buildBoxGeometry(options());const auto g0=mesh.poissonCellEdgeCoefficients(0),g1=mesh.poissonCellEdgeCoefficients(1);
    const double expected=(.3*.3-.25)/(2*.3)+(.9*.9-.25)/(2*.9);
    REQUIRE(g0[0]==0.);REQUIRE(g1[0]==Catch::Approx(expected).epsilon(1e-13));
    REQUIRE(mesh.lastGeometryBuildReport().transferredCellBoxEdges==1);
    std::array<double,4> total{};
    for(Index c=0;c<2;++c){const auto M=detail::cellBoxNodeMeasures(mesh,c);for(int k=0;k<3;++k)total[mesh.getCell(c).node_ids[k]]+=M[k];}
    for(int i=0;i<4;++i)REQUIRE(total[i]==Catch::Approx(raw[i]).epsilon(1e-13));
    const auto M=detail::cellBoxNodeMeasures(mesh,0);
    REQUIRE(M[0]+M[1]+M[2]>1.5e-13); // transfer changes cell partitions, not global nodal sums
    REQUIRE(detail::computeNodeVolumes(mesh)==originalVolumes);
}

TEST_CASE("Delaunay generation rejects unsupported edges and conflicting inputs", "[element_box]")
{
    auto nondel=pairMesh(.3);REQUIRE_THROWS(nondel.buildBoxGeometry(options()));
    auto boundary=pairMesh(.9,true);REQUIRE_THROWS(boundary.buildBoxGeometry(options()));
    auto material=pairMesh(.9,false,true);REQUIRE_THROWS(material.buildBoxGeometry(options()));
    auto mesh=pairMesh();auto o=options();o.poissonPermittivityPolicy=BoxGeometryBuilder::PoissonPermittivityPolicy::LegacyAverage;
    REQUIRE_THROWS(mesh.buildBoxGeometry(o));o=options();o.poissonCellEdgeCoefficients.resize(2);REQUIRE_THROWS(mesh.buildBoxGeometry(o));
    REQUIRE_THROWS(parseBoxGeometryOptions(nlohmann::json{{"mesh_geometry",{{"cell_box_policy","unknown"}}}}));
}

TEST_CASE("Element box mobility weights nodal Masetti values with actual cell measures", "[element_box]")
{
    auto mesh=pairMesh();mesh.buildBoxGeometry(options());MaterialDatabase db;DopingModel doping(4);
    doping.setNodeDoping(0,1e21,0);doping.setNodeDoping(1,2e24,0);doping.setNodeDoping(2,1e23,0);doping.setNodeDoping(3,1e22,0);
    const auto cfg=mobilityConfig();const auto model=makeMobilityModel(cfg);const auto mats=detail::buildCellMaterials(mesh,db,300);
    const auto ec=detail::buildEdgeCellMap(mesh);
    for(auto carrier:{CarrierType::Electron,CarrierType::Hole}) {
        const auto mat=db.getMaterial("Si");
        auto mu=[&](double n){return carrier==CarrierType::Electron?model->electronMobility(mat,n,0,0):model->holeMobility(mat,n,0,0);};
        // Upper cell g01=0 and equal oblique edges give weights (1/4,1/4,1/2).
        const double expected=.25*mu(1e21)+.25*mu(2e24)+.5*mu(1e23);
        REQUIRE(detail::elementBoxCellMobility(mesh,doping,*model,mat,0,carrier)==Catch::Approx(expected).epsilon(1e-13));
        Index edge=0;while(!(mesh.getEdge(edge).n0==0&&mesh.getEdge(edge).n1==2))++edge;
        REQUIRE(detail::elementBoxEdgeMobility(ec,mesh,doping,*model,mats,edge,carrier,cfg)==Catch::Approx(expected).epsilon(1e-13));
        REQUIRE(std::abs(expected/mu(.25e21+.5e24+.5e23)-1)>1e-3);
    }
    auto bad=cfg;bad.model="masetti_field";REQUIRE_THROWS(detail::validateElementBoxMobilityContext(mesh,bad,scaling().regionResolvedInterfaceAssembly));
    REQUIRE_THROWS(detail::validateElementBoxMobilityContext(mesh,cfg,{}));
    auto s=scaling();s.regionResolvedInterfaceAssembly.transportEdgeCouplingRatios=std::vector<Real>(mesh.numEdges(),1.);
    REQUIRE_THROWS(detail::computeEffectiveTransportEdgeCouplings(mesh,ec,mats,s.regionResolvedInterfaceAssembly));
    REQUIRE_THROWS(mobilityModelConfigFromJson(nlohmann::json{{"model","phumob"},{"edge_averaging","element_box"}}));
}

TEST_CASE("Element box Newton derivatives retain all state-independent mobility effects", "[element_box]")
{
    auto mesh=pairMesh();mesh.buildBoxGeometry(options());MaterialDatabase db;DopingModel doping(4);
    for(Index i=0;i<4;++i)doping.setNodeDoping(i,1e22*(i+1),0);
    RecombinationModelConfig recombination;
    CoupledDDAssembler a(mesh,db,doping,constants::kb*300/constants::q,mobilityConfig(),recombination,{}, {},{}, {},scaling());
    CoupledDDState state{VectorXd::Zero(4),VectorXd::Zero(4),VectorXd::Zero(4)};
    state.psi<<.01,.015,-.01,.008;state.phin<<.001,-.002,.003,.004;state.phip<<.002,.004,-.002,.003;
    const auto x=a.pack(state);const auto J=a.assembleJacobian(x,{});VectorXd v(12);for(int i=0;i<12;++i)v[i]=std::sin(.4*(i+1));
    const VectorXd analytic=J*v;
    for(double h:{1e-5,5e-6,2.5e-6}) {
        const VectorXd fd=(a.residual(x+h*v,{})-a.residual(x-h*v,{}))/(2*h);
        for(int b=0;b<3;++b)REQUIRE((fd.segment(4*b,4)-analytic.segment(4*b,4)).norm()/analytic.segment(4*b,4).norm()<1e-6);
    }
    REQUIRE_THROWS(DDAssembler(mesh,db,doping,.02585,mobilityConfig(),recombination));
}

TEST_CASE("PhuMob box candidate requires explicit geometry and complete live populations", "[element_box][phumob_box]")
{
    auto mesh=pairMesh();mesh.buildBoxGeometry(options());MaterialDatabase db;DopingModel doping(4);
    auto cfg=mobilityModelConfigFromJson(nlohmann::json{{"model","phumob"},
        {"edge_averaging","element_box_phumob"},{"doping_concentration_basis","total_impurity"}});
    REQUIRE_NOTHROW(detail::validateElementBoxMobilityContext(mesh,cfg,scaling().regionResolvedInterfaceAssembly));
    REQUIRE_THROWS(mobilityModelConfigFromJson(nlohmann::json{{"model","phumob_lombardi"},{"edge_averaging","element_box_phumob"}}));
    REQUIRE_THROWS(detail::validateElementBoxMobilityContext(mesh,cfg,{}));
    for(Index i=0;i<4;++i)doping.setNodeDoping(i,1e23,2e21);
    const auto model=makeMobilityModel(cfg);const auto mats=detail::buildCellMaterials(mesh,db,300);
    const auto ec=detail::buildEdgeCellMap(mesh);
    VectorXd n=VectorXd::Constant(4,8e22),p=VectorXd::Constant(4,3e20);
    detail::EdgeMobilityCarrierState populations{n(0),n(2),p(0),p(2),&n,&p};
    Index edge=0;while(!(mesh.getEdge(edge).n0==0&&mesh.getEdge(edge).n1==2))++edge;
    REQUIRE_THROWS(detail::elementBoxEdgeMobility(ec,mesh,doping,*model,mats,edge,CarrierType::Electron,cfg));
    const Real first=detail::elementBoxEdgeMobility(ec,mesh,doping,*model,mats,edge,CarrierType::Electron,cfg,&populations);
    // Change only the third vertex, keeping the legacy endpoint state fixed.
    n(1)*=4.;p(1)*=2.;
    const Real second=detail::elementBoxEdgeMobility(ec,mesh,doping,*model,mats,edge,CarrierType::Electron,cfg,&populations);
    const auto scalar=[&](Index i){return evaluatePhuMobScalar(CarrierType::Electron,
        {doping.donors(i),doping.acceptors(i),n(i),p(i),300.},cfg.phuMob).mobility;};
    // The upper cell has independently known box fractions (1/4,1/4,1/2).
    const Real expected=.25*scalar(0)+.25*scalar(1)+.5*scalar(2);
    REQUIRE(std::abs(second/expected-1)<2e-13);
    REQUIRE(std::abs(second/first-1)>1e-3);
    REQUIRE_THROWS(DDAssembler(mesh,db,doping,constants::Vt_300,cfg,{}));
}

TEST_CASE("PhuMob box assembly retains a weak third-vertex column and port conservation", "[element_box][phumob_box]")
{
    auto mesh=pairMesh();mesh.buildBoxGeometry(options());
    Contact contact;contact.name="test";contact.node_ids={0};mesh.addContact(contact);
    MaterialDatabase db;DopingModel doping(4);
    for(Index i=0;i<4;++i)doping.setNodeDoping(i,5e27,2e23);
    auto cfg=mobilityConfig();cfg.model="phumob";cfg.edgeAveraging="element_box_phumob";
    const auto rec=recombinationModelConfig({"none"});
    CoupledDDAssembler assembler(mesh,db,doping,constants::Vt_300,cfg,rec,{}, {},{}, {},scaling());
    const Real ni=db.getMaterial("Si",300.).ni;
    CoupledDDState state;
    state.psi=VectorXd(4);state.psi<<0.,.01,-.02,0.;
    state.phin=state.psi.array()-constants::Vt_300*std::log(5e27/ni);
    state.phip=state.psi.array()+constants::Vt_300*std::log(1e3/ni);
    const auto x=assembler.pack(state);const auto J=assembler.assembleJacobian(x,{});
    const auto edges=assembler.sgEdgeFluxDiagnostics(x,{});
    const auto edge=std::find_if(edges.begin(),edges.end(),[](const auto& e){return e.node0==0&&e.node1==1;});
    REQUIRE(edge!=edges.end());
    // Edge 0-3 is at flat quasi-Fermi potential. Thus column p_3 in row n_0
    // only receives mobility coupling through the opposite edge 0-1.
    // Independent triangle cotangents and squared lengths give this fraction.
    const Real diagonalG=(.9*.9-.25)/(2*.9)+(.3*.3-.25)/(2*.3);
    const Real outerG=.5/(2*.9);
    const Real thirdMeasure=2*outerG*(.25+.9*.9)/4;
    const Real volume=2*(diagonalG/4+outerG*(.25+.9*.9)/4)+thirdMeasure;
    const Real hpDerivative=1.8252671486063443e-28;
    const Real expected=edge->electronFlux/edge->electronMobility_m2_V_s*
        hpDerivative*(thirdMeasure/volume)/constants::Vt_300;
    REQUIRE(expected!=0.);
    REQUIRE(std::abs(J.coeff(4,11)/expected-1)<1e-10);
    for(int block:{1,2})for(int col=0;col<12;++col) {
        Real sum=0,scale=0;
        for(int row=4*block;row<4*(block+1);++row){sum+=J.coeff(row,col);scale+=std::abs(J.coeff(row,col));}
        REQUIRE(std::abs(sum)<=5e-15*scale);
    }
    DDSolution solution;solution.psi=state.psi;solution.phin=state.phin;solution.phip=state.phip;
    solution.n=assembler.electronDensity(x);solution.p=assembler.holeDensity(x);
    ContactCurrent port(mesh,db,doping,cfg,300.,scaling());
    const auto actual=port.compute(solution,"test"),reference=port.computeFromResidual(assembler,x,"test");
    REQUIRE(std::abs(actual.totalCurrent/reference.totalCurrent-1)<2e-12);
    REQUIRE_THROWS(assembler.transportEdgeJacobianDiagnostics(x,{},1e-6));
}

TEST_CASE("PhuMob box full Jacobian columns converge under independent state differences", "[element_box][phumob_box]")
{
    auto mesh=pairMesh();mesh.buildBoxGeometry(options());MaterialDatabase db;DopingModel doping(4);
    for(Index i=0;i<4;++i)doping.setNodeDoping(i,1e23*(i+1),2e22*(4-i));
    auto cfg=mobilityConfig();cfg.model="phumob";cfg.edgeAveraging="element_box_phumob";
    CoupledDDAssembler assembler(mesh,db,doping,constants::Vt_300,cfg,recombinationModelConfig({"none"}),{}, {},{}, {},scaling());
    const Real ni=db.getMaterial("Si",300.).ni;
    CoupledDDState state;state.psi=VectorXd(4);state.psi<<.01,.015,-.01,.008;
    state.phin=VectorXd(4);state.phip=VectorXd(4);
    for(int i=0;i<4;++i){
        state.phin(i)=state.psi(i)-constants::Vt_300*std::log(3e22*(i+1)/ni);
        state.phip(i)=state.psi(i)+constants::Vt_300*std::log(2e22*(4-i)/ni);
    }
    const auto x=assembler.pack(state);const auto J=assembler.assembleJacobian(x,{});
    for(int column=0;column<12;++column) {
        VectorXd direction=VectorXd::Zero(12);direction(column)=1.;const VectorXd analytic=J*direction;
        for(Real h:{1e-5,3e-6,1e-6}) {
            const VectorXd fd=(assembler.residual(x+h*direction,{})-assembler.residual(x-h*direction,{}))/(2*h);
            for(int block=0;block<3;++block) {
                const Real scale=analytic.segment(4*block,4).norm();
                if(scale==0.)continue;
                CAPTURE(column,block,h);
                REQUIRE((fd.segment(4*block,4)-analytic.segment(4*block,4)).norm()/scale<2e-6);
            }
        }
    }
}

TEST_CASE("Independent SRH area preserves Poisson and uses one residual Jacobian measure", "[srh_volume]") {
    auto mesh=pairMesh();mesh.buildBoxGeometry(options());MaterialDatabase db;DopingModel doping(4);
    auto cfg=mobilityConfig();cfg.model="phumob";cfg.edgeAveraging="element_box_phumob";
    auto rec=recombinationModelConfig({"srh"},1e-12,3e-12);
    auto baseScaling=scaling();baseScaling.regionResolvedInterfaceAssembly.poissonChargeNodeVolume="signed_transport";
    CoupledDDAssembler base(mesh,db,doping,constants::Vt_300,cfg,rec,{}, {},{}, {},baseScaling);
    CoupledDDState state;state.psi=VectorXd::Constant(4,.01);state.phin=VectorXd::Constant(4,-.04);state.phip=VectorXd::Constant(4,.04);
    const auto x=base.pack(state);const auto mats=detail::buildCellMaterials(mesh,db,300.);
    const auto signedVolumes=detail::computeTransportSignedAverageBoxNodeVolumes(mesh,mats);
    for(double fraction:{-.01,0.,.01,1.}) {
        auto spec=baseScaling;spec.regionResolvedInterfaceAssembly.srhSignedTransportVolumeFraction=fraction;
        CoupledDDAssembler candidate(mesh,db,doping,constants::Vt_300,cfg,rec,{}, {},{}, {},spec);
        for(bool split:{false,true}) {
            base.enableSplitDDState(split);candidate.enableSplitDDState(split);
            const auto r0=base.residual(x,{}),r=candidate.residual(x,{});
            REQUIRE((r.head(4)-r0.head(4)).norm()==0.);
            for(int i=0;i<4;++i) {
                const double v=mesh.getNode(i).volume;
                const double ratio=fraction==1.?signedVolumes[i]/v:1.+fraction*(signedVolumes[i]/v-1.);
                for(int block:{1,2})REQUIRE(r(4*block+i)/r0(4*block+i)==Catch::Approx(ratio).epsilon(1e-12));
            }
            const auto J=candidate.assembleJacobian(x,{});
            for(int col=0;col<12;++col) {
                VectorXd d=VectorXd::Zero(12);d(col)=1e-7;
                VectorXd fd;
                if(split){VectorXd a,b;fd=candidate.splitRuntime()->symmetricDifference(x,d,{},a,b)/1e-7;}
                else fd=(candidate.residual(x+d,{})-candidate.residual(x-d,{}))/(2e-7);
                const VectorXd analytic=J.col(col);
                REQUIRE((fd-analytic).norm()/analytic.norm()<2e-6);
            }
        }
    }
    auto invalid=baseScaling;invalid.regionResolvedInterfaceAssembly.srhSignedTransportVolumeFraction=1.;
    REQUIRE_THROWS(CoupledDDAssembler(mesh,db,doping,constants::Vt_300,cfg,recombinationModelConfig({"srh","auger"}),{}, {},{}, {},invalid));
    invalid.regionResolvedInterfaceAssembly.transportNodeVolume=true;
    REQUIRE_THROWS(CoupledDDAssembler(mesh,db,doping,constants::Vt_300,cfg,rec,{}, {},{}, {},invalid));
    REQUIRE_THROWS(newtonConfigFromJson(nlohmann::json{{"region_resolved_interface_assembly",{{"srh_signed_transport_volume_fraction",1.01}}}}));
}
