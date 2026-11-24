#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <nlohmann/json.hpp>
#include "vela/mesh/DelaunayBox.h"
#include "vela/mesh/AverageBox.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/BVProcessProbe.h"
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

TEST_CASE("Average intersection recovers Delaunay transfer and removes a non-Delaunay face", "[element_box][average_box]")
{
    auto o=options();o.cellBoxPolicy=BoxGeometryBuilder::CellBoxPolicy::AverageIntersection;
    auto mesh=pairMesh();const auto reference=buildDelaunayCellBox(mesh);
    mesh.buildBoxGeometry(o);
    for(Index c=0;c<mesh.numCells();++c)
        for(int k=0;k<3;++k)
            REQUIRE(mesh.poissonCellEdgeCoefficients(c)[k]==
                Catch::Approx(reference.coefficients[c][k]).margin(1e-14));
    auto nonDelaunay=pairMesh(.3);
    nonDelaunay.buildBoxGeometry(o);
    REQUIRE(nonDelaunay.poissonCellEdgeCoefficients(0)[0]==0.);
    REQUIRE(nonDelaunay.poissonCellEdgeCoefficients(1)[0]==0.);
    // Symmetric diamond: endpoint face lengths / edge length are .3 and
    // 5/6, giving 17/30 after averaging. The diagonal must not conduct.
    for(Index c=0;c<2;++c) {
        for(int k=1;k<3;++k)
            REQUIRE(nonDelaunay.poissonCellEdgeCoefficients(c)[k]==Catch::Approx(17./30.));
        for(auto v:detail::cellBoxNodeMeasures(nonDelaunay,c)) REQUIRE(v>0.);
    }
    // Cocircular limit: both diagonal faces have zero length without an
    // arbitrary coefficient cutoff; oblique edges remain finite.
    // Move the upper vertex too, making the four vertices cocircular.
    // Recreate because DeviceMesh exposes const node access.
    DeviceMesh q;
    for(Index i=0;i<4;++i){Node n;n.id=i;n.x=(i==0?0:i==1?1:.5)*1e-6;n.y=(i==2?.5:i==3?-.5:0)*1e-6;q.addNode(n);}
    for(Index c=0;c<2;++c){Cell t;t.id=c;t.type=CellType::Tri3;t.region_id=0;t.node_ids=c==0?std::vector<Index>{0,1,2}:std::vector<Index>{1,0,3};q.addCell(t);}
    Region reg;reg.id=0;reg.name="Si";reg.material="Si";reg.cell_ids={0,1};q.addRegion(reg);
    q.buildEdges();q.buildBoxGeometry(o);
    REQUIRE(q.poissonCellEdgeCoefficients(0)[0]==0.);
    REQUIRE(q.poissonCellEdgeCoefficients(1)[0]==0.);
    const auto mats=detail::buildCellMaterials(q,MaterialDatabase{},300);
    const auto charge=detail::computePoissonChargeVolumes(q,mats,PoissonChargeVolumePolicy::MaterialLocalBox);
    Real total=0;for(auto v:charge)total+=v;
    REQUIRE(total==Catch::Approx(.5e-12).epsilon(1e-13));
    REQUIRE(detail::recombinationNodeVolumes(q,mats,"material_local_box")==charge);
}

TEST_CASE("Average intersection is invariant to winding and uniform length units", "[element_box][average_box]")
{
    const auto mesh=pairMesh(.3);const auto expected=buildAverageCellBox(mesh);
    DeviceMesh transformed;
    for(const auto& p:mesh.nodes()){auto n=p;n.x=p.x*1e6+3;n.y=p.y*1e6-7;transformed.addNode(n);}
    for(const auto& t:mesh.cells()){auto c=t;std::reverse(c.node_ids.begin(),c.node_ids.end());transformed.addCell(c);}
    const auto actual=buildAverageCellBox(transformed);
    for(Index c=0;c<2;++c)for(int k=0;k<3;++k)
        REQUIRE(actual.coefficients[c][k]==Catch::Approx(expected.coefficients[c][(1-k+3)%3]).margin(1e-13));
    auto o=options();o.cellBoxPolicy=BoxGeometryBuilder::CellBoxPolicy::AverageIntersection;
    o.poissonPermittivityPolicy=BoxGeometryBuilder::PoissonPermittivityPolicy::LegacyAverage;
    REQUIRE_THROWS_AS(transformed.buildBoxGeometry(o),std::invalid_argument);
    o=options();o.cellBoxPolicy=BoxGeometryBuilder::CellBoxPolicy::AverageIntersection;o.poissonCellEdgeCoefficients.resize(2);
    REQUIRE_THROWS_AS(transformed.buildBoxGeometry(o),std::invalid_argument);
    REQUIRE(parseBoxGeometryOptions(nlohmann::json{{"mesh_geometry",{{"cell_box_policy","average_intersection"},{"poisson_permittivity_policy","cell_material"}}}}).cellBoxPolicy==BoxGeometryBuilder::CellBoxPolicy::AverageIntersection);
}

TEST_CASE("Average box charge and recombination exclude dielectric ownership", "[element_box][average_box]")
{
    auto mesh=pairMesh(.3,false,true);auto o=options();
    o.cellBoxPolicy=BoxGeometryBuilder::CellBoxPolicy::AverageIntersection;
    mesh.buildBoxGeometry(o);
    const auto mats=detail::buildCellMaterials(mesh,MaterialDatabase{},300.);
    const auto charge=detail::computePoissonChargeVolumes(mesh,mats,PoissonChargeVolumePolicy::MaterialLocalBox);
    REQUIRE(charge[3]==0.); // exclusively oxide
    REQUIRE(charge[0]==Catch::Approx((17./30.)*.34e-12/4).epsilon(1e-13));
    REQUIRE(charge[1]==Catch::Approx(charge[0]).epsilon(1e-13));
    REQUIRE(charge[2]==Catch::Approx(2*charge[0]).epsilon(1e-13));
    REQUIRE(detail::recombinationNodeVolumes(mesh,mats,"material_local_box")==charge);
    REQUIRE(parsePoissonChargeVolumePolicy({{"discretization",{{"poisson_charge_volume_policy","material_local_box"}}}})==PoissonChargeVolumePolicy::MaterialLocalBox);
    auto invalid=pairMesh();auto cell=invalid.getCell(0);cell.id=2;invalid.addCell(cell);
    REQUIRE_THROWS_AS(buildAverageCellBox(invalid),std::invalid_argument);
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

TEST_CASE("Fermi PhuMob box includes degenerate and opposite-vertex density derivatives", "[element_box][phumob_box][fermi]") {
    const std::string transportMode = GENERATE("generalized_einstein", "effective_band_sg");
    CAPTURE(transportMode);
    const bool average=GENERATE(false,true);
    auto mesh=pairMesh(average?.3:.9);auto geometry=options();
    if(average)geometry.cellBoxPolicy=BoxGeometryBuilder::CellBoxPolicy::AverageIntersection;
    mesh.buildBoxGeometry(geometry);MaterialDatabase db;DopingModel doping(4);
    for(Index i=0;i<4;++i)doping.setNodeDoping(i,1e26*(i+1),2e25);
    auto cfg=mobilityConfig();cfg.model="phumob";cfg.edgeAveraging="element_box_phumob";
    auto scales=scaling();auto sources=recombinationModelConfig({"none"});
    if(average){scales.poissonChargeVolumePolicy=PoissonChargeVolumePolicy::MaterialLocalBox;
        sources=recombinationModelConfig({"srh","auger"});sources.volumePolicy="material_local_box";}
    CoupledDDAssembler assembler(mesh,db,doping,constants::Vt_300,cfg,sources,
        {},{}, {},{},scales,{},CarrierStatisticsConfig{"fermi_dirac", transportMode});
    const auto mat=db.getMaterial("Si",300.);
    CoupledDDState state{VectorXd(4),VectorXd(4),VectorXd(4)};
    for(int i=0;i<4;++i) {
        state.psi(i)=.02*i;
        state.phin(i)=electronQuasiFermiPotential(mat.ni,*mat.Nc_m3,state.psi(i),1e26*(i+1),constants::Vt_300,CarrierStatisticsConfig{"fermi_dirac", transportMode});
        state.phip(i)=holeQuasiFermiPotential(mat.ni,*mat.Nv_m3,state.psi(i),2e25*(4-i),constants::Vt_300,CarrierStatisticsConfig{"fermi_dirac", transportMode});
    }
    const auto x=assembler.pack(state);const auto J=assembler.assembleJacobian(x,{});
    for(int col=0;col<12;++col) {
        VectorXd d=VectorXd::Zero(12);d(col)=1.;
        const VectorXd analytic=J*d;
        for(Real h:{3e-6,1e-6}) {
            const VectorXd fd=(assembler.residual(x+h*d,{})-assembler.residual(x-h*d,{}))/(2*h);
            for(int block=0;block<3;++block) {
                const Real scale=std::max(analytic.segment(4*block,4).norm(),fd.segment(4*block,4).norm());
                if(scale==0.)continue;
                CAPTURE(col,block,h);
                REQUIRE((fd.segment(4*block,4)-analytic.segment(4*block,4)).norm()/scale<2e-5);
            }
        }
    }
}

TEST_CASE("Element PhuMob constant material bypasses density and field scattering", "[element_box][phumob_box]") {
    auto mesh=pairMesh();mesh.buildBoxGeometry(options());MaterialDatabase db;DopingModel doping(4);
    auto cfg=mobilityConfig();cfg.model="phumob";cfg.edgeAveraging="element_box_phumob";cfg.constantMaterials={"Si"};
    const auto model=makeMobilityModel(cfg);const auto mats=detail::buildCellMaterials(mesh,db,300.);
    const auto ec=detail::buildEdgeCellMap(mesh);VectorXd n=VectorXd::Constant(4,1e27),p=VectorXd::Constant(4,1e22);
    detail::EdgeMobilityCarrierState population{n(0),n(1),p(0),p(1),&n,&p};
    for(CarrierType carrier:{CarrierType::Electron,CarrierType::Hole}) {
        const Real mu=detail::elementBoxEdgeMobility(ec,mesh,doping,*model,mats,0,carrier,cfg,&population);
        REQUIRE(mu==Catch::Approx(carrier==CarrierType::Electron?mats[0].mun:mats[0].mup).epsilon(1e-14));
        REQUIRE(detail::elementBoxPhuMobDerivatives(ec,mesh,doping,mats,0,carrier,cfg,n,p).empty());
    }
}

TEST_CASE("PhuMob SG source signals population domain errors without regularizing density", "[element_box][phumob_avalanche]") {
    for (Real p:{0.,-1.,std::numeric_limits<Real>::infinity(),std::numeric_limits<Real>::quiet_NaN()}) {
        REQUIRE_THROWS_AS(detail::phumobSourceLog(p),PopulationDomainError);
        REQUIRE_THROWS_AS(detail::phumobSourceLog(detail::Tri3LocalForwardDual(p)),PopulationDomainError);
    }
    const Real tiny=std::numeric_limits<Real>::min();
    REQUIRE(detail::phumobSourceLog(tiny)==std::log(tiny));
    REQUIRE(detail::phumobSourceLog(1.)==0.);
}

TEST_CASE("Fermi PhuMob element avalanche preserves pair generation and its full Jacobian", "[element_box][fermi][phumob_avalanche]") {
    const std::string transportMode = GENERATE("generalized_einstein", "effective_band_sg");
    CAPTURE(transportMode);
    auto mesh=surfaceMesh();MaterialDatabase db;DopingModel doping(5);
    for(Index i=0;i<5;++i)doping.setNodeDoping(i,2e25*(i+1),1e24);
    auto cfg=elementLombardiConfig();cfg.model="phumob_field_lombardi";
    cfg.highFieldDrivingForce="quasi_fermi_gradient";
    cfg.highFieldGradientDiscretization="element_vertex_partial_layer";
    ImpactIonizationModelConfig impact;impact.model="van_overstraeten";
    impact.contactElectricFieldFallback=GENERATE(false,true);
    CAPTURE(impact.contactElectricFieldFallback);
    impact.generation="current_density";impact.drivingForce="quasi_fermi_gradient";
    impact.currentApproximation=GENERATE("phumob_element_density_grad_qf","phumob_element_sg_best_vertex");
    CAPTURE(impact.currentApproximation);
    impact.quasiFermiGradientDiscretization="cell_gradient";
    impact.sourceMappingMode="element_vertex_box_measure";impact.sourceJacobianMode="local_ad";
    REQUIRE_NOTHROW(detail::validateImpactIonizationDrivingForce(impact,"test"));
    auto invalid=impact;invalid.sourceJacobianMode="finite_difference";
    REQUIRE_THROWS(detail::validateImpactIonizationDrivingForce(invalid,"test"));
    CoupledDDAssembler assembler(mesh,db,doping,constants::Vt_300,cfg,recombinationModelConfig({"none"}),
        {},impact, {},{},scaling(),{},CarrierStatisticsConfig{"fermi_dirac", transportMode});
    ImpactIonizationModelConfig noImpact;
    CoupledDDAssembler transport(mesh,db,doping,constants::Vt_300,cfg,recombinationModelConfig({"none"}),
        {},noImpact, {},{},scaling(),{},CarrierStatisticsConfig{"fermi_dirac", transportMode});
    const auto mat=db.getMaterial("Si",300.);
    CoupledDDState state{VectorXd(5),VectorXd(5),VectorXd(5)};
    for(int i=0;i<5;++i) {
        state.psi(i)=.15*i;
        state.phin(i)=electronQuasiFermiPotential(mat.ni,*mat.Nc_m3,state.psi(i),2e25*(i+1),constants::Vt_300,CarrierStatisticsConfig{"fermi_dirac", transportMode});
        state.phip(i)=holeQuasiFermiPotential(mat.ni,*mat.Nv_m3,state.psi(i),1e24*(6-i),constants::Vt_300,CarrierStatisticsConfig{"fermi_dirac", transportMode});
    }
    const auto x=assembler.pack(state),y=transport.pack(state);
    REQUIRE_FALSE(assembler.impactIonizationActiveBranchFingerprint(x).empty());
    const VectorXd source=assembler.residual(x,{})-transport.residual(y,{});
    REQUIRE(source.head(5).norm()==0.);
    REQUIRE(source.segment(5,5).norm()>0.);
    // Assembled carrier equations both use particle continuity, so pair
    // generation enters the two residual blocks with the same sign.
    REQUIRE((source.segment(5,5)-source.tail(5)).norm()/source.segment(5,5).norm()<1e-12);
    REQUIRE(source(8)==0.); REQUIRE(source(13)==0.); // oxide-only vertex
    DDSolution solution;solution.psi=state.psi;solution.phin=state.phin;solution.phip=state.phip;
    solution.n=assembler.electronDensity(x);solution.p=assembler.holeDensity(x);
    const auto probe=evaluateBVProcessProbe(mesh,doping,solution,cfg,impact,{},db,300.,1.);
    auto uncachedProbeConfig=cfg;
    uncachedProbeConfig.highFieldCells.clear();
    const auto uncachedProbe=evaluateBVProcessProbe(mesh,doping,solution,uncachedProbeConfig,impact,{},db,300.,1.);
    REQUIRE(uncachedProbe.records.size()==probe.records.size());
    REQUIRE(uncachedProbe.totalSourceIntegral==Catch::Approx(probe.totalSourceIntegral).epsilon(2e-12));
    REQUIRE(probe.totalSourceIntegral==Catch::Approx(-source.segment(5,5).sum()).epsilon(2e-12));
    for(const auto& record:probe.records) {
        REQUIRE(record.supportKind=="phumob_element_vertex");
        REQUIRE(record.activeBranches.find("flux_sign=not_applicable")!=std::string::npos);
        REQUIRE(std::isfinite(record.lowFieldMobility));
        REQUIRE(record.finalMobility<=record.lowFieldMobility);
        const auto expectedField = impact.contactElectricFieldFallback && record.contactAdjacent
            ? record.electricFieldVector.norm() : (impact.currentApproximation=="phumob_element_sg_best_vertex"
                ? (record.contactAdjacent?record.quasiFermiGradientVector.norm():record.highFieldDrive) : record.quasiFermiGradientVector.norm());
        REQUIRE(record.impactField==Catch::Approx(expectedField).epsilon(2e-12));
        if(impact.currentApproximation=="phumob_element_sg_best_vertex") {
            REQUIRE(record.currentVectorProvenance=="sg_best_vertex_effective_band");
            REQUIRE(record.generationRate==Catch::Approx(record.alpha*record.currentVector.norm()).epsilon(2e-12));
        }
    }
    const auto J=assembler.assembleJacobian(x,{});
    auto sourceJ=J;sourceJ-=transport.assembleJacobian(y,{});
    for(int col=0;col<15;++col) {
        VectorXd d=VectorXd::Zero(15);d(col)=1.;
        for(Real h:{3e-6,1e-6}) {
            const VectorXd fd=(assembler.residual(x+h*d,{})-assembler.residual(x-h*d,{}))/(2*h);
            const VectorXd exact=J*d;
            const VectorXd sourceFd=(assembler.residual(x+h*d,{})-transport.residual(y+h*d,{})
                -assembler.residual(x-h*d,{})+transport.residual(y-h*d,{}))/(2*h);
            const VectorXd sourceExact=sourceJ*d;
            const Real sourceScale=std::max(sourceFd.norm(),sourceExact.norm());
            CAPTURE(impact.currentApproximation,col,h);
            if(sourceScale>1e-30)
                REQUIRE((sourceFd-sourceExact).norm()/sourceScale<5e-5);
            for(int block=0;block<3;++block) {
                const Real size=std::max(fd.segment(5*block,5).norm(),exact.segment(5*block,5).norm());
                if(size==0.)continue;
                CAPTURE(col,block,h);
                REQUIRE((fd.segment(5*block,5)-exact.segment(5*block,5)).norm()/size<3e-5);
            }
        }
    }
}

TEST_CASE("PhuMob material routing leaves PolySilicon at its material mobility with no avalanche", "[element_box][phumob_avalanche]") {
    DeviceMesh mesh;
    for(int cell=0;cell<2;++cell) {
        for(int k=0;k<3;++k) {
            Node node;node.id=mesh.numNodes();node.x=(3*cell+(k==1?1.:k==2?.5:0.))*1e-6;
            node.y=(k==2?.8660254037844386:0.)*1e-6;mesh.addNode(node);
        }
        Cell tri;tri.id=cell;tri.type=CellType::Tri3;tri.region_id=cell;
        tri.node_ids={Index(3*cell),Index(3*cell+1),Index(3*cell+2)};mesh.addCell(tri);
        Region region;region.id=cell;region.name=region.material=cell?"PolySilicon":"Si";
        region.cell_ids={Index(cell)};mesh.addRegion(region);
    }
    mesh.buildEdges();mesh.buildBoxGeometry(options());MaterialDatabase db;
    auto poly=db.getMaterial("Si",300.);poly.name="PolySilicon";poly.mun=.1417;poly.mup=.04705;db.addMaterial(poly);
    auto mats=detail::buildCellMaterials(mesh,db,300.);DopingModel doping(6);
    for(Index i=0;i<6;++i)doping.setNodeDoping(i,3e26,0.);
    auto cfg=mobilityConfig();cfg.model="phumob";cfg.edgeAveraging="element_box_phumob";
    cfg.constantMaterials={"PolySilicon"};const auto model=makeMobilityModel(cfg);
    VectorXd n=VectorXd::Constant(6,3e26),p=VectorXd::Constant(6,1e5);
    detail::EdgeMobilityCarrierState population{n(0),n(1),p(0),p(1),&n,&p};
    for(CarrierType carrier:{CarrierType::Electron,CarrierType::Hole}) {
        const Real constant=carrier==CarrierType::Electron?poly.mun:poly.mup;
        REQUIRE(detail::elementBoxCellMobility(mesh,doping,*model,mats[1],1,carrier,&cfg,&population)==constant);
        REQUIRE(detail::elementBoxCellMobility(mesh,doping,*model,mats[0],0,carrier,&cfg,&population)<constant);
    }
    std::array<Real,3> z{0.,.1,.2},density{3e26,3e26,3e26};
    const auto source=detail::phumobElementAvalancheSourceIntegralsLocal(
        ImpactIonizationModelConfig{},cfg,mesh,doping,mats,1,z,z,z,density,density);
    for(int k=0;k<3;++k) REQUIRE(source.combined[k]==0.);
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

TEST_CASE("Diagnostic transport QF linearization preserves residuals and source derivatives", "[element_box][diagnostic_qf]") {
    DeviceMesh mesh;
    for(int y=-1;y<4;++y)for(int x=0;x<4;++x){Node n;n.id=mesh.numNodes();n.x=x*1e-8;n.y=y*1e-8;mesh.addNode(n);}
    Region si;si.id=0;si.name="Si";si.material="Si";Region ox;ox.id=1;ox.name="oxide";ox.material="SiO2";
    for(int y=0;y<4;++y)for(int x=0;x<3;++x)for(int t=0;t<2;++t) {
        Index a=y*4+x,b=a+1,c=a+4,d=c+1;Cell cell;cell.id=mesh.numCells();cell.type=CellType::Tri3;cell.region_id=y==0?1:0;
        cell.node_ids=t?std::vector<Index>{a,d,c}:std::vector<Index>{a,b,d};
        (y==0?ox:si).cell_ids.push_back(cell.id);mesh.addCell(cell);
    }
    mesh.addRegion(si);mesh.addRegion(ox);mesh.buildEdges();mesh.buildBoxGeometry(options());Contact contact;contact.name="test";contact.node_ids={19};mesh.addContact(contact);
    MaterialDatabase db;DopingModel doping(20);
    for(Index i=0;i<20;++i)doping.setNodeDoping(i,2e25*(i+1),1e24);
    auto cfg=elementLombardiConfig();cfg.model="phumob_field_lombardi";
    cfg.highFieldDrivingForce="quasi_fermi_gradient";
    cfg.highFieldGradientDiscretization="element_vertex_partial_layer";
    ImpactIonizationModelConfig impact;impact.model="van_overstraeten";
    impact.contactElectricFieldFallback=GENERATE(false,true);
    CAPTURE(impact.contactElectricFieldFallback);
    impact.generation="current_density";impact.drivingForce="quasi_fermi_gradient";
    impact.currentApproximation=GENERATE(std::string("phumob_element_density_grad_qf"),
                                         std::string("phumob_element_sg_best_vertex"));
    impact.quasiFermiGradientDiscretization="cell_gradient";
    impact.sourceMappingMode="element_vertex_box_measure";impact.sourceJacobianMode="local_ad";
    REQUIRE_NOTHROW(detail::validateImpactIonizationDrivingForce(impact,"test"));
    auto laggedImpact=impact;laggedImpact.sourceJacobianMode="frozen";
    const bool sgSource=impact.currentApproximation=="phumob_element_sg_best_vertex";
    if (sgSource)
        REQUIRE_NOTHROW(detail::validateImpactIonizationDrivingForce(laggedImpact,"test"));
    else
        REQUIRE_THROWS(detail::validateImpactIonizationDrivingForce(laggedImpact,"test"));
    CoupledDDAssembler assembler(mesh,db,doping,constants::Vt_300,cfg,recombinationModelConfig({"none"}),
        {},impact, {},{},scaling(),{},CarrierStatisticsConfig{"fermi_dirac"});
    ImpactIonizationModelConfig noImpact;
    CoupledDDAssembler transport(mesh,db,doping,constants::Vt_300,cfg,recombinationModelConfig({"none"}),
        {},noImpact, {},{},scaling(),{},CarrierStatisticsConfig{"fermi_dirac"});
    const auto mat=db.getMaterial("Si",300.);
    CoupledDDState state{VectorXd(20),VectorXd(20),VectorXd(20)};
    for(int i=0;i<20;++i) {
        state.psi(i)=.15*i;
        state.phin(i)=electronQuasiFermiPotential(mat.ni,*mat.Nc_m3,state.psi(i),2e25*(i+1),constants::Vt_300,CarrierStatisticsConfig{"fermi_dirac"});
        state.phip(i)=holeQuasiFermiPotential(mat.ni,*mat.Nv_m3,state.psi(i),1e24*(21-i),constants::Vt_300,CarrierStatisticsConfig{"fermi_dirac"});
    }
    const auto x=assembler.pack(state),y=transport.pack(state);
    REQUIRE_FALSE(assembler.impactIonizationActiveBranchFingerprint(x).empty());
    const VectorXd residual=assembler.residual(x,{});
    const SparseMatrixd exact=assembler.assembleJacobian(x,{}), exactTransport=transport.assembleJacobian(y,{});
    if (sgSource) {
        CoupledDDAssembler laggedAssembler(mesh,db,doping,constants::Vt_300,cfg,
            recombinationModelConfig({"none"}),{},laggedImpact,{},{},scaling(),{},
            CarrierStatisticsConfig{"fermi_dirac"});
        REQUIRE(laggedAssembler.residual(x,{})==residual);
        const SparseMatrixd laggedJacobian=laggedAssembler.assembleJacobian(x,{});
        const SparseMatrixd expected=assembler.assembleSourceDiagnosticJacobian(x,{},2);
        REQUIRE((laggedJacobian-expected).norm()==0.);
        REQUIRE((laggedJacobian-exact).norm()>0.);
    }
    assembler.setDiagnosticFreezeTransportQfMobility(true);
    transport.setDiagnosticFreezeTransportQfMobility(true);
    REQUIRE((assembler.residual(x,{}).array()==residual.array()).all());
    const SparseMatrixd frozen=assembler.assembleJacobian(x,{}), frozenTransport=transport.assembleJacobian(y,{});
    const SparseMatrixd matched=assembler.assembleSourceDiagnosticJacobian(x,{},1);
    const SparseMatrixd lagged=assembler.assembleSourceDiagnosticJacobian(x,{},2);
    REQUIRE((assembler.residual(x,{}).array()==residual.array()).all());
    REQUIRE((lagged-frozenTransport).norm()<=1e-12*std::max(1.,frozenTransport.norm()));
    const SparseMatrixd sourceMobility=frozen-matched;
    REQUIRE(sourceMobility.norm()>0.);
    REQUIRE(sourceMobility.topRows(20).norm()==0.);
    REQUIRE(sourceMobility.leftCols(20).norm()==0.);
    REQUIRE((SparseMatrixd(sourceMobility.middleRows(20,20))-SparseMatrixd(sourceMobility.bottomRows(20))).norm()
            <=1e-12*sourceMobility.norm());
    CoupledDDBoundaryConditions dbc;dbc.psi[19]=0.;dbc.phin[19]=0.;dbc.phip[19]=0.;
    for(int mode:{0,1,2}) {
        const auto j=assembler.assembleSourceDiagnosticJacobian(x,dbc,mode);
        for(int row:{19,39,59}) {
            REQUIRE(j.coeff(row,row)==1.);
            REQUIRE(VectorXd(j.row(row).transpose()).squaredNorm()==1.);
        }
    }
    REQUIRE_THROWS(assembler.assembleSourceDiagnosticJacobian(x,{},3));
    REQUIRE_THROWS(transport.assembleSourceDiagnosticJacobian(y,{},1));
    const SparseMatrixd removed=exact-frozen, removedTransport=exactTransport-frozenTransport;
    REQUIRE(removed.norm()>0.);
    REQUIRE((removed-removedTransport).norm()/removed.norm()<1e-12);
    REQUIRE(removed.leftCols(20).norm()==0.);
    REQUIRE(removed.topRows(20).norm()==0.);
    for(int col=0;col<60;++col) {
        VectorXd v=removed.col(col);
        REQUIRE(std::abs(v.segment(20,20).sum())<=1e-12*v.norm());
        REQUIRE(std::abs(v.tail(20).sum())<=1e-12*v.norm());
    }
    assembler.setDiagnosticFreezeTransportQfMobility(false);
    REQUIRE((assembler.assembleJacobian(x,{})-exact).norm()==0.);
    auto invalidMobility=cfg;invalidMobility.jacobianFieldDerivatives=false;
    REQUIRE_THROWS(CoupledDDAssembler(mesh,db,doping,constants::Vt_300,invalidMobility,
        recombinationModelConfig({"none"}),{},impact,{}, {},scaling(),{},CarrierStatisticsConfig{"fermi_dirac"}));
}

TEST_CASE("SG source reconstruction retains equilibrium response and affine currents", "[phumob_avalanche][sg_source]") {
    const auto mesh=surfaceMesh();
    auto cell=mesh.getCell(0);
    for(bool electron:{false,true}) for(bool reverse:{false,true}) {
        if(reverse)std::swap(cell.node_ids[1],cell.node_ids[2]);
        std::array<Real,3> density{2e23,2e23,2e23},qf;
        for(int k=0;k<3;++k) {
            const auto& node=mesh.getNode(cell.node_ids[k]);
            qf[k]=2e6*node.x-3e6*node.y;
        }
        for(int best=0;best<3;++best) {
            const auto j=detail::phumobSourceSgVector(mesh,cell,best,density,qf,constants::Vt_300,electron);
            REQUIRE(j[0]==Catch::Approx(-4e29).epsilon(2e-13));
            REQUIRE(j[1]==Catch::Approx(6e29).epsilon(2e-13));
        }
        using D=detail::Tri3LocalForwardDual;
        std::array<D,3> dn{D(2e23),D(2e23),D(2e23)},df;
        for(int k=0;k<3;++k)df[k]=D::variable(0.,k);
        const auto j=detail::phumobSourceSgVector(mesh,cell,0,dn,df,constants::Vt_300,electron);
        const auto gradient=detail::localAdTri3Gradient(mesh,cell,df);
        for(int d=0;d<2;++d) {
            REQUIRE(j[d].value==0.);
            for(int k=0;k<3;++k)
                REQUIRE(j[d].derivative[k]==Catch::Approx(-2e23*gradient[d].derivative[k]).epsilon(2e-13));
        }
    }
    // Extreme QF drop must not overflow expm1; reversal preserves orientation.
    const auto f=detail::phumobSourceSgDrop(1e10,1e20,800.,constants::Vt_300);
    const auto r=detail::phumobSourceSgDrop(1e20,1e10,-800.,constants::Vt_300);
    REQUIRE(std::isfinite(f));REQUIRE(f==-r);
}

TEST_CASE("SG source applies partial layer projection to current and force", "[phumob_avalanche][sg_source]") {
    auto mesh=surfaceMesh();MaterialDatabase db;DopingModel doping(5);
    auto cfg=elementLombardiConfig();cfg.model="phumob_field_lombardi";
    cfg.highFieldDrivingForce="quasi_fermi_gradient";
    cfg.highFieldGradientDiscretization="element_vertex_partial_layer";
    const auto mats=detail::buildCellMaterials(mesh,db,300.);
    const auto ec=detail::buildEdgeCellMap(mesh);
    VectorXd potential=VectorXd::Zero(5);
    detail::updateElementHighFieldGeometry(cfg,mesh,ec,mats);
    detail::updateSurfaceMobilityCellGeometry(cfg,mesh,ec,potential,1.,&mats);
    // Supply a rank-one x-tangent partial layer. The geometry routine is tested
    // separately; here source current must use that same map as alpha and HFS.
    cfg.highFieldCells[0].contact=false;
    cfg.highFieldCells[0].stencil={1e8,0.,0.,0.};
    ImpactIonizationModelConfig impact;impact.model="van_overstraeten";
    impact.currentApproximation="phumob_element_sg_best_vertex";
    std::array<Real,3> psi{0.,0.,0.},fn{0.,.2,-.1},fp{0.,-.1,.2};
    std::array<Real,3> n{2e25,2e25,2e25},p{3e20,3e20,3e20};
    detail::PhuMobElementAvalancheDiagnostics diag;
    const auto result=detail::phumobElementAvalancheSourceIntegralsLocal(
        impact,cfg,mesh,doping,mats,0,psi,fn,fp,n,p,&diag);
    REQUIRE(diag.electronDrive==Catch::Approx(2e7));
    REQUIRE(diag.holeDrive==Catch::Approx(1e7));
    REQUIRE(diag.electronSgCurrent[1]==0.);REQUIRE(diag.holeSgCurrent[1]==0.);
    REQUIRE(diag.electronSgCurrent[0]==Catch::Approx(-2e25*2e7*diag.electronMobility[0]).epsilon(2e-13));
    REQUIRE(diag.holeSgCurrent[0]==Catch::Approx(3e20*1e7*diag.holeMobility[0]).epsilon(2e-13));
    const auto volume=detail::cellBoxNodeMeasures(mesh,0);
    for(int k=0;k<3;++k) {
        REQUIRE(result.electron[k]/volume[k]==Catch::Approx(result.electron[0]/volume[0]).epsilon(2e-13));
        REQUIRE(result.hole[k]/volume[k]==Catch::Approx(result.hole[0]/volume[0]).epsilon(2e-13));
    }
}

TEST_CASE("Effective band configuration preserves port and residual consistency", "[element_box][fermi][effective_band]")
{
    const std::string mode=GENERATE("generalized_einstein","effective_band_sg");
    const CarrierStatisticsConfig statistics{"fermi_dirac",mode};
    const nlohmann::json j={{"carrier_statistics",{{"model","fermi_dirac"},{"transport",mode}}}};
    REQUIRE(newtonConfigFromJson(j).carrierStatistics.transport==mode);
    auto bad=j;bad["carrier_statistics"]["transport"]="unknown";
    REQUIRE_THROWS(newtonConfigFromJson(bad));
    bad=j;bad["carrier_statistics"]["model"]="boltzmann";
    if(mode=="effective_band_sg") REQUIRE_THROWS(newtonConfigFromJson(bad));
    auto mesh=pairMesh();mesh.buildBoxGeometry(options());
    Contact contact;contact.name="test";contact.node_ids={0};mesh.addContact(contact);
    MaterialDatabase db;DopingModel doping(4);
    auto cfg=mobilityConfig();cfg.model="phumob";cfg.edgeAveraging="element_box_phumob";
    for(Index i=0;i<4;++i)doping.setNodeDoping(i,1e26,1e24);
    CoupledDDAssembler a(mesh,db,doping,constants::Vt_300,cfg,recombinationModelConfig({"none"}),{},{},{},{},scaling(),{},statistics);
    const auto mat=db.getMaterial("Si",300.);
    CoupledDDState s{VectorXd(4),VectorXd(4),VectorXd(4)};
    for(int i=0;i<4;++i){
        s.psi(i)=.02*i;
        s.phin(i)=electronQuasiFermiPotential(mat.ni,*mat.Nc_m3,s.psi(i),1e26*(i+1),constants::Vt_300,statistics);
        s.phip(i)=holeQuasiFermiPotential(mat.ni,*mat.Nv_m3,s.psi(i),2e25*(4-i),constants::Vt_300,statistics);
    }
    const auto x=a.pack(s);
    DDSolution sol;sol.psi=s.psi;sol.phin=s.phin;sol.phip=s.phip;
    sol.n=a.electronDensity(x);sol.p=a.holeDensity(x);
    ContactCurrent port(mesh,db,doping,cfg,300.,scaling(),{},statistics);
    const auto direct=port.compute(sol,"test"), residual=port.computeFromResidual(a,x,"test");
    REQUIRE(direct.electronCurrent==Catch::Approx(residual.electronCurrent).epsilon(2e-12));
    REQUIRE(direct.holeCurrent==Catch::Approx(residual.holeCurrent).epsilon(2e-12));
    REQUIRE(direct.electronDriftCurrent+direct.electronDiffusionCurrent==Catch::Approx(direct.electronCurrent).epsilon(2e-12));
    REQUIRE(direct.holeDriftCurrent+direct.holeDiffusionCurrent==Catch::Approx(direct.holeCurrent).epsilon(2e-12));
    if(mode=="effective_band_sg"){
        DensityGradientQuantumPotentialConfig quantum;quantum.enabled=true;
        REQUIRE_THROWS(ContactCurrent(mesh,db,doping,cfg,300.,scaling(),{},statistics,quantum));
        REQUIRE_THROWS(CoupledDDAssembler(mesh,db,doping,constants::Vt_300,cfg,recombinationModelConfig({"none"}),{},{},{},{},scaling(),{},statistics,quantum));
    }
}
