#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>
#include "vela/mesh/DelaunayBox.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/simulation/ConfigParsing.h"
#include "vela/solver/NewtonSolver.h"

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
