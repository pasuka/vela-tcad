#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>
#include "vela/equation/PoissonAssembler.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/simulation/ConfigParsing.h"
#include <limits>

using namespace vela;
namespace {
using Policy = BoxGeometryBuilder::PoissonPermittivityPolicy;
DeviceMesh patch(bool homogeneous = false, bool reverse = false, double scale = 1e-6, double height = 1.)
{
    DeviceMesh mesh;
    const std::array<std::array<double,2>,4> xy = {{{0,0},{1,0},{.5,height},{.5,-.6}}};
    for (Index i=0;i<4;++i) { Node n; n.id=i; n.x=xy[i][0]*scale; n.y=xy[i][1]*scale; mesh.addNode(n); }
    for (Index i=0;i<2;++i) {
        Cell c; c.id=i; c.type=CellType::Tri3; c.region_id=i;
        const Index side=reverse?1-i:i;
        c.node_ids=side==0?std::vector<Index>{0,1,2}:std::vector<Index>{1,0,3};
        if(reverse) std::reverse(c.node_ids.begin(),c.node_ids.end());
        mesh.addCell(c);
        Region r; r.id=i; r.name="r"+std::to_string(i); r.cell_ids={i};
        r.material=(homogeneous || side==0)?"Si":"SiO2"; mesh.addRegion(r);
    }
    mesh.buildEdges(); return mesh;
}
Eigen::MatrixXd matrix(const DeviceMesh& mesh)
{
    MaterialDatabase db; DopingModel doping(mesh.numNodes());
    PoissonAssembler a(mesh,db,doping); a.assemble(); return Eigen::MatrixXd(a.matrix());
}
BoxGeometryBuilder::Options cellOptions()
{
    BoxGeometryBuilder::Options o; o.poissonPermittivityPolicy=Policy::CellMaterial; return o;
}
}

TEST_CASE("Permittivity uses material-weighted unequal interface box measures", "[permittivity]")
{
    auto mesh=patch(); const auto old=matrix(mesh); mesh.buildBoxGeometry(cellOptions());
    MaterialDatabase db;
    const double si=db.getMaterial("Si").eps_r*constants::eps0;
    const double ox=db.getMaterial("SiO2").eps_r*constants::eps0;
    // Independent analytic cotangents for two acute triangles on a horizontal edge.
    const double gSi=(1.-.25)/2., gOx=(.6*.6-.25)/(2.*.6);
    const auto K=matrix(mesh);
    REQUIRE(-K(0,1)==Catch::Approx(si*gSi+ox*gOx).epsilon(2e-14));
    REQUIRE(-old(0,1)==Catch::Approx(.5*(si+ox)*(gSi+gOx)).epsilon(2e-14));
    REQUIRE(std::abs(K(0,1)/old(0,1)-1.)>.1);
    REQUIRE((K-K.transpose()).norm()==0.);
    REQUIRE((K*VectorXd::Ones(4)).norm()<1e-24);
    VectorXd v(4);v<<.1,-.3,.7,.2;
    REQUIRE(v.dot(K*v)>0.);
    auto reordered=patch(false,true);reordered.buildBoxGeometry(cellOptions());
    REQUIRE((matrix(reordered)-K).norm()<1e-24);
}

TEST_CASE("Cell permittivity preserves homogeneous geometry and 2D coordinate scaling", "[permittivity]")
{
    for(bool fallback:{false,true}) for(double height:{1.,.2}) {
        auto mesh=patch(true,false,1e-6,height);auto o=cellOptions();o.fallbackNegativeCotangent=fallback;
        auto legacy=o;legacy.poissonPermittivityPolicy=Policy::LegacyAverage;
        mesh.buildBoxGeometry(legacy);auto old=matrix(mesh);mesh.buildBoxGeometry(o);
        REQUIRE((matrix(mesh)-old).norm()/old.norm()<3e-16);
    }
    auto a=patch(false,false,1e-6),b=patch(false,false,1.);
    a.buildBoxGeometry(cellOptions());b.buildBoxGeometry(cellOptions());
    REQUIRE((matrix(a)-matrix(b)).norm()/matrix(a).norm()<1e-14);
}

TEST_CASE("Layered dielectric patch preserves normal displacement and series capacitance", "[permittivity]")
{
    DeviceMesh mesh;
    for(Index i=0;i<6;++i) { Node n;n.id=i;n.x=(i/2==2?3.:double(i/2))*1e-6;n.y=(i%2)*1e-6;mesh.addNode(n); }
    const std::array<std::array<Index,3>,4> triangles={{{0,2,3},{0,3,1},{2,4,5},{2,5,3}}};
    for(Index i=0;i<4;++i) { Cell c;c.id=i;c.type=CellType::Tri3;c.region_id=i/2;c.node_ids.assign(triangles[i].begin(),triangles[i].end());mesh.addCell(c); }
    for(Index i=0;i<2;++i) { Region r;r.id=i;r.name="r"+std::to_string(i);r.material=i==0?"Si":"SiO2";r.cell_ids={2*i,2*i+1};mesh.addRegion(r); }
    mesh.buildEdges();mesh.buildBoxGeometry(cellOptions());
    MaterialDatabase db;DopingModel doping(6);PoissonAssembler assembler(mesh,db,doping);
    assembler.assemble();const Eigen::MatrixXd K=assembler.matrix();assembler.applyDirichlet({{0,0.},{1,0.},{4,1.},{5,1.}});
    const VectorXd psi=Eigen::MatrixXd(assembler.matrix()).fullPivLu().solve(assembler.rhs());
    const double si=db.getMaterial("Si").eps_r*constants::eps0,ox=db.getMaterial("SiO2").eps_r*constants::eps0;
    const double expected=ox/(2*si+ox);REQUIRE(psi[2]==Catch::Approx(expected).epsilon(1e-13));REQUIRE(psi[3]==Catch::Approx(expected).epsilon(1e-13));
    REQUIRE(si*psi[2]/1e-6==Catch::Approx(ox*(1-psi[2])/2e-6).epsilon(1e-13));
    const VectorXd reaction=K*psi;const double capacitancePerDepth=1e-6/(1e-6/si+2e-6/ox);
    REQUIRE(reaction[4]+reaction[5]==Catch::Approx(capacitancePerDepth).epsilon(1e-13));
    REQUIRE(std::abs(reaction.sum())<1e-24);
}

TEST_CASE("Poisson Gummel Newton and electrode reaction share cell permittivity", "[permittivity]")
{
    for(bool supplied:{false,true}) {
        auto mesh=patch();auto o=cellOptions();
        if(supplied) o.poissonCellEdgeCoefficients={{{0,1,2},{.2,.3,.4}},{{1,0,3},{.5,.6,.7}}};
        mesh.buildBoxGeometry(o);const auto K=matrix(mesh);
        MaterialDatabase db;DopingModel doping(4);
        DDAssembler dd(mesh,db,doping,.02585,1e-6,1e-6);
        dd.assemblePoissonWithCarriers(VectorXd::Zero(4),VectorXd::Zero(4),VectorXd::Zero(4));
        REQUIRE((Eigen::MatrixXd(dd.matrix())-K).norm()<1e-24);
        CoupledDDAssembler coupled(mesh,db,doping,.02585,1e-6,1e-6);
        CoupledDDState s{VectorXd::Zero(4),VectorXd::Zero(4),VectorXd::Zero(4)};
        const auto x=coupled.pack(s);VectorXd direction=VectorXd::Zero(12),v(4);v<<.11,-.2,.3,-.07;
        for(int k=0;k<3;++k)direction.segment(4*k,4)=v;
        // Moving psi and both quasi-Fermi potentials together fixes n and p.
        const VectorXd Jv=coupled.assembleJacobian(x,{})*direction;
        REQUIRE((Jv.head(4)-K*v).norm()<1e-23);
        for(double h:{1e-3,5e-4,2.5e-4}) {
            const VectorXd delta=(coupled.residual(x+h*direction,{})-coupled.residual(x-h*direction,{}))/(2*h);
            REQUIRE((delta.head(4)-K*v).norm()/(K*v).norm()<1e-10);
        }
        const VectorXd reaction=coupled.poissonDirichletReactionChargePerMeter(x+direction)-
            coupled.poissonDirichletReactionChargePerMeter(x);
        REQUIRE((reaction-K*v).norm()/(K*v).norm()<1e-12);
    }
}

TEST_CASE("Explicit cell geometry changes only dielectric flux and resets on rebuild", "[permittivity]")
{
    auto mesh=patch();const auto old=matrix(mesh);const auto edges=detail::computeEdgeCouplings(mesh);
    const auto volumes=detail::computeNodeVolumes(mesh);auto o=cellOptions();
    o.poissonCellEdgeCoefficients={{{0,1,2},{.2,.3,.4}},{{1,0,3},{.5,.6,.7}}};
    mesh.buildBoxGeometry(o);MaterialDatabase db;
    const double expected=constants::eps0*(.2*db.getMaterial("Si").eps_r+.5*db.getMaterial("SiO2").eps_r);
    REQUIRE(-matrix(mesh)(0,1)==Catch::Approx(expected).epsilon(1e-14));
    REQUIRE(detail::computeEdgeCouplings(mesh)==edges);
    REQUIRE(detail::computeNodeVolumes(mesh)==volumes);
    REQUIRE_THROWS_AS(detail::poissonEdgeCoefficient(mesh,db,detail::buildEdgeCellMap(mesh),0,true),std::invalid_argument);
    mesh.buildBoxGeometry();REQUIRE((matrix(mesh)-old).norm()==0.);
}

TEST_CASE("Cell coefficient input rejects topology size value and policy mismatches", "[permittivity]")
{
    auto mesh=patch();auto o=cellOptions();
    const std::vector<BoxGeometryBuilder::PoissonCellEdgeCoefficients> valid={{{0,1,2},{.2,.3,.4}},{{1,0,3},{.5,.6,.7}}};
    o.poissonCellEdgeCoefficients=valid;o.poissonCellEdgeCoefficients.pop_back();REQUIRE_THROWS(mesh.buildBoxGeometry(o));
    o.poissonCellEdgeCoefficients=valid;o.poissonCellEdgeCoefficients[0].nodeIds[0]=3;REQUIRE_THROWS(mesh.buildBoxGeometry(o));
    for(double bad:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        o.poissonCellEdgeCoefficients=valid;o.poissonCellEdgeCoefficients[1].coefficients[0]=bad;REQUIRE_THROWS(mesh.buildBoxGeometry(o));
    }
    o.poissonCellEdgeCoefficients=valid;o.poissonPermittivityPolicy=Policy::LegacyAverage;REQUIRE_THROWS(mesh.buildBoxGeometry(o));
    REQUIRE(parseBoxGeometryOptions(nlohmann::json::object()).poissonPermittivityPolicy==Policy::LegacyAverage);
    auto deck=nlohmann::json::parse(R"({"mesh_geometry":{"poisson_permittivity_policy":"cell_material","poisson_cell_edge_coefficients":[{"cell_id":0,"node_ids":[0,1,2],"coefficients":[0.2,0.3,0.4]},{"cell_id":1,"node_ids":[1,0,3],"coefficients":[0.5,0.6,0.7]}]}})");
    REQUIRE_NOTHROW(mesh.buildBoxGeometry(parseBoxGeometryOptions(deck)));
    deck["mesh_geometry"]["poisson_cell_edge_coefficients"][0]["cell_id"]=1;REQUIRE_THROWS(parseBoxGeometryOptions(deck));
    deck["mesh_geometry"]["poisson_cell_edge_coefficients"]=nlohmann::json::array();REQUIRE_THROWS(parseBoxGeometryOptions(deck));
    deck["mesh_geometry"]={{"poisson_permittivity_policy","unknown"}};REQUIRE_THROWS(parseBoxGeometryOptions(deck));
}
