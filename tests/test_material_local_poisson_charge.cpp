#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/DDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/simulation/ConfigParsing.h"

#include <limits>
#include <nlohmann/json.hpp>

using namespace vela;

namespace {

DeviceMesh makeHorizontalInterfaceMesh()
{
    DeviceMesh mesh;
    const Real length = 1.0e-6;

    Node n0; n0.id = 0; n0.x = 0.0;          n0.y = 0.0;     mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = length;       n1.y = 0.0;     mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.5 * length; n2.y = -length; mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.5 * length; n3.y = length;  mesh.addNode(n3);

    Cell siliconCell;
    siliconCell.id = 0;
    siliconCell.type = CellType::Tri3;
    siliconCell.region_id = 0;
    siliconCell.node_ids = {0, 1, 2};
    mesh.addCell(siliconCell);

    Cell oxideCell;
    oxideCell.id = 1;
    oxideCell.type = CellType::Tri3;
    oxideCell.region_id = 1;
    oxideCell.node_ids = {0, 3, 1};
    mesh.addCell(oxideCell);

    Region silicon;
    silicon.id = 0;
    silicon.name = "channel";
    silicon.material = "Si";
    silicon.cell_ids = {0};
    mesh.addRegion(silicon);

    Region oxide;
    oxide.id = 1;
    oxide.name = "gate_oxide";
    oxide.material = "SiO2";
    oxide.cell_ids = {1};
    mesh.addRegion(oxide);

    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeSiliconTriangleMesh()
{
    DeviceMesh mesh;
    Node n0; n0.id = 0; n0.x = 0.0;    n0.y = 0.0;    mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0e-6; n1.y = 0.0;    mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.0;    n2.y = 1.0e-6; mesh.addNode(n2);
    Cell cell;
    cell.id = 0;
    cell.type = CellType::Tri3;
    cell.region_id = 0;
    cell.node_ids = {0, 1, 2};
    mesh.addCell(cell);
    Region silicon;
    silicon.id = 0;
    silicon.name = "silicon";
    silicon.material = "Si";
    silicon.cell_ids = {0};
    mesh.addRegion(silicon);
    mesh.buildEdges();
    return mesh;
}

} // namespace

TEST_CASE("material-local Poisson charge volume excludes oxide cell area",
          "[poisson][charge_volume]")
{
    const DeviceMesh mesh = makeHorizontalInterfaceMesh();
    const MaterialDatabase matdb;
    const auto materials = detail::buildCellMaterials(mesh, matdb, constants::T0);
    const auto global = detail::computePoissonChargeVolumes(
        mesh, materials, PoissonChargeVolumePolicy::Global);
    const auto local = detail::computePoissonChargeVolumes(
        mesh, materials, PoissonChargeVolumePolicy::MaterialLocalBarycentric);

    REQUIRE(local.at(0) == Catch::Approx(0.5 * global.at(0)).epsilon(1.0e-14));
    REQUIRE(local.at(1) == Catch::Approx(0.5 * global.at(1)).epsilon(1.0e-14));
    REQUIRE(local.at(2) == global.at(2));
    REQUIRE(local.at(3) == 0.0);
}

TEST_CASE("material-local Poisson charge volume preserves silicon-only mesh",
          "[poisson][charge_volume]")
{
    const DeviceMesh mesh = makeSiliconTriangleMesh();
    const MaterialDatabase matdb;
    const auto materials = detail::buildCellMaterials(mesh, matdb, constants::T0);
    const auto global = detail::computePoissonChargeVolumes(
        mesh, materials, PoissonChargeVolumePolicy::Global);
    const auto local = detail::computePoissonChargeVolumes(
        mesh, materials, PoissonChargeVolumePolicy::MaterialLocalBarycentric);
    REQUIRE(local == global);
}

TEST_CASE("material-local Poisson charge volume changes residual and Jacobian together",
          "[poisson][charge_volume][jacobian]")
{
    const DeviceMesh mesh = makeHorizontalInterfaceMesh();
    const MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 1.0e21, 0.0}, {"gate_oxide", 0.0, 0.0}});

    auto build = [&](PoissonChargeVolumePolicy policy) {
        DDScalingSpec scaling;
        scaling.poissonChargeVolumePolicy = policy;
        return CoupledDDAssembler(
            mesh, matdb, doping, constants::Vt_300,
            MobilityModelConfig{},
            recombinationModelConfig({"srh"}, 1.0e-6, 1.0e-6),
            BandgapNarrowingConfig{}, ImpactIonizationModelConfig{},
            {}, {}, scaling);
    };

    CoupledDDAssembler baseline = build(PoissonChargeVolumePolicy::Global);
    CoupledDDAssembler candidate = build(
        PoissonChargeVolumePolicy::MaterialLocalBarycentric);
    CoupledDDState state{
        VectorXd::Zero(4), VectorXd::Constant(4, 0.01),
        VectorXd::Constant(4, -0.02)};
    const VectorXd baselineState = baseline.pack(state);
    const VectorXd candidateState = candidate.pack(state);
    const CoupledDDBoundaryConditions bcs;

    const VectorXd baselineResidual = baseline.residual(baselineState, bcs);
    const VectorXd candidateResidual = candidate.residual(candidateState, bcs);
    REQUIRE((candidateResidual.segment(4, 8) -
             baselineResidual.segment(4, 8)).norm() == 0.0);
    const auto baselineTerms = baseline.poissonTermDiagnostics(baselineState, bcs);
    const auto candidateTerms = candidate.poissonTermDiagnostics(candidateState, bcs);
    REQUIRE(candidateTerms.at(0).dopingCharge ==
            Catch::Approx(0.5 * baselineTerms.at(0).dopingCharge).epsilon(1.0e-12));
    REQUIRE(candidateTerms.at(0).electronCharge ==
            Catch::Approx(0.5 * baselineTerms.at(0).electronCharge).epsilon(1.0e-12));

    const SparseMatrixd baselineJacobian = baseline.assembleJacobian(baselineState, bcs);
    const SparseMatrixd candidateJacobian = candidate.assembleJacobian(candidateState, bcs);
    REQUIRE(candidateJacobian.coeff(0, 4) ==
            Catch::Approx(0.5 * baselineJacobian.coeff(0, 4)).epsilon(1.0e-12));
    REQUIRE(candidateJacobian.coeff(0, 8) ==
            Catch::Approx(0.5 * baselineJacobian.coeff(0, 8)).epsilon(1.0e-12));
}

TEST_CASE("Poisson charge volume policy parsing is explicit and barycentric-only",
          "[poisson][charge_volume][config]")
{
    REQUIRE(parsePoissonChargeVolumePolicy(nlohmann::json::object()) ==
            PoissonChargeVolumePolicy::Global);
    REQUIRE(parsePoissonChargeVolumePolicy(nlohmann::json{
                {"discretization", {{"poisson_charge_volume_policy", "material_local"}}}}) ==
            PoissonChargeVolumePolicy::MaterialLocalBarycentric);
    REQUIRE_THROWS_WITH(
        parsePoissonChargeVolumePolicy(nlohmann::json{
            {"mesh_geometry", {{"node_volume_policy", "mixed_voronoi"}}},
            {"discretization", {{"poisson_charge_volume_policy", "material_local"}}}}),
        Catch::Matchers::ContainsSubstring("qualified only"));
}

TEST_CASE("Mixed material-local Poisson preserves volume and continuity with full charge derivatives",
          "[poisson][charge_volume][jacobian]")
{
    auto mesh=makeHorizontalInterfaceMesh();
    BoxGeometryBuilder::Options geometry;
    geometry.nodeVolumePolicy=BoxGeometryBuilder::NodeVolumePolicy::MixedVoronoi;
    mesh.buildBoxGeometry(geometry);
    MaterialDatabase materials;
    const auto cellMaterials=detail::buildCellMaterials(mesh,materials,300.);
    const auto local=detail::computePoissonChargeVolumes(mesh,cellMaterials,PoissonChargeVolumePolicy::MaterialLocalMixedVoronoi);
    const auto global=detail::computeNodeVolumes(mesh);
    CHECK(local[0]==Catch::Approx(.5*global[0]).epsilon(1e-14));
    CHECK(local[1]==Catch::Approx(.5*global[1]).epsilon(1e-14));
    CHECK(local[2]==global[2]);CHECK(local[3]==0.);
    CHECK(local[0]+local[1]+local[2]==Catch::Approx(.5e-12).epsilon(1e-14));
    const auto doping=DopingModel::fromMeshAndRegions(mesh,{{"channel",1e21,0.}});
    DDScalingSpec scaling;scaling.poissonChargeVolumePolicy=PoissonChargeVolumePolicy::MaterialLocalMixedVoronoi;
    CoupledDDAssembler candidate(mesh,materials,doping,constants::Vt_300,{}, {}, {}, {}, {}, {},scaling);
    scaling.poissonChargeVolumePolicy=PoissonChargeVolumePolicy::Global;
    CoupledDDAssembler baseline(mesh,materials,doping,constants::Vt_300,{}, {}, {}, {}, {}, {},scaling);
    CoupledDDState state{VectorXd::Constant(4,.03),VectorXd::Constant(4,.01),VectorXd::Constant(4,-.02)};
    const auto x=candidate.pack(state);CoupledDDBoundaryConditions bc;
    CHECK(candidate.residual(x,bc).tail(8)==baseline.residual(x,bc).tail(8));
    const auto jac=candidate.assembleJacobian(x,bc);
    for(int col=0;col<12;++col){
        auto a=x,b=x;a[col]+=1e-6;b[col]-=1e-6;
        const VectorXd fd=(candidate.residual(a,bc)-candidate.residual(b,bc))/(2e-6);
        const VectorXd analytic=VectorXd(jac.col(col));
        CHECK((analytic-fd).norm()/std::max({analytic.norm(),fd.norm(),1e-30})<1e-6);
    }
    const nlohmann::json cfg{{"mesh_geometry",{{"node_volume_policy","mixed_voronoi"}}},
        {"discretization",{{"poisson_charge_volume_policy","material_local_mixed_voronoi"}}}};
    CHECK(parsePoissonChargeVolumePolicy(cfg)==PoissonChargeVolumePolicy::MaterialLocalMixedVoronoi);
    auto bad=cfg;bad["mesh_geometry"]["node_volume_policy"]="barycentric";
    CHECK_THROWS(parsePoissonChargeVolumePolicy(bad));
    auto single=makeSiliconTriangleMesh();single.buildBoxGeometry(geometry);
    CHECK(detail::computePoissonChargeVolumes(single,detail::buildCellMaterials(single,materials,300.),
        PoissonChargeVolumePolicy::MaterialLocalMixedVoronoi)==detail::computeNodeVolumes(single));
}

TEST_CASE("Poisson carrier assembly reports a non-finite carrier state",
          "[poisson][diagnostic]")
{
    const DeviceMesh mesh = makeHorizontalInterfaceMesh();
    const MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 1.0e21, 0.0}, {"gate_oxide", 0.0, 0.0}});
    const VectorXd psi = VectorXd::Zero(4);
    const VectorXd electrons = VectorXd::Constant(4, 1.0e16);
    VectorXd holes = VectorXd::Constant(4, 1.0e16);
    holes(0) = std::numeric_limits<Real>::infinity();

    RecombinationModelConfig recombination;
    recombination.mechanisms = {"none"};
    DDAssembler assembler(mesh, matdb, doping, constants::Vt_300,
                          MobilityModelConfig{}, recombination,
                          BandgapNarrowingConfig{},
                          ImpactIonizationModelConfig{}, {}, {});

    REQUIRE_THROWS_WITH(
        assembler.assemblePoissonWithCarriers(electrons, holes, psi),
        Catch::Matchers::ContainsSubstring("non-finite carrier diagonal at node 0"));
}
