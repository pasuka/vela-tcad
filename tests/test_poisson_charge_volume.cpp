#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/MobilityModel.h"

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

    Cell silicon;
    silicon.id = 0;
    silicon.type = CellType::Tri3;
    silicon.region_id = 0;
    silicon.node_ids = {0, 1, 2};
    mesh.addCell(silicon);

    Cell oxide;
    oxide.id = 1;
    oxide.type = CellType::Tri3;
    oxide.region_id = 1;
    oxide.node_ids = {0, 3, 1};
    mesh.addCell(oxide);

    Region channel;
    channel.id = 0;
    channel.name = "channel";
    channel.material = "Si";
    channel.cell_ids = {0};
    mesh.addRegion(channel);

    Region gateOxide;
    gateOxide.id = 1;
    gateOxide.name = "gate_oxide";
    gateOxide.material = "SiO2";
    gateOxide.cell_ids = {1};
    mesh.addRegion(gateOxide);

    mesh.buildEdges();
    return mesh;
}

} // namespace

TEST_CASE("electron-only Poisson volume split preserves every other charge term",
          "[poisson][charge_volume][interface_geometry]")
{
    const DeviceMesh mesh = makeHorizontalInterfaceMesh();
    const MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 1.0e21, 0.0}, {"gate_oxide", 0.0, 0.0}});
    const auto cellMaterials = detail::buildCellMaterials(
        mesh, matdb, constants::T0);
    const auto allVolumes = detail::computeNodeVolumes(mesh);
    const auto transportVolumes = detail::computeTransportNodeVolumes(
        mesh, cellMaterials);

    MobilityModelConfig mobility;
    mobility.model = "constant";
    RecombinationModelConfig recombination;
    recombination.mechanisms = {"none"};
    DDScalingSpec baselineScaling;
    DDScalingSpec splitScaling;
    splitScaling.regionResolvedInterfaceAssembly
        .poissonElectronTransportNodeVolume = true;
    CoupledDDAssembler baseline(
        mesh, matdb, doping, constants::Vt_300, mobility, recombination,
        {}, {}, {}, {}, baselineScaling);
    CoupledDDAssembler split(
        mesh, matdb, doping, constants::Vt_300, mobility, recombination,
        {}, {}, {}, {}, splitScaling);

    CoupledDDState state;
    state.psi = VectorXd::Zero(4);
    state.phin = VectorXd::Zero(4);
    state.phip = VectorXd::Zero(4);
    const VectorXd baselineX = baseline.pack(state);
    const VectorXd splitX = split.pack(state);
    const VectorXd baselineReaction =
        baseline.poissonDirichletReactionChargePerMeter(baselineX);
    const VectorXd splitReaction =
        split.poissonDirichletReactionChargePerMeter(splitX);
    const Real electronDensity = matdb.getMaterial("Si").ni;
    const Real expectedReactionChange = constants::q * electronDensity *
        (transportVolumes.at(0) - allVolumes.at(0));
    REQUIRE(splitReaction(0) - baselineReaction(0) ==
            Catch::Approx(expectedReactionChange));
    REQUIRE(splitReaction(2) - baselineReaction(2) == Catch::Approx(0.0));
    REQUIRE(splitReaction(3) - baselineReaction(3) == Catch::Approx(0.0));

    const CoupledDDBoundaryConditions bcs;
    const SparseMatrixd baselineJacobian =
        baseline.assembleJacobian(baselineX, bcs);
    const SparseMatrixd splitJacobian = split.assembleJacobian(splitX, bcs);
    const Real expectedJacobianChange = expectedReactionChange /
        constants::Vt_300;
    REQUIRE(splitJacobian.coeff(0, 0) - baselineJacobian.coeff(0, 0) ==
            Catch::Approx(expectedJacobianChange));
}

TEST_CASE("hole and dopant Poisson volume splits are independent and additive",
          "[poisson][charge_volume][interface_geometry]")
{
    const DeviceMesh mesh = makeHorizontalInterfaceMesh();
    const MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 1.0e21, 0.0}, {"gate_oxide", 0.0, 0.0}});
    const auto cellMaterials = detail::buildCellMaterials(
        mesh, matdb, constants::T0);
    const auto allVolumes = detail::computeNodeVolumes(mesh);
    const auto transportVolumes = detail::computeTransportNodeVolumes(
        mesh, cellMaterials);

    MobilityModelConfig mobility;
    mobility.model = "constant";
    RecombinationModelConfig recombination;
    recombination.mechanisms = {"none"};
    DDScalingSpec baselineScaling;
    DDScalingSpec holeScaling;
    holeScaling.regionResolvedInterfaceAssembly
        .poissonHoleTransportNodeVolume = true;
    DDScalingSpec dopantScaling;
    dopantScaling.regionResolvedInterfaceAssembly
        .poissonDopantTransportNodeVolume = true;
    DDScalingSpec combinedScaling;
    combinedScaling.regionResolvedInterfaceAssembly
        .poissonElectronTransportNodeVolume = true;
    combinedScaling.regionResolvedInterfaceAssembly
        .poissonHoleTransportNodeVolume = true;
    combinedScaling.regionResolvedInterfaceAssembly
        .poissonDopantTransportNodeVolume = true;

    auto makeAssembler = [&](const DDScalingSpec& scaling) {
        return CoupledDDAssembler(
            mesh, matdb, doping, constants::Vt_300, mobility, recombination,
            {}, {}, {}, {}, scaling);
    };
    auto baseline = makeAssembler(baselineScaling);
    auto hole = makeAssembler(holeScaling);
    auto dopant = makeAssembler(dopantScaling);
    auto combined = makeAssembler(combinedScaling);

    CoupledDDState state;
    state.psi = VectorXd::Zero(4);
    state.phin = VectorXd::Zero(4);
    state.phip = VectorXd::Zero(4);
    const VectorXd baselineX = baseline.pack(state);
    const Real volumeDelta = transportVolumes.at(0) - allVolumes.at(0);
    const Real intrinsicDensity = matdb.getMaterial("Si").ni;
    const Real holeExpected = -constants::q * intrinsicDensity * volumeDelta;
    const Real dopantExpected = -constants::q *
        doping.netDoping(0) * volumeDelta;
    const auto reaction = [&](CoupledDDAssembler& assembler) {
        return assembler.poissonDirichletReactionChargePerMeter(
            assembler.pack(state));
    };
    const VectorXd baselineReaction = reaction(baseline);
    REQUIRE(reaction(hole)(0) - baselineReaction(0) ==
            Catch::Approx(holeExpected));
    REQUIRE(reaction(dopant)(0) - baselineReaction(0) ==
            Catch::Approx(dopantExpected));
    REQUIRE(reaction(combined)(0) - baselineReaction(0) ==
            Catch::Approx(dopantExpected));

    const CoupledDDBoundaryConditions bcs;
    const SparseMatrixd baselineJacobian =
        baseline.assembleJacobian(baselineX, bcs);
    const SparseMatrixd holeJacobian =
        hole.assembleJacobian(hole.pack(state), bcs);
    REQUIRE(holeJacobian.coeff(0, 0) - baselineJacobian.coeff(0, 0) ==
            Catch::Approx(-holeExpected / constants::Vt_300));
    REQUIRE(holeJacobian.coeff(0, 8) - baselineJacobian.coeff(0, 8) ==
            Catch::Approx(holeExpected / constants::Vt_300));
}
