#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DensityGradientQuantumPotential.h"
#include "vela/solver/NewtonSolver.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

using namespace vela;

namespace {

DeviceMesh stripMesh()
{
    DeviceMesh mesh;
    mesh.addNode({0, 0.0, 0.0});
    mesh.addNode({1, 1.0e-9, 0.0});
    mesh.addNode({2, 2.0e-9, 0.0});
    mesh.addNode({3, 0.0, 1.0e-9});
    mesh.addNode({4, 1.0e-9, 1.0e-9});
    mesh.addNode({5, 2.0e-9, 1.0e-9});
    mesh.addRegion({0, "silicon", "Silicon", {0, 1, 2, 3}});
    mesh.addCell({0, CellType::Tri3, 0, {0, 1, 4}});
    mesh.addCell({1, CellType::Tri3, 0, {0, 4, 3}});
    mesh.addCell({2, CellType::Tri3, 0, {1, 2, 5}});
    mesh.addCell({3, CellType::Tri3, 0, {1, 5, 4}});
    mesh.buildEdges();
    return mesh;
}

DeviceMesh twoMaterialStripMesh()
{
    DeviceMesh mesh;
    mesh.addNode({0, 0.0, 0.0});
    mesh.addNode({1, 1.0e-9, 0.0});
    mesh.addNode({2, 2.0e-9, 0.0});
    mesh.addNode({3, 0.0, 1.0e-9});
    mesh.addNode({4, 1.0e-9, 1.0e-9});
    mesh.addNode({5, 2.0e-9, 1.0e-9});
    mesh.addRegion({0, "silicon", "Silicon", {0, 1}});
    mesh.addRegion({1, "oxide", "SiO2", {2, 3}});
    mesh.addCell({0, CellType::Tri3, 0, {0, 1, 4}});
    mesh.addCell({1, CellType::Tri3, 0, {0, 4, 3}});
    mesh.addCell({2, CellType::Tri3, 1, {1, 2, 5}});
    mesh.addCell({3, CellType::Tri3, 1, {1, 5, 4}});
    mesh.buildEdges();
    return mesh;
}

DeviceMesh twoInsulatorReentrantCornerMesh()
{
    DeviceMesh mesh;
    mesh.addNode({0, 0.0, 0.0});
    mesh.addNode({1, 1.0e-9, 0.0});
    mesh.addNode({2, 0.5e-9, 0.8660254037844386e-9});
    mesh.addNode({3, -0.5e-9, 0.8660254037844386e-9});
    mesh.addNode({4, -1.0e-9, 0.0});
    mesh.addNode({5, -0.5e-9, -0.8660254037844386e-9});
    mesh.addNode({6, 0.5e-9, -0.8660254037844386e-9});
    mesh.addRegion({0, "oxide", "SiO2", {0, 1}});
    mesh.addRegion({1, "nitride", "Nitride", {2, 3, 4, 5}});
    mesh.addCell({0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.addCell({1, CellType::Tri3, 0, {0, 2, 3}});
    mesh.addCell({2, CellType::Tri3, 1, {0, 3, 4}});
    mesh.addCell({3, CellType::Tri3, 1, {0, 4, 5}});
    mesh.addCell({4, CellType::Tri3, 1, {0, 5, 6}});
    mesh.addCell({5, CellType::Tri3, 1, {0, 6, 1}});
    mesh.buildEdges();
    return mesh;
}

DeviceMesh obtuseTriangleMesh()
{
    DeviceMesh mesh;
    mesh.addNode({0, 0.0, 0.0});
    mesh.addNode({1, 2.0e-9, 0.0});
    mesh.addNode({2, 0.5e-9, 0.2e-9});
    mesh.addRegion({0, "silicon", "Silicon", {0}});
    mesh.addCell({0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.buildEdges();
    return mesh;
}

} // namespace

TEST_CASE("Density-gradient coefficient follows Ancona-Tiersten form",
          "[density_gradient][quantum]")
{
    DensityGradientQuantumPotentialConfig config;
    config.gamma = 3.6;
    config.effectiveMassRatio = 1.0618016171622988;
    const Real hbar = constants::h / (2.0 * std::acos(-1.0));
    const Real expected = config.gamma * hbar * hbar /
        (6.0 * config.effectiveMassRatio * constants::m0 * constants::q);
    REQUIRE(densityGradientCoefficientVm2(config) ==
            Catch::Approx(expected).epsilon(1.0e-14));
}

TEST_CASE("DG population Poisson preserves a sub-ULP electrostatic correction",
          "[density_gradient][population][compensation]")
{
    const auto mesh = stripMesh();
    MaterialDatabase materials;
    Material silicon = materials.getMaterial("Si");
    silicon.name = "Silicon";
    materials.addMaterial(silicon);
    const DopingModel doping(mesh.numNodes());
    const int nodes = static_cast<int>(mesh.numNodes());
    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = 0.03125;
    scaling.stateVariables = "potential_density";
    const CoupledDDAssembler assembler(
        mesh, materials, doping, constants::Vt_300, {}, {}, {}, {}, {}, {}, scaling);

    VectorXd state = VectorXd::Constant(3 * nodes, 4.);
    VectorXd low = VectorXd::Zero(nodes);
    low[0] = std::ldexp(1., -54);
    REQUIRE(state[0] + low[0] == state[0]);
    const auto residual = assembler.compensatedPopulationPoissonResidual(state, low, {});
    REQUIRE(residual[0] > 0.);

    // A flat high part cancels and the independently resolvable perturbation
    // must give the same linear response after rescaling.
    VectorXd resolved = state;
    resolved[0] += 0.25;
    const auto reference = assembler.compensatedPopulationPoissonResidual(
        resolved, VectorXd::Zero(nodes), {});
    REQUIRE((residual - (low[0] / 0.25) * reference).isZero(1.e-28));
    REQUIRE(std::abs(residual.sum()) < 1.e-28);
}

TEST_CASE("Fitted DG box exposes both material-side drive derivatives",
          "[density_gradient][quantum][coupling]")
{
    const auto mesh = twoMaterialStripMesh();
    DensityGradientQuantumPotentialConfig config;
    config.globalDiscretization = "sentaurus_box";
    config.sentaurusInterfaceSiliconReactionWeight = 1.3;
    config.sentaurusInterfaceInsulatorAtSiliconReactionWeight = 2.1;
    constexpr Real vt = 0.026;
    std::vector<DensityGradientCellMaterial> cells(mesh.numCells());
    for (Index c = 0; c < cells.size(); ++c) {
        auto& data = cells[c];
        data.cellId = c;
        data.isTransport = c < 2;
        data.coefficientVm2 = densityGradientCoefficientVm2(data.isTransport ? 3.6 : 1.0, 0.5);
        for (int k = 0; k < 3; ++k) {
            const Index node = mesh.getCell(c).node_ids[k];
            data.materialBandDrive_V[k] = (data.isTransport ? 0.12 : 0.08) + 0.01 * node;
            data.dynamicDrivingPotential_V[k] = (data.isTransport ? -0.02 : 0.013) * node;
        }
    }
    const VectorXd band = VectorXd::Constant(6, 0.1);
    const VectorXd state = VectorXd::LinSpaced(6, -0.16, -0.07);
    const std::vector<bool> active{true, true, true, false, true, true};
    const std::unordered_map<Index, Real> dirichlet{{0, 0.01}, {5, 0.0}};
    const auto evaluate = [&](const auto& input, const VectorXd& drive, const VectorXd& phi) {
        return linearizeElectronDensityGradientBox(
            mesh, input, drive, active, dirichlet, vt, PhysicalUnitSystem::legacySI(), config, phi);
    };
    const auto jacobian = evaluate(cells, band, state);
    const VectorXd phiDirection = VectorXd::LinSpaced(6, 0.3, -0.2);
    const VectorXd materialDirection = VectorXd::LinSpaced(6, -0.7, 0.5);
    const VectorXd transportDirection = VectorXd::LinSpaced(6, -0.2, -0.8);
    const VectorXd insulatorDirection = VectorXd::LinSpaced(6, 0.7, -0.35);
    const auto perturb = [&](Real step) {
        auto shifted = cells;
        for (auto& data : shifted) {
            const auto& nodes = mesh.getCell(data.cellId).node_ids;
            for (int k = 0; k < 3; ++k) {
                data.materialBandDrive_V[k] += step * materialDirection[nodes[k]];
                data.dynamicDrivingPotential_V[k] += step *
                    (data.isTransport ? transportDirection[nodes[k]] : insulatorDirection[nodes[k]]);
            }
        }
        return evaluate(shifted, band + step * materialDirection, state + step * phiDirection).residual;
    };
    constexpr Real delta = 1.0e-7;
    const VectorXd finiteDifference = (perturb(delta) - perturb(-delta)) / (2.0 * delta);
    const VectorXd assembled = jacobian.potentialLike * phiDirection
        + jacobian.materialDriveDiagonal.cwiseProduct(materialDirection)
        + jacobian.transportDrive * transportDirection
        + jacobian.insulatorDrive * insulatorDirection;
    for (int row = 0; row < 6; ++row)
        REQUIRE(assembled[row] == Catch::Approx(finiteDifference[row]).epsilon(1.0e-7));
    REQUIRE(jacobian.residual[3] == state[3]);
    REQUIRE(jacobian.materialDriveDiagonal[0] == 1.0);
    REQUIRE(jacobian.materialDriveDiagonal[3] == 0.0);
    REQUIRE(jacobian.transportDrive.coeff(1, 0) != 0.0);
    REQUIRE(jacobian.insulatorDrive.coeff(1, 2) != 0.0);
}

TEST_CASE("Coupled DG Newton preserves referenced quasi-Fermi contact voltages",
          "[density_gradient][quantum][coupling]")
{
    auto mesh = stripMesh();
    mesh.addContact({0, "source", 0, {0, 3}, {}});
    mesh.addContact({1, "drain", 0, {2, 5}, {}});
    MaterialDatabase materials;
    auto silicon = materials.getMaterial("Si");
    silicon.name = "Silicon";
    materials.addMaterial(silicon);
    const auto doping = DopingModel::fromMeshAndRegions(
        mesh, std::vector<RegionDopingSpec>{{"silicon", 1.e21, 0.}});
    NewtonConfig cfg;
    cfg.maxIter = 40;
    cfg.warmStart = true;
    cfg.quasiFermiReference = "contact_basin";
    cfg.electronQuantumPotential.enabled = true;
    cfg.electronQuantumPotential.includeInsulators = true;
    cfg.electronQuantumPotential.globalDiscretization = "sentaurus_box";
    cfg.electronQuantumPotential.couplingMode = "monolithic_box";
    const NewtonSolver equilibrium(mesh, materials, doping,
        {{"source", 0.}, {"drain", 0.}}, cfg);
    const auto initial = equilibrium.solve();
    CAPTURE(initial.failureDiagnostics.failureReason);
    REQUIRE(initial.converged);

    constexpr Real bias = 1.e-3;
    const NewtonSolver biased(mesh, materials, doping, {{"source", 0.}, {"drain", bias}}, cfg);
    const auto result = biased.solve(initial.solution);
    CAPTURE(result.failureDiagnostics.failureReason);
    REQUIRE(result.converged);
    REQUIRE(result.convergenceReason == "four_field_dg_newton+electron_density_gradient");
    REQUIRE(result.solution.hasReferencedElectronQuasiFermi());
    REQUIRE(result.solution.hasReferencedHoleQuasiFermi());
    REQUIRE(result.solution.electronQuantumPotentialLike.size() == mesh.numNodes());
    REQUIRE(result.solution.electronQuantumPotentialLike.allFinite());
    for (int node : {0, 3})
        REQUIRE(result.solution.electronQuasiFermiAt(node) == Catch::Approx(0.).margin(1.e-14));
    for (int node : {2, 5})
        REQUIRE(std::abs(result.solution.electronQuasiFermiAt(node)) == Catch::Approx(bias));
}

TEST_CASE("Local DG transport tangent matches the live quantum field response",
          "[density_gradient][quantum][coupling]")
{
    const auto mesh = stripMesh();
    MaterialDatabase materials;
    auto silicon = materials.getMaterial("Si");
    silicon.name = "Silicon";
    materials.addMaterial(silicon);
    const DopingModel doping(mesh.numNodes());
    MobilityModelConfig mobility;
    mobility.jacobianFieldDerivatives = false;
    CarrierStatisticsConfig statistics;
    statistics.model = "fermi_dirac";
    DensityGradientQuantumPotentialConfig quantum;
    quantum.enabled = true;
    RecombinationModelConfig recombination;
    recombination.srhDopingDependence.densityCoupling = "quantum_band_edge";
    CoupledDDAssembler assembler(mesh, materials, doping, constants::Vt_300,
        mobility, recombination, {}, {}, {}, {}, {}, {}, statistics, quantum);
    REQUIRE(assembler.supportsLocalDgTransportTangent());

    constexpr int nodes = 6;
    const VectorXd state = assembler.pack({
        VectorXd::LinSpaced(nodes, 0.15, 0.2),
        VectorXd::LinSpaced(nodes, -0.03, 0.05),
        VectorXd::Zero(nodes)});
    const VectorXd shift = VectorXd::LinSpaced(nodes, -0.01, 0.02);
    const VectorXd direction = VectorXd::LinSpaced(nodes, -0.4, 0.3);
    CoupledDDBoundaryConditions boundary;
    boundary.phin.emplace(0, 0.);
    boundary.phip.emplace(5, 0.);
    assembler.setElectronQuantumPotential(shift);
    SparseMatrixd tangent;
    assembler.assembleJacobianWithDgTransportTangent(state, boundary, tangent);
    const VectorXd action = tangent * direction;

    constexpr Real step = 1.e-6;
    assembler.setElectronQuantumPotential(shift + step * direction);
    const VectorXd plus = assembler.residual(state, boundary);
    assembler.setElectronQuantumPotential(shift - step * direction);
    const VectorXd reference = (plus - assembler.residual(state, boundary)) / (2. * step);
    for (int offset : {nodes, 2 * nodes}) {
        CAPTURE(offset);
        const auto expected = reference.segment(offset, nodes);
        REQUIRE(expected.cwiseAbs().maxCoeff() > 0.);
        REQUIRE((action.segment(offset, nodes) - expected).cwiseAbs().maxCoeff()
            <= 2.e-6 * expected.cwiseAbs().maxCoeff());
    }
    REQUIRE(action[nodes] == 0.);
    REQUIRE(action[3 * nodes - 1] == 0.);

    const VectorXd equilibrium = assembler.pack({
        VectorXd::Constant(nodes, 0.15), VectorXd::Zero(nodes), VectorXd::Zero(nodes)});
    assembler.setElectronQuantumPotential(VectorXd::Zero(nodes));
    SparseMatrixd equilibriumTangent;
    assembler.assembleJacobianWithDgTransportTangent(equilibrium, boundary, equilibriumTangent);
    REQUIRE(equilibriumTangent.nonZeros() == tangent.nonZeros());
    REQUIRE(std::equal(tangent.outerIndexPtr(), tangent.outerIndexPtr() + tangent.outerSize() + 1,
                       equilibriumTangent.outerIndexPtr()));
    REQUIRE(std::equal(tangent.innerIndexPtr(), tangent.innerIndexPtr() + tangent.nonZeros(),
                       equilibriumTangent.innerIndexPtr()));
    REQUIRE(std::count(equilibriumTangent.valuePtr(),
                       equilibriumTangent.valuePtr() + equilibriumTangent.nonZeros(), 0.) > 0);
}

TEST_CASE("DEVSIM oxide WKB penetration depth matches the analytic barrier",
          "[density_gradient][quantum][oxide_oracle]")
{
    const Real depth = densityGradientOxidePenetrationDepthM(0.4, 3.15);
    const Real hbar = constants::h / (2.0 * std::acos(-1.0));
    const Real expected = hbar / std::sqrt(
        2.0 * 0.4 * constants::m0 * constants::q * 3.15);
    REQUIRE(depth == Catch::Approx(expected).epsilon(1.0e-14));
    REQUIRE(depth * 1.0e9 == Catch::Approx(0.1739).epsilon(2.0e-3));
}

TEST_CASE("Global Eq 231 rejects unsupported experimental controls",
          "[density_gradient][quantum][oxide_oracle]")
{
    const DeviceMesh mesh = twoMaterialStripMesh();
    const VectorXd drive = VectorXd::Zero(6);
    const std::vector<bool> active(6, true);
    DensityGradientQuantumPotentialConfig config;
    config.globalDiscretization = "unknown";
    REQUIRE_THROWS_AS(
        solveElectronDensityGradientPotentialLikeGlobal(
            mesh, {}, drive, active, {}, constants::Vt_300,
            PhysicalUnitSystem::legacySI(), config),
        std::invalid_argument);

    config.globalDiscretization = "p1_direct";
    config.oxideBoundary = "devsim_wkb";
    config.oxideBarrierHeight_V = 0.0;
    REQUIRE_THROWS_AS(
        solveElectronDensityGradientPotentialLikeGlobal(
            mesh, {}, drive, active, {}, constants::Vt_300,
            PhysicalUnitSystem::legacySI(), config),
        std::invalid_argument);

    config.globalDiscretization = "gss_density_fitted";
    config.oxideBarrierHeight_V = 3.15;
    REQUIRE_THROWS_AS(
        solveElectronDensityGradientPotentialLikeGlobal(
            mesh, {}, drive, active, {}, constants::Vt_300,
            PhysicalUnitSystem::legacySI(), config),
        std::invalid_argument);
}

TEST_CASE("GSS fitted operator keeps Lambda continuous across a band step",
          "[density_gradient][quantum][gss][material_interface][manufactured]")
{
    const DeviceMesh mesh = twoMaterialStripMesh();
    const Real vt = constants::Vt_300;
    std::vector<DensityGradientCellMaterial> materials;
    materials.reserve(mesh.numCells());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        DensityGradientCellMaterial material;
        material.cellId = cellId;
        material.isTransport = cellId < 2;
        material.coefficientVm2 = densityGradientCoefficientVm2(
            material.isTransport ? 3.6 : 1.0,
            material.isTransport ? 1.0 : 0.42);
        const Real bandDrive = material.isTransport ? 0.25 : 3.4;
        const Real sideW = material.isTransport ? -0.3 : -8.0;
        for (int local = 0; local < 3; ++local) {
            material.materialBandDrive_V[local] = bandDrive;
            // Each material has a constant but discontinuous sqrt(n) trace.
            // The fitted normal flux is therefore zero on both sides while
            // the common exact Lambda remains zero.
            material.dynamicDrivingPotential_V[local] =
                sideW * vt - bandDrive;
        }
        materials.push_back(material);
    }
    DensityGradientQuantumPotentialConfig config;
    config.globalDiscretization = "gss_density_fitted";
    config.maxIterations = 20;
    config.damping = 1.0;
    config.relativeTolerance = 0.0;
    config.absoluteTolerance_V = 1.0e-12;
    const auto result = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, VectorXd::Zero(6), std::vector<bool>(6, true),
        {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}}, vt,
        PhysicalUnitSystem::legacySI(), config,
        VectorXd::Constant(6, 1.0e-3));
    REQUIRE(result.converged);
    REQUIRE(result.residualInfinityNorm < 1.0e-9);
    REQUIRE(result.potential_V.lpNorm<Eigen::Infinity>() < 1.0e-10);
    REQUIRE(result.potentialLike_V.size() == 0);

    // Keep the same manufactured solution with both inactive and prescribed
    // rows, including an inactive node that also carries a Dirichlet value.
    std::vector<bool> active(6, true);
    active[0] = active[3] = false;
    const auto constrained = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, VectorXd::Zero(6), active,
        {{0, 0.0}, {2, 0.0}, {5, 0.0}}, vt, PhysicalUnitSystem::legacySI(), config,
        VectorXd::Constant(6, 1.0e-3));
    REQUIRE(constrained.converged);
    REQUIRE(constrained.potential_V.lpNorm<Eigen::Infinity>() < 1.0e-10);

    config.globalDiscretization = "p1_lambda_direct";
    const auto expanded = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, VectorXd::Zero(6), std::vector<bool>(6, true),
        {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}}, vt,
        PhysicalUnitSystem::legacySI(), config,
        VectorXd::Constant(6, 1.0e-3));
    REQUIRE(expanded.converged);
    REQUIRE(expanded.residualInfinityNorm < 1.0e-9);
    REQUIRE(expanded.potential_V.lpNorm<Eigen::Infinity>() < 1.0e-10);
    REQUIRE(expanded.potentialLike_V.size() == 0);
}

TEST_CASE("GSS fitted operator preserves the linear sqrt-density limit",
          "[density_gradient][quantum][gss][manufactured]")
{
    const DeviceMesh mesh = stripMesh();
    const Real vt = constants::Vt_300;
    const Real u[6] = {
        1.0, 1.0001, 1.0002, 1.0, 1.0001, 1.0002};
    std::vector<DensityGradientCellMaterial> materials;
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        DensityGradientCellMaterial material;
        material.cellId = cellId;
        material.coefficientVm2 = densityGradientCoefficientVm2(3.6, 1.0);
        const Cell& cell = mesh.getCell(cellId);
        for (int local = 0; local < 3; ++local) {
            const int node = static_cast<int>(cell.node_ids[local]);
            material.dynamicDrivingPotential_V[local] =
                2.0 * vt * std::log(u[node]);
        }
        materials.push_back(material);
    }
    DensityGradientQuantumPotentialConfig config;
    config.globalDiscretization = "gss_density_fitted";
    config.maxIterations = 30;
    config.damping = 1.0;
    config.relativeTolerance = 0.0;
    config.absoluteTolerance_V = 1.0e-12;
    const auto result = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, VectorXd::Zero(6), std::vector<bool>(6, true),
        {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}}, vt,
        PhysicalUnitSystem::legacySI(), config,
        VectorXd::Constant(6, 1.0e-4));
    REQUIRE(result.converged);
    REQUIRE(result.residualInfinityNorm < 1.0e-9);
    REQUIRE(result.potential_V.lpNorm<Eigen::Infinity>() < 1.0e-9);
}

TEST_CASE("Sentaurus box uses the maximum material-side reaction trace",
          "[density_gradient][quantum][sentaurus_box][material_interface][manufactured]")
{
    const DeviceMesh mesh = twoMaterialStripMesh();
    const Real vt = constants::Vt_300;
    std::vector<DensityGradientCellMaterial> materials;
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        DensityGradientCellMaterial material;
        material.cellId = cellId;
        material.isTransport = cellId < 2;
        material.coefficientVm2 = densityGradientCoefficientVm2(
            material.isTransport ? 3.6 : 1.0,
            material.isTransport ? 1.0 : 0.42);
        const Real bandDrive = material.isTransport ? 0.2 : 0.5;
        for (int local = 0; local < 3; ++local) {
            material.materialBandDrive_V[local] = bandDrive;
            material.dynamicDrivingPotential_V[local] = 0.0;
        }
        materials.push_back(material);
    }

    DensityGradientQuantumPotentialConfig config;
    config.globalDiscretization = "sentaurus_box";
    config.maxIterations = 3;
    config.relativeTolerance = 0.0;
    config.absoluteTolerance_V = 1.0e-12;
    // Deliberately choose the Silicon-side output owner. sentaurus_box must
    // still use the maximum material-side drive (0.5 V) at shared vertices.
    const VectorXd outputBandDrive = VectorXd::Constant(6, 0.2);
    const VectorXd exactPotentialLike = VectorXd::Constant(6, -0.5);
    const auto result = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, outputBandDrive, std::vector<bool>(6, true),
        {{0, -0.3}, {2, 0.0}, {3, -0.3}, {5, 0.0}}, vt,
        PhysicalUnitSystem::legacySI(), config, VectorXd::Zero(6),
        exactPotentialLike);

    REQUIRE(result.converged);
    REQUIRE(result.residualInfinityNorm < 1.0e-12);
    REQUIRE(std::abs(result.potential_V(1)) < 1.0e-12);
    REQUIRE(std::abs(result.potential_V(4)) < 1.0e-12);
    REQUIRE((result.potentialLike_V - exactPotentialLike)
                .lpNorm<Eigen::Infinity>() < 1.0e-12);

    SECTION("explicit interface trace offset has a restart-consistent root") {
        config.sentaurusInterfaceInsulatorHalfJumpOffset = 0.02012;
        config.sentaurusInterfaceSiliconHalfJumpOffset = -4.6008840569922854e-5;
        config.sentaurusInterfacePolysiliconHalfJumpOffset =
            0.0026674992132365016;
        config.sentaurusInterfaceSiliconReactionWeight = 0.3370294035925238;
        config.sentaurusInterfacePolysiliconReactionWeight = 1.044775636206889;
        config.sentaurusInterfaceInsulatorAtSiliconReactionWeight =
            2.5714393713226125;
        config.sentaurusInterfaceInsulatorAtPolysiliconReactionWeight =
            2.591446049803917;
        config.sentaurusInterfaceSiliconReactionOffset_V = -2.0e-4;
        config.sentaurusInterfaceInsulatorAtSiliconReactionOffset_V = -1.5e-3;
        config.maxIterations = 100;
        config.damping = 0.5;
        const auto shifted = solveElectronDensityGradientPotentialLikeGlobal(
            mesh, materials, outputBandDrive, std::vector<bool>(6, true),
            {{0, -0.3}, {2, 0.0}, {3, -0.3}, {5, 0.0}}, vt,
            PhysicalUnitSystem::legacySI(), config, VectorXd::Zero(6),
            exactPotentialLike);
        REQUIRE(shifted.converged);
        REQUIRE(shifted.residualInfinityNorm < 1.0e-9);
        REQUIRE((shifted.potentialLike_V - exactPotentialLike)
                    .lpNorm<Eigen::Infinity>() > 1.0e-8);

        const auto restarted = solveElectronDensityGradientPotentialLikeGlobal(
            mesh, materials, outputBandDrive, std::vector<bool>(6, true),
            {{0, -0.3}, {2, 0.0}, {3, -0.3}, {5, 0.0}}, vt,
            PhysicalUnitSystem::legacySI(), config, shifted.potential_V,
            shifted.potentialLike_V);
        REQUIRE(restarted.converged);
        REQUIRE((restarted.potentialLike_V - shifted.potentialLike_V)
                    .lpNorm<Eigen::Infinity>() < 1.0e-9);
    }
}

TEST_CASE("Sentaurus box keeps the signed AverageBox measure on an obtuse cell",
          "[density_gradient][quantum][sentaurus_box][average_box][obtuse]")
{
    const DeviceMesh mesh = obtuseTriangleMesh();
    const Real x[3] = {0.0, 2.0e-9, 0.5e-9};
    const Real y[3] = {0.0, 0.0, 0.2e-9};
    const Real area = 0.5 * std::abs(
        (x[1] - x[0]) * (y[2] - y[0]) -
        (x[2] - x[0]) * (y[1] - y[0]));
    const Real b[3] = {y[1] - y[2], y[2] - y[0], y[0] - y[1]};
    const Real c[3] = {x[2] - x[1], x[0] - x[2], x[1] - x[0]};
    const auto distanceSquared = [&](int first, int second) {
        return (x[second] - x[first]) * (x[second] - x[first]) +
            (y[second] - y[first]) * (y[second] - y[first]);
    };
    const auto cotangent = [&](int vertex, int first, int second) {
        const Real firstX = x[first] - x[vertex];
        const Real firstY = y[first] - y[vertex];
        const Real secondX = x[second] - x[vertex];
        const Real secondY = y[second] - y[vertex];
        return (firstX * secondX + firstY * secondY) /
            std::abs(firstX * secondY - firstY * secondX);
    };
    const Real averageBoxMeasure = 0.125 * (
        distanceSquared(2, 1) * cotangent(0, 2, 1) +
        distanceSquared(2, 0) * cotangent(1, 2, 0));
    REQUIRE(averageBoxMeasure > 0.5 * area);

    const Real w[3] = {0.05, 0.05, 0.0};
    Real fittedFlux = 0.0;
    for (int neighbour : {0, 1}) {
        const Real stiffness = area *
            (b[2] * b[neighbour] + c[2] * c[neighbour]) /
            (4.0 * area * area);
        const Real halfJump = 0.5 * (w[neighbour] - w[2]);
        fittedFlux += -2.0 * stiffness *
            (halfJump + 0.5 * halfJump * halfJump);
    }
    const Real coefficient = densityGradientCoefficientVm2(3.6, 1.0);
    const Real exactLambda =
        -fittedFlux * coefficient / (2.0 * averageBoxMeasure);

    DensityGradientCellMaterial material;
    material.cellId = 0;
    material.coefficientVm2 = coefficient;
    for (int local = 0; local < 3; ++local) {
        material.dynamicDrivingPotential_V[local] =
            exactLambda + constants::Vt_300 * w[local];
    }
    DensityGradientQuantumPotentialConfig config;
    config.globalDiscretization = "sentaurus_box";
    config.maxIterations = 4;
    config.damping = 1.0;
    config.relativeTolerance = 0.0;
    config.absoluteTolerance_V = 1.0e-12;
    const auto result = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, {material}, VectorXd::Zero(3), std::vector<bool>(3, true),
        {{0, exactLambda}, {1, exactLambda}}, constants::Vt_300,
        PhysicalUnitSystem::legacySI(), config,
        VectorXd::Constant(3, exactLambda));

    REQUIRE(result.converged);
    REQUIRE(result.residualInfinityNorm < 1.0e-9);
    REQUIRE(result.potential_V(2) ==
            Catch::Approx(exactLambda).margin(1.0e-11));
}

TEST_CASE("Sentaurus box closes a two-four insulator re-entrant corner",
          "[density_gradient][quantum][sentaurus_box][insulator_corner]")
{
    const DeviceMesh mesh = twoInsulatorReentrantCornerMesh();
    const Real vt = constants::Vt_300;
    std::vector<DensityGradientCellMaterial> materials;
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        DensityGradientCellMaterial material;
        material.cellId = cellId;
        material.isTransport = false;
        const bool oxide = cellId < 2;
        material.coefficientVm2 = densityGradientCoefficientVm2(
            1.0, oxide ? 0.5 : 0.42);
        const Real bandDrive = oxide ? 0.2 : 0.5;
        for (int local = 0; local < 3; ++local) {
            material.materialBandDrive_V[local] = bandDrive;
            material.dynamicDrivingPotential_V[local] = 0.0;
        }
        materials.push_back(material);
    }

    DensityGradientQuantumPotentialConfig config;
    config.globalDiscretization = "sentaurus_box";
    config.maxIterations = 80;
    config.damping = 0.5;
    config.relativeTolerance = 0.0;
    config.absoluteTolerance_V = 1.0e-12;
    const std::unordered_map<Index, Real> boundary = {
        {1, 0.1}, {2, 0.0}, {3, 0.0},
        {4, 0.0}, {5, 0.0}, {6, 0.0}};
    const VectorXd outputBandDrive = VectorXd::Constant(7, 0.5);
    const auto generic = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, outputBandDrive, std::vector<bool>(7, true),
        boundary, vt, PhysicalUnitSystem::legacySI(), config,
        VectorXd::Zero(7), VectorXd::Zero(7));
    REQUIRE(generic.converged);

    config.sentaurusInsulatorReentrantCornerReactionWeight = 0.5;
    const auto closed = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, outputBandDrive, std::vector<bool>(7, true),
        boundary, vt, PhysicalUnitSystem::legacySI(), config,
        generic.potential_V, generic.potentialLike_V);
    REQUIRE(closed.converged);
    REQUIRE(closed.residualInfinityNorm < 1.0e-9);
    REQUIRE(std::abs(closed.potentialLike_V(0) -
                     generic.potentialLike_V(0)) > 1.0e-6);

    const auto restarted = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, outputBandDrive, std::vector<bool>(7, true),
        boundary, vt, PhysicalUnitSystem::legacySI(), config,
        closed.potential_V, closed.potentialLike_V);
    REQUIRE(restarted.converged);
    REQUIRE((restarted.potentialLike_V - closed.potentialLike_V)
                .lpNorm<Eigen::Infinity>() < 1.0e-9);
}

TEST_CASE("Sentaurus Eq 233 step function is stable at zero",
          "[density_gradient][quantum][step_boundary]")
{
    REQUIRE(densityGradientStepFunction(0.0) ==
            Catch::Approx(0.5).margin(1.0e-15));
    REQUIRE(densityGradientStepFunctionDerivative(0.0) ==
            Catch::Approx(1.0 / 6.0).margin(1.0e-15));
    const Real x = 0.37;
    const Real h = 1.0e-6;
    const Real finiteDifference = (
        densityGradientStepFunction(x + h) -
        densityGradientStepFunction(x - h)) / (2.0 * h);
    REQUIRE(densityGradientStepFunctionDerivative(x) ==
            Catch::Approx(finiteDifference).epsilon(1.0e-7));
    REQUIRE(densityGradientStepFunction(-100.0) > 0.0);
}

TEST_CASE("Sentaurus Eq 233 step boundary changes semiconductor quantum state",
          "[density_gradient][quantum][step_boundary]")
{
    const DeviceMesh mesh = stripMesh();
    VectorXd density(6);
    density << 1.0e20, 1.0e22, 1.0e22, 1.0e20, 1.0e22, 1.0e22;
    const VectorXd coefficient = VectorXd::Constant(
        6, densityGradientCoefficientVm2(3.6, 1.0));
    DensityGradientQuantumPotentialConfig config;
    config.maxIterations = 120;
    config.damping = 0.4;
    const std::vector<bool> active = {true, true, false, true, true, false};
    const auto homogeneous = solveElectronDensityGradientPotential(
        mesh, density, coefficient, active, {{0, 0.0}, {3, 0.0}},
        constants::Vt_300, PhysicalUnitSystem::legacySI(), config);
    const Index interfaceEdge = [&] {
        for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
            const Edge& edge = mesh.getEdge(edgeId);
            if ((edge.n0 == 1 && edge.n1 == 4) ||
                (edge.n0 == 4 && edge.n1 == 1))
                return edgeId;
        }
        return mesh.numEdges();
    }();
    REQUIRE(interfaceEdge < mesh.numEdges());
    DensityGradientStepBoundary boundary;
    boundary.edgeId = interfaceEdge;
    boundary.barrierHeightN0_V = 3.15;
    boundary.barrierHeightN1_V = 3.15;
    boundary.barrierEffectiveMassRatio = 0.42;
    boundary.solvedGamma = 3.6;
    boundary.theta = 0.5;
    const auto stepped = solveElectronDensityGradientPotential(
        mesh, density, coefficient, active, {{0, 0.0}, {3, 0.0}},
        {boundary}, constants::Vt_300, PhysicalUnitSystem::legacySI(), config,
        homogeneous.potential_V);
    REQUIRE(stepped.potential_V.allFinite());
    REQUIRE(stepped.converged);
    REQUIRE((stepped.potential_V - homogeneous.potential_V)
                .lpNorm<Eigen::Infinity>() > 1.0e-6);
}

TEST_CASE("Global potential-like formulation solves one-material Eq 231",
          "[density_gradient][quantum][material_interface]")
{
    const DeviceMesh mesh = stripMesh();
    VectorXd density(6);
    density << 1.0e20, 1.0e22, 1.0e20, 1.0e20, 1.0e22, 1.0e20;
    DensityGradientQuantumPotentialConfig config;
    config.maxIterations = 120;
    config.damping = 0.5;
    const auto diagnosticDirectory =
        std::filesystem::temp_directory_path() / "vela_eq231_residual_test";
    std::filesystem::remove_all(diagnosticDirectory);
    config.residualDiagnosticPrefix =
        (diagnosticDirectory / "audit").string();
    const std::vector<bool> active(6, true);
    const auto reference = solveElectronDensityGradientPotential(
        mesh, density, active, {{0, 0.0}, {3, 0.0}}, constants::Vt_300,
        PhysicalUnitSystem::legacySI(), config);
    std::vector<DensityGradientCellMaterial> materials;
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        DensityGradientCellMaterial material;
        material.cellId = cellId;
        material.coefficientVm2 = densityGradientCoefficientVm2(3.6, 1.0);
        const Cell& cell = mesh.getCell(cellId);
        for (int local = 0; local < 3; ++local) {
            const int node = static_cast<int>(cell.node_ids[local]);
            material.dynamicDrivingPotential_V[local] =
                constants::Vt_300 * std::log(density(node) / 1.0e22);
            material.initialLambda_V[local] = reference.potential_V(node);
        }
        materials.push_back(material);
    }
    const auto global = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, VectorXd::Zero(6), active,
        {{0, 0.0}, {3, 0.0}}, constants::Vt_300,
        PhysicalUnitSystem::legacySI(), config, reference.potential_V);
    REQUIRE(global.converged);
    REQUIRE(global.potential_V.allFinite());
    REQUIRE(global.potentialLike_V.size() == 6);
    REQUIRE(global.potential_V(1) > 0.0);
    REQUIRE(global.potential_V(4) > 0.0);
    REQUIRE(std::filesystem::exists(diagnosticDirectory / "audit_cells.csv"));
    REQUIRE(std::filesystem::exists(diagnosticDirectory / "audit_nodes.csv"));
    REQUIRE(std::filesystem::exists(diagnosticDirectory / "audit_regions.csv"));
    REQUIRE(std::filesystem::exists(diagnosticDirectory / "audit_summary.txt"));
    std::ifstream nodeDiagnostic(diagnosticDirectory / "audit_nodes.csv");
    std::string header;
    std::getline(nodeDiagnostic, header);
    REQUIRE(header.find("stiffness,gradient_squared,reaction,raw_total") !=
            std::string::npos);
    nodeDiagnostic.close();
    config.residualDiagnosticPrefix.clear();
    const auto restarted = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, VectorXd::Zero(6), active,
        {{0, 0.0}, {3, 0.0}}, constants::Vt_300,
        PhysicalUnitSystem::legacySI(), config, VectorXd::Zero(6),
        global.potentialLike_V);
    REQUIRE(restarted.converged);
    REQUIRE((restarted.potential_V - global.potential_V)
                .lpNorm<Eigen::Infinity>() < 1.0e-7);
    std::filesystem::remove_all(diagnosticDirectory);
}

TEST_CASE("Material-resolved density-gradient coefficient scales material rows",
          "[density_gradient][quantum][material_interface]")
{
    const DeviceMesh mesh = stripMesh();
    VectorXd density(6);
    density << 1.0e20, 1.0e22, 1.0e20, 1.0e20, 1.0e22, 1.0e20;
    VectorXd coefficient = VectorXd::Constant(
        6, densityGradientCoefficientVm2(3.6, 1.0));
    coefficient(2) = densityGradientCoefficientVm2(1.0, 0.42);
    coefficient(5) = coefficient(2);
    DensityGradientQuantumPotentialConfig config;
    config.maxIterations = 80;
    const auto result = solveElectronDensityGradientPotential(
        mesh, density, coefficient, std::vector<bool>(6, true),
        {{0, 0.0}, {3, 0.0}}, constants::Vt_300,
        PhysicalUnitSystem::legacySI(), config);
    REQUIRE(result.potential_V.allFinite());
    REQUIRE(result.potential_V(1) > 0.0);
}

TEST_CASE("Global Eq 231 cancels electrostatic drive in an insulator",
          "[density_gradient][quantum][material_interface][manufactured]")
{
    const DeviceMesh mesh = twoMaterialStripMesh();
    constexpr Real electricSlope_V_per_m = 2.0e7;
    const Real vt = constants::Vt_300;
    std::vector<DensityGradientCellMaterial> materials;
    materials.reserve(mesh.numCells());
    VectorXd outputShift(6);
    VectorXd potentialLike(6);
    for (int node = 0; node < 6; ++node) {
        const Real psi = electricSlope_V_per_m * mesh.getNode(node).x;
        outputShift(node) = psi;
        potentialLike(node) = -psi;
    }
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        DensityGradientCellMaterial material;
        material.cellId = cellId;
        material.isTransport = cellId < 2;
        material.coefficientVm2 = densityGradientCoefficientVm2(
            material.isTransport ? 3.6 : 1.0,
            material.isTransport ? 1.0 : 0.42);
        const Cell& cell = mesh.getCell(cellId);
        for (int local = 0; local < 3; ++local) {
            const int node = static_cast<int>(cell.node_ids[local]);
            const Real psi = electricSlope_V_per_m * mesh.getNode(node).x;
            if (material.isTransport) {
                // xi=eta=1: choosing phin=psi makes A=0 and Lambda=0.
                material.materialBandDrive_V[local] = psi;
                material.dynamicDrivingPotential_V[local] = -psi;
            } else {
                // xi=eta=0: Phi/q=-psi makes
                // A=grad((-psi-Phi/q)/Vt)=0. The region-side material shift
                // reconstructs Lambda=0 even at the shared interface nodes.
                material.materialBandDrive_V[local] = psi;
                material.dynamicDrivingPotential_V[local] = -psi;
            }
        }
        materials.push_back(material);
    }
    DensityGradientQuantumPotentialConfig config;
    config.maxIterations = 5;
    config.absoluteTolerance_V = 1.0e-12;
    config.globalDiscretization = "exponential_fitted";
    const auto result = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, outputShift, std::vector<bool>(6, true),
        {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}}, vt,
        PhysicalUnitSystem::legacySI(), config, VectorXd::Zero(6),
        potentialLike);
    REQUIRE(result.potential_V.allFinite());
    REQUIRE(result.converged);
    REQUIRE(result.residualInfinityNorm < 1.0e-9);
    REQUIRE(result.potential_V.lpNorm<Eigen::Infinity>() < 1.0e-10);
}

TEST_CASE("DEVSIM MOSCAP oxide interface oracle applies the WKB face source",
          "[density_gradient][quantum][oxide_oracle]")
{
    const DeviceMesh mesh = twoMaterialStripMesh();
    std::vector<DensityGradientCellMaterial> materials;
    materials.reserve(mesh.numCells());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        DensityGradientCellMaterial material;
        material.cellId = cellId;
        material.isTransport = cellId < 2;
        material.coefficientVm2 = densityGradientCoefficientVm2(
            material.isTransport ? 3.6 : 1.0,
            material.isTransport ? 1.0 : 0.42);
        materials.push_back(material);
    }
    DensityGradientQuantumPotentialConfig config;
    config.maxIterations = 80;
    config.damping = 0.5;
    config.absoluteTolerance_V = 1.0e-11;
    config.globalDiscretization = "cvfem_full";
    const auto noBoundary = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, VectorXd::Zero(6), std::vector<bool>(6, true),
        {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}}, constants::Vt_300,
        PhysicalUnitSystem::legacySI(), config);
    REQUIRE(noBoundary.converged);
    REQUIRE(noBoundary.potential_V.lpNorm<Eigen::Infinity>() < 1.0e-12);

    config.oxideBoundary = "devsim_wkb";
    config.oxideQuantumMassRatio = 0.14;
    config.oxideBarrierMassRatio = 0.4;
    config.oxideBarrierHeight_V = 3.15;
    const auto withBoundary = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, VectorXd::Zero(6), std::vector<bool>(6, true),
        {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}}, constants::Vt_300,
        PhysicalUnitSystem::legacySI(), config);
    REQUIRE(withBoundary.converged);
    REQUIRE(withBoundary.potential_V.allFinite());
    REQUIRE(withBoundary.potential_V(1) > 0.0);
    REQUIRE(withBoundary.potential_V(4) > 0.0);
    // The diagonal split makes the two interface control volumes slightly
    // different; the WKB response should nevertheless agree to a few percent.
    REQUIRE(withBoundary.potential_V(1) ==
            Catch::Approx(withBoundary.potential_V(4)).epsilon(3.0e-2));
}

TEST_CASE("CVFEM full Eq 231 preserves the zero manufactured state",
          "[density_gradient][quantum][cvfem][manufactured]")
{
    const DeviceMesh mesh = twoMaterialStripMesh();
    std::vector<DensityGradientCellMaterial> materials;
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        DensityGradientCellMaterial material;
        material.cellId = cellId;
        material.isTransport = cellId < 2;
        material.coefficientVm2 = densityGradientCoefficientVm2(
            material.isTransport ? 3.6 : 1.0,
            material.isTransport ? 1.0 : 0.42);
        materials.push_back(material);
    }
    DensityGradientQuantumPotentialConfig config;
    config.globalDiscretization = "cvfem_full";
    config.maxIterations = 3;
    config.absoluteTolerance_V = 1.0e-13;
    const auto result = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, VectorXd::Zero(6), std::vector<bool>(6, true),
        {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}}, constants::Vt_300,
        PhysicalUnitSystem::legacySI(), config);
    REQUIRE(result.converged);
    REQUIRE(result.residualInfinityNorm < 1.0e-13);
    REQUIRE(result.potential_V.lpNorm<Eigen::Infinity>() < 1.0e-13);
}

TEST_CASE("Exponential-fitted Eq 231 preserves a linear auxiliary field",
          "[density_gradient][quantum][material_interface][manufactured]")
{
    const DeviceMesh mesh = stripMesh();
    const Real vt = constants::Vt_300;
    std::vector<DensityGradientCellMaterial> materials;
    const Real u[6] = {1.0, 1.2, 1.4, 1.0, 1.2, 1.4};
    VectorXd exactPotential(6);
    for (int node = 0; node < 6; ++node)
        exactPotential(node) = -2.0 * vt * std::log(u[node]);
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        DensityGradientCellMaterial material;
        material.cellId = cellId;
        material.coefficientVm2 = densityGradientCoefficientVm2(3.6, 1.0);
        const Cell& cell = mesh.getCell(cellId);
        for (int local = 0; local < 3; ++local) {
            const int node = static_cast<int>(cell.node_ids[local]);
            material.dynamicDrivingPotential_V[local] = 0.0;
            // Cancel the reaction term while retaining the manufactured w.
            material.materialBandDrive_V[local] = -exactPotential(node);
        }
        materials.push_back(material);
    }
    DensityGradientQuantumPotentialConfig config;
    config.maxIterations = 5;
    config.absoluteTolerance_V = 1.0e-12;
    config.globalDiscretization = "exponential_fitted";
    const auto result = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, -exactPotential, std::vector<bool>(6, true),
        {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}}, vt,
        PhysicalUnitSystem::legacySI(), config, VectorXd::Zero(6),
        exactPotential);
    REQUIRE(result.converged);
    REQUIRE(result.residualInfinityNorm < 1.0e-12);
    REQUIRE((result.potentialLike_V - exactPotential)
                .lpNorm<Eigen::Infinity>() < 1.0e-12);

    config.globalDiscretization = "conservative_sqrt_fitted";
    const auto conservative = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, -exactPotential, std::vector<bool>(6, true),
        {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}}, vt,
        PhysicalUnitSystem::legacySI(), config, VectorXd::Zero(6),
        exactPotential);
    REQUIRE(conservative.converged);
    REQUIRE(conservative.residualInfinityNorm < 1.0e-12);
    REQUIRE((conservative.potentialLike_V - exactPotential)
                .lpNorm<Eigen::Infinity>() < 1.0e-12);

    config.globalDiscretization = "p1_direct";
    config.maxIterations = 1;
    const auto expanded = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, -exactPotential, std::vector<bool>(6, true),
        {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}}, vt,
        PhysicalUnitSystem::legacySI(), config, VectorXd::Zero(6),
        exactPotential);
    REQUIRE(expanded.residualInfinityNorm > 1.0e-4);
}

TEST_CASE("Global Eq 231 does not report convergence after line-search stagnation",
          "[density_gradient][quantum][material_interface][manufactured]")
{
    const DeviceMesh mesh = stripMesh();
    std::vector<DensityGradientCellMaterial> materials;
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        DensityGradientCellMaterial material;
        material.cellId = cellId;
        material.coefficientVm2 = densityGradientCoefficientVm2(3.6, 1.0);
        const Cell& cell = mesh.getCell(cellId);
        for (int local = 0; local < 3; ++local) {
            const int node = static_cast<int>(cell.node_ids[local]);
            // Force an extreme exponential ratio and a reaction target that
            // cannot be accepted by the first capped Newton direction.
            material.dynamicDrivingPotential_V[local] =
                node == 1 || node == 4 ? 4.0 : 0.0;
            material.materialBandDrive_V[local] = 1.0;
        }
        materials.push_back(material);
    }
    DensityGradientQuantumPotentialConfig config;
    config.maxIterations = 5;
    config.maxUpdate_V = 0.1;
    config.damping = 0.5;
    const auto result = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, VectorXd::Constant(6, 1.0),
        std::vector<bool>(6, true), {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}},
        constants::Vt_300, PhysicalUnitSystem::legacySI(), config,
        VectorXd::Zero(6), VectorXd::Zero(6));
    REQUIRE_FALSE(result.converged);
    REQUIRE(result.residualInfinityNorm > 1.0);
}

TEST_CASE("Global Eq 231 requires residual convergence as well as a small step",
          "[density_gradient][quantum][material_interface][manufactured]")
{
    const DeviceMesh mesh = stripMesh();
    std::vector<DensityGradientCellMaterial> materials;
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        DensityGradientCellMaterial material;
        material.cellId = cellId;
        material.coefficientVm2 = densityGradientCoefficientVm2(3.6, 1.0);
        const Cell& cell = mesh.getCell(cellId);
        for (int local = 0; local < 3; ++local) {
            const int node = static_cast<int>(cell.node_ids[local]);
            material.materialBandDrive_V[local] = 0.5;
            material.dynamicDrivingPotential_V[local] =
                (node == 1 || node == 4) ? 0.25 : 0.0;
        }
        materials.push_back(material);
    }
    DensityGradientQuantumPotentialConfig config;
    config.maxIterations = 2;
    config.maxUpdate_V = 1.0e-12;
    config.absoluteTolerance_V = 1.0e-10;
    config.relativeTolerance = 0.0;
    const auto result = solveElectronDensityGradientPotentialLikeGlobal(
        mesh, materials, VectorXd::Constant(6, 0.5),
        std::vector<bool>(6, true), {{0, 0.0}, {2, 0.0}, {3, 0.0}, {5, 0.0}},
        constants::Vt_300, PhysicalUnitSystem::legacySI(), config);
    REQUIRE_FALSE(result.converged);
    REQUIRE(result.residualInfinityNorm > config.absoluteTolerance_V);
}

TEST_CASE("Uniform density has zero electron quantum potential",
          "[density_gradient][quantum]")
{
    const DeviceMesh mesh = stripMesh();
    const VectorXd density = VectorXd::Constant(6, 1.0e22);
    const auto result = solveElectronDensityGradientPotential(
        mesh, density, std::vector<bool>(6, true), {{0, 0.0}, {3, 0.0}},
        constants::Vt_300);
    REQUIRE(result.converged);
    REQUIRE(result.potential_V.lpNorm<Eigen::Infinity>() < 1.0e-12);
}

TEST_CASE("Carrier confinement raises the electron quantum potential",
          "[density_gradient][quantum]")
{
    const DeviceMesh mesh = stripMesh();
    VectorXd density(6);
    density << 1.0e20, 1.0e22, 1.0e20, 1.0e20, 1.0e22, 1.0e20;
    DensityGradientQuantumPotentialConfig config;
    config.damping = 0.5;
    config.maxIterations = 80;
    const auto result = solveElectronDensityGradientPotential(
        mesh, density, std::vector<bool>(6, true), {{0, 0.0}, {3, 0.0}},
        constants::Vt_300, PhysicalUnitSystem::legacySI(), config);
    REQUIRE(result.converged);
    REQUIRE(result.potential_V(1) > 0.0);
    REQUIRE(result.potential_V(4) > 0.0);
    REQUIRE(result.potential_V(0) == 0.0);
    REQUIRE(result.potential_V(3) == 0.0);
}

TEST_CASE("Inactive insulator nodes are pinned and decoupled",
          "[density_gradient][quantum]")
{
    const DeviceMesh mesh = stripMesh();
    VectorXd density = VectorXd::Constant(6, 1.0e22);
    density(2) = 0.0;
    density(5) = 0.0;
    const std::vector<bool> active = {true, true, false, true, true, false};
    const auto result = solveElectronDensityGradientPotential(
        mesh, density, active, {{0, 0.0}, {3, 0.0}}, constants::Vt_300);
    REQUIRE(result.converged);
    REQUIRE(result.potential_V(2) == 0.0);
    REQUIRE(result.potential_V(5) == 0.0);
}

TEST_CASE("Material coefficient validation rejects an active nonpositive mass term",
          "[density_gradient][quantum][material_interface]")
{
    const DeviceMesh mesh = stripMesh();
    const VectorXd density = VectorXd::Constant(6, 1.0e22);
    VectorXd coefficient = VectorXd::Constant(
        6, densityGradientCoefficientVm2(3.6, 1.0));
    coefficient(2) = 0.0;
    REQUIRE_THROWS_AS(
        solveElectronDensityGradientPotential(
            mesh, density, coefficient, std::vector<bool>(6, true),
            {{0, 0.0}, {3, 0.0}}, constants::Vt_300),
        std::invalid_argument);
}

TEST_CASE("Density-gradient update cap stabilizes a steep imported potential",
          "[density_gradient][quantum]")
{
    const DeviceMesh mesh = stripMesh();
    VectorXd density(6);
    density << 1.0e5, 1.0e27, 1.0e5, 1.0e5, 1.0e27, 1.0e5;
    VectorXd initial(6);
    initial << 0.0, 0.31, -0.31, 0.0, 0.31, -0.31;
    DensityGradientQuantumPotentialConfig config;
    config.maxIterations = 80;
    config.maxUpdate_V = 0.05;
    const auto result = solveElectronDensityGradientPotential(
        mesh, density, std::vector<bool>(6, true), {{0, 0.0}, {3, 0.0}},
        constants::Vt_300, PhysicalUnitSystem::legacySI(), config, initial);
    REQUIRE(result.potential_V.allFinite());
    REQUIRE(std::isfinite(result.residualInfinityNorm));
}
