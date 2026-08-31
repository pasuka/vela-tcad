#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <Eigen/SparseLU>
#include <nlohmann/json.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/core/PerformanceProfiler.h"
#include "vela/core/UnitScaling.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/DDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/numerics/ResidualNorm.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/post/ContactCurrent.h"
#include "vela/solver/GummelSolver.h"
#include "vela/solver/NewtonSolver.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using namespace vela;

static DeviceMesh makePNMesh();

TEST_CASE("Newton JSON parses Sentaurus electron density-gradient controls",
          "[newton][density_gradient]")
{
    const NewtonConfig cfg = newtonConfigFromJson(nlohmann::json{
        {"electron_quantum_potential", {
            {"enabled", true},
            {"coupling_mode", "frozen"},
            {"transport_coupling", "sentaurus_exponential"},
            {"transport_coupling_weight", 0.75},
            {"formulation", "density_based"},
            {"gamma", 3.6},
            {"effective_mass_ratio", 1.0},
            {"coefficient_mass_ratio", 1.09},
            {"include_insulators", true},
            {"insulator_gamma", 1.0},
            {"insulator_effective_mass_ratio", 0.42},
            {"insulator_coefficient_mass_ratio", 0.40},
            {"interface_boundary", "sentaurus_step"},
            {"theta", 0.5},
            {"conduction_band_narrowing_fraction", 0.5},
            {"max_iterations", 12},
            {"outer_max_iterations", 9},
            {"relative_tolerance", 2.0e-7},
            {"absolute_tolerance_V", 3.0e-10},
            {"outer_absolute_tolerance_V", 5.0e-4},
            {"damping", 0.4},
            {"max_update_V", 0.08},
            {"outer_acceleration", "aitken"},
            {"outer_relaxation", 0.6},
            {"outer_relaxation_min", 0.2},
            {"outer_relaxation_max", 1.2},
            {"residual_diagnostic_prefix", "eq231/audit"},
            {"residual_diagnostic_use_initial_state", true},
            {"global_discretization", "gss_density_fitted"},
            {"sentaurus_interface_insulator_half_jump_offset", 0.02012},
            {"sentaurus_interface_silicon_half_jump_offset", -4.6e-5},
            {"sentaurus_interface_polysilicon_half_jump_offset", 0.00267},
            {"sentaurus_interface_silicon_reaction_weight", 0.34},
            {"sentaurus_interface_polysilicon_reaction_weight", 1.04},
            {"sentaurus_interface_insulator_at_silicon_reaction_weight", 2.57},
            {"sentaurus_interface_insulator_at_polysilicon_reaction_weight", 2.59},
            {"sentaurus_interface_silicon_reaction_offset_V", -2.0e-4},
            {"sentaurus_interface_polysilicon_reaction_offset_V", 3.0e-4},
            {"sentaurus_interface_insulator_at_silicon_reaction_offset_V", -1.5e-3},
            {"sentaurus_interface_insulator_at_polysilicon_reaction_offset_V", 1.7e-3},
            {"sentaurus_insulator_reentrant_corner_reaction_weight", 0.97},
            {"oxide_boundary", "none"},
            {"oxide_quantum_mass_ratio", 0.14},
            {"oxide_barrier_mass_ratio", 0.4},
            {"oxide_barrier_height_V", 3.15},
        }},
    });
    REQUIRE(cfg.electronQuantumPotential.enabled);
    REQUIRE(cfg.electronQuantumPotential.couplingMode == "frozen");
    REQUIRE(cfg.electronQuantumPotential.transportCoupling ==
            "sentaurus_exponential");
    REQUIRE(cfg.electronQuantumPotential.transportCouplingWeight ==
            Catch::Approx(0.75));
    REQUIRE(cfg.electronQuantumPotential.formulation == "density_based");
    REQUIRE(cfg.electronQuantumPotential.gamma == Catch::Approx(3.6));
    REQUIRE(cfg.electronQuantumPotential.effectiveMassRatio ==
            Catch::Approx(1.0));
    REQUIRE(cfg.electronQuantumPotential.coefficientMassRatio ==
            Catch::Approx(1.09));
    REQUIRE(cfg.electronQuantumPotential.includeInsulators);
    REQUIRE(cfg.electronQuantumPotential.insulatorGamma == Catch::Approx(1.0));
    REQUIRE(cfg.electronQuantumPotential.insulatorEffectiveMassRatio ==
            Catch::Approx(0.42));
    REQUIRE(cfg.electronQuantumPotential.insulatorCoefficientMassRatio ==
            Catch::Approx(0.40));
    REQUIRE(cfg.electronQuantumPotential.interfaceBoundary == "sentaurus_step");
    REQUIRE(cfg.electronQuantumPotential.theta == Catch::Approx(0.5));
    REQUIRE(cfg.electronQuantumPotential.conductionBandNarrowingFraction ==
            Catch::Approx(0.5));
    REQUIRE(cfg.electronQuantumPotential.maxIterations == 12);
    REQUIRE(cfg.electronQuantumPotential.outerMaxIterations == 9);
    REQUIRE(cfg.electronQuantumPotential.relativeTolerance == Catch::Approx(2.0e-7));
    REQUIRE(cfg.electronQuantumPotential.absoluteTolerance_V == Catch::Approx(3.0e-10));
    REQUIRE(cfg.electronQuantumPotential.outerAbsoluteTolerance_V.has_value());
    REQUIRE(cfg.electronQuantumPotential.outerAbsoluteTolerance_V.value() ==
            Catch::Approx(5.0e-4));
    REQUIRE(cfg.electronQuantumPotential.damping == Catch::Approx(0.4));
    REQUIRE(cfg.electronQuantumPotential.maxUpdate_V == Catch::Approx(0.08));
    REQUIRE(cfg.electronQuantumPotential.outerAcceleration == "aitken");
    REQUIRE(cfg.electronQuantumPotential.outerRelaxation == Catch::Approx(0.6));
    REQUIRE(cfg.electronQuantumPotential.outerRelaxationMin == Catch::Approx(0.2));
    REQUIRE(cfg.electronQuantumPotential.outerRelaxationMax == Catch::Approx(1.2));
    REQUIRE(cfg.electronQuantumPotential.residualDiagnosticPrefix ==
            "eq231/audit");
    REQUIRE(cfg.electronQuantumPotential.residualDiagnosticUseInitialState);
    REQUIRE(cfg.electronQuantumPotential.globalDiscretization ==
            "gss_density_fitted");
    REQUIRE(cfg.electronQuantumPotential.
                sentaurusInterfaceInsulatorHalfJumpOffset ==
            Catch::Approx(0.02012));
    REQUIRE(cfg.electronQuantumPotential.sentaurusInterfaceSiliconHalfJumpOffset ==
            Catch::Approx(-4.6e-5));
    REQUIRE(cfg.electronQuantumPotential.
                sentaurusInterfacePolysiliconHalfJumpOffset ==
            Catch::Approx(0.00267));
    REQUIRE(cfg.electronQuantumPotential.sentaurusInterfaceSiliconReactionWeight ==
            Catch::Approx(0.34));
    REQUIRE(cfg.electronQuantumPotential.
                sentaurusInterfacePolysiliconReactionWeight ==
            Catch::Approx(1.04));
    REQUIRE(cfg.electronQuantumPotential.
                sentaurusInterfaceInsulatorAtSiliconReactionWeight ==
            Catch::Approx(2.57));
    REQUIRE(cfg.electronQuantumPotential.
                sentaurusInterfaceInsulatorAtPolysiliconReactionWeight ==
            Catch::Approx(2.59));
    REQUIRE(cfg.electronQuantumPotential.
                sentaurusInterfaceSiliconReactionOffset_V ==
            Catch::Approx(-2.0e-4));
    REQUIRE(cfg.electronQuantumPotential.
                sentaurusInterfacePolysiliconReactionOffset_V ==
            Catch::Approx(3.0e-4));
    REQUIRE(cfg.electronQuantumPotential.
                sentaurusInterfaceInsulatorAtSiliconReactionOffset_V ==
            Catch::Approx(-1.5e-3));
    REQUIRE(cfg.electronQuantumPotential.
                sentaurusInterfaceInsulatorAtPolysiliconReactionOffset_V ==
            Catch::Approx(1.7e-3));
    REQUIRE(cfg.electronQuantumPotential.
                sentaurusInsulatorReentrantCornerReactionWeight ==
            Catch::Approx(0.97));
    REQUIRE(cfg.electronQuantumPotential.oxideBoundary == "none");
    REQUIRE(cfg.electronQuantumPotential.oxideQuantumMassRatio ==
            Catch::Approx(0.14));
    REQUIRE(cfg.electronQuantumPotential.oxideBarrierMassRatio ==
            Catch::Approx(0.4));
    REQUIRE(cfg.electronQuantumPotential.oxideBarrierHeight_V ==
            Catch::Approx(3.15));
}

TEST_CASE("Newton JSON rejects unknown density-gradient transport coupling",
          "[newton][density_gradient]")
{
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"electron_quantum_potential", {
                {"enabled", true},
                {"transport_coupling", "unknown"},
            }},
        }),
        std::invalid_argument);
}

TEST_CASE("Newton JSON parses and validates Sentaurus-default SRH density coupling",
          "[newton][density_gradient][srh][json]")
{
    const NewtonConfig cfg = newtonConfigFromJson(nlohmann::json{
        {"srh_density_coupling", "sentaurus_default"},
    });
    REQUIRE(cfg.srhDopingDependence.densityCoupling == "sentaurus_default");
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"srh_density_coupling", "unknown"},
        }),
        std::invalid_argument);
}

TEST_CASE("Newton JSON rejects explicit oxide plus WKB quantum closure",
          "[newton][density_gradient]")
{
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"electron_quantum_potential", {
                {"enabled", true},
                {"include_insulators", true},
                {"oxide_boundary", "devsim_wkb"},
            }},
        }),
        std::invalid_argument);
}

TEST_CASE("Frozen electron quantum potential shifts density and preserves flat-QF flux",
          "[newton][density_gradient][jacobian]")
{
    const DeviceMesh mesh = makePNMesh();
    MaterialDatabase materials;
    DopingModel doping(mesh.numNodes());
    const Real thermalVoltage = constants::Vt_300;
    const MobilityModelConfig mobility;
    const RecombinationModelConfig recombination =
        recombinationModelConfig({"none"});
    DensityGradientQuantumPotentialConfig quantumConfig;
    quantumConfig.enabled = true;
    CoupledDDAssembler assembler(
        mesh, materials, doping, thermalVoltage, mobility, recombination,
        {}, {}, {}, {}, {}, {}, {}, quantumConfig);
    VectorXd quantum = VectorXd::Constant(static_cast<int>(mesh.numNodes()), 0.04);
    assembler.setElectronQuantumPotential(quantum);
    CoupledDDState state;
    state.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phip = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    const VectorXd x = assembler.pack(state);
    const VectorXd density = assembler.electronDensity(x);
    REQUIRE(density(0) / assembler.intrinsicDensity().at(0) ==
            Catch::Approx(std::exp(-0.04 / thermalVoltage)).epsilon(1.0e-12));

    CoupledDDBoundaryConditions boundaries;
    const auto diagnostics = assembler.sgEdgeFluxDiagnostics(x, boundaries);
    for (const auto& edge : diagnostics) {
        REQUIRE(std::abs(edge.electronFlux) < 1.0e-12);
        REQUIRE(edge.electronDensity0_m3 == Catch::Approx(density(edge.node0)));
        REQUIRE(edge.electronDensity1_m3 == Catch::Approx(density(edge.node1)));
    }

    const SparseMatrixd analytic = assembler.assembleJacobian(x, boundaries);
    const SparseMatrixd finiteDifference =
        assembler.finiteDifferenceJacobian(x, boundaries, 1.0e-7);
    const Real denominator = std::max<Real>(finiteDifference.norm(), 1.0e-30);
    REQUIRE((analytic - finiteDifference).norm() / denominator < 5.0e-5);
}

TEST_CASE("Carrier transport couple override changes SG residual and Jacobian but not Poisson",
          "[newton][sg][carrier_transport_profile][jacobian]")
{
    DeviceMesh baselineMesh = makePNMesh();
    DeviceMesh candidateMesh = makePNMesh();
    const Edge target = candidateMesh.getEdge(0);
    candidateMesh.setTransportCouple(0, 0.5 * target.couple);

    MaterialDatabase materials;
    DopingModel baselineDoping(baselineMesh.numNodes());
    DopingModel candidateDoping(candidateMesh.numNodes());
    const RecombinationModelConfig noRecombination =
        recombinationModelConfig({"none"});
    CoupledDDAssembler baseline(
        baselineMesh, materials, baselineDoping, constants::Vt_300,
        MobilityModelConfig{}, noRecombination);
    CoupledDDAssembler candidate(
        candidateMesh, materials, candidateDoping, constants::Vt_300,
        MobilityModelConfig{}, noRecombination);

    const int n = static_cast<int>(candidateMesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(n, -0.02, 0.03);
    state.phin = VectorXd::LinSpaced(n, -0.01, 0.02);
    state.phip = VectorXd::LinSpaced(n, 0.015, -0.01);
    CoupledDDBoundaryConditions boundaries;
    const VectorXd baselineX = baseline.pack(state);
    const VectorXd candidateX = candidate.pack(state);

    const auto baselinePoisson = baseline.poissonTermDiagnostics(
        baselineX, boundaries);
    const auto candidatePoisson = candidate.poissonTermDiagnostics(
        candidateX, boundaries);
    for (Index node = 0; node < candidateMesh.numNodes(); ++node) {
        REQUIRE(candidatePoisson[node].dielectricFlux ==
                Catch::Approx(baselinePoisson[node].dielectricFlux));
    }

    const auto baselineEdges = baseline.sgEdgeFluxDiagnostics(
        baselineX, boundaries);
    const auto candidateEdges = candidate.sgEdgeFluxDiagnostics(
        candidateX, boundaries);
    const auto matchesTarget = [&](const CoupledDDEdgeFluxDiagnostic& edge) {
        return edge.node0 == target.n0 && edge.node1 == target.n1;
    };
    const auto baselineTarget = std::find_if(
        baselineEdges.begin(), baselineEdges.end(), matchesTarget);
    const auto candidateTarget = std::find_if(
        candidateEdges.begin(), candidateEdges.end(), matchesTarget);
    REQUIRE(baselineTarget != baselineEdges.end());
    REQUIRE(candidateTarget != candidateEdges.end());
    REQUIRE(candidateTarget->couple_m == Catch::Approx(0.5 * baselineTarget->couple_m));
    REQUIRE(candidateTarget->electronFlux ==
            Catch::Approx(0.5 * baselineTarget->electronFlux));
    REQUIRE(candidateTarget->holeFlux ==
            Catch::Approx(0.5 * baselineTarget->holeFlux));

    const SparseMatrixd analytic = candidate.assembleJacobian(
        candidateX, boundaries);
    const SparseMatrixd finiteDifference = candidate.finiteDifferenceJacobian(
        candidateX, boundaries, 1.0e-7);
    REQUIRE((analytic - finiteDifference).norm() /
            std::max<Real>(finiteDifference.norm(), 1.0e-30) < 5.0e-5);
}

TEST_CASE("Sentaurus exponential DG coupling multiplies classical Fermi density",
          "[newton][density_gradient][fermi]")
{
    const DeviceMesh mesh = makePNMesh();
    MaterialDatabase materials;
    DopingModel doping(mesh.numNodes());
    const Real thermalVoltage = constants::Vt_300;
    const MobilityModelConfig mobility;
    const RecombinationModelConfig recombination =
        recombinationModelConfig({"none"});
    const CarrierStatisticsConfig statistics{"fermi_dirac"};

    DensityGradientQuantumPotentialConfig directConfig;
    directConfig.enabled = true;
    directConfig.transportCoupling = "direct_band_edge";
    DensityGradientQuantumPotentialConfig exponentialConfig = directConfig;
    exponentialConfig.transportCoupling = "sentaurus_exponential";
    DensityGradientQuantumPotentialConfig classicalConfig;

    auto makeAssembler = [&](const DensityGradientQuantumPotentialConfig& config) {
        return CoupledDDAssembler(
            mesh, materials, doping, thermalVoltage, mobility, recombination,
            {}, {}, {}, {}, {}, {}, statistics, config);
    };
    CoupledDDAssembler direct = makeAssembler(directConfig);
    CoupledDDAssembler exponential = makeAssembler(exponentialConfig);
    CoupledDDAssembler classical = makeAssembler(classicalConfig);
    const VectorXd quantum =
        VectorXd::Constant(static_cast<int>(mesh.numNodes()), 0.04);
    direct.setElectronQuantumPotential(quantum);
    exponential.setElectronQuantumPotential(quantum);

    CoupledDDState state;
    state.psi = VectorXd::Constant(static_cast<int>(mesh.numNodes()), 0.8);
    state.phin = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phip = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    const VectorXd directDensity = direct.electronDensity(direct.pack(state));
    const VectorXd exponentialDensity =
        exponential.electronDensity(exponential.pack(state));
    const VectorXd classicalDensity =
        classical.electronDensity(classical.pack(state));

    REQUIRE(exponentialDensity(0) / classicalDensity(0) ==
            Catch::Approx(std::exp(-0.04 / thermalVoltage)).epsilon(1.0e-12));
    REQUIRE(std::abs(exponentialDensity(0) / directDensity(0) - 1.0) > 1.0e-3);

    CoupledDDBoundaryConditions boundaries;
    const auto diagnostics = exponential.sgEdgeFluxDiagnostics(
        exponential.pack(state), boundaries);
    for (const auto& edge : diagnostics)
        REQUIRE(std::abs(edge.electronFlux) < 1.0e-12);
}

TEST_CASE("SG edge diagnostic converts native line flux to particles per metre",
          "[newton][sg][unit_scaling]")
{
    const DeviceMesh mesh = makePNMesh();
    MaterialDatabase materials;
    DopingModel doping(mesh.numNodes());
    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = 1.0;
    scaling.C0 = 1.0;
    scaling.D0 = 1.0;
    scaling.L0 = 1.0;
    scaling.mu0 = 1.0;
    scaling.unitSystem = PhysicalUnitSystem::tcadInternal();
    scaling.fieldFromCoordinateDeltaFactor =
        scaling.unitSystem.fieldFromCoordinateDeltaFactor();
    scaling.currentDensityLineIntegralFactor =
        scaling.unitSystem.currentDensityAM2PerInternal() *
        scaling.unitSystem.lengthMPerInternal();
    CoupledDDAssembler assembler(
        mesh, materials, doping, constants::Vt_300, {},
        recombinationModelConfig({"none"}), {}, {}, {}, {}, scaling);
    CoupledDDState state;
    state.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phip = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin(1) = 0.01;
    CoupledDDBoundaryConditions boundaries;
    const auto diagnostics = assembler.sgEdgeFluxDiagnostics(
        assembler.pack(state), boundaries);
    bool foundNonzero = false;
    bool foundZeroCoupleResponse = false;
    for (const auto& edge : diagnostics) {
        if (edge.couple_m == 0.0 &&
            std::abs(edge.electronFluxPerInternalCouple) > 1.0e-30) {
            foundZeroCoupleResponse = true;
            REQUIRE(edge.electronFlux == 0.0);
            REQUIRE(edge.electronParticleLineFlux_per_m_s == 0.0);
        }
        if (std::abs(edge.electronFlux) < 1.0e-30)
            continue;
        foundNonzero = true;
        REQUIRE(edge.electronFluxPerInternalCouple ==
                Catch::Approx(edge.electronFlux / edge.couple_m));
        REQUIRE(edge.electronParticleLineFlux_per_m_s ==
                Catch::Approx(
                    edge.electronFlux * scaling.currentDensityLineIntegralFactor));
        REQUIRE(edge.electronParticleLineFluxPerInternalCouple_per_m_s ==
                Catch::Approx(
                    edge.electronParticleLineFlux_per_m_s / edge.couple_m));
    }
    REQUIRE(foundNonzero);
    REQUIRE(foundZeroCoupleResponse);
}

static DeviceMesh makePNMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;

    Node n0; n0.id = 0; n0.x = 0.0;     n0.y = 0.0;     mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;       n1.y = 0.0;     mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = L;       n2.y = L;       mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0;     n3.y = L;       mesh.addNode(n3);
    Node n4; n4.id = 4; n4.x = 0.5 * L; n4.y = 0.5 * L; mesh.addNode(n4);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0; c0.node_ids = {0, 1, 4}; mesh.addCell(c0);
    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 0; c1.node_ids = {1, 2, 4}; mesh.addCell(c1);
    Cell c2; c2.id = 2; c2.type = CellType::Tri3; c2.region_id = 1; c2.node_ids = {2, 3, 4}; mesh.addCell(c2);
    Cell c3; c3.id = 3; c3.type = CellType::Tri3; c3.region_id = 1; c3.node_ids = {3, 0, 4}; mesh.addCell(c3);

    Region r0; r0.id = 0; r0.name = "n_region"; r0.material = "Si"; r0.cell_ids = {0, 1}; mesh.addRegion(r0);
    Region r1; r1.id = 1; r1.name = "p_region"; r1.material = "Si"; r1.cell_ids = {2, 3}; mesh.addRegion(r1);

    Contact cathode; cathode.id = 0; cathode.name = "cathode"; cathode.region_id = 0; cathode.node_ids = {1, 2}; mesh.addContact(cathode);
    Contact anode; anode.id = 1; anode.name = "anode"; anode.region_id = 1; anode.node_ids = {0, 3}; mesh.addContact(anode);

    mesh.buildEdges();
    return mesh;
}

static DopingModel makePNDoping(const DeviceMesh& mesh)
{
    std::vector<RegionDopingSpec> specs = {
        {"n_region", 1.0e21, 0.0},
        {"p_region", 0.0, 1.0e21},
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

static DeviceMesh makeOxideMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;

    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;   n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = L;   n2.y = L;   mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0; n3.y = L;   mesh.addNode(n3);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0; c0.node_ids = {0, 1, 2}; mesh.addCell(c0);
    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 0; c1.node_ids = {0, 2, 3}; mesh.addCell(c1);

    Region r0; r0.id = 0; r0.name = "oxide"; r0.material = "SiO2"; r0.cell_ids = {0, 1}; mesh.addRegion(r0);

    mesh.buildEdges();
    return mesh;
}

static DeviceMesh makeContactedOxideMesh()
{
    DeviceMesh mesh = makeOxideMesh();
    Contact gate;
    gate.id = 0;
    gate.name = "gate";
    gate.region_id = 0;
    gate.node_ids = {0, 1, 2, 3};
    mesh.addContact(gate);
    return mesh;
}

static DeviceMesh makePartiallyContactedOxideMesh()
{
    DeviceMesh mesh = makeOxideMesh();
    Contact gate;
    gate.id = 0;
    gate.name = "gate";
    gate.region_id = 0;
    gate.node_ids = {0};
    mesh.addContact(gate);
    return mesh;
}

static std::unordered_map<std::string, Real> zeroBias()
{
    return {{"anode", 0.0}, {"cathode", 0.0}};
}

static NewtonConfig newtonConfig()
{
    NewtonConfig cfg;
    cfg.maxIter = 10;
    cfg.reltol = 1.0e-7;
    cfg.abstol = 1.0e-20;
    cfg.dampingFactor = 1.0;
    cfg.lineSearch = true;
    cfg.verbose = false;
    return cfg;
}

static void requireFiniteNewtonSolution(const NewtonResult& result, Index nodeCount);

TEST_CASE("NewtonSolver: PN diode equilibrium converges", "[newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), newtonConfig());

    REQUIRE(result.converged);
    REQUIRE(result.iters >= 0);
    REQUIRE(result.finalResidualNorm <= result.initialResidualNorm);
}

TEST_CASE("NewtonSolver: accepts a qualified initial numerical-floor restart",
          "[newton][restart][stall_floor]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    const NewtonResult equilibrium =
        runNewton(mesh, matdb, doping, zeroBias(), newtonConfig());
    REQUIRE(equilibrium.converged);

    DDSolution restart = equilibrium.solution;
    restart.psi(4) += 1.0e-10;

    NewtonConfig cfg = newtonConfig();
    cfg.warmStart = true;
    cfg.abstol = 1.0e-30;
    NewtonSolver auditSolver(mesh, matdb, doping, zeroBias(), cfg);
    const NewtonResidualEvaluation audit = auditSolver.evaluateResidual(restart);
    REQUIRE(audit.blockNorms.combined > cfg.abstol);
    cfg.stallResidualFloor = 2.0 * audit.blockNorms.combined;

    const NewtonResult reclosed =
        runNewton(mesh, matdb, doping, zeroBias(), restart, cfg);
    REQUIRE(reclosed.converged);
    REQUIRE(reclosed.iters == 0);
    REQUIRE(reclosed.convergenceReason == "initial_stall_residual_floor");

    DDSolution unsafeRestart = restart;
    unsafeRestart.phin(4) += 1.0e-3;
    NewtonSolver unsafeAuditSolver(mesh, matdb, doping, zeroBias(), cfg);
    const NewtonResidualEvaluation unsafeAudit =
        unsafeAuditSolver.evaluateResidual(unsafeRestart);
    NewtonConfig guardedCfg = cfg;
    guardedCfg.maxIter = 1;
    guardedCfg.stallResidualFloor = 2.0 * unsafeAudit.blockNorms.combined;
    const NewtonResult guarded =
        runNewton(mesh, matdb, doping, zeroBias(), unsafeRestart, guardedCfg);
    REQUIRE(guarded.convergenceReason != "initial_stall_residual_floor");
    REQUIRE(guarded.convergenceReason != "max_iter_stall_residual_floor");
}

TEST_CASE("NewtonSolver: enforced block absolute convergence prevents scalar tolerance masking",
          "[newton][convergence][block_absolute]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const NewtonResult equilibrium =
        runNewton(mesh, matdb, doping, zeroBias(), newtonConfig());
    REQUIRE(equilibrium.converged);

    DDSolution restart = equilibrium.solution;
    restart.phin(4) += 1.0e-7;
    NewtonConfig cfg = newtonConfig();
    cfg.warmStart = true;
    cfg.maxIter = 0;
    cfg.stallResidualFloor = 0.0;
    const NewtonResidualEvaluation audit =
        NewtonSolver(mesh, matdb, doping, zeroBias(), cfg).evaluateResidual(restart);
    cfg.abstol = 1.0e100;

    const NewtonResult scalarAccepted =
        runNewton(mesh, matdb, doping, zeroBias(), restart, cfg);
    REQUIRE(scalarAccepted.converged);
    REQUIRE(scalarAccepted.convergenceReason == "initial_abstol");

    cfg.blockAbsoluteConvergence.mode = "enforce";
    cfg.blockAbsoluteConvergence.psiResidualCeiling = 1.0e100;
    cfg.blockAbsoluteConvergence.electronResidualCeiling =
        std::max(0.5 * audit.blockNorms.phin, 1.0e-30);
    cfg.blockAbsoluteConvergence.holeResidualCeiling = 1.0e100;
    const NewtonResult blockRejected =
        runNewton(mesh, matdb, doping, zeroBias(), restart, cfg);
    REQUIRE_FALSE(blockRejected.converged);
    REQUIRE(blockRejected.failureDiagnostics.failureReason ==
            "block_absolute_convergence");

    cfg.blockAbsoluteConvergence.electronResidualCeiling = 1.0e100;
    const NewtonResult blockAccepted =
        runNewton(mesh, matdb, doping, zeroBias(), restart, cfg);
    REQUIRE(blockAccepted.converged);
    REQUIRE(blockAccepted.convergenceReason == "initial_block_abstol");
}

TEST_CASE("NewtonSolver: preserves the best accepted iterate on failure",
          "[newton][diagnostics][best_iterate]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const NewtonResult equilibrium =
        runNewton(mesh, matdb, doping, zeroBias(), newtonConfig());
    REQUIRE(equilibrium.converged);

    DDSolution restart = equilibrium.solution;
    restart.psi(4) += 1.0e-7;
    NewtonConfig cfg = newtonConfig();
    cfg.warmStart = true;
    cfg.maxIter = 0;
    cfg.abstol = 0.0;
    cfg.stallResidualFloor = 0.0;
    const NewtonResult failed =
        runNewton(mesh, matdb, doping, zeroBias(), restart, cfg);
    REQUIRE_FALSE(failed.converged);
    REQUIRE(failed.hasBestSolution);
    REQUIRE(failed.bestIteration == 0);
    REQUIRE(failed.bestResidualNorm == Catch::Approx(failed.initialResidualNorm));
    REQUIRE((failed.bestSolution.psi - failed.solution.psi).norm() ==
            Catch::Approx(0.0).margin(1.0e-20));
}

TEST_CASE("NewtonSolver: contact-majority branch guard rejects unsafe convergence",
          "[newton][contact][branch_guard]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    for (Index node = 0; node < mesh.numNodes(); ++node)
        doping.setNodeDoping(node, 1.0e21, 0.0);

    NewtonConfig cfg = newtonConfig();
    cfg.warmStart = true;
    cfg.maxIter = 0;
    cfg.abstol = 1.0e100;
    cfg.stallResidualFloor = 0.0;
    cfg.contactMajorityQfBranchDropLimit_V = 1.0e-6;

    const int n = static_cast<int>(mesh.numNodes());
    DDSolution unsafe;
    unsafe.psi = VectorXd::Zero(n);
    unsafe.phin = VectorXd::Zero(n);
    unsafe.phip = VectorXd::Zero(n);
    unsafe.n = VectorXd::Constant(n, 1.0e21);
    unsafe.p = VectorXd::Constant(n, 1.0e21);
    unsafe.phin(4) = 1.0e-3;

    const NewtonResult rejected =
        runNewton(mesh, matdb, doping, zeroBias(), unsafe, cfg);
    REQUIRE_FALSE(rejected.converged);
    REQUIRE(rejected.failureDiagnostics.failureReason ==
            "contact_majority_qf_branch_guard");
    REQUIRE_FALSE(rejected.hasBestSolution);
}

TEST_CASE("NewtonSolver: ohmic contact BC resists compensated-node polarity flips",
          "[newton][contact_bc]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    // Inject an opposite-sign outlier on one anode node to mimic imported
    // compensated/tie ownership artifacts.
    doping.setNodeDoping(0, 8.0e20, 0.0);

    NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), newtonConfig());
    REQUIRE(result.converged);

    // Both anode nodes should keep p-side built-in sign (negative psi at 0 V).
    REQUIRE(result.solution.psi(0) < 0.0);
    REQUIRE(result.solution.psi(3) < 0.0);
}

TEST_CASE("NewtonSolver: PN diode equilibrium converges with unit_scaling state",
          "[newton][scaling]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};

    NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), cfg);

    REQUIRE(result.converged);
    REQUIRE(result.iters >= 0);
    REQUIRE(result.finalResidualNorm <= result.initialResidualNorm);
    requireFiniteNewtonSolution(result, mesh.numNodes());
}

TEST_CASE("NewtonSolver: pseudo-arclength corrector advances a converged device point",
          "[newton][arclength]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};

    NewtonSolver solver(mesh, matdb, doping, zeroBias(), cfg);
    const NewtonResult equilibrium = solver.solve();
    REQUIRE(equilibrium.converged);

    ArclengthSystem system = solver.makeArclengthSystem("anode");

    ArclengthState point;
    point.x = solver.packArclengthState(equilibrium.solution);
    point.lambda = 0.0;

    // The packed equilibrium state already satisfies F(x, lambda = 0) ~ 0, so the
    // bordered corrector reuses the same scaled residual as the Newton solve.
    const VectorXd f0 = system.residual(point.x, point.lambda);
    REQUIRE(f0.lpNorm<Eigen::Infinity>() < 1.0e-6);

    PseudoArclengthConfig arcCfg;
    arcCfg.enabled = true;
    arcCfg.initialStep = 0.02;
    arcCfg.minStep = 1.0e-5;
    arcCfg.maxStep = 0.1;
    arcCfg.growthFactor = 1.1;
    arcCfg.shrinkFactor = 0.5;
    arcCfg.maxCorrectorIterations = 40;
    arcCfg.correctorTolerance = 1.0e-7;
    arcCfg.maxStepRetries = 8;
    arcCfg.parameterScale = 1.0;

    PseudoArclengthContinuation continuation(system, arcCfg);
    const ArclengthTangent tangent = continuation.computeTangent(point, +1.0);

    const ArclengthStepResult result =
        continuation.step(point, tangent, arcCfg.initialStep);
    REQUIRE(result.converged);
    REQUIRE(result.residualNorm <= arcCfg.correctorTolerance);

    // The continuation parameter (anode bias) advanced off equilibrium and the
    // corrected device state still solves the coupled drift-diffusion residual.
    REQUIRE(std::abs(result.state.lambda) > 1.0e-6);
    const VectorXd f1 = system.residual(result.state.x, result.state.lambda);
    REQUIRE(f1.lpNorm<Eigen::Infinity>() < 1.0e-6);
}

TEST_CASE("NewtonSolver: arclength quasi-Fermi cap uniformly scales bordered update",
          "[newton][arclength][scaling]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const int interiorNode = 4;
    doping.setNodeDoping(interiorNode, 0.0, 1.0e21);

    NewtonConfig uncappedCfg = newtonConfig();
    uncappedCfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    NewtonSolver uncappedSolver(mesh, matdb, doping, zeroBias(), uncappedCfg);
    ArclengthSystem uncappedSystem = uncappedSolver.makeArclengthSystem("anode");
    REQUIRE_FALSE(static_cast<bool>(uncappedSystem.limitUpdate));

    NewtonConfig cappedCfg = uncappedCfg;
    cappedCfg.quasiFermiUpdateLimit_V = 0.05;
    cappedCfg.quasiFermiUpdateLimitMinority_V = 0.01;
    NewtonSolver cappedSolver(mesh, matdb, doping, zeroBias(), cappedCfg);
    ArclengthSystem cappedSystem = cappedSolver.makeArclengthSystem("anode");
    REQUIRE(static_cast<bool>(cappedSystem.limitUpdate));

    const NewtonResult equilibrium = cappedSolver.solve();
    REQUIRE(equilibrium.converged);
    const Real potentialScale = cappedSolver.evaluateResidual(equilibrium.solution).potentialScale;
    const int N = static_cast<int>(mesh.numNodes());

    const VectorXd x = cappedSolver.packArclengthState(equilibrium.solution);
    VectorXd deltaX = VectorXd::Zero(3 * N);
    deltaX(interiorNode) = 2.0;
    deltaX(N + interiorNode) = 2.0;
    deltaX(2 * N + interiorNode) = -2.0;
    Real deltaLambda = 0.25;

    cappedSystem.limitUpdate(x, deltaX, deltaLambda);

    const Real expectedScale =
        cappedCfg.quasiFermiUpdateLimitMinority_V / (2.0 * potentialScale);
    REQUIRE(deltaX(interiorNode) == Catch::Approx(2.0 * expectedScale));
    REQUIRE(deltaLambda == Catch::Approx(0.25 * expectedScale));
    REQUIRE(std::abs(deltaX(N + interiorNode) * potentialScale) ==
            Catch::Approx(cappedCfg.quasiFermiUpdateLimitMinority_V).margin(1.0e-12));
    REQUIRE(std::abs(deltaX(2 * N + interiorNode) * potentialScale) ==
            Catch::Approx(cappedCfg.quasiFermiUpdateLimitMinority_V).margin(1.0e-12));
}
TEST_CASE("NewtonSolver: high-doping unit-scaled PN cold start reaches near-zero 0V current",
          "[newton][scaling][bgn]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e23, 0.0},
        {"p_region", 0.0, 1.0e23},
    });

    NewtonConfig cfg;
    cfg.maxIter = 30;
    cfg.reltol = 1.0e-6;
    cfg.abstol = 1.0e-18;
    cfg.dampingFactor = 1.0;
    cfg.lineSearch = true;
    cfg.verbose = false;
    cfg.maxUpdate = 5.0;
    cfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    cfg.mobility = mobilityModelConfig("caughey_thomas_field");
    cfg.recombination = {"srh", "auger"};
    cfg.bandgapNarrowing = bandgapNarrowingConfig("slotboom");

    const NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), cfg);

    REQUIRE(result.converged);
    requireFiniteNewtonSolution(result, mesh.numNodes());
    ContactCurrent current(mesh, matdb, doping, cfg.mobility, cfg.temperature_K);
    const ContactCurrentResult anode = current.compute(result.solution, "anode");
    REQUIRE(std::abs(anode.totalCurrent) < 1.0e-9);
}

TEST_CASE("ContactCurrent: preserved edge hole QF drop changes reporting only",
          "[contact_current][qf_floor]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e23, 0.0},
        {"p_region", 0.0, 1.0e23},
    });

    DDSolution solution;
    const int N = static_cast<int>(mesh.numNodes());
    solution.psi = VectorXd::Constant(N, -13.2);
    solution.phin = VectorXd::Constant(N, -12.8);
    solution.phip = VectorXd::Constant(N, -12.8);
    solution.n = VectorXd::Constant(N, 1.0e10);
    solution.p = VectorXd::Constant(N, 1.0e23);

    MobilityModelConfig mobility = mobilityModelConfig("masetti_field");
    mobility.highFieldDrivingForce = "quasi_fermi_gradient";
    ContactCurrent current(mesh, matdb, doping, mobility, 300.0);

    const ContactCurrentDetailedResult baseline =
        current.computeDetailed(solution, "anode");
    REQUIRE(std::abs(baseline.totals.holeCurrent) < 1.0e-30);

    REQUIRE_FALSE(baseline.edges.empty());
    const Index edgeId = baseline.edges.front().edgeId;
    ContactCurrentEdgeOverrides overrides;
    const Real ulpAtBias = std::abs(std::nextafter(-12.8, -std::numeric_limits<Real>::infinity()) + 12.8);
    overrides.holeQuasiFermiDropByEdge[edgeId] = -5.0 * ulpAtBias;

    const ContactCurrentDetailedResult reported =
        current.computeDetailed(solution, "anode", overrides);

    REQUIRE(reported.totals.electronCurrent ==
            Catch::Approx(baseline.totals.electronCurrent));
    REQUIRE(std::abs(reported.totals.holeCurrent) > 1.0e-30);
    REQUIRE(reported.totals.totalCurrent ==
            Catch::Approx(reported.totals.electronCurrent - reported.totals.holeCurrent));
    REQUIRE(reported.totals.totalCurrent !=
            Catch::Approx(baseline.totals.totalCurrent));

    const auto changed = std::find_if(
        reported.edges.begin(),
        reported.edges.end(),
        [edgeId](const ContactCurrentEdgeDiagnostic& edge) {
            return edge.edgeId == edgeId;
        });
    REQUIRE(changed != reported.edges.end());
    REQUIRE(changed->holeQfDropOverrideApplied);
    REQUIRE(changed->phip1 - changed->phip0 ==
            Catch::Approx(-5.0 * ulpAtBias));
}

TEST_CASE("ContactCurrent: raw residual method matches SG flux without avalanche",
          "[contact_current][residual]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    MobilityModelConfig mobility = mobilityModelConfig("constant");
    const RecombinationModelConfig noRecombination = recombinationModelConfig({"none"});
    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        mobility,
        noRecombination);

    CoupledDDState state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.03);
    state.phin = VectorXd::LinSpaced(N, -0.04, 0.01);
    state.phip = VectorXd::LinSpaced(N, 0.02, -0.03);

    const VectorXd x = assembler.pack(state);
    DDSolution solution;
    solution.psi = state.psi;
    solution.phin = state.phin;
    solution.phip = state.phip;
    solution.n = assembler.electronDensity(x);
    solution.p = assembler.holeDensity(x);

    ContactCurrent current(mesh, matdb, doping, mobility, constants::T0);
    const ContactCurrentResult sgFlux = current.compute(solution, "anode");
    const ContactCurrentResult residual =
        current.computeFromResidual(assembler, x, "anode");

    REQUIRE(residual.electronCurrent == Catch::Approx(sgFlux.electronCurrent));
    REQUIRE(residual.holeCurrent == Catch::Approx(sgFlux.holeCurrent));
    REQUIRE(residual.totalCurrent == Catch::Approx(sgFlux.totalCurrent));
}

TEST_CASE("ContactCurrent: Fermi-Dirac generalized SG matches coupled residual",
          "[contact_current][residual][fermi_dirac]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    MobilityModelConfig mobility = mobilityModelConfig("constant");
    const RecombinationModelConfig noRecombination = recombinationModelConfig({"none"});
    const CarrierStatisticsConfig statistics{"fermi_dirac"};
    CoupledDDAssembler assembler(
        mesh, matdb, doping, constants::Vt_300, mobility, noRecombination,
        {}, {}, {}, {}, {}, {}, statistics);

    CoupledDDState state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.01, 0.72);
    state.phin = VectorXd::LinSpaced(N, -0.025, 0.015);
    state.phip = VectorXd::LinSpaced(N, 0.01, -0.02);

    const VectorXd x = assembler.pack(state);
    DDSolution solution;
    solution.psi = state.psi;
    solution.phin = state.phin;
    solution.phip = state.phip;
    solution.n = assembler.electronDensity(x);
    solution.p = assembler.holeDensity(x);

    ContactCurrent current(
        mesh, matdb, doping, mobility, constants::T0, {}, {}, statistics);
    const ContactCurrentResult sgFlux = current.compute(solution, "anode");
    const ContactCurrentResult residual =
        current.computeFromResidual(assembler, x, "anode");

    REQUIRE(residual.electronCurrent == Catch::Approx(sgFlux.electronCurrent));
    REQUIRE(residual.holeCurrent == Catch::Approx(sgFlux.holeCurrent));
    REQUIRE(residual.totalCurrent == Catch::Approx(sgFlux.totalCurrent));
    REQUIRE(sgFlux.electronCurrent == Catch::Approx(
        sgFlux.electronDriftCurrent + sgFlux.electronDiffusionCurrent)
        .epsilon(1.0e-12).margin(1.0e-18));
    REQUIRE(sgFlux.holeCurrent == Catch::Approx(
        sgFlux.holeDriftCurrent + sgFlux.holeDiffusionCurrent)
        .epsilon(1.0e-12).margin(1.0e-18));
}

TEST_CASE("ContactCurrent: unit-scaled residual method matches SG flux without avalanche",
          "[contact_current][residual][scaling]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = constants::Vt_300;
    scaling.C0 = 1.0;
    scaling.mu0 = 1.0;
    scaling.D0 = 1.0;
    scaling.L0 = 1.0;
    scaling.permittivityReference_F_per_m = constants::eps0 * 11.7;
    scaling.fieldFromCoordinateDeltaFactor = 1.0e4;
    scaling.currentDensityLineIntegralFactor = 1.0e-2;

    MobilityModelConfig mobility = mobilityModelConfig("constant");
    const RecombinationModelConfig noRecombination = recombinationModelConfig({"none"});
    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        mobility,
        noRecombination,
        BandgapNarrowingConfig{},
        ImpactIonizationModelConfig{},
        {},
        {},
        scaling);

    CoupledDDState state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.03) / scaling.V0;
    state.phin = VectorXd::LinSpaced(N, -0.04, 0.01) / scaling.V0;
    state.phip = VectorXd::LinSpaced(N, 0.02, -0.03) / scaling.V0;

    const VectorXd x = assembler.pack(state);
    DDSolution solution;
    solution.psi = state.psi * scaling.V0;
    solution.phin = state.phin * scaling.V0;
    solution.phip = state.phip * scaling.V0;
    solution.n = assembler.electronDensity(x);
    solution.p = assembler.holeDensity(x);

    ContactCurrent current(mesh, matdb, doping, mobility, constants::T0, scaling);
    const ContactCurrentResult sgFlux = current.compute(solution, "anode");
    const ContactCurrentResult residual =
        current.computeFromResidual(assembler, x, "anode");

    REQUIRE(residual.electronCurrent == Catch::Approx(sgFlux.electronCurrent));
    REQUIRE(residual.holeCurrent == Catch::Approx(sgFlux.holeCurrent));
    REQUIRE(residual.totalCurrent == Catch::Approx(sgFlux.totalCurrent));
}

TEST_CASE("ContactCurrent: referenced quasi-Fermi increments preserve residual current",
          "[contact_current][residual][qf-reference]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    MobilityModelConfig mobility = mobilityModelConfig("constant");
    CoupledDDAssembler assembler(
        mesh, matdb, doping, constants::Vt_300, mobility,
        recombinationModelConfig({"none"}));
    assembler.setQuasiFermiReferences(0.0, -20.0);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -20.02, 0.03);
    state.phin = VectorXd::LinSpaced(N, -0.04, 0.01);
    state.phip = VectorXd::Constant(N, -20.0);
    state.phip(N / 2) += 2.0e-15;
    const VectorXd x = assembler.pack(state);

    DDSolution solution;
    solution.psi = state.psi;
    solution.phin = state.phin;
    solution.phip = state.phip;
    solution.phinIncrement = x.segment(N, N);
    solution.phipIncrement = x.segment(2 * N, N);
    solution.electronQfReference_V = 0.0;
    solution.holeQfReference_V = -20.0;
    solution.n = assembler.electronDensity(x);
    solution.p = assembler.holeDensity(x);

    ContactCurrent current(mesh, matdb, doping, mobility, constants::T0);
    const ContactCurrentResult sgFlux = current.compute(solution, "anode");
    const ContactCurrentResult residual =
        current.computeFromResidual(assembler, x, "anode");
    REQUIRE(residual.electronCurrent ==
            Catch::Approx(sgFlux.electronCurrent).epsilon(1.0e-12));
    REQUIRE(residual.holeCurrent ==
            Catch::Approx(sgFlux.holeCurrent).epsilon(1.0e-12));
    REQUIRE(residual.totalCurrent ==
            Catch::Approx(sgFlux.totalCurrent).epsilon(1.0e-12));
}

TEST_CASE("ContactCurrent: electron reference preserves a sub-ULP contact drop",
          "[contact_current][residual][qf-reference][deep-off]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    MobilityModelConfig mobility = mobilityModelConfig("constant");
    BandgapNarrowingConfig bgn;
    bgn.model = "old_slotboom";
    CoupledDDAssembler assembler(
        mesh, matdb, doping, constants::Vt_300, mobility,
        recombinationModelConfig({"none"}), bgn);
    assembler.setQuasiFermiReferences(1.1, 0.0);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.03);
    state.phin = VectorXd::Constant(N, 1.1);
    state.phip = VectorXd::Zero(N);
    VectorXd x = assembler.pack(state);
    x(N + 4) = 5.0e-18;

    DDSolution solution;
    const CoupledDDState physical = assembler.unpack(x);
    solution.psi = physical.psi;
    solution.phin = physical.phin;
    solution.phip = physical.phip;
    solution.phinIncrement = x.segment(N, N);
    solution.phipIncrement = x.segment(2 * N, N);
    solution.electronQfReference_V = 1.1;
    solution.holeQfReference_V = 0.0;
    solution.n = assembler.electronDensity(x);
    solution.p = assembler.holeDensity(x);

    // The physical output cannot represent this drop around 1.1 V, while the
    // reference-coordinate state still retains it for residual and terminal
    // current evaluation.
    REQUIRE(solution.phin(0) == solution.phin(4));
    REQUIRE(solution.phinIncrement(0) != solution.phinIncrement(4));

    ContactCurrent current(
        mesh, matdb, doping, mobility, constants::T0, {}, bgn);
    const ContactCurrentResult sgFlux = current.compute(solution, "anode");
    const ContactCurrentResult residual =
        current.computeFromResidual(assembler, x, "anode");
    CAPTURE(sgFlux.electronCurrent, residual.electronCurrent,
            sgFlux.totalCurrent, residual.totalCurrent,
            residual.electronCurrent / sgFlux.electronCurrent);
    REQUIRE(sgFlux.electronCurrent != 0.0);
    REQUIRE(residual.electronCurrent ==
            Catch::Approx(sgFlux.electronCurrent).epsilon(1.0e-12));
    REQUIRE(residual.totalCurrent ==
            Catch::Approx(sgFlux.totalCurrent).epsilon(1.0e-12));
}

TEST_CASE("ContactCurrent: contact-basin references preserve local sub-ULP flux",
          "[contact_current][qf-reference][contact-basin][precision]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    MobilityModelConfig mobility = mobilityModelConfig("constant");
    BandgapNarrowingConfig bgn;
    bgn.model = "old_slotboom";
    CoupledDDAssembler assembler(
        mesh, matdb, doping, constants::Vt_300, mobility,
        recombinationModelConfig({"none"}), bgn);
    assembler.setQuasiFermiReferences(1.1, 0.0);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.03);
    state.phin = VectorXd::Constant(N, 1.1);
    state.phip = VectorXd::Zero(N);
    VectorXd x = assembler.pack(state);
    x(N + 4) = 5.0e-18;

    DDSolution global;
    const CoupledDDState physical = assembler.unpack(x);
    global.psi = physical.psi;
    global.phin = physical.phin;
    global.phip = physical.phip;
    global.phinIncrement = x.segment(N, N);
    global.phipIncrement = x.segment(2 * N, N);
    global.electronQfReference_V = 1.1;
    global.holeQfReference_V = 0.0;
    global.n = assembler.electronDensity(x);
    global.p = assembler.holeDensity(x);

    DDSolution partitioned = global;
    partitioned.electronQfReference.resize(N);
    partitioned.electronQfReference << 1.1, 0.0, 0.0, 1.1, 1.1;
    partitioned.holeQfReference = VectorXd::Zero(N);
    for (int i = 0; i < N; ++i) {
        if (partitioned.electronQfReference(i) ==
            global.electronQfReference_V) {
            partitioned.phinIncrement(i) = global.phinIncrement(i);
            continue;
        }
        const long double physicalQf =
            static_cast<long double>(global.electronQfReference_V) +
            static_cast<long double>(global.phinIncrement(i));
        partitioned.phinIncrement(i) = static_cast<Real>(
            physicalQf - static_cast<long double>(
                partitioned.electronQfReference(i)));
    }

    ContactCurrent current(
        mesh, matdb, doping, mobility, constants::T0, {}, bgn,
        CarrierStatisticsConfig{"fermi_dirac"});
    const ContactCurrentDetailedResult globalResult =
        current.computeDetailed(global, "anode");
    const ContactCurrentDetailedResult partitionedResult =
        current.computeDetailed(partitioned, "anode");

    REQUIRE(globalResult.totals.electronCurrent != 0.0);
    int sameBasinEdges = 0;
    for (const ContactCurrentEdgeDiagnostic& partitionEdge :
         partitionedResult.edges) {
        const int node0 = static_cast<int>(partitionEdge.node0);
        const int node1 = static_cast<int>(partitionEdge.node1);
        if (partitioned.electronQfReference(node0) !=
            partitioned.electronQfReference(node1)) {
            continue;
        }
        const auto globalEdge = std::find_if(
            globalResult.edges.begin(), globalResult.edges.end(),
            [&](const ContactCurrentEdgeDiagnostic& candidate) {
                return candidate.edgeId == partitionEdge.edgeId;
            });
        REQUIRE(globalEdge != globalResult.edges.end());
        const Real edgeScale = std::max<Real>(
            1.0e-300, std::abs(globalEdge->electronCurrent));
        REQUIRE(std::abs(
                    partitionEdge.electronCurrent -
                    globalEdge->electronCurrent) <=
                1.0e-13 * edgeScale);
        ++sameBasinEdges;
    }
    REQUIRE(sameBasinEdges > 0);
    const Real comparisonScale = std::max<Real>(
        1.0e-300, std::abs(partitionedResult.totals.electronCurrent));
    REQUIRE(std::abs(
                partitionedResult.precision.electronCurrentCompensated -
                partitionedResult.totals.electronCurrent) <=
            1.0e-15 * comparisonScale);
    REQUIRE(std::isfinite(
        partitionedResult.precision.electronCurrentLongDoubleReference));
    REQUIRE(partitionedResult.precision.electronCurrentLongDoubleReference != 0.0);
}

TEST_CASE("CoupledDDAssembler: nodal quasi-Fermi repartition preserves equations",
          "[newton][qf-reference][contact-basin][invariance]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    MobilityModelConfig mobility = mobilityModelConfig("constant");
    BandgapNarrowingConfig bgn;
    bgn.model = "old_slotboom";
    CarrierStatisticsConfig statistics;
    statistics.model = "fermi_dirac";
    CoupledDDAssembler global(
        mesh, matdb, doping, constants::Vt_300, mobility,
        recombinationModelConfig({"none"}), bgn, {}, {}, {}, {}, {}, statistics);
    CoupledDDAssembler partitioned(
        mesh, matdb, doping, constants::Vt_300, mobility,
        recombinationModelConfig({"none"}), bgn, {}, {}, {}, {}, {}, statistics);
    global.setQuasiFermiReferences(1.1, 0.0);
    VectorXd electronReference(5);
    electronReference << 0.0, 1.1, 1.1, 0.0, 0.0;
    const VectorXd holeReference = VectorXd::Zero(5);
    partitioned.setQuasiFermiReferenceFields(
        electronReference, holeReference);

    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(5, -0.03, 0.08);
    state.phin.resize(5);
    state.phin << 0.0, 1.1, 1.1, 0.0, 0.27;
    state.phip = VectorXd::LinSpaced(5, -0.01, 0.02);
    const VectorXd globalX = global.pack(state);
    const VectorXd partitionedX = partitioned.pack(state);
    const CoupledDDBoundaryConditions bcs;

    const VectorXd globalResidual = global.residual(globalX, bcs);
    const VectorXd partitionedResidual = partitioned.residual(partitionedX, bcs);
    const Real residualScale = std::max<Real>(1.0, globalResidual.norm());
    REQUIRE((partitionedResidual - globalResidual).norm() <=
            2.0e-12 * residualScale);
    REQUIRE((partitioned.electronDensity(partitionedX) -
             global.electronDensity(globalX)).norm() <=
            2.0e-13 * global.electronDensity(globalX).norm());
    REQUIRE((partitioned.holeDensity(partitionedX) -
             global.holeDensity(globalX)).norm() <=
            2.0e-13 * global.holeDensity(globalX).norm());

    const SparseMatrixd globalJacobian = global.assembleJacobian(globalX, bcs);
    const SparseMatrixd partitionedJacobian =
        partitioned.assembleJacobian(partitionedX, bcs);
    const Real jacobianScale = std::max<Real>(1.0, globalJacobian.norm());
    REQUIRE((partitionedJacobian - globalJacobian).norm() <=
            2.0e-9 * jacobianScale);
}

TEST_CASE("CoupledDDAssembler: source-basin coordinate retains sub-ULP Newton update",
          "[newton][qf-reference][contact-basin][precision]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    MobilityModelConfig mobility = mobilityModelConfig("constant");
    CoupledDDAssembler global(
        mesh, matdb, doping, constants::Vt_300, mobility,
        recombinationModelConfig({"none"}));
    CoupledDDAssembler partitioned(
        mesh, matdb, doping, constants::Vt_300, mobility,
        recombinationModelConfig({"none"}));
    global.setQuasiFermiReferences(1.1, 0.0);
    VectorXd electronReference(5);
    electronReference << 0.0, 1.1, 1.1, 0.0, 0.0;
    partitioned.setQuasiFermiReferenceFields(
        electronReference, VectorXd::Zero(5));

    CoupledDDState state;
    state.psi = VectorXd::Zero(5);
    state.phin.resize(5);
    state.phin << 0.0, 1.1, 1.1, 0.0, 0.0;
    state.phip = VectorXd::Zero(5);
    VectorXd globalX = global.pack(state);
    VectorXd partitionedX = partitioned.pack(state);
    const int sourceInteriorRow = static_cast<int>(mesh.numNodes()) + 4;
    const Real delta = 5.0e-18;
    const Real globalBefore = globalX(sourceInteriorRow);
    const Real partitionedBefore = partitionedX(sourceInteriorRow);
    globalX(sourceInteriorRow) += delta;
    partitionedX(sourceInteriorRow) += delta;

    REQUIRE(globalX(sourceInteriorRow) == globalBefore);
    REQUIRE(partitionedX(sourceInteriorRow) != partitionedBefore);
    REQUIRE(partitioned.unpack(partitionedX).phin(4) == delta);
}

TEST_CASE("NewtonSolver: Gummel initial guess reduces Newton iterations", "[newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    NewtonResult fromDefault = runNewton(mesh, matdb, doping, zeroBias(), cfg);

    GummelConfig gcfg;
    gcfg.maxIter = 8;
    gcfg.reltol = 1.0e-6;
    gcfg.dampingPsi = 0.5;
    DDSolution gummel = runGummel(mesh, matdb, doping, zeroBias(), gcfg);
    NewtonResult fromGummel = runNewton(mesh, matdb, doping, zeroBias(), gummel, cfg);

    REQUIRE(fromDefault.converged);
    REQUIRE(fromGummel.converged);
    REQUIRE(fromGummel.iters <= fromDefault.iters);
}

TEST_CASE("NewtonSolver: unit_scaling accepts a physical Gummel warm initial guess",
          "[newton][scaling][warm_start]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    GummelConfig gcfg;
    gcfg.maxIter = 8;
    gcfg.reltol = 1.0e-6;
    gcfg.dampingPsi = 0.5;
    gcfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    DDSolution gummel = runGummel(mesh, matdb, doping, zeroBias(), gcfg);

    NewtonConfig cfg = newtonConfig();
    cfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    cfg.warmStart = true;

    NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), gummel, cfg);

    REQUIRE(result.converged);
    REQUIRE(result.finalResidualNorm <= result.initialResidualNorm);
    requireFiniteNewtonSolution(result, mesh.numNodes());
}

TEST_CASE("NewtonSolver: no NaN or Inf and carriers stay positive", "[newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), newtonConfig());
    REQUIRE(result.converged);

    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(std::isfinite(result.solution.psi(i)));
        REQUIRE(std::isfinite(result.solution.phin(i)));
        REQUIRE(std::isfinite(result.solution.phip(i)));
        REQUIRE(std::isfinite(result.solution.n(i)));
        REQUIRE(std::isfinite(result.solution.p(i)));
        REQUIRE(result.solution.n(i) > 0.0);
        REQUIRE(result.solution.p(i) > 0.0);
    }
}

TEST_CASE("CoupledDDAssembler: zero-mobility continuity rows are pinned", "[newton][coupled]")
{
    DeviceMesh mesh = makeOxideMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());
    CoupledDDAssembler assembler(mesh, matdb, doping, constants::Vt_300, 1.0e-6, 1.0e-6);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::Zero(N);
    state.phin = VectorXd::LinSpaced(N, 0.1, 0.4);
    state.phip = VectorXd::LinSpaced(N, -0.4, -0.1);
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    for (Index i = 0; i < mesh.numNodes(); ++i)
        bcs.psi[i] = 0.0;

    const VectorXd r = assembler.residual(x, bcs);
    for (int i = 0; i < N; ++i) {
        REQUIRE(r(N + i) == Catch::Approx(state.phin(i)));
        REQUIRE(r(2 * N + i) == Catch::Approx(state.phip(i)));
    }

    const SparseMatrixd J = assembler.finiteDifferenceJacobian(x, bcs);
    for (int i = 0; i < N; ++i) {
        REQUIRE(J.coeff(N + i, N + i) == Catch::Approx(1.0));
        REQUIRE(J.coeff(2 * N + i, 2 * N + i) == Catch::Approx(1.0));
    }

    Eigen::SparseLU<SparseMatrixd> lu;
    lu.compute(J);
    REQUIRE(lu.info() == Eigen::Success);
}

TEST_CASE("CoupledDDAssembler: profiling separates Jacobian assembly phases",
          "[newton][coupled][performance_profiler]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    CoupledDDAssembler assembler(
        mesh, matdb, doping, constants::Vt_300, 1.0e-6, 1.0e-6);

    const int nodeCount = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::Zero(nodeCount);
    state.phin = VectorXd::Zero(nodeCount);
    state.phip = VectorXd::Zero(nodeCount);
    const VectorXd x = assembler.pack(state);
    CoupledDDBoundaryConditions bcs;
    bcs.psi[0] = 0.0;
    bcs.phin[0] = 0.0;
    bcs.phip[0] = 0.0;

    PerformanceProfiler profiler({true, "unused.json"});
    VectorXd perturbed = x;
    perturbed(1) += 0.01;
    perturbed(nodeCount + 1) -= 0.002;
    perturbed(2 * nodeCount + 1) += 0.003;
    {
        ActivePerformanceProfilerScope active(&profiler);
        (void)assembler.assembleJacobian(x, bcs);
        (void)assembler.assembleJacobian(perturbed, bcs);
    }

    const nlohmann::json profile = profiler.toJson();
    std::unordered_map<std::string, std::uint64_t> stageCalls;
    for (const auto& stage : profile.at("stages"))
        stageCalls[stage.at("name").get<std::string>()] =
            stage.at("calls").get<std::uint64_t>();
    for (const std::string name : {
             "jacobian.topology", "jacobian.edge_physics",
             "jacobian.cell_physics", "jacobian.node_sources",
             "jacobian.boundary_rows", "jacobian.initialize_values",
             "jacobian.finalize_triplets"}) {
        REQUIRE(stageCalls.at(name) == 2);
    }
    REQUIRE(profile.at("counters").at("jacobian.pattern_observations") == 2);
    REQUIRE(profile.at("counters").value("jacobian.pattern_change_count", 0) == 0);
    REQUIRE(profile.at("counters").at("jacobian.pattern_build_calls") == 1);
    REQUIRE(profile.at("observations").at("jacobian.triplet_count").at("min") > 0.0);
    REQUIRE(profile.at("observations").at("jacobian.triplet_capacity").at("max") == 0.0);
    REQUIRE(profile.at("observations").at("jacobian.nonzero_count").at("min") > 0.0);
    REQUIRE(profile.at("observations")
                .at("jacobian.structural_nonzero_count").at("min") >=
            profile.at("observations").at("jacobian.nonzero_count").at("max"));
}

TEST_CASE("CoupledDDAssembler: edge avalanche profiling counters preserve exact accounting",
          "[newton][coupled][performance_profiler][impact]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    MobilityModelConfig mobility = mobilityModelConfig("constant");
    mobility.highFieldDrivingForce = "quasi_fermi_gradient";

    ImpactIonizationModelConfig impact;
    impact.model = "selberherr";
    impact.drivingForce = "eparallel";
    impact.generation = "current_density";
    impact.currentApproximation = "nodal_vector_current_reconstructed";
    impact.sourceVolumePolicy = "edge_box";
    impact.electronA = 1.0e8;
    impact.electronB = 1.0e-30;
    impact.holeA = 1.0e8;
    impact.holeB = 1.0e-30;

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        mobility,
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        impact);

    const int nodeCount = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(nodeCount, 0.0, -0.8);
    state.phin = VectorXd::LinSpaced(nodeCount, 0.0, -1.4);
    state.phip = VectorXd::LinSpaced(nodeCount, -0.2, -0.9);
    const VectorXd x = assembler.pack(state);

    PerformanceProfiler profiler({true, "unused.json"});
    {
        ActivePerformanceProfilerScope active(&profiler);
        (void)assembler.assembleJacobian(x, {});
    }

    const nlohmann::json counters = profiler.toJson().at("counters");
    const std::uint64_t evaluations =
        counters.at("jacobian.edge_avalanche_source_evaluations");
    const std::uint64_t baseEvaluations =
        counters.at("jacobian.edge_avalanche_base_evaluations");
    const std::uint64_t perturbedEvaluations =
        counters.at("jacobian.edge_avalanche_perturbed_evaluations");
    const std::uint64_t fluxRequests =
        counters.at("jacobian.edge_avalanche_neighbor_flux_requests");
    const std::uint64_t fluxHits =
        counters.at("jacobian.avalanche_flux_cache_hits");
    const std::uint64_t fluxMisses =
        counters.at("jacobian.avalanche_flux_cache_misses");
    const std::uint64_t baseFluxReuses =
        counters.at("jacobian.edge_avalanche_base_flux_reuses");
    const std::uint64_t perturbedFluxRecomputations =
        counters.at("jacobian.edge_avalanche_perturbed_flux_recomputations");
    const std::uint64_t carrierSideEvaluations =
        counters.at("jacobian.edge_avalanche_carrier_side_evaluations");
    const std::uint64_t carrierSideReuses =
        counters.at("jacobian.edge_avalanche_carrier_side_reuses");
    const std::uint64_t electricVectorReuses =
        counters.at("jacobian.edge_avalanche_electric_vector_reuses");
    const std::uint64_t electricVectorRecomputations =
        counters.at("jacobian.edge_avalanche_electric_vector_recomputations");

    REQUIRE(baseEvaluations > 0);
    REQUIRE(perturbedEvaluations > baseEvaluations);
    REQUIRE(evaluations == baseEvaluations + perturbedEvaluations);
    REQUIRE(fluxRequests == fluxHits + fluxMisses + baseFluxReuses);
    REQUIRE(perturbedFluxRecomputations == fluxMisses);
    REQUIRE(carrierSideEvaluations + carrierSideReuses == 2 * evaluations);
    REQUIRE(carrierSideReuses > 0);
    REQUIRE(electricVectorReuses + electricVectorRecomputations <= evaluations);
    REQUIRE(electricVectorReuses > 0);
    REQUIRE(electricVectorRecomputations > 0);
    REQUIRE(counters.at(
        "jacobian.edge_avalanche_precomputed_scatter_adds") > 0);
    REQUIRE(counters.at(
        "jacobian.edge_avalanche_preparsed_config_evaluations") ==
        evaluations);
    REQUIRE(counters.at("jacobian.edge_avalanche_nodal_reconstruction_calls") > 0);
    REQUIRE(counters.at("jacobian.edge_avalanche_direct_flux_evaluations") == 0);
    REQUIRE(counters.at("jacobian.edge_avalanche_direct_flux_skips") > 0);
}

TEST_CASE("CoupledDDAssembler: nodal Eparallel avalanche Jacobian matches residual finite differences",
          "[newton][coupled][impact][edge_flux_audit]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    MobilityModelConfig mobility = mobilityModelConfig("masetti_field");
    mobility.highFieldDrivingForce = "quasi_fermi_gradient";
    mobility.highFieldGradientDiscretization = "transport_cell_vector";
    mobility.jacobianFieldDerivatives = true;

    ImpactIonizationModelConfig impact =
        impactIonizationModelConfig("van_overstraeten");
    impact.couplingMode = "self_consistent";
    impact.drivingForce = "eparallel";
    impact.generation = "current_density";
    impact.currentApproximation = "nodal_vector_current_reconstructed";
    impact.currentMagnitudeMode = "edge_scalar_abs";
    impact.sourceVolumePolicy = "genius_truncated";
    impact.sourceMappingMode = "node_F_node_alpha_node_G";

    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = constants::Vt_300;
    scaling.C0 = 1.0;
    scaling.mu0 = 1.0;
    scaling.D0 = 1.0;
    scaling.L0 = 1.0;
    scaling.permittivityReference_F_per_m = constants::eps0 * 11.7;
    scaling.fieldFromCoordinateDeltaFactor = 1.0e4;
    scaling.currentDensityLineIntegralFactor = 1.0e-2;

    const CarrierStatisticsConfig statistics{"fermi_dirac"};
    const auto makeAssembler = [&](const ImpactIonizationModelConfig& config) {
        return std::make_unique<CoupledDDAssembler>(
            mesh, matdb, doping, constants::Vt_300, mobility,
            recombinationModelConfig({"none"}), BandgapNarrowingConfig{},
            config, std::vector<RegionFixedChargeSpec>{},
            std::vector<InterfaceSheetChargeSpec>{}, scaling,
            CarrierDiagonalFloorRegularizationConfig{}, statistics);
    };
    const auto baseline = makeAssembler(ImpactIonizationModelConfig{});
    const auto coupled = makeAssembler(impact);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.10, 0.85) / scaling.V0;
    state.phin = VectorXd::LinSpaced(N, -0.15, 0.20) / scaling.V0;
    state.phip = VectorXd::LinSpaced(N, 0.10, -0.35) / scaling.V0;
    const VectorXd x = coupled->pack(state);

    const Eigen::MatrixXd analyticImpact = Eigen::MatrixXd(
        coupled->assembleJacobian(x, {}) - baseline->assembleJacobian(x, {}));
    const Eigen::MatrixXd residualFiniteDifferenceImpact = Eigen::MatrixXd(
        coupled->finiteDifferenceJacobian(x, {}, 1.0e-7) -
        baseline->finiteDifferenceJacobian(x, {}, 1.0e-7));
    const Real relativeError =
        (analyticImpact - residualFiniteDifferenceImpact).norm() /
        std::max<Real>(1.0, residualFiniteDifferenceImpact.norm());
    CAPTURE(relativeError);
    REQUIRE(relativeError < 5.0e-4);
}

TEST_CASE("CoupledDDAssembler: contact-referenced quasi-Fermi coordinates preserve physical state",
          "[newton][coupled][qf-reference]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    CoupledDDAssembler absolute(
        mesh, matdb, doping, constants::Vt_300, 1.0e-6, 1.0e-6);
    CoupledDDAssembler referenced(
        mesh, matdb, doping, constants::Vt_300, 1.0e-6, 1.0e-6);
    referenced.setQuasiFermiReferences(0.3, -0.2);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.05, 0.05);
    state.phin = VectorXd::LinSpaced(N, 0.299, 0.301);
    state.phip = VectorXd::LinSpaced(N, -0.201, -0.199);

    const VectorXd absoluteX = absolute.pack(state);
    const VectorXd referencedX = referenced.pack(state);
    const CoupledDDState roundTrip = referenced.unpack(referencedX);
    REQUIRE((roundTrip.psi - state.psi).norm() == Catch::Approx(0.0));
    REQUIRE((roundTrip.phin - state.phin).norm() == Catch::Approx(0.0).margin(1.0e-15));
    REQUIRE((roundTrip.phip - state.phip).norm() == Catch::Approx(0.0).margin(1.0e-15));

    CoupledDDBoundaryConditions bcs;
    bcs.phin[0] = state.phin(0);
    bcs.phip[0] = state.phip(0);
    const VectorXd absoluteResidual = absolute.residual(absoluteX, bcs);
    const VectorXd referencedResidual = referenced.residual(referencedX, bcs);
    REQUIRE((referencedResidual - absoluteResidual).norm()
            <= 1.0e-12 * std::max<Real>(1.0, absoluteResidual.norm()));
}


TEST_CASE("CoupledDDAssembler: analytic pinned rows suppress zero-rate recombination derivatives", "[newton][coupled]")
{
    DeviceMesh mesh = makeOxideMesh();
    MaterialDatabase matdb;
    Material zeroMobilitySemiconductor;
    zeroMobilitySemiconductor.name = "SiO2";
    zeroMobilitySemiconductor.eps_r = 11.7;
    zeroMobilitySemiconductor.ni = 1.0e16;
    zeroMobilitySemiconductor.mun = 0.0;
    zeroMobilitySemiconductor.mup = 0.0;
    matdb.addMaterial(zeroMobilitySemiconductor);

    DopingModel doping(mesh.numNodes());
    CoupledDDAssembler assembler(mesh, matdb, doping, constants::Vt_300, 1.0e-7, 1.0e-7);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::Zero(N);
    state.phin = VectorXd::LinSpaced(N, 0.05, 0.2);
    state.phip = state.phin;
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    const SparseMatrixd J = assembler.assembleJacobian(x, bcs);
    const Eigen::MatrixXd dense = Eigen::MatrixXd(J);

    for (int i = 0; i < N; ++i) {
        const int electronRow = N + i;
        const int holeRow = 2 * N + i;
        for (int col = 0; col < 3 * N; ++col) {
            REQUIRE(dense(electronRow, col) == Catch::Approx(col == electronRow ? 1.0 : 0.0));
            REQUIRE(dense(holeRow, col) == Catch::Approx(col == holeRow ? 1.0 : 0.0));
        }
    }
}

static Real triangleArea(const DeviceMesh& mesh, const Cell& cell)
{
    const Node& a = mesh.getNode(cell.node_ids[0]);
    const Node& b = mesh.getNode(cell.node_ids[1]);
    const Node& c = mesh.getNode(cell.node_ids[2]);
    return 0.5 * std::abs(
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

static Real nodeVolume(const DeviceMesh& mesh, Index node)
{
    Real volume = 0.0;
    for (const Cell& cell : mesh.cells()) {
        if (std::find(cell.node_ids.begin(), cell.node_ids.end(), node) != cell.node_ids.end())
            volume += triangleArea(mesh, cell) / static_cast<Real>(cell.node_ids.size());
    }
    return volume;
}

TEST_CASE("CoupledDDAssembler: TCAD SRH source scaling covers residual diagnostics and Jacobian",
          "[newton][coupled][scaling][source_units][srh]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    DDScalingSpec legacyScaling;
    legacyScaling.enabled = true;
    legacyScaling.V0 = 1.0;
    legacyScaling.C0 = 1.0;
    legacyScaling.mu0 = 1.0;
    legacyScaling.D0 = 1.0;
    legacyScaling.L0 = 1.0;
    legacyScaling.permittivityReference_F_per_m = constants::eps0 * 11.7;
    legacyScaling.unitSystem = PhysicalUnitSystem::legacySI();

    DDScalingSpec tcadScaling = legacyScaling;
    tcadScaling.unitSystem = PhysicalUnitSystem::tcadInternal();

    const RecombinationModelConfig noRecombination =
        recombinationModelConfig({"none"});
    const RecombinationModelConfig srh =
        recombinationModelConfig({"srh"}, 1.0e-7, 2.0e-7);

    const auto makeAssembler = [&](const RecombinationModelConfig& recombination,
                                   const DDScalingSpec& scaling) {
        return CoupledDDAssembler(
            mesh,
            matdb,
            doping,
            constants::Vt_300,
            mobilityModelConfig("constant"),
            recombination,
            BandgapNarrowingConfig{},
            ImpactIonizationModelConfig{},
            {},
            {},
            scaling);
    };

    CoupledDDAssembler baseline = makeAssembler(noRecombination, legacyScaling);
    CoupledDDAssembler legacy = makeAssembler(srh, legacyScaling);
    CoupledDDAssembler tcad = makeAssembler(srh, tcadScaling);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::Zero(N);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);
    const VectorXd x = legacy.pack(state);
    const CoupledDDBoundaryConditions bcs;

    const auto legacyTerms = legacy.carrierContinuityTermDiagnostics(x, bcs);
    const auto tcadTerms = tcad.carrierContinuityTermDiagnostics(x, bcs);
    const Index node = 4;
    REQUIRE(std::abs(legacyTerms[node].electronRecombination) > 0.0);
    REQUIRE(tcadTerms[node].electronRecombination ==
            Catch::Approx(legacyTerms[node].electronRecombination * 1.0e-4)
                .epsilon(1.0e-12));
    REQUIRE(tcadTerms[node].holeRecombination ==
            Catch::Approx(legacyTerms[node].holeRecombination * 1.0e-4)
                .epsilon(1.0e-12));

    const VectorXd legacySourceResidual =
        legacy.residual(x, bcs).segment(N, 2 * N);
    const VectorXd tcadSourceResidual =
        tcad.residual(x, bcs).segment(N, 2 * N);
    REQUIRE(tcadSourceResidual.norm() ==
            Catch::Approx(legacySourceResidual.norm() * 1.0e-4)
                .epsilon(1.0e-12));

    const Eigen::MatrixXd baselineJacobian =
        Eigen::MatrixXd(baseline.assembleJacobian(x, bcs));
    const Eigen::MatrixXd legacySourceJacobian =
        Eigen::MatrixXd(legacy.assembleJacobian(x, bcs)) - baselineJacobian;
    const Eigen::MatrixXd tcadSourceJacobian =
        Eigen::MatrixXd(tcad.assembleJacobian(x, bcs)) - baselineJacobian;
    REQUIRE(legacySourceJacobian.norm() > 0.0);
    REQUIRE(tcadSourceJacobian.norm() ==
            Catch::Approx(legacySourceJacobian.norm() * 1.0e-4)
                .epsilon(1.0e-5));
}

TEST_CASE("CoupledDDAssembler: TCAD impact source scaling covers residual diagnostics and Jacobian",
          "[newton][coupled][scaling][source_units][impact]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    DDScalingSpec legacyScaling;
    legacyScaling.enabled = true;
    legacyScaling.V0 = 1.0;
    legacyScaling.C0 = 1.0;
    legacyScaling.mu0 = 1.0;
    legacyScaling.D0 = 1.0;
    legacyScaling.L0 = 1.0;
    legacyScaling.permittivityReference_F_per_m = constants::eps0 * 11.7;
    legacyScaling.unitSystem = PhysicalUnitSystem::legacySI();

    DDScalingSpec tcadScaling = legacyScaling;
    tcadScaling.unitSystem = PhysicalUnitSystem::tcadInternal();

    MobilityModelConfig mobility = mobilityModelConfig("constant");
    mobility.highFieldDrivingForce = "quasi_fermi_gradient";

    ImpactIonizationModelConfig impact;
    impact.model = "selberherr";
    impact.drivingForce = "quasi_fermi_gradient";
    impact.generation = "current_density";
    impact.currentApproximation = "cell_reconstructed";
    impact.sourceVolumePolicy = "edge_box";
    impact.electronA = 1.0e8;
    impact.electronB = 1.0e-30;
    impact.holeA = 1.0e8;
    impact.holeB = 1.0e-30;

    const auto makeAssembler = [&](const ImpactIonizationModelConfig& impactConfig,
                                   const DDScalingSpec& scaling) {
        return CoupledDDAssembler(
            mesh,
            matdb,
            doping,
            constants::Vt_300,
            mobility,
            recombinationModelConfig({"none"}),
            BandgapNarrowingConfig{},
            impactConfig,
            {},
            {},
            scaling);
    };

    CoupledDDAssembler baseline =
        makeAssembler(ImpactIonizationModelConfig{}, legacyScaling);
    CoupledDDAssembler legacy = makeAssembler(impact, legacyScaling);
    CoupledDDAssembler tcad = makeAssembler(impact, tcadScaling);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, 0.0, -0.8);
    state.phin = VectorXd::LinSpaced(N, 0.0, -1.4);
    state.phip = VectorXd::LinSpaced(N, -0.2, -0.9);
    const VectorXd x = legacy.pack(state);
    const CoupledDDBoundaryConditions bcs;

    const auto legacyTerms = legacy.carrierContinuityTermDiagnostics(x, bcs);
    const auto tcadTerms = tcad.carrierContinuityTermDiagnostics(x, bcs);
    Real legacyImpactNorm = 0.0;
    Real tcadImpactNorm = 0.0;
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        legacyImpactNorm += std::abs(legacyTerms[node].electronImpact);
        tcadImpactNorm += std::abs(tcadTerms[node].electronImpact);
    }
    REQUIRE(legacyImpactNorm > 0.0);
    REQUIRE(tcadImpactNorm ==
            Catch::Approx(legacyImpactNorm * 1.0e-4).epsilon(1.0e-10));

    const VectorXd legacyResidual = legacy.residual(x, bcs);
    const VectorXd tcadResidual = tcad.residual(x, bcs);
    for (int node = 0; node < N; ++node) {
        REQUIRE(legacyResidual(N + node) - tcadResidual(N + node) ==
                Catch::Approx(
                    legacyTerms[static_cast<Index>(node)].electronImpact -
                    tcadTerms[static_cast<Index>(node)].electronImpact)
                    .epsilon(1.0e-8));
        REQUIRE(legacyResidual(2 * N + node) - tcadResidual(2 * N + node) ==
                Catch::Approx(
                    legacyTerms[static_cast<Index>(node)].holeImpact -
                    tcadTerms[static_cast<Index>(node)].holeImpact)
                    .epsilon(1.0e-8));
    }

    const Eigen::MatrixXd baselineJacobian =
        Eigen::MatrixXd(baseline.assembleJacobian(x, bcs));
    const Eigen::MatrixXd legacySourceJacobian =
        Eigen::MatrixXd(legacy.assembleJacobian(x, bcs)) - baselineJacobian;
    const Eigen::MatrixXd tcadSourceJacobian =
        Eigen::MatrixXd(tcad.assembleJacobian(x, bcs)) - baselineJacobian;
    REQUIRE(legacySourceJacobian.norm() > 0.0);
    REQUIRE(tcadSourceJacobian.norm() ==
            Catch::Approx(legacySourceJacobian.norm() * 1.0e-4)
                .epsilon(1.0e-4));
}

TEST_CASE("CoupledDDAssembler: carrier diagonal floor anchors depleted minority row", "[newton][coupled]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    RecombinationModelConfig recombinationConfig = recombinationModelConfig({"srh"}, 1.0e-7, 1.0e-7);
    CarrierDiagonalFloorRegularizationConfig floorConfig;
    floorConfig.enabled = true;

    CoupledDDAssembler baseline(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        MobilityModelConfig{},
        recombinationConfig);
    CoupledDDAssembler regularized(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        MobilityModelConfig{},
        recombinationConfig,
        BandgapNarrowingConfig{},
        ImpactIonizationModelConfig{},
        {},
        {},
        {},
        floorConfig);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::Zero(N);
    state.phin = VectorXd::Zero(N);
    state.phip = VectorXd::Zero(N);
    const Index depletedNode = 4;
    state.phin(static_cast<int>(depletedNode)) = 600.0 * constants::Vt_300;
    const VectorXd x = baseline.pack(state);

    CoupledDDBoundaryConditions bcs;
    const SparseMatrixd J0 = baseline.assembleJacobian(x, bcs);
    const SparseMatrixd J1 = regularized.assembleJacobian(x, bcs);

    const int electronRow = N + static_cast<int>(depletedNode);
    const Real ni = baseline.intrinsicDensity().at(static_cast<std::size_t>(depletedNode));
    const Real expectedFloor =
        nodeVolume(mesh, depletedNode) * ni / ((recombinationConfig.taun + recombinationConfig.taup) * constants::Vt_300);

    REQUIRE(std::abs(J0.coeff(electronRow, electronRow)) < expectedFloor);
    REQUIRE(std::abs(J1.coeff(electronRow, electronRow)) ==
            Catch::Approx(expectedFloor).epsilon(1.0e-12));

    const int normalElectronRow = N;
    REQUIRE(J1.coeff(normalElectronRow, normalElectronRow) ==
            Catch::Approx(J0.coeff(normalElectronRow, normalElectronRow)));
}
TEST_CASE("CoupledDDAssembler: analytic Jacobian matches finite differences on small mesh", "[newton][coupled]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    CoupledDDAssembler assembler(mesh, matdb, doping, constants::Vt_300, 1.0e-7, 1.0e-7);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.04, 0.05);
    state.phin = VectorXd::LinSpaced(N, 0.01, -0.015);
    state.phip = VectorXd::LinSpaced(N, -0.02, 0.012);
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    bcs.psi[0] = state.psi(0);
    bcs.phin[0] = state.phin(0);
    bcs.phip[0] = state.phip(0);
    bcs.psi[2] = state.psi(2);
    bcs.phin[2] = state.phin(2);
    bcs.phip[2] = state.phip(2);

    const SparseMatrixd Ja = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd Jfd = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd diff = Eigen::MatrixXd(Ja - Jfd);
    const Eigen::MatrixXd ref = Eigen::MatrixXd(Jfd);
    const Real rel = diff.norm() / std::max<Real>(1.0, ref.norm());

    REQUIRE(rel < 5.0e-5);
}

TEST_CASE("CoupledDDAssembler: Sentaurus-default Fermi SRH source has analytic DG Jacobian",
          "[newton][coupled][density_gradient][srh][fermi][jacobian]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    CarrierStatisticsConfig statistics{"fermi_dirac"};
    DensityGradientQuantumPotentialConfig quantum;
    quantum.enabled = true;
    quantum.couplingMode = "frozen";

    RecombinationModelConfig srh =
        recombinationModelConfig({"srh"}, 1.1e-7, 2.3e-7);
    srh.srhDopingDependence.densityCoupling = "sentaurus_default";
    RecombinationModelConfig none = recombinationModelConfig({"none"});
    none.srhDopingDependence.densityCoupling = "sentaurus_default";
    CoupledDDAssembler withSrh(
        mesh, matdb, doping, constants::Vt_300, {}, srh, {}, {}, {}, {}, {}, {},
        statistics, quantum);
    CoupledDDAssembler withoutSrh(
        mesh, matdb, doping, constants::Vt_300, {}, none, {}, {}, {}, {}, {}, {},
        statistics, quantum);

    const int N = static_cast<int>(mesh.numNodes());
    VectorXd quantumPotential = VectorXd::LinSpaced(N, 0.018, 0.047);
    withSrh.setElectronQuantumPotential(quantumPotential);
    withoutSrh.setElectronQuantumPotential(quantumPotential);
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.035, 0.052);
    state.phin = VectorXd::LinSpaced(N, 0.012, -0.019);
    state.phip = VectorXd::LinSpaced(N, -0.016, 0.014);
    const VectorXd x = withSrh.pack(state);
    const CoupledDDBoundaryConditions bcs;

    const Eigen::MatrixXd analyticSource = Eigen::MatrixXd(
        withSrh.assembleJacobian(x, bcs) -
        withoutSrh.assembleJacobian(x, bcs));
    const Eigen::MatrixXd finiteDifferenceSource = Eigen::MatrixXd(
        withSrh.finiteDifferenceJacobian(x, bcs, 1.0e-8) -
        withoutSrh.finiteDifferenceJacobian(x, bcs, 1.0e-8));
    const Real relativeError =
        (analyticSource - finiteDifferenceSource).norm() /
        std::max<Real>(1.0e-30, finiteDifferenceSource.norm());
    REQUIRE(finiteDifferenceSource.norm() > 0.0);
    REQUIRE(relativeError < 2.0e-4);

    RecombinationModelConfig quantumSrh = srh;
    quantumSrh.srhDopingDependence.densityCoupling = "quantum";
    CoupledDDAssembler quantumCoupled(
        mesh, matdb, doping, constants::Vt_300, {}, quantumSrh, {}, {}, {}, {},
        {}, {}, statistics, quantum);
    quantumCoupled.setElectronQuantumPotential(quantumPotential);
    REQUIRE((withSrh.residual(x, bcs) - quantumCoupled.residual(x, bcs)).norm() >
            0.0);

    const VectorXd zeroQuantum = VectorXd::Zero(N);
    withSrh.setElectronQuantumPotential(zeroQuantum);
    quantumCoupled.setElectronQuantumPotential(zeroQuantum);
    REQUIRE((withSrh.residual(x, bcs) - quantumCoupled.residual(x, bcs)).norm() ==
            Catch::Approx(0.0).margin(1.0e-18));
}

TEST_CASE("CoupledDDAssembler: classical Fermi Auger source has analytic Jacobian",
          "[newton][coupled][auger][fermi][jacobian]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const CarrierStatisticsConfig statistics{"fermi_dirac"};

    RecombinationModelConfig classical = recombinationModelConfig({"auger"});
    classical.augerExcessProduct = "classical_np";
    // Amplify the source above transport subtraction round-off in this
    // source-only analytic-versus-finite-difference comparison.
    classical.augerCn = 2.9e-28;
    classical.augerCp = 1.028e-28;
    RecombinationModelConfig generalized = classical;
    generalized.augerExcessProduct = "generalized_fermi";
    const RecombinationModelConfig none = recombinationModelConfig({"none"});

    CoupledDDAssembler withClassical(
        mesh, matdb, doping, constants::Vt_300, {}, classical, {}, {}, {}, {},
        {}, {}, statistics);
    CoupledDDAssembler withGeneralized(
        mesh, matdb, doping, constants::Vt_300, {}, generalized, {}, {}, {}, {},
        {}, {}, statistics);
    CoupledDDAssembler withoutAuger(
        mesh, matdb, doping, constants::Vt_300, {}, none, {}, {}, {}, {}, {}, {},
        statistics);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.06, 0.08);
    state.phin = VectorXd::LinSpaced(N, 0.03, -0.025);
    state.phip = VectorXd::LinSpaced(N, -0.04, 0.02);
    const VectorXd x = withClassical.pack(state);
    const CoupledDDBoundaryConditions bcs;

    const Eigen::MatrixXd analyticSource = Eigen::MatrixXd(
        withClassical.assembleJacobian(x, bcs) -
        withoutAuger.assembleJacobian(x, bcs));
    const Eigen::MatrixXd finiteDifferenceSource = Eigen::MatrixXd(
        withClassical.finiteDifferenceJacobian(x, bcs, 1.0e-8) -
        withoutAuger.finiteDifferenceJacobian(x, bcs, 1.0e-8));
    const Real relativeError =
        (analyticSource - finiteDifferenceSource).norm() /
        std::max<Real>(1.0e-30, finiteDifferenceSource.norm());

    REQUIRE(finiteDifferenceSource.norm() > 0.0);
    REQUIRE(relativeError < 2.0e-4);
    REQUIRE((withClassical.residual(x, bcs) -
             withGeneralized.residual(x, bcs)).norm() > 0.0);
}

TEST_CASE("CoupledDDAssembler: analytic Jacobian matches finite differences with varying intrinsic density",
          "[newton][coupled][bgn]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e24, 0.0},
        {"p_region", 0.0, 1.0e21},
    });

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        MobilityModelConfig{},
        recombinationModelConfig({"none"}),
        bandgapNarrowingConfig("slotboom"));

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.025);
    state.phin = VectorXd::LinSpaced(N, 0.006, -0.008);
    state.phip = VectorXd::LinSpaced(N, -0.007, 0.005);
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    bcs.psi[0] = state.psi(0);
    bcs.phin[0] = state.phin(0);
    bcs.phip[0] = state.phip(0);
    bcs.psi[2] = state.psi(2);
    bcs.phin[2] = state.phin(2);
    bcs.phip[2] = state.phip(2);

    const SparseMatrixd Ja = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd Jfd = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd diff = Eigen::MatrixXd(Ja - Jfd);
    const Eigen::MatrixXd ref = Eigen::MatrixXd(Jfd);
    const Real rel = diff.norm() / std::max<Real>(1.0, ref.norm());

    REQUIRE(rel < 1.0e-4);
}

TEST_CASE("CoupledDDAssembler: analytic Jacobian matches finite differences at BV absolute potential scale",
          "[newton][coupled][bgn][bv]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e24, 0.0},
        {"p_region", 0.0, 1.0e21},
    });

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        MobilityModelConfig{},
        recombinationModelConfig({"none"}),
        bandgapNarrowingConfig("slotboom"));

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi.resize(N);
    state.phin.resize(N);
    state.phip.resize(N);
    state.psi << -13.20365, -13.20340, -13.20315, -13.20390, -13.20355;
    state.phin << -12.79890, -12.79930, -12.79970, -12.80010, -12.79910;
    state.phip << -12.80020, -12.79980, -12.79940, -12.80000, -12.79960;
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    bcs.psi[0] = state.psi(0);
    bcs.phin[0] = state.phin(0);
    bcs.phip[0] = state.phip(0);
    bcs.psi[1] = state.psi(1);
    bcs.phin[1] = state.phin(1);
    bcs.phip[1] = state.phip(1);
    bcs.psi[2] = state.psi(2);
    bcs.phin[2] = state.phin(2);
    bcs.phip[2] = state.phip(2);
    bcs.psi[3] = state.psi(3);
    bcs.phin[3] = state.phin(3);
    bcs.phip[3] = state.phip(3);

    const SparseMatrixd Ja = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd Jfd = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-8);
    const Eigen::MatrixXd diff = Eigen::MatrixXd(Ja - Jfd);
    const Eigen::MatrixXd ref = Eigen::MatrixXd(Jfd);
    const Real rel = diff.norm() / std::max<Real>(1.0, ref.norm());

    REQUIRE(rel < 1.0e-5);
}

TEST_CASE("CoupledDDAssembler: transport Jacobian captures quasi-Fermi high-field mobility",
          "[newton][coupled][mobility][field]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e23, 0.0},
        {"p_region", 0.0, 1.0e23},
    });

    MobilityModelConfig mobility = mobilityModelConfig("masetti_field");
    mobility.highFieldDrivingForce = "quasi_fermi_gradient";

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        mobility,
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{});

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.025);
    state.phin = VectorXd::LinSpaced(N, 0.7, -0.7);
    state.phip = VectorXd::LinSpaced(N, -0.65, 0.65);
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    bcs.psi[0] = state.psi(0);
    bcs.phin[0] = state.phin(0);
    bcs.phip[0] = state.phip(0);
    bcs.psi[2] = state.psi(2);
    bcs.phin[2] = state.phin(2);
    bcs.phip[2] = state.phip(2);

    const SparseMatrixd Ja = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd Jfd = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd diff = Eigen::MatrixXd(Ja - Jfd);
    const Eigen::MatrixXd ref = Eigen::MatrixXd(Jfd);
    const Real rel = diff.norm() / std::max<Real>(1.0, ref.norm());

    REQUIRE(rel < 1.0e-4);

    mobility.jacobianFieldDerivatives = false;
    CoupledDDAssembler frozenJacobianAssembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        mobility,
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{});
    const SparseMatrixd frozenJa = frozenJacobianAssembler.assembleJacobian(x, bcs);
    const Real frozenRel = Eigen::MatrixXd(frozenJa - Jfd).norm() /
        std::max<Real>(1.0, ref.norm());

    REQUIRE(frozenRel > rel * 10.0);
}
TEST_CASE("CoupledDDAssembler: Slotboom BGN uses total impurity at compensated nodes",
          "[newton][coupled][bgn][doping]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());
    doping.setNodeDoping(4, 1.0e23, 1.0e23);

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        MobilityModelConfig{},
        recombinationModelConfig({"none"}),
        bandgapNarrowingConfig("slotboom"));

    const Material& si = matdb.getMaterial("Si");
    const SlotboomBandgapNarrowing bgn(bandgapNarrowingConfig("slotboom"));
    const Real expected = effectiveIntrinsicDensity(
        si.ni,
        constants::Vt_300,
        bgn.deltaEg(doping.totalImpurity(4), 0.0, 0.0));

    REQUIRE(doping.netDoping(4) == Catch::Approx(0.0));
    REQUIRE(expected > si.ni);
    REQUIRE(assembler.intrinsicDensity().at(4) == Catch::Approx(expected));
}

TEST_CASE("CoupledDDAssembler: scaled state residual and Jacobian are consistent",
          "[newton][coupled][scaling]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    UnitScalingSystem::AutoInputs inputs =
        UnitScalingSystem::autoInputsFrom(mesh, doping, matdb, 1.0e16);
    const UnitScalingSystem sc = UnitScalingSystem::fromInputs(
        300.0, constants::eps0 * 11.7, inputs);

    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = sc.V0();
    scaling.C0 = sc.C0();
    scaling.mu0 = sc.mu0();
    scaling.D0 = sc.D0();
    scaling.L0 = sc.L0();
    scaling.permittivityReference_F_per_m = constants::eps0 * 11.7;

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        MobilityModelConfig{},
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        ImpactIonizationModelConfig{},
        {},
        {},
        scaling);

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.04, 0.05) / scaling.V0;
    state.phin = VectorXd::LinSpaced(N, 0.01, -0.015) / scaling.V0;
    state.phip = VectorXd::LinSpaced(N, -0.02, 0.012) / scaling.V0;
    const VectorXd x = assembler.pack(state);

    CoupledDDBoundaryConditions bcs;
    bcs.psi[0] = state.psi(0);
    bcs.phin[0] = state.phin(0);
    bcs.phip[0] = state.phip(0);
    bcs.psi[2] = state.psi(2);
    bcs.phin[2] = state.phin(2);
    bcs.phip[2] = state.phip(2);

    const VectorXd r = assembler.residual(x, bcs);
    REQUIRE(r.allFinite());
    REQUIRE(r.norm() < 1.0e8);

    const SparseMatrixd Ja = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd Jfd = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd diff = Eigen::MatrixXd(Ja - Jfd);
    const Eigen::MatrixXd ref = Eigen::MatrixXd(Jfd);
    const Real rel = diff.norm() / std::max<Real>(1.0, ref.norm());
    REQUIRE(rel < 1.0e-8);
}

TEST_CASE("NewtonSolver: evaluateStep reports one physical Newton correction", "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonStepEvaluation step = solver.evaluateStep(state);

    REQUIRE(step.residual.raw.size() == 3 * N);
    REQUIRE(step.deltaPsi.size() == N);
    REQUIRE(step.deltaPhin.size() == N);
    REQUIRE(step.deltaPhip.size() == N);
    REQUIRE(step.trialSolution.psi.size() == N);
    REQUIRE(step.stepNorm > 0.0);
    REQUIRE(step.rawStepNorm > 0.0);
    REQUIRE(step.trialSolution.psi(0) ==
            Catch::Approx(state.psi(0) + step.deltaPsi(0)));
    REQUIRE(step.trialSolution.phin(1) ==
            Catch::Approx(state.phin(1) + step.deltaPhin(1)));
    REQUIRE(step.trialSolution.phip(2) ==
            Catch::Approx(state.phip(2) + step.deltaPhip(2)));
    REQUIRE(step.trialResidual.blockNorms.combined < step.residual.blockNorms.combined);
}

TEST_CASE("NewtonSolver: feedback substitutions share one Jacobian and preserve boundary rows",
          "[newton][diagnostics][feedback-substitution]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    NewtonConfig cfg = newtonConfig();
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.warmStart = true;

    NewtonSolver solver(mesh, matdb, doping, zeroBias(), cfg);
    const NewtonResult equilibrium = solver.solve();
    REQUIRE(equilibrium.converged);

    DDSolution replacement = equilibrium.solution;
    replacement.phin(4) += 0.01;
    replacement.phip(4) -= 0.015;
    replacement.n(4) *= 1.5;
    replacement.p(4) *= 0.75;

    const auto variants =
        solver.evaluateFeedbackSubstitutions(equilibrium.solution, replacement);
    REQUIRE(variants.size() == 8);
    REQUIRE(variants[0].variant == "baseline");
    REQUIRE(variants[1].variant == "electron_density_only");
    REQUIRE(variants[2].variant == "hole_density_only");
    REQUIRE(variants[3].variant == "density_only");
    REQUIRE(variants[4].variant == "electron_qfp_only");
    REQUIRE(variants[5].variant == "hole_qfp_only");
    REQUIRE(variants[6].variant == "qfp_only");
    REQUIRE(variants[7].variant == "density_qfp");
    REQUIRE_FALSE(variants[0].replacesDensity);
    REQUIRE_FALSE(variants[0].replacesQuasiFermi);
    REQUIRE(variants[1].replacesDensity);
    REQUIRE_FALSE(variants[1].replacesQuasiFermi);
    REQUIRE(variants[2].replacesDensity);
    REQUIRE_FALSE(variants[2].replacesQuasiFermi);
    REQUIRE(variants[4].replacesQuasiFermi);
    REQUIRE(variants[5].replacesQuasiFermi);
    REQUIRE(variants[7].replacesDensity);
    REQUIRE(variants[7].replacesQuasiFermi);

    const int N = static_cast<int>(mesh.numNodes());
    const std::vector<int> contactNodes = {0, 1, 2, 3};
    for (const auto& variant : variants) {
        REQUIRE(variant.residual.raw.size() == 3 * N);
        REQUIRE(variant.desiredResidual.size() == 3 * N);
        REQUIRE(variant.carrierTerms.size() == static_cast<std::size_t>(N));
        REQUIRE(variant.deltaPsi.allFinite());
        REQUIRE(variant.deltaPhin.allFinite());
        REQUIRE(variant.deltaPhip.allFinite());
        REQUIRE(variant.carrierOnlyDeltaPhin.allFinite());
        REQUIRE(variant.carrierOnlyDeltaPhip.allFinite());
        REQUIRE((variant.desiredResidual - variants[0].desiredResidual).norm() ==
                Catch::Approx(0.0).margin(0.0));
        for (const int node : contactNodes) {
            REQUIRE(variant.residual.raw(node) == variants[0].residual.raw(node));
            REQUIRE(variant.residual.raw(N + node) ==
                    variants[0].residual.raw(N + node));
            REQUIRE(variant.residual.raw(2 * N + node) ==
                    variants[0].residual.raw(2 * N + node));
        }
        const auto& center = variant.carrierTerms[4];
        const Real electronSum = center.electronFlux
            + center.electronRecombination
            + center.electronImpact
            + center.electronGauge
            + center.electronBoundary;
        const Real holeSum = center.holeFlux
            + center.holeRecombination
            + center.holeImpact
            + center.holeGauge
            + center.holeBoundary;
        REQUIRE(electronSum ==
                Catch::Approx(variant.residual.raw(N + 4)).margin(1.0e-12));
        REQUIRE(holeSum ==
                Catch::Approx(variant.residual.raw(2 * N + 4)).margin(1.0e-12));
    }
    REQUIRE((variants[1].residual.raw - variants[0].residual.raw).norm() > 0.0);
    REQUIRE((variants[4].residual.raw - variants[0].residual.raw).norm() > 0.0);
}

TEST_CASE("NewtonSolver: Poisson-QFP cross-block decomposition closes to full Newton step",
          "[newton][diagnostics][poisson-qfp-cross-block]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    NewtonConfig cfg = newtonConfig();
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.warmStart = true;

    NewtonSolver solver(mesh, matdb, doping, zeroBias(), cfg);
    const NewtonResult equilibrium = solver.solve();
    REQUIRE(equilibrium.converged);

    DDSolution replacement = equilibrium.solution;
    replacement.phin(4) += 0.01;
    replacement.phip(4) -= 0.015;
    const auto evaluation =
        solver.evaluatePoissonQfpCrossBlockDecomposition(
            equilibrium.solution, replacement);
    const int N = static_cast<int>(mesh.numNodes());

    REQUIRE(evaluation.residual.raw.size() == 3 * N);
    REQUIRE(evaluation.jacobianPsiPsi.rows() == N);
    REQUIRE(evaluation.jacobianPsiPsi.cols() == N);
    REQUIRE(evaluation.jacobianPsiQfp.rows() == N);
    REQUIRE(evaluation.jacobianPsiQfp.cols() == 2 * N);
    REQUIRE(evaluation.jacobianQfpPsi.rows() == 2 * N);
    REQUIRE(evaluation.jacobianQfpPsi.cols() == N);
    REQUIRE(evaluation.jacobianQfpQfp.rows() == 2 * N);
    REQUIRE(evaluation.jacobianQfpQfp.cols() == 2 * N);
    REQUIRE(evaluation.effectiveSchurLoop.rows() == 2 * N);
    REQUIRE(evaluation.effectiveSchurLoop.cols() == 2 * N);
    REQUIRE(evaluation.loopComponents.size() == 3);
    REQUIRE(evaluation.jacobianPsiQfpNorm > 0.0);
    REQUIRE(evaluation.jacobianQfpPsiNorm > 0.0);

    REQUIRE(
        (evaluation.noPsiQfpDeltaPsi - evaluation.independentDeltaPsi).norm()
        == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(
        (evaluation.noQfpPsiDeltaPhin - evaluation.independentDeltaPhin).norm()
        == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(
        (evaluation.noQfpPsiDeltaPhip - evaluation.independentDeltaPhip).norm()
        == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(
        (evaluation.schurDeltaPsi - evaluation.fullRawDeltaPsi).norm()
        == Catch::Approx(0.0).margin(1.0e-10));
    REQUIRE(
        (evaluation.schurDeltaPhin - evaluation.fullRawDeltaPhin).norm()
        == Catch::Approx(0.0).margin(1.0e-10));
    REQUIRE(
        (evaluation.schurDeltaPhip - evaluation.fullRawDeltaPhip).norm()
        == Catch::Approx(0.0).margin(1.0e-10));
    REQUIRE(evaluation.fullLinearClosureNorm < 1.0e-8);
    REQUIRE(evaluation.schurRelativeClosure < 1.0e-10);
    REQUIRE(evaluation.loopComponentClosureNorm < 1.0e-8);
    REQUIRE(evaluation.psiQfpProduct.size() == N);
    REQUIRE(evaluation.qfpPsiProduct.size() == 2 * N);
    REQUIRE(evaluation.qfpFiniteDifferenceDirectionPhin.size() == N);
    REQUIRE(evaluation.psiFiniteDifferenceDirection.size() == N);
    REQUIRE(
        evaluation.psiQfpDirectionalDerivativeRelativeError < 1.0e-4);
    INFO("C analytic norm: "
         << evaluation.analyticQfpPsiDirectionalDerivative.norm());
    INFO("C finite-difference norm: "
         << evaluation.finiteDifferenceQfpPsiDirectionalDerivative.norm());
    INFO("C relative error: "
         << evaluation.qfpPsiDirectionalDerivativeRelativeError);
    REQUIRE(
        evaluation.qfpPsiDirectionalDerivativeRelativeError < 1.0e-4);
    REQUIRE(evaluation.jacobianPsiPsiCondition.numericalRank > 0);
    REQUIRE(evaluation.jacobianQfpQfpCondition.numericalRank > 0);
    REQUIRE(evaluation.schurCondition.numericalRank > 0);
}

TEST_CASE("NewtonSolver: evaluateDirectionalDerivative compares analytic and finite-difference Jv",
          "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    DDSolution perturbation;
    perturbation.psi = VectorXd::Zero(N);
    perturbation.phin = VectorXd::Zero(N);
    perturbation.phip = VectorXd::Zero(N);
    perturbation.psi(4) = 0.5e-6;
    perturbation.phin(4) = -0.5e-6;

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonDirectionalDerivativeEvaluation jvp =
        solver.evaluateDirectionalDerivative(state, perturbation);

    REQUIRE(jvp.residual.raw.size() == 3 * N);
    REQUIRE(jvp.analyticJv.size() == 3 * N);
    REQUIRE(jvp.finiteDifferenceJv.size() == 3 * N);
    REQUIRE(jvp.perturbationPsi.size() == N);
    REQUIRE(jvp.perturbationPhin.size() == N);
    REQUIRE(jvp.perturbationPhip.size() == N);
    REQUIRE(jvp.analyticJv.norm() > 0.0);
    REQUIRE(jvp.finiteDifferenceJv.norm() > 0.0);
    REQUIRE(jvp.relativeError < 1.0e-6);
    REQUIRE(jvp.perturbationPsi(4) == Catch::Approx(0.5e-6));
    REQUIRE(jvp.perturbationPhin(4) == Catch::Approx(-0.5e-6));

    cfg.quasiFermiReference = "contact_majority";
    NewtonSolver referencedSolver(mesh, matdb, doping, biases, cfg);
    const NewtonDirectionalDerivativeEvaluation referencedJvp =
        referencedSolver.evaluateDirectionalDerivative(state, perturbation);
    REQUIRE(referencedJvp.perturbationNorm ==
            Catch::Approx(jvp.perturbationNorm).epsilon(1.0e-12));
    REQUIRE((referencedJvp.analyticJv - jvp.analyticJv).norm() <
            1.0e-11 * std::max<Real>(1.0, jvp.analyticJv.norm()));
    REQUIRE((referencedJvp.finiteDifferenceJv - jvp.finiteDifferenceJv).norm() <
            1.0e-8 * std::max<Real>(1.0, jvp.finiteDifferenceJv.norm()));
}

TEST_CASE("NewtonSolver: constant-field contact HFS Jacobian uses limited mobility",
          "[newton][diagnostics][mobility][contact]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.mobility.model = "constant_field";
    cfg.mobility.highFieldDrivingForce = "quasi_fermi_gradient";
    cfg.mobility.contactElectricFieldFallback = true;
    cfg.mobility.electronField.saturationVelocity = 1.0e4;
    cfg.mobility.holeField.saturationVelocity = 1.0e4;

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution state;
    state.psi = VectorXd::LinSpaced(N, -0.2, 0.2);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    DDSolution perturbation;
    perturbation.psi = VectorXd::Zero(N);
    perturbation.phin = VectorXd::Zero(N);
    perturbation.phip = VectorXd::Zero(N);
    perturbation.phin(4) = 1.0e-7;

    NewtonSolver solver(
        mesh, matdb, doping,
        {{"anode", -0.1}, {"cathode", 0.0}}, cfg);
    const NewtonDirectionalDerivativeEvaluation jvp =
        solver.evaluateDirectionalDerivative(state, perturbation);

    INFO("analytic norm: " << jvp.analyticNorm);
    INFO("finite-difference norm: " << jvp.finiteDifferenceNorm);
    INFO("relative error: " << jvp.relativeError);
    REQUIRE(jvp.finiteDifferenceNorm > 0.0);
    REQUIRE(jvp.relativeError < 1.0e-6);
}

TEST_CASE("NewtonSolver: evaluateJacobianBlockAudit reports finite block rows",
          "[newton][diagnostics][coupled]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"srh"};
    cfg.warmStart = true;
    cfg.impactIonization.model = "van_overstraeten";
    cfg.impactIonization.drivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.generation = "current_density";
    cfg.impactIonization.currentApproximation = "density_gradient";

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::LinSpaced(N, -0.015, 0.015);
    state.phip = VectorXd::LinSpaced(N, 0.012, -0.012);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const auto rows = solver.evaluateJacobianBlockAudit(state, 1.0e-7);
    const auto hasBlock = [&](const std::string& name) {
        return std::any_of(rows.begin(), rows.end(), [&](const auto& row) {
            return row.block == name &&
                   std::isfinite(row.analyticNorm) &&
                   std::isfinite(row.fdNorm) &&
                   std::isfinite(row.diffNorm) &&
                   std::isfinite(row.relDiff) &&
                   std::isfinite(row.analyticPsiColumnNorm) &&
                   std::isfinite(row.fdPsiColumnNorm) &&
                   std::isfinite(row.diffPsiColumnNorm) &&
                   std::isfinite(row.relPsiColumnDiff) &&
                   std::isfinite(row.analyticPhinColumnNorm) &&
                   std::isfinite(row.fdPhinColumnNorm) &&
                   std::isfinite(row.diffPhinColumnNorm) &&
                   std::isfinite(row.relPhinColumnDiff) &&
                   std::isfinite(row.analyticPhipColumnNorm) &&
                   std::isfinite(row.fdPhipColumnNorm) &&
                   std::isfinite(row.diffPhipColumnNorm) &&
                   std::isfinite(row.relPhipColumnDiff);
        });
    };

    REQUIRE(hasBlock("poisson"));
    REQUIRE(hasBlock("transport"));
    REQUIRE(hasBlock("srh_auger"));
    REQUIRE(hasBlock("sg_avalanche"));
    REQUIRE(hasBlock("dirichlet_or_gauge"));
}

TEST_CASE("NewtonSolver: evaluateJacobianBlockAudit can restrict expensive block rows",
          "[newton][diagnostics][coupled]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"srh"};
    cfg.warmStart = true;
    cfg.impactIonization.model = "van_overstraeten";
    cfg.impactIonization.drivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.generation = "current_density";
    cfg.impactIonization.currentApproximation = "density_gradient";

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::LinSpaced(N, -0.015, 0.015);
    state.phip = VectorXd::LinSpaced(N, 0.012, -0.012);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const auto rows = solver.evaluateJacobianBlockAudit(
        state, 1.0e-7, std::vector<std::string>{"sg_avalanche"});

    REQUIRE(rows.size() == 1);
    REQUIRE(rows.front().block == "sg_avalanche");
    REQUIRE(std::isfinite(rows.front().analyticNorm));
    REQUIRE(std::isfinite(rows.front().fdNorm));
    REQUIRE(std::isfinite(rows.front().diffNorm));
    REQUIRE(std::isfinite(rows.front().relDiff));

    cfg.recombination = {"auger"};
    cfg.augerExcessProduct = "classical_np";
    cfg.carrierStatistics.model = "fermi_dirac";
    NewtonSolver sourceSolver(mesh, matdb, doping, biases, cfg);
    const auto sourceRows = sourceSolver.evaluateJacobianBlockAudit(
        state,
        1.0e-7,
        std::vector<std::string>{"srh_auger"},
        "double_symmetric",
        std::vector<Index>{1});
    REQUIRE(sourceRows.size() == 1);
    REQUIRE(sourceRows.front().block == "srh_auger");
    REQUIRE(sourceRows.front().relDiff < 1.0e-5);
}


TEST_CASE("NewtonSolver: element-edge avalanche block audit is independently source-only",
          "[newton][diagnostics][impact][element_edge]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -1.0},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.mobility.model = "constant";
    cfg.mobility.highFieldDrivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.model = "selberherr";
    cfg.impactIonization.drivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.generation = "current_density";
    cfg.impactIonization.currentApproximation = "element_edge_sg_gss_laux";
    cfg.impactIonization.sourceMappingMode = "element_vertex_box_measure";
    cfg.impactIonization.electronA = 1.0;
    cfg.impactIonization.electronB = 1.0e-30;
    cfg.impactIonization.holeA = 1.0;
    cfg.impactIonization.holeB = 1.0e-30;

    DDSolution state;
    state.psi = (VectorXd(5) << 0.05, -0.30, -0.75, -0.10, -0.42).finished();
    state.phin = (VectorXd(5) << 0.02, -0.65, -1.10, -0.08, -0.27).finished();
    state.phip = (VectorXd(5) << -0.15, -0.35, -0.82, -0.60, -0.24).finished();

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const auto sourceRows = solver.evaluateJacobianBlockAudit(
        state, 1.0e-7, std::vector<std::string>{"sg_avalanche"});
    const auto branchResolvedRows = solver.evaluateJacobianBlockAudit(
        state,
        1.0e-15,
        std::vector<std::string>{"sg_avalanche"},
        "multiprecision_branch_resolved");
    const auto transportRows = solver.evaluateJacobianBlockAudit(
        state, 1.0e-7, std::vector<std::string>{"transport"});

    REQUIRE(sourceRows.size() == 1);
    REQUIRE(branchResolvedRows.size() == 1);
    REQUIRE(transportRows.size() == 1);
    REQUIRE(sourceRows.front().analyticNorm > 0.0);
    REQUIRE(sourceRows.front().fdNorm > 0.0);
    REQUIRE(sourceRows.front().relDiff <= 1.0e-8);
    REQUIRE(sourceRows.front().analyticElectronPhinNorm > 0.0);
    REQUIRE(sourceRows.front().analyticElectronPhipNorm > 0.0);
    REQUIRE(sourceRows.front().analyticHolePhinNorm > 0.0);
    REQUIRE(sourceRows.front().analyticHolePhipNorm > 0.0);
    REQUIRE(sourceRows.front().relElectronPhinDiff <= 5.0e-5);
    REQUIRE(sourceRows.front().relElectronPhipDiff <= 5.0e-5);
    REQUIRE(sourceRows.front().relHolePhinDiff <= 5.0e-5);
    REQUIRE(sourceRows.front().relHolePhipDiff <= 5.0e-5);
    REQUIRE(branchResolvedRows.front().relDiff <= 1.0e-8);
    REQUIRE_FALSE(sourceRows.front().configurationFingerprint.empty());
    REQUIRE_FALSE(sourceRows.front().activeBranchFingerprint.empty());
    REQUIRE(
        branchResolvedRows.front().configurationFingerprint ==
        sourceRows.front().configurationFingerprint);
    REQUIRE(
        branchResolvedRows.front().activeBranchFingerprint ==
        sourceRows.front().activeBranchFingerprint);
    REQUIRE(sourceRows.front().analyticNorm != Catch::Approx(transportRows.front().analyticNorm));
    REQUIRE(sourceRows.front().fdNorm != Catch::Approx(transportRows.front().fdNorm));
}
TEST_CASE("NewtonSolver: cell-reconstructed SG avalanche Jacobian matches midpoint residual",
          "[newton][diagnostics][impact]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -1.0},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.mobility.model = "constant";
    cfg.mobility.highFieldDrivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.model = "selberherr";
    cfg.impactIonization.drivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.generation = "current_density";
    cfg.impactIonization.currentApproximation = "cell_reconstructed";
    cfg.impactIonization.sourceVolumePolicy = "edge_box";
    cfg.impactIonization.electronA = 1.0;
    cfg.impactIonization.electronB = 1.0e-30;
    cfg.impactIonization.holeA = 1.0;
    cfg.impactIonization.holeB = 1.0e-30;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, 0.0, -0.8);
    state.phin = VectorXd::LinSpaced(N, 0.0, -1.4);
    state.phip = VectorXd::LinSpaced(N, -0.2, -0.9);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const auto rows = solver.evaluateJacobianBlockAudit(
        state, 1.0e-7, std::vector<std::string>{"sg_avalanche"});

    REQUIRE(rows.size() == 1);
    REQUIRE(rows.front().block == "sg_avalanche");
    REQUIRE(rows.front().fdNorm > 0.0);
    REQUIRE(rows.front().analyticNorm > 0.0);
    REQUIRE(rows.front().relDiff < 5.0e-3);
}

TEST_CASE("NewtonSolver: triangle GSS local-AD Jacobian matches the nine-column cell residual",
          "[newton][diagnostics][impact][triangle_gss]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -1.0},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.mobility.model = "constant";
    cfg.mobility.highFieldDrivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.model = "selberherr";
    cfg.impactIonization.drivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.generation = "current_density";
    cfg.impactIonization.currentApproximation = "cell_reconstructed";
    cfg.impactIonization.currentMagnitudeMode = "edge_scalar_abs";
    cfg.impactIonization.cellReconstructedMidpointDensity = "gss_logistic";
    cfg.impactIonization.quasiFermiGradientDiscretization = "cell_gradient";
    cfg.impactIonization.sourceVolumePolicy = "genius_truncated";
    cfg.impactIonization.sourceVolumeFactor = 0.0;
    cfg.impactIonization.sourceGeometryScale = 1.0;
    cfg.impactIonization.edgeSourcePartition = "symmetric";
    cfg.impactIonization.sourceMappingMode = "triangle_gss_gradqf_truncated";
    cfg.impactIonization.sourceJacobianMode = "local_ad";
    cfg.impactIonization.electronA = 1.0;
    cfg.impactIonization.electronB = 1.0e-30;
    cfg.impactIonization.holeA = 1.0;
    cfg.impactIonization.holeB = 1.0e-30;

    DDSolution state;
    state.psi = (VectorXd(5) << 0.05, -0.30, -0.75, -0.10, -0.42).finished();
    state.phin = (VectorXd(5) << 0.02, -0.65, -1.10, -0.08, -0.27).finished();
    state.phip = (VectorXd(5) << -0.15, -0.35, -0.82, -0.60, -0.24).finished();

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const auto rows = solver.evaluateJacobianBlockAudit(
        state, 1.0e-7, std::vector<std::string>{"sg_avalanche"});

    REQUIRE(rows.size() == 1);
    REQUIRE(rows.front().block == "sg_avalanche");
    REQUIRE(rows.front().fdNorm > 0.0);
    REQUIRE(rows.front().analyticNorm > 0.0);
    REQUIRE(rows.front().relDiff < 2.0e-5);
}

TEST_CASE("NewtonSolver: frozen triangle GSS avalanche Jacobian keeps only the residual",
          "[newton][diagnostics][impact][triangle_gss]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -1.0},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.mobility.model = "constant";
    cfg.mobility.highFieldDrivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.model = "selberherr";
    cfg.impactIonization.drivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.generation = "current_density";
    cfg.impactIonization.currentApproximation = "cell_reconstructed";
    cfg.impactIonization.currentMagnitudeMode = "edge_scalar_abs";
    cfg.impactIonization.cellReconstructedMidpointDensity = "gss_logistic";
    cfg.impactIonization.quasiFermiGradientDiscretization = "cell_gradient";
    cfg.impactIonization.sourceVolumePolicy = "genius_truncated";
    cfg.impactIonization.sourceVolumeFactor = 0.0;
    cfg.impactIonization.sourceGeometryScale = 1.0;
    cfg.impactIonization.edgeSourcePartition = "symmetric";
    cfg.impactIonization.sourceMappingMode = "triangle_gss_gradqf_truncated";
    cfg.impactIonization.sourceJacobianMode = "frozen";
    cfg.impactIonization.electronA = 1.0;
    cfg.impactIonization.electronB = 1.0e-30;
    cfg.impactIonization.holeA = 1.0;
    cfg.impactIonization.holeB = 1.0e-30;

    DDSolution state;
    state.psi = (VectorXd(5) << 0.05, -0.30, -0.75, -0.10, -0.42).finished();
    state.phin = (VectorXd(5) << 0.02, -0.65, -1.10, -0.08, -0.27).finished();
    state.phip = (VectorXd(5) << -0.15, -0.35, -0.82, -0.60, -0.24).finished();

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const auto rows = solver.evaluateJacobianBlockAudit(
        state, 1.0e-7, std::vector<std::string>{"sg_avalanche"});

    REQUIRE(rows.size() == 1);
    REQUIRE(rows.front().fdNorm > 0.0);
    REQUIRE(rows.front().analyticNorm == Catch::Approx(0.0).margin(1.0e-30));
}

TEST_CASE("NewtonSolver: arithmetic cell-reconstructed SG avalanche Jacobian matches residual",
          "[newton][diagnostics][impact]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -1.0},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.mobility.model = "constant";
    cfg.mobility.highFieldDrivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.model = "selberherr";
    cfg.impactIonization.drivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.generation = "current_density";
    cfg.impactIonization.currentApproximation = "cell_reconstructed";
    cfg.impactIonization.cellReconstructedMidpointDensity = "arithmetic";
    cfg.impactIonization.sourceVolumePolicy = "edge_box";
    cfg.impactIonization.electronA = 1.0;
    cfg.impactIonization.electronB = 1.0e-30;
    cfg.impactIonization.holeA = 1.0;
    cfg.impactIonization.holeB = 1.0e-30;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, 0.0, -0.8);
    state.phin = VectorXd::LinSpaced(N, 0.0, -1.4);
    state.phip = VectorXd::LinSpaced(N, -0.2, -0.9);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const auto rows = solver.evaluateJacobianBlockAudit(
        state, 1.0e-7, std::vector<std::string>{"sg_avalanche"});

    REQUIRE(rows.size() == 1);
    REQUIRE(rows.front().block == "sg_avalanche");
    REQUIRE(rows.front().fdNorm > 0.0);
    REQUIRE(rows.front().analyticNorm > 0.0);
    REQUIRE(rows.front().relDiff < 5.0e-3);
}

TEST_CASE("NewtonSolver: conserved-total-current SG avalanche Jacobian matches residual",
          "[newton][diagnostics][impact]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -1.0},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.mobility.model = "constant";
    cfg.mobility.highFieldDrivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.model = "selberherr";
    cfg.impactIonization.drivingForce = "quasi_fermi_gradient";
    cfg.impactIonization.generation = "current_density";
    cfg.impactIonization.currentApproximation = "conserved_total_current";
    cfg.impactIonization.sourceVolumePolicy = "edge_box";
    cfg.impactIonization.electronA = 1.0;
    cfg.impactIonization.electronB = 1.0e-30;
    cfg.impactIonization.holeA = 1.0;
    cfg.impactIonization.holeB = 1.0e-30;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, 0.0, -0.8);
    state.phin = VectorXd::LinSpaced(N, 0.0, -1.4);
    state.phip = VectorXd::LinSpaced(N, -0.2, -0.9);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    PerformanceProfiler profiler({true, "unused.json"});
    const auto rows = [&] {
        ActivePerformanceProfilerScope active(&profiler);
        return solver.evaluateJacobianBlockAudit(
            state, 1.0e-7, std::vector<std::string>{"sg_avalanche"});
    }();

    REQUIRE(rows.size() == 1);
    REQUIRE(rows.front().block == "sg_avalanche");
    REQUIRE(rows.front().fdNorm > 0.0);
    REQUIRE(rows.front().analyticNorm > 0.0);
    REQUIRE(rows.front().relDiff < 5.0e-3);
    REQUIRE(profiler.toJson().at("counters").at(
        "jacobian.edge_avalanche_carrier_side_reuses") == 0);
}

TEST_CASE("NewtonSolver: evaluateBlockStep freezes complementary unknown blocks",
          "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonBlockStepEvaluation poisson =
        solver.evaluateBlockStep(state, "poisson_only");
    const NewtonBlockStepEvaluation carriers =
        solver.evaluateBlockStep(state, "carrier_only");

    REQUIRE(poisson.mode == "poisson_only");
    REQUIRE(carriers.mode == "carrier_only");
    REQUIRE(poisson.residual.raw.size() == 3 * N);
    REQUIRE(carriers.residual.raw.size() == 3 * N);
    REQUIRE(poisson.deltaPsi.size() == N);
    REQUIRE(carriers.deltaPhin.size() == N);
    REQUIRE(poisson.deltaPsi.norm() > 0.0);
    REQUIRE(poisson.deltaPhin.norm() == Catch::Approx(0.0));
    REQUIRE(poisson.deltaPhip.norm() == Catch::Approx(0.0));
    REQUIRE(carriers.deltaPsi.norm() == Catch::Approx(0.0));
    REQUIRE(carriers.deltaPhin.norm() + carriers.deltaPhip.norm() > 0.0);
    REQUIRE(poisson.trialSolution.psi(0) ==
            Catch::Approx(state.psi(0) + poisson.deltaPsi(0)));
    REQUIRE(carriers.trialSolution.phin(1) ==
            Catch::Approx(state.phin(1) + carriers.deltaPhin(1)));
}

TEST_CASE("NewtonSolver: builds a Poisson-block initialized cold state",
          "[newton][poisson_block_initialization]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    std::unordered_map<std::string, Real> biases = {
        {"anode", -1.0e-14},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.maxUpdate = 0.0;

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);

    const auto init = solver.buildPoissonBlockInitialization();
    const auto coldResidual = solver.evaluateResidual(init.coldInitial);
    const auto poissonResidual = solver.evaluateResidual(init.poissonBlockInitial);

    REQUIRE(init.coldInitial.psi.size() == init.poissonBlockInitial.psi.size());
    REQUIRE(init.rawStepNorm > 0.0);
    REQUIRE(poissonResidual.blockNorms.psi < coldResidual.blockNorms.psi);
    REQUIRE(poissonResidual.blockNorms.phin ==
            Catch::Approx(coldResidual.blockNorms.phin));
    REQUIRE(poissonResidual.blockNorms.phip ==
            Catch::Approx(coldResidual.blockNorms.phip));
}

TEST_CASE("NewtonSolver: ABA Poisson solve reconstructs contact-basin QFs",
          "[newton][poisson_only][aba]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.maxIter = 100;
    cfg.reltol = 1.0e-8;
    cfg.abstol = 1.0e-8;

    const std::unordered_map<std::string, Real> equilibriumBiases = {
        {"anode", 0.0}, {"cathode", 0.0}};
    const NewtonResult equilibrium = runNewton(
        mesh, matdb, doping, equilibriumBiases, cfg);
    REQUIRE(equilibrium.converged);

    const std::unordered_map<std::string, Real> biasedContacts = {
        {"anode", -0.05}, {"cathode", 0.0}};
    const NewtonResult aba = runNewtonPoissonOnly(
        mesh, matdb, doping, biasedContacts, equilibrium.solution, cfg);
    REQUIRE(aba.converged);
    REQUIRE(aba.convergenceReason.find("poisson_only") != std::string::npos);

    // Node 4 is the sole interior node. Its carrier state is recomputed from
    // the nearest majority-carrier contact basins, while both contact QFs are
    // still projected exactly.
    REQUIRE(aba.solution.n(4) > 0.0);
    REQUIRE(aba.solution.p(4) > 0.0);
    REQUIRE(aba.solution.phin(0) == Catch::Approx(-0.05));
    REQUIRE(aba.solution.phip(0) == Catch::Approx(-0.05));
    REQUIRE((aba.solution.psi - equilibrium.solution.psi).norm() > 0.0);
}

TEST_CASE("NewtonSolver: Poisson-only accepts a qualified numerical-floor stall",
          "[newton][poisson_only][stall_floor]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig equilibriumCfg = newtonConfig();
    equilibriumCfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    equilibriumCfg.recombination = {"none"};
    const NewtonResult equilibrium =
        runNewton(mesh, matdb, doping, zeroBias(), equilibriumCfg);
    REQUIRE(equilibrium.converged);

    DDSolution restart = equilibrium.solution;
    restart.psi(4) += 1.0e-5;

    NewtonConfig cfg = equilibriumCfg;
    cfg.warmStart = true;
    cfg.maxIter = 20;
    cfg.reltol = 1.0e-30;
    cfg.abstol = 1.0e-30;
    cfg.poissonLineSearchStallResidualFloor = 1.0e-6;
    cfg.poissonLineSearchStallCarrierResidualFloor = 1.0e-6;

    const NewtonResult result =
        runNewtonPoissonOnly(mesh, matdb, doping, zeroBias(), restart, cfg);
    REQUIRE(result.converged);
    REQUIRE(result.finalBlockNorms.psi <=
            cfg.poissonLineSearchStallResidualFloor);
    REQUIRE(result.finalBlockNorms.phin <=
            cfg.poissonLineSearchStallCarrierResidualFloor);
    REQUIRE(result.finalBlockNorms.phip <=
            cfg.poissonLineSearchStallCarrierResidualFloor);
    REQUIRE(result.convergenceReason == "poisson_only_line_search_stall_floor");
}

TEST_CASE("NewtonSolver: evaluateRegularizedCarrierStep damps carrier-only correction",
          "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.maxUpdate = 0.0;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonBlockStepEvaluation baseline =
        solver.evaluateBlockStep(state, "carrier_only");
    const NewtonRegularizedCarrierStepEvaluation zero =
        solver.evaluateRegularizedCarrierStep(state, 0.0);
    const NewtonRegularizedCarrierStepEvaluation damped =
        solver.evaluateRegularizedCarrierStep(state, 10.0);

    REQUIRE(zero.regularizationScale == Catch::Approx(0.0));
    REQUIRE(damped.regularizationScale == Catch::Approx(10.0));
    REQUIRE(zero.deltaPsi.norm() == Catch::Approx(0.0));
    REQUIRE(damped.deltaPsi.norm() == Catch::Approx(0.0));
    REQUIRE((zero.deltaPhin - baseline.deltaPhin).norm() ==
            Catch::Approx(0.0).margin(1.0e-12));
    REQUIRE((zero.deltaPhip - baseline.deltaPhip).norm() ==
            Catch::Approx(0.0).margin(1.0e-12));
    REQUIRE(damped.rawStepNorm < baseline.rawStepNorm);
}

TEST_CASE("NewtonSolver: evaluateCarrierRowDiagnostics reports carrier row stiffness",
          "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;
    cfg.maxUpdate = 1.0;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonCarrierRowDiagnosticsEvaluation rows =
        solver.evaluateCarrierRowDiagnostics(state);

    REQUIRE(rows.rows.size() == static_cast<std::size_t>(N));
    REQUIRE(rows.rawCarrierStepNorm > 0.0);
    REQUIRE(rows.cappedCarrierStepNorm > 0.0);
    REQUIRE(rows.rows[0].nodeId == 0);
    REQUIRE(rows.rows[0].electronRowAbsSum >= std::abs(rows.rows[0].electronDiagonal));
    REQUIRE(rows.rows[0].holeRowAbsSum >= std::abs(rows.rows[0].holeDiagonal));
    REQUIRE(rows.rows[0].electronRowL2Norm >= 0.0);
    REQUIRE(rows.rows[0].holeRowL2Norm >= 0.0);
    REQUIRE(rows.rows[0].rawDeltaPhin_V != Catch::Approx(0.0));
    REQUIRE(std::abs(rows.rows[0].cappedDeltaPhin_V) <=
            cfg.maxUpdate * rows.potentialScale + 1.0e-12);
}

TEST_CASE("NewtonSolver: Poisson linear diagnostics compare equivalent equilibrated solve",
          "[newton][diagnostics][linear_equilibration]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1}, {"cathode", 0.0}};

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"none"};
    cfg.warmStart = true;

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution state;
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const auto diagnostics =
        solver.evaluatePoissonLinearDiagnostics(state, 4);

    REQUIRE(diagnostics.rows.size() == static_cast<std::size_t>(N));
    REQUIRE(diagnostics.focusNode == 4);
    REQUIRE(std::find(
        diagnostics.focusPatchNodes.begin(),
        diagnostics.focusPatchNodes.end(), 4) !=
        diagnostics.focusPatchNodes.end());
    REQUIRE(diagnostics.focusPatchRawCondition.numericalRank > 0);
    REQUIRE(diagnostics.focusPatchEquilibratedCondition.numericalRank > 0);
    REQUIRE(diagnostics.rawStepNorm > 0.0);
    REQUIRE(diagnostics.equilibratedStepNorm > 0.0);
    REQUIRE(diagnostics.equilibratedRelativeLinearClosure <=
            diagnostics.rawRelativeLinearClosure * 1.01 + 1.0e-12);
    REQUIRE(diagnostics.rows[4].poissonRowL2Norm > 0.0);
    REQUIRE(diagnostics.rows[4].poissonColumnL2Norm > 0.0);
}

TEST_CASE("NewtonSolver: carrier block decomposition separates scale and coupling",
          "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"srh"};
    cfg.warmStart = true;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonCarrierBlockDecompositionEvaluation decomposition =
        solver.evaluateCarrierBlockDecomposition(state);

    const Index freeUnknowns =
        decomposition.freeElectronUnknowns + decomposition.freeHoleUnknowns;
    REQUIRE(freeUnknowns > 0);
    REQUIRE(decomposition.columns.size() ==
            static_cast<std::size_t>(freeUnknowns));
    REQUIRE(decomposition.singularModes.size() ==
            static_cast<std::size_t>(freeUnknowns));
    REQUIRE(decomposition.rawCondition.rows == freeUnknowns);
    REQUIRE(decomposition.rawCondition.columns == freeUnknowns);
    REQUIRE(decomposition.rawCondition.largestSingularValue > 0.0);
    REQUIRE(decomposition.freeColumnNormSpread >= 1.0);
    REQUIRE(decomposition.freeRowNormSpread >= 1.0);
    REQUIRE(decomposition.solveVariants.size() == 6);
    REQUIRE(decomposition.solveVariants[0].name == "full");
    REQUIRE(decomposition.solveVariants[1].name == "row_scaled");
    REQUIRE(decomposition.solveVariants[1].relativeDifferenceFromFull ==
            Catch::Approx(0.0).margin(1.0e-9));
    REQUIRE(decomposition.solveVariants[1].cosineWithFull ==
            Catch::Approx(1.0).margin(1.0e-9));
    Real stepEnergySum = 0.0;
    Real rhsEnergySum = 0.0;
    for (const auto& mode : decomposition.singularModes) {
        stepEnergySum += mode.stepEnergyFraction;
        rhsEnergySum += mode.rhsEnergyFraction;
        REQUIRE(std::abs(mode.jacobianProjectionClosure) <=
                std::max(1.0e-12, std::abs(mode.singularValue) * 1.0e-10));
        REQUIRE(std::abs(mode.rhsProjectionClosure) <=
                std::max(1.0e-12, std::abs(mode.rhsProjection) * 1.0e-10));
    }
    REQUIRE(stepEnergySum == Catch::Approx(1.0).margin(1.0e-10));
    REQUIRE(rhsEnergySum == Catch::Approx(1.0).margin(1.0e-10));
    REQUIRE(decomposition.transportCrossNorm == Catch::Approx(0.0).margin(1.0e-18));
}

TEST_CASE("NewtonSolver: evaluateCarrierTermDiagnostics decomposes continuity residual",
          "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    std::unordered_map<std::string, Real> biases = {
        {"anode", -0.1},
        {"cathode", 0.0},
    };

    NewtonConfig cfg;
    cfg.inputScaling.mode = UnitScalingMode::UnitScaling;
    cfg.recombination = {"srh"};
    cfg.warmStart = true;

    DDSolution state;
    const int N = static_cast<int>(mesh.numNodes());
    state.psi = VectorXd::LinSpaced(N, -0.02, 0.02);
    state.phin = VectorXd::Constant(N, -0.01);
    state.phip = VectorXd::Constant(N, 0.01);

    NewtonSolver solver(mesh, matdb, doping, biases, cfg);
    const NewtonCarrierTermDiagnosticsEvaluation terms =
        solver.evaluateCarrierTermDiagnostics(state);

    REQUIRE(terms.rows.size() == static_cast<std::size_t>(N));
    const auto& center = terms.rows[4];
    REQUIRE(center.nodeId == 4);
    REQUIRE(center.electronBoundary == Catch::Approx(0.0));
    const Real electronSum = center.electronFlux
        + center.electronRecombination
        + center.electronImpact
        + center.electronGauge
        + center.electronBoundary;
    const Real holeSum = center.holeFlux
        + center.holeRecombination
        + center.holeImpact
        + center.holeGauge
        + center.holeBoundary;
    REQUIRE(electronSum == Catch::Approx(center.electronResidual).margin(1.0e-18));
    REQUIRE(holeSum == Catch::Approx(center.holeResidual).margin(1.0e-18));

    const NewtonCarrierTermDiagnosticsEvaluation solvedTerms =
        solver.evaluateCarrierTermDiagnostics(state, true);
    REQUIRE(solvedTerms.rows.size() == static_cast<std::size_t>(N));
    for (int node = 0; node < N; ++node) {
        REQUIRE(solvedTerms.rows[static_cast<std::size_t>(node)].electronResidual ==
                Catch::Approx(solvedTerms.residual.raw(N + node)).margin(1.0e-18));
        REQUIRE(solvedTerms.rows[static_cast<std::size_t>(node)].holeResidual ==
                Catch::Approx(solvedTerms.residual.raw(2 * N + node)).margin(1.0e-18));
    }
}

TEST_CASE("NewtonSolver: defaults to analytic Jacobian", "[newton]")
{
    const NewtonConfig cfg;
    REQUIRE(cfg.jacobian == "analytic");
    REQUIRE(cfg.quasiFermiReference == "none");
    REQUIRE_FALSE(cfg.warmStart);
    REQUIRE(cfg.quasiFermiUpdateLimit_V == Catch::Approx(0.0));
    REQUIRE(cfg.quasiFermiUpdateLimitMinority_V == Catch::Approx(0.0));
    REQUIRE(cfg.quasiFermiUpdateLimitMode == "componentwise_poisson_recorrect");
    REQUIRE(cfg.quasiFermiTrustRegionGrowthFactor == Catch::Approx(1.0));
    REQUIRE(cfg.quasiFermiTrustRegionMaxMultiplier == Catch::Approx(1.0));
    REQUIRE(cfg.quasiFermiTrustRegionExpansionThreshold == Catch::Approx(0.75));
    REQUIRE(cfg.quasiFermiTrustRegionShrinkFactor == Catch::Approx(1.0));
    REQUIRE(cfg.quasiFermiTrustRegionMinMultiplier == Catch::Approx(1.0));
    REQUIRE(cfg.carrierRegularizationScale == Catch::Approx(0.0));
    REQUIRE_FALSE(cfg.carrierDiagonalFloor.enabled);
    REQUIRE(cfg.contactBoundaryReconstruction == "dominant_signed_contact_mean");
        REQUIRE(cfg.contactBoundaryMinorityElectronRelaxation);
        REQUIRE(cfg.contactBoundaryMinorityElectronRelaxationBiasThreshold_V ==
            Catch::Approx(0.1));
        REQUIRE(cfg.contactBoundaryMinorityElectronRelaxationTwoTerminalOnly);
            REQUIRE(cfg.contactBoundaryMinorityElectronRelaxationContactSide == "p_contact_only");

    const NewtonConfig debugCfg = newtonConfigFromJson(nlohmann::json{
        {"jacobian", "finite_difference"},
        {"quasi_fermi_reference", "contact_majority"},
        {"warm_start", true},
        {"quasi_fermi_update_limit_V", 0.0259},
        {"quasi_fermi_update_limit_minority_V", 0.01},
        {"quasi_fermi_update_limit_mode", "uniform_trust_region"},
        {"quasi_fermi_trust_region_growth_factor", 2.0},
        {"quasi_fermi_trust_region_max_multiplier", 8.0},
        {"quasi_fermi_trust_region_expansion_threshold", 0.8},
        {"quasi_fermi_trust_region_shrink_factor", 0.5},
        {"quasi_fermi_trust_region_min_multiplier", 1.0},
        {"carrier_regularization_scale", 3.0},
        {"carrier_diagonal_floor", nlohmann::json{{"enabled", true}, {"scale", 2.0}, {"minority_density_ratio", 0.25}}},
        {"local_update_diagnostics", nlohmann::json{
            {"enabled", true}, {"csv_file", "local_updates.csv"},
            {"nodes", nlohmann::json::array({1, 3})},
            {"first_iterations", 4}, {"every_iterations", 7}}},
        {"contact_boundary_reconstruction", "legacy_node_local"},
    });
    REQUIRE(debugCfg.jacobian == "finite_difference");
    REQUIRE(debugCfg.quasiFermiReference == "contact_majority");
    REQUIRE(debugCfg.warmStart);
    REQUIRE(debugCfg.quasiFermiUpdateLimit_V == Catch::Approx(0.0259));
    REQUIRE(debugCfg.quasiFermiUpdateLimitMinority_V == Catch::Approx(0.01));
    REQUIRE(debugCfg.quasiFermiUpdateLimitMode == "uniform_trust_region");
    REQUIRE(debugCfg.quasiFermiTrustRegionGrowthFactor == Catch::Approx(2.0));
    REQUIRE(debugCfg.quasiFermiTrustRegionMaxMultiplier == Catch::Approx(8.0));
    REQUIRE(debugCfg.quasiFermiTrustRegionExpansionThreshold == Catch::Approx(0.8));
    REQUIRE(debugCfg.quasiFermiTrustRegionShrinkFactor == Catch::Approx(0.5));
    REQUIRE(debugCfg.quasiFermiTrustRegionMinMultiplier == Catch::Approx(1.0));
    REQUIRE(debugCfg.carrierRegularizationScale == Catch::Approx(3.0));
    REQUIRE(debugCfg.carrierDiagonalFloor.enabled);
    REQUIRE(debugCfg.carrierDiagonalFloor.scale == Catch::Approx(2.0));
    REQUIRE(debugCfg.carrierDiagonalFloor.minorityDensityRatio == Catch::Approx(0.25));
    REQUIRE(debugCfg.localUpdateDiagnostics.enabled);
    REQUIRE(debugCfg.localUpdateDiagnostics.csvFile == "local_updates.csv");
    REQUIRE(debugCfg.localUpdateDiagnostics.nodes == std::vector<Index>{1, 3});
    REQUIRE(debugCfg.localUpdateDiagnostics.firstIterations == 4);
    REQUIRE(debugCfg.localUpdateDiagnostics.everyIterations == 7);
    REQUIRE(debugCfg.contactBoundaryReconstruction == "legacy_node_local");

    const NewtonConfig basinCfg = newtonConfigFromJson(
        nlohmann::json{{"quasi_fermi_reference", "contact_basin"}});
    REQUIRE(basinCfg.quasiFermiReference == "contact_basin");
}

TEST_CASE("NewtonSolver: unit_scaling config records scaled mode and preserves analytic Jacobian",
          "[newton][scaling][config]")
{
    const NewtonConfig cfg = newtonConfigFromJson(
        nlohmann::json{{"jacobian", "analytic"}},
        UnitScalingConfig{UnitScalingMode::UnitScaling});

    REQUIRE(cfg.inputScaling.isUnitScaling());
    REQUIRE(cfg.jacobian == "analytic");
}

TEST_CASE("NewtonSolver: parses Fermi-Dirac carrier statistics", "[newton][config][fermi_dirac]")
{
    const NewtonConfig stringConfig = newtonConfigFromJson(
        nlohmann::json{{"carrier_statistics", "fermi_dirac"}});
    REQUIRE(stringConfig.carrierStatistics.model == "fermi_dirac");

    const NewtonConfig objectConfig = newtonConfigFromJson(nlohmann::json{
        {"carrier_statistics", {{"model", "fermi_dirac"}}}
    });
    REQUIRE(objectConfig.carrierStatistics.model == "fermi_dirac");
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"carrier_statistics", "unknown"}}),
        std::invalid_argument);
}

TEST_CASE("NewtonSolver: warm start preserves supplied quasi-Fermi guess", "[newton][warm_start]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution initial;
    initial.psi = VectorXd::Zero(N);
    initial.phin = VectorXd::Zero(N);
    initial.phip = VectorXd::Zero(N);
    const int interiorNode = 4;
    initial.phin(interiorNode) = 0.02;
    initial.phip(interiorNode) = -0.015;

    NewtonConfig coldCfg = newtonConfig();
    coldCfg.maxIter = 0;
    coldCfg.reltol = 0.0;
    coldCfg.abstol = 0.0;
    coldCfg.warmStart = false;

    NewtonConfig warmCfg = coldCfg;
    warmCfg.warmStart = true;

    const NewtonResult cold = runNewton(mesh, matdb, doping, zeroBias(), initial, coldCfg);
    const NewtonResult warm = runNewton(mesh, matdb, doping, zeroBias(), initial, warmCfg);

    REQUIRE_FALSE(warm.converged);

    REQUIRE(cold.solution.phin(interiorNode) == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(cold.solution.phip(interiorNode) == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(warm.solution.phin(interiorNode) ==
            Catch::Approx(initial.phin(interiorNode)).margin(1.0e-14));
    REQUIRE(warm.solution.phip(interiorNode) ==
            Catch::Approx(initial.phip(interiorNode)).margin(1.0e-14));
    REQUIRE(warm.initialResidualNorm != Catch::Approx(cold.initialResidualNorm));
}

TEST_CASE("NewtonSolver: warm start projects contact nodes to the current bias",
          "[newton][warm_start][contacts]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution stale;
    stale.psi = VectorXd::Zero(N);
    stale.phin = VectorXd::Zero(N);
    stale.phip = VectorXd::Zero(N);
    stale.n = VectorXd::Ones(N);
    stale.p = VectorXd::Ones(N);

    const int interiorNode = 4;
    stale.phin(interiorNode) = 0.02;
    stale.phip(interiorNode) = -0.015;

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.warmStart = true;

    const Real anodeBias = -0.25;
    const NewtonResult result = runNewton(
        mesh,
        matdb,
        doping,
        {{"anode", anodeBias}, {"cathode", 0.0}},
        stale,
        cfg);

    REQUIRE_FALSE(result.converged);
    REQUIRE(result.solution.phin(interiorNode) ==
            Catch::Approx(stale.phin(interiorNode)).margin(1.0e-14));
    REQUIRE(result.solution.phip(interiorNode) ==
            Catch::Approx(stale.phip(interiorNode)).margin(1.0e-14));
    for (Index nid : mesh.getContact(1).node_ids) {
        const int ii = static_cast<int>(nid);
        REQUIRE(result.solution.phip(ii) == Catch::Approx(anodeBias).margin(1.0e-14));
    }
}



TEST_CASE("NewtonSolver: block residual norm balances mixed equation blocks", "[newton][residual]")
{
    const int N = 2;
    VectorXd initial(3 * N);
    initial << 1.0e-18, -1.0e-18,
               10.0, -10.0,
               5.0, -5.0;
    VectorXd current(3 * N);
    current << 1.0e-18, 0.0,
               1.0, -1.0,
               0.5, -0.5;

    const ResidualBlockNormValue initialBlocks = ResidualNorm::computeBlocks(initial, N);
    const ResidualBlockNormValue currentBlocks = ResidualNorm::computeBlocks(current, N);

    REQUIRE(initialBlocks.psi == Catch::Approx(std::sqrt(2.0) * 1.0e-18));
    REQUIRE(initialBlocks.phin == Catch::Approx(std::sqrt(200.0)));
    REQUIRE(initialBlocks.phip == Catch::Approx(std::sqrt(50.0)));

    const Real balanced = ResidualNorm::normalizedBlockL2(currentBlocks, initialBlocks);
    REQUIRE(balanced == Catch::Approx(std::sqrt(0.5 + 0.01 + 0.01)));

    ResidualBlockWeights continuityOnly;
    continuityOnly.psi = 0.0;
    continuityOnly.phin = 1.0;
    continuityOnly.phip = 4.0;
    const Real weighted = ResidualNorm::normalizedBlockL2(
        currentBlocks, initialBlocks, continuityOnly);
    REQUIRE(weighted == Catch::Approx(std::sqrt(0.01 + 4.0 * 0.01)));
}

TEST_CASE("NewtonSolver: evaluates residual for an externally supplied state",
          "[newton][residual]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    NewtonSolver solver(mesh, matdb, doping, zeroBias(), cfg);

    const int n = static_cast<int>(mesh.numNodes());
    DDSolution state;
    state.psi = VectorXd::LinSpaced(n, -0.02, 0.03);
    state.phin = VectorXd::Constant(n, 0.004);
    state.phip = VectorXd::Constant(n, -0.003);
    state.n = VectorXd::Zero(n);
    state.p = VectorXd::Zero(n);

    const NewtonResidualEvaluation residual = solver.evaluateResidual(state);

    REQUIRE(residual.raw.size() == 3 * n);
    REQUIRE(residual.blockNorms.psi > 0.0);
    REQUIRE(residual.blockNorms.phin > 0.0);
    REQUIRE(residual.blockNorms.phip > 0.0);
    REQUIRE(residual.blockNorms.combined == Catch::Approx(residual.raw.norm()));
    REQUIRE(residual.intrinsicDensity.size() == static_cast<std::size_t>(n));
    REQUIRE(residual.potentialScale > 0.0);
}

TEST_CASE("NewtonSolver: diagnostic probes preserve referenced sub-ULP increments",
          "[newton][diagnostics][qf-reference][precision]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    NewtonConfig cfg = newtonConfig();
    cfg.quasiFermiReference = "contact_majority";
    cfg.warmStart = true;
    const std::unordered_map<std::string, Real> biases = {
        {"anode", 0.0}, {"cathode", 1.1}};
    NewtonSolver solver(mesh, matdb, doping, biases, cfg);

    const int n = static_cast<int>(mesh.numNodes());
    DDSolution zero;
    zero.psi = VectorXd::Constant(n, 1.1);
    zero.phin = VectorXd::Constant(n, 1.1);
    zero.phip = VectorXd::Zero(n);
    zero.phinIncrement = VectorXd::Zero(n);
    zero.phipIncrement = VectorXd::Zero(n);
    zero.electronQfReference = VectorXd::Constant(n, 1.1);
    zero.holeQfReference = VectorXd::Zero(n);
    zero.electronQfReference_V = 1.1;
    zero.holeQfReference_V = 0.0;
    zero.n = VectorXd::Zero(n);
    zero.p = VectorXd::Zero(n);

    DDSolution perturbed = zero;
    const Real subUlpIncrement = 6.0e-17;
    REQUIRE(1.1 + subUlpIncrement == 1.1);
    perturbed.phinIncrement(4) = subUlpIncrement;

    const NewtonResidualEvaluation zeroResidual = solver.evaluateResidual(zero);
    const NewtonResidualEvaluation perturbedResidual =
        solver.evaluateResidual(perturbed);
    REQUIRE((perturbedResidual.raw - zeroResidual.raw).norm() > 0.0);
}

TEST_CASE("NewtonSolver: reports maximum contact majority quasi-Fermi drop", "[newton][contact][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    for (Index node = 0; node < mesh.numNodes(); ++node)
        doping.setNodeDoping(node, 1.0e21, 0.0);

    NewtonSolver solver(mesh, matdb, doping, zeroBias(), newtonConfig());

    const int n = static_cast<int>(mesh.numNodes());
    DDSolution state;
    state.psi = VectorXd::Zero(n);
    state.phin = VectorXd::Zero(n);
    state.phip = VectorXd::Zero(n);
    state.n = VectorXd::Constant(n, 1.0e21);
    state.p = VectorXd::Constant(n, 1.0e21);

    state.phin(1) = 1.2e-10;
    state.phip(1) = 1.0e-5;
    REQUIRE(solver.maxContactMajorityQuasiFermiDrop(state) == Catch::Approx(1.2e-10));
    REQUIRE(solver.maxContactMajorityQuasiFermiDrop(state, {"cathode"}) ==
            Catch::Approx(1.2e-10));
    REQUIRE(solver.maxContactMajorityQuasiFermiDrop(state, {"anode"}) ==
            Catch::Approx(0.0));

    // A contact may share mesh nodes/edges with an insulator.  A quasi-Fermi
    // value on a zero-density edge is not a transport contact drop and must not
    // block numerical-floor convergence.
    state.n(1) = 0.0;
    state.phin(1) = 2.0;
    REQUIRE(solver.maxContactMajorityQuasiFermiDrop(state) == Catch::Approx(0.0));
}
TEST_CASE("NewtonSolver: parses block residual norm controls", "[newton][config]")
{
    const NewtonConfig cfg = newtonConfigFromJson(nlohmann::json{
        {"residual_norm", "block"},
        {"line_search_mode", "block_filter"},
        {"residual_filter_gamma", 2.0e-4},
        {"residual_filter_envelope_factor", 1.75},
        {"max_update", 0.25},
        {"stall_residual_floor", 2.0e-8},
        {"poisson_line_search_stall_residual_floor", 3.0e-7},
        {"poisson_line_search_stall_relative_increase", 2.0e-5},
        {"poisson_line_search_stall_carrier_residual_floor", 4.0e-8},
        {"poisson_line_search_stall_contact_majority_qf_drop_limit_V", 7.0e-11},
        {"block_absolute_convergence", {
            {"mode", "enforce"},
            {"psi_residual_ceiling", 8.0e-8},
            {"electron_residual_ceiling", 9.0e-9},
            {"hole_residual_ceiling", 1.0e-9}
        }},
        {"contact_majority_qf_branch_drop_limit_V", 6.0e-10},
        {"contact_majority_qf_branch_guard_contacts", {"anode"}},
        {"linear_equilibration", {{"mode", "l2_row_column"}}},
        {"carrier_row_qualified_stall_acceptance", true},
        {"auger_cn_m6_per_s", 4.0e-43},
        {"auger_cp_m6_per_s", 2.0e-43},
        {"auger_excess_product", "classical_np"},
        {"residual_weights", {{"psi", 0.25}, {"phin", 2.0}, {"phip", 3.0}}},
        {"residual_scales", {{"psi", 1.0e-18}, {"phin", 2.0e4}, {"phip", 3.0e4}}}
    });

    REQUIRE(cfg.residualNorm == "block");
    REQUIRE(cfg.lineSearchMode == "block_filter");
    REQUIRE(cfg.residualFilterGamma == Catch::Approx(2.0e-4));
    REQUIRE(cfg.residualFilterEnvelopeFactor == Catch::Approx(1.75));
    REQUIRE(cfg.maxUpdate == Catch::Approx(0.25));
    REQUIRE(cfg.stallResidualFloor == Catch::Approx(2.0e-8));
    REQUIRE(cfg.poissonLineSearchStallResidualFloor == Catch::Approx(3.0e-7));
    REQUIRE(cfg.poissonLineSearchStallRelativeIncrease == Catch::Approx(2.0e-5));
    REQUIRE(cfg.poissonLineSearchStallCarrierResidualFloor == Catch::Approx(4.0e-8));
    REQUIRE(cfg.poissonLineSearchStallContactMajorityQfDropLimit_V == Catch::Approx(7.0e-11));
    REQUIRE(cfg.blockAbsoluteConvergence.mode == "enforce");
    REQUIRE(cfg.blockAbsoluteConvergence.psiResidualCeiling == Catch::Approx(8.0e-8));
    REQUIRE(cfg.blockAbsoluteConvergence.electronResidualCeiling == Catch::Approx(9.0e-9));
    REQUIRE(cfg.blockAbsoluteConvergence.holeResidualCeiling == Catch::Approx(1.0e-9));
    REQUIRE(cfg.contactMajorityQfBranchDropLimit_V == Catch::Approx(6.0e-10));
    REQUIRE(cfg.contactMajorityQfBranchGuardContacts ==
            std::vector<std::string>{"anode"});
    REQUIRE(cfg.linearEquilibration.mode == "l2_row_column");
    REQUIRE(cfg.carrierRowQualifiedStallAcceptance);
    REQUIRE(cfg.residualWeightPsi == Catch::Approx(0.25));
    REQUIRE(cfg.residualWeightPhin == Catch::Approx(2.0));
    REQUIRE(cfg.residualWeightPhip == Catch::Approx(3.0));
    REQUIRE(cfg.residualScalePsi == Catch::Approx(1.0e-18));
    REQUIRE(cfg.residualScalePhin == Catch::Approx(2.0e4));
    REQUIRE(cfg.residualScalePhip == Catch::Approx(3.0e4));
    REQUIRE(cfg.augerCn == Catch::Approx(4.0e-43));
    REQUIRE(cfg.augerCp == Catch::Approx(2.0e-43));
    REQUIRE(cfg.augerExcessProduct == "classical_np");

    const NewtonConfig boundaryCfg = newtonConfigFromJson(nlohmann::json{
        {"contact_boundary_minority_electron_relaxation", false},
        {"contact_boundary_minority_electron_relaxation_bias_threshold_V", 0.2},
        {"contact_boundary_minority_electron_relaxation_two_terminal_only", false},
        {"contact_boundary_minority_electron_relaxation_contact_side", "both_contacts"},
        {"contact_boundary_minority_electron_relaxation_strength", 0.5},
    });
    REQUIRE_FALSE(boundaryCfg.contactBoundaryMinorityElectronRelaxation);
    REQUIRE(boundaryCfg.contactBoundaryMinorityElectronRelaxationBiasThreshold_V ==
            Catch::Approx(0.2));
    REQUIRE_FALSE(boundaryCfg.contactBoundaryMinorityElectronRelaxationTwoTerminalOnly);
    REQUIRE(boundaryCfg.contactBoundaryMinorityElectronRelaxationContactSide ==
            "both_contacts");
    REQUIRE(boundaryCfg.contactBoundaryMinorityElectronRelaxationStrength ==
            Catch::Approx(0.5));

    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"residual_norm", "unknown"}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"line_search_mode", "unknown"}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"linear_equilibration", "unknown"}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"residual_filter_gamma", 0.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"residual_filter_envelope_factor", 0.5}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"max_update", -1.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{
            "electron_quantum_potential",
            {{"outer_absolute_tolerance_V", -1.0}}
        }}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"quasi_fermi_update_limit_V", -1.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"quasi_fermi_update_limit_minority_V", -1.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"quasi_fermi_update_limit_mode", "unknown"}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"quasi_fermi_trust_region_growth_factor", 0.5}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"quasi_fermi_trust_region_max_multiplier", 0.5}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"quasi_fermi_trust_region_expansion_threshold", 0.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"quasi_fermi_trust_region_shrink_factor", 0.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"quasi_fermi_trust_region_max_multiplier", 2.0},
            {"quasi_fermi_trust_region_min_multiplier", 3.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"quasi_fermi_trust_region_min_multiplier", 0.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"carrier_regularization_scale", -1.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"poisson_line_search_stall_residual_floor", -1.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{
            "block_absolute_convergence", {{"mode", "unknown"}}
        }}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{
            "block_absolute_convergence", {
                {"mode", "enforce"},
                {"psi_residual_ceiling", 1.0e-8},
                {"electron_residual_ceiling", 0.0},
                {"hole_residual_ceiling", 1.0e-9}
            }
        }}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{
            "contact_majority_qf_branch_drop_limit_V", -1.0
        }}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"poisson_line_search_stall_relative_increase", -1.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"poisson_line_search_stall_carrier_residual_floor", -1.0}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"carrier_diagonal_floor", nlohmann::json{{"enabled", true}, {"scale", -1.0}}}
        }),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{
            "quasi_fermi_update_limit_V",
            std::numeric_limits<Real>::infinity()
        }}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"contact_boundary_minority_electron_relaxation_bias_threshold_V", -1.0}
        }),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{"contact_boundary_reconstruction", "unexpected_mode"}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{
            "contact_boundary_minority_electron_relaxation_contact_side",
            "unexpected_side"
        }}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{{
            "contact_boundary_minority_electron_relaxation_strength",
            1.5
        }}),
        std::invalid_argument);
}

TEST_CASE("NewtonSolver: rejects disabled residual weights", "[newton][config]")
{
    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"residual_weights", {{"psi", 0.0}, {"phin", 0.0}, {"phip", 0.0}}}
        }),
        std::invalid_argument);

    REQUIRE_THROWS_AS(
        newtonConfigFromJson(nlohmann::json{
            {"residual_weights", {{"psi", -1.0}, {"phin", 0.0}, {"phip", 1.0}}}
        }),
        std::invalid_argument);

    const NewtonConfig cfg = newtonConfigFromJson(nlohmann::json{
        {"residual_weights", {{"psi", 0.0}, {"phin", 0.0}, {"phip", 1.0}}}
    });
    REQUIRE(cfg.residualWeightPsi == Catch::Approx(0.0));
    REQUIRE(cfg.residualWeightPhin == Catch::Approx(0.0));
    REQUIRE(cfg.residualWeightPhip == Catch::Approx(1.0));
}

TEST_CASE("NewtonSolver: verbose false suppresses failure diagnostics", "[newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution initial;
    initial.psi = VectorXd::LinSpaced(N, -0.03, 0.04);
    initial.phin = VectorXd::Constant(N, 0.01);
    initial.phip = VectorXd::Constant(N, -0.01);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.warmStart = true;
    cfg.verbose = false;

    std::ostringstream capturedStderr;
    std::streambuf* previousStderr = std::cerr.rdbuf(capturedStderr.rdbuf());
    const NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), initial, cfg);
    std::cerr.rdbuf(previousStderr);

    REQUIRE_FALSE(result.converged);
    REQUIRE(capturedStderr.str().empty());
}

TEST_CASE("NewtonSolver: pure-insulator zero carriers are valid during line search",
          "[newton][line_search][oxide]")
{
    DeviceMesh mesh = makePartiallyContactedOxideMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution initial;
    initial.psi = VectorXd::Zero(N);
    initial.phin = VectorXd::Constant(N, 0.1);
    initial.phip = VectorXd::Constant(N, -0.1);
    initial.phin(0) = 0.0;
    initial.phip(0) = 0.0;

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 3;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.verbose = false;
    cfg.lineSearch = true;
    cfg.warmStart = true;

    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"gate", 0.0}}, initial, cfg);

    REQUIRE(result.converged);
    REQUIRE(result.iters >= 1);
    REQUIRE_FALSE(result.history.empty());
    REQUIRE(result.trace.size() >= 2);
    REQUIRE(result.trace.front().event == "initial");
    REQUIRE(result.trace.back().event == "accepted_iteration");
    REQUIRE(result.trace.back().lineSearchAccepted);
    REQUIRE_FALSE(
        result.trace.back().sourceJacobianActiveBranchFingerprint.empty());
    REQUIRE(result.finalResidualNorm < result.initialResidualNorm);
    REQUIRE((result.solution.psi - initial.psi).norm() == Catch::Approx(0.0));
    REQUIRE(result.solution.phin.norm() == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(result.solution.phip.norm() == Catch::Approx(0.0).margin(1.0e-14));
}

TEST_CASE("NewtonSolver: metal gate constrains only electrostatic potential",
          "[newton][contact_bc][metal_gate][oxide]")
{
    DeviceMesh mesh = makePartiallyContactedOxideMesh();
    MaterialDatabase matdb;
    DopingModel doping(mesh.numNodes());

    const int N = static_cast<int>(mesh.numNodes());
    DDSolution initial;
    initial.psi = VectorXd::Constant(N, 0.55);
    initial.phin = VectorXd::Constant(N, 0.12);
    initial.phip = VectorXd::Constant(N, -0.08);

    ContactBoundarySpec gateSpec;
    gateSpec.name = "gate";
    gateSpec.type = ContactType::MetalGate;
    gateSpec.bias = 0.0;
    gateSpec.flatbandVoltage = -0.55;
    ContactSpecsMap specs{{"gate", gateSpec}};

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 4;
    cfg.reltol = 0.0;
    cfg.abstol = 1.0e-14;
    cfg.verbose = false;
    cfg.warmStart = true;

    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"gate", 0.55}}, initial, cfg, {}, {}, specs);

    REQUIRE(result.converged);
    for (int i = 0; i < N; ++i) {
        REQUIRE(result.solution.psi(i) == Catch::Approx(0.55).margin(1.0e-12));
        REQUIRE(result.solution.phin(i) == Catch::Approx(0.0).margin(1.0e-12));
        REQUIRE(result.solution.phip(i) == Catch::Approx(0.0).margin(1.0e-12));
    }
}

TEST_CASE("NewtonSolver: carrier regularization damps coupled Newton carrier mode",
          "[newton][line_search][regularization]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const std::unordered_map<std::string, Real> biases = {
        {"anode", -0.05},
        {"cathode", 0.0},
    };

    NewtonConfig seedCfg = newtonConfig();
    seedCfg.maxIter = 0;
    seedCfg.reltol = 0.0;
    seedCfg.abstol = 0.0;
    seedCfg.warmStart = true;
    DDSolution initial = runNewton(mesh, matdb, doping, biases, seedCfg).solution;
    const int interiorNode = 4;
    initial.phin(interiorNode) = 0.18;
    initial.phip(interiorNode) = -0.16;

    NewtonConfig baselineCfg = newtonConfig();
    baselineCfg.maxIter = 1;
    baselineCfg.reltol = 0.0;
    baselineCfg.abstol = 0.0;
    baselineCfg.lineSearch = false;
    baselineCfg.warmStart = true;

    NewtonConfig regularizedCfg = baselineCfg;
    regularizedCfg.carrierRegularizationScale = 10.0;

    const NewtonResult baseline = runNewton(
        mesh, matdb, doping, biases, initial, baselineCfg);
    const NewtonResult regularized = runNewton(
        mesh, matdb, doping, biases, initial, regularizedCfg);

    REQUIRE(baseline.iters == 1);
    REQUIRE(regularized.iters == 1);
    REQUIRE_FALSE(baseline.history.empty());
    REQUIRE_FALSE(regularized.history.empty());

    REQUIRE(regularized.history.front().rawStepNorm < baseline.history.front().rawStepNorm);
    REQUIRE((regularized.solution.psi - initial.psi).norm() > 0.0);
    REQUIRE((regularized.solution.phin - initial.phin).norm() <
            (baseline.solution.phin - initial.phin).norm());
}

TEST_CASE("NewtonSolver: max_update limits a large Newton step before line search",
          "[newton][line_search]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 1;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.verbose = false;
    cfg.lineSearch = false;
    cfg.maxUpdate = 0.05;

    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"anode", 0.05}, {"cathode", 0.0}}, cfg);

    REQUIRE(result.iters > 0);
    REQUIRE_FALSE(result.history.empty());
    REQUIRE(result.history.front().rawStepNorm == Catch::Approx(result.history.front().stepNorm));
    REQUIRE(result.history.front().stepNorm <=
            std::sqrt(static_cast<Real>(3 * mesh.numNodes())) * cfg.maxUpdate);
}

TEST_CASE("NewtonSolver: quasi-Fermi update limit caps updates in unit_scaling",
          "[newton][line_search][scaling]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const std::unordered_map<std::string, Real> biases = {
        {"anode", -0.05},
        {"cathode", 0.0},
    };

    NewtonConfig seedCfg = newtonConfig();
    seedCfg.maxIter = 0;
    seedCfg.reltol = 0.0;
    seedCfg.abstol = 0.0;
    seedCfg.warmStart = true;
    seedCfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    DDSolution initial = runNewton(mesh, matdb, doping, biases, seedCfg).solution;
    const int interiorNode = 4;
    initial.phin(interiorNode) = 0.18;
    initial.phip(interiorNode) = -0.16;

    NewtonConfig unclampedCfg = newtonConfig();
    unclampedCfg.maxIter = 1;
    unclampedCfg.reltol = 0.0;
    unclampedCfg.abstol = 0.0;
    unclampedCfg.dampingFactor = 1.0;
    unclampedCfg.lineSearch = false;
    unclampedCfg.warmStart = true;
    unclampedCfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};

    NewtonConfig clampedCfg = unclampedCfg;
    clampedCfg.quasiFermiUpdateLimit_V = 0.01;

    const NewtonResult unclamped = runNewton(
        mesh, matdb, doping, biases, initial, unclampedCfg);
    const NewtonResult clamped = runNewton(
        mesh, matdb, doping, biases, initial, clampedCfg);

    REQUIRE(unclamped.iters == 1);
    REQUIRE(clamped.iters == 1);

    const VectorXd unclampedPhinDelta = unclamped.solution.phin - initial.phin;
    const VectorXd unclampedPhipDelta = unclamped.solution.phip - initial.phip;
    const VectorXd clampedPhinDelta = clamped.solution.phin - initial.phin;
    const VectorXd clampedPhipDelta = clamped.solution.phip - initial.phip;

    const Real limit = clampedCfg.quasiFermiUpdateLimit_V;
    const Real maxUnclampedQfDelta = std::max(
        unclampedPhinDelta.cwiseAbs().maxCoeff(),
        unclampedPhipDelta.cwiseAbs().maxCoeff());

    REQUIRE(maxUnclampedQfDelta > limit);
    REQUIRE(clampedPhinDelta.cwiseAbs().maxCoeff() <= limit + 1.0e-12);
    REQUIRE(clampedPhipDelta.cwiseAbs().maxCoeff() <= limit + 1.0e-12);
    REQUIRE(clampedPhinDelta(interiorNode) ==
            Catch::Approx((unclampedPhinDelta(interiorNode) > 0.0 ? 1.0 : -1.0) * limit)
                .margin(1.0e-12));
    REQUIRE(clampedPhipDelta(interiorNode) ==
            Catch::Approx((unclampedPhipDelta(interiorNode) > 0.0 ? 1.0 : -1.0) * limit)
                .margin(1.0e-12));
}

TEST_CASE("NewtonSolver: minority quasi-Fermi update limit caps only the minority carrier per node",
          "[newton][line_search][scaling]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    // Force the single interior node to be strongly p-type so that electrons are
    // its minority carrier (phin minority) and holes are its majority carrier.
    const int interiorNode = 4;
    doping.setNodeDoping(interiorNode, 0.0, 1.0e21);

    const std::unordered_map<std::string, Real> biases = {
        {"anode", -0.05},
        {"cathode", 0.0},
    };

    NewtonConfig seedCfg = newtonConfig();
    seedCfg.maxIter = 0;
    seedCfg.reltol = 0.0;
    seedCfg.abstol = 0.0;
    seedCfg.warmStart = true;
    seedCfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    DDSolution initial = runNewton(mesh, matdb, doping, biases, seedCfg).solution;
    // Inject a large quasi-Fermi perturbation on both carriers at the interior node
    // so a single Newton step drives both updates well past the global cap.
    initial.phin(interiorNode) = 0.30;
    initial.phip(interiorNode) = -0.30;

    NewtonConfig globalCfg = newtonConfig();
    globalCfg.maxIter = 1;
    globalCfg.reltol = 0.0;
    globalCfg.abstol = 0.0;
    globalCfg.dampingFactor = 1.0;
    globalCfg.lineSearch = false;
    globalCfg.warmStart = true;
    globalCfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    globalCfg.quasiFermiUpdateLimit_V = 0.05;

    NewtonConfig minorityCfg = globalCfg;
    minorityCfg.quasiFermiUpdateLimitMinority_V = 0.01;

    const NewtonResult globalResult = runNewton(
        mesh, matdb, doping, biases, initial, globalCfg);
    const NewtonResult minorityResult = runNewton(
        mesh, matdb, doping, biases, initial, minorityCfg);

    REQUIRE(globalResult.iters == 1);
    REQUIRE(minorityResult.iters == 1);

    const Real globalLimit = globalCfg.quasiFermiUpdateLimit_V;
    const Real minorityLimit = minorityCfg.quasiFermiUpdateLimitMinority_V;

    const Real globalPhinDelta =
        globalResult.solution.phin(interiorNode) - initial.phin(interiorNode);
    const Real globalPhipDelta =
        globalResult.solution.phip(interiorNode) - initial.phip(interiorNode);
    const Real minorityPhinDelta =
        minorityResult.solution.phin(interiorNode) - initial.phin(interiorNode);
    const Real minorityPhipDelta =
        minorityResult.solution.phip(interiorNode) - initial.phip(interiorNode);

    // Global cap clips both carriers identically at the p-type interior node.
    REQUIRE(std::abs(globalPhinDelta) == Catch::Approx(globalLimit).margin(1.0e-12));
    REQUIRE(std::abs(globalPhipDelta) == Catch::Approx(globalLimit).margin(1.0e-12));

    // Minority cap tightens only the minority electron quasi-Fermi update; the
    // majority hole quasi-Fermi update keeps the looser global cap.
    REQUIRE(std::abs(minorityPhinDelta) == Catch::Approx(minorityLimit).margin(1.0e-12));
    REQUIRE(std::abs(minorityPhipDelta) == Catch::Approx(globalLimit).margin(1.0e-12));
    REQUIRE(std::abs(minorityPhinDelta) < std::abs(globalPhinDelta));
}

TEST_CASE("NewtonSolver: uniform quasi-Fermi trust region preserves the complete Newton direction",
          "[newton][line_search][scaling][trust-region]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const std::unordered_map<std::string, Real> biases = {
        {"anode", -0.05},
        {"cathode", 0.0},
    };

    NewtonConfig seedCfg = newtonConfig();
    seedCfg.maxIter = 0;
    seedCfg.reltol = 0.0;
    seedCfg.abstol = 0.0;
    seedCfg.warmStart = true;
    seedCfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    DDSolution initial = runNewton(mesh, matdb, doping, biases, seedCfg).solution;
    initial.phin(4) = 0.18;
    initial.phip(4) = -0.16;

    NewtonConfig rawCfg = newtonConfig();
    rawCfg.warmStart = true;
    rawCfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    const NewtonStepEvaluation raw =
        NewtonSolver(mesh, matdb, doping, biases, rawCfg).evaluateStep(initial);

    NewtonConfig trustCfg = rawCfg;
    trustCfg.quasiFermiUpdateLimit_V = 0.01;
    trustCfg.quasiFermiUpdateLimitMode = "uniform_trust_region";
    const NewtonStepEvaluation trust =
        NewtonSolver(mesh, matdb, doping, biases, trustCfg).evaluateStep(initial);

    const Real maxRawQf = std::max(
        raw.deltaPhin.cwiseAbs().maxCoeff(),
        raw.deltaPhip.cwiseAbs().maxCoeff());
    REQUIRE(maxRawQf > trustCfg.quasiFermiUpdateLimit_V);
    const Real expectedScale = trustCfg.quasiFermiUpdateLimit_V / maxRawQf;

    REQUIRE(trust.deltaPhin.cwiseAbs().maxCoeff() <=
            trustCfg.quasiFermiUpdateLimit_V + 1.0e-12);
    REQUIRE(trust.deltaPhip.cwiseAbs().maxCoeff() <=
            trustCfg.quasiFermiUpdateLimit_V + 1.0e-12);
    REQUIRE((trust.deltaPsi - expectedScale * raw.deltaPsi).norm() < 1.0e-11);
    REQUIRE((trust.deltaPhin - expectedScale * raw.deltaPhin).norm() < 1.0e-11);
    REQUIRE((trust.deltaPhip - expectedScale * raw.deltaPhip).norm() < 1.0e-11);
    REQUIRE(trust.stepNorm == Catch::Approx(expectedScale * raw.stepNorm).epsilon(1.0e-10));
}

TEST_CASE("NewtonSolver: max-iteration exit honors stall residual floor", "[newton][line_search]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.stallResidualFloor = 1.0e9;

    const NewtonResult guarded = runNewton(
        mesh, matdb, doping, {{"anode", 0.05}, {"cathode", 0.0}}, cfg);

    REQUIRE_FALSE(guarded.converged);

    // The floor remains explicitly available for diagnostic workflows that
    // opt out of the contact-majority QF safety guard.
    cfg.poissonLineSearchStallContactMajorityQfDropLimit_V = 0.0;
    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"anode", 0.05}, {"cathode", 0.0}}, cfg);

    REQUIRE(result.converged);
    REQUIRE(result.finalResidualNorm <= cfg.stallResidualFloor);
}


TEST_CASE("NewtonSolver: optionally records line-search diagnostics in history", "[newton][diagnostics]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.diagnostics = true;

    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"anode", 0.05}, {"cathode", 0.0}}, cfg);

    REQUIRE(result.converged);
    REQUIRE_FALSE(result.history.empty());
    REQUIRE(result.trace.size() == result.history.size() + 1);
    REQUIRE(result.trace.front().iter == 0);
    REQUIRE(result.trace.front().event == "initial");
    REQUIRE_FALSE(
        result.trace.front().sourceJacobianActiveBranchFingerprint.empty());
    const NewtonIterationInfo& first = result.history.front();
    REQUIRE(first.iter == 1);
    REQUIRE(first.rawStepNorm >= first.stepNorm);
    REQUIRE(first.stepNorm == Catch::Approx(first.dampingFactor * first.rawStepNorm));
    REQUIRE(first.relativeResidualNorm == Catch::Approx(
        ResidualNorm::relative(first.residualNorm, result.initialResidualNorm)));
    REQUIRE(first.lineSearchAccepted);
    REQUIRE(first.lineSearchAttempts >= 1);
    REQUIRE(first.lineSearchHistory.size() == static_cast<std::size_t>(first.lineSearchAttempts));
    REQUIRE(first.lineSearchHistory.back().accepted);
    REQUIRE(first.lineSearchHistory.back().damping == Catch::Approx(first.dampingFactor));
    REQUIRE(first.lineSearchHistory.back().residualNorm == Catch::Approx(first.residualNorm));
    REQUIRE(first.rowScaledBlockResiduals.psi >= 0.0);
    REQUIRE(first.rowScaledBlockResiduals.phin >= 0.0);
    REQUIRE(first.rowScaledBlockResiduals.phip >= 0.0);
    REQUIRE(first.topPoissonResidual.nodeId >= 0);
    REQUIRE(first.topElectronResidual.nodeId >= 0);
    REQUIRE(first.topHoleResidual.nodeId >= 0);
    REQUIRE(first.topPoissonResidual.absoluteResidual ==
            Catch::Approx(std::abs(first.topPoissonResidual.signedResidual)));
    REQUIRE_FALSE(first.sourceJacobianActiveBranchFingerprint.empty());
    REQUIRE(first.event == "accepted_iteration");

    const NewtonConfig parsed = newtonConfigFromJson(nlohmann::json{{"diagnostic_history", true}});
    REQUIRE(parsed.diagnostics);
}

TEST_CASE("NewtonSolver: configured temperature is parsed and passed to initial Gummel guess", "[newton][temperature]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg300 = newtonConfig();
    cfg300.maxIter = 0;
    cfg300.temperature_K = 300.0;
    const NewtonResult result300 = runNewton(mesh, matdb, doping, zeroBias(), cfg300);

    NewtonConfig cfg600 = cfg300;
    cfg600.temperature_K = 600.0;
    const NewtonResult result600 = runNewton(mesh, matdb, doping, zeroBias(), cfg600);

    const Real builtIn300 = result300.solution.psi(1) - result300.solution.psi(0);
    const Real builtIn600 = result600.solution.psi(1) - result600.solution.psi(0);
    REQUIRE(builtIn600 > 0.0);
    REQUIRE(builtIn600 < builtIn300);

    const NewtonConfig parsed = newtonConfigFromJson(nlohmann::json{{"temperature_K", 325.0}});
    REQUIRE(parsed.temperature_K == Catch::Approx(325.0));
    REQUIRE_THROWS_AS(newtonConfigFromJson(nlohmann::json{{"temperature_K", -1.0}}),
                      std::invalid_argument);
}

TEST_CASE("Newton solver parses Sentaurus E2 band-to-band parameters",
          "[newton][btbt][json]")
{
    const NewtonConfig cfg = newtonConfigFromJson(nlohmann::json{
        {"band_to_band", {
            {"model", "e2"},
            {"A_cm_inv_s_inv_V_inv2", 3.4e21},
            {"B_V_per_cm", 22.6e6},
            {"source_integration", "transport_node_lumped"},
            {"jacobian", "frozen_field"},
        }},
    });
    REQUIRE(cfg.bandToBand.model == "e2");
    REQUIRE(cfg.bandToBand.prefactorA_SI == Catch::Approx(3.4e23));
    REQUIRE(cfg.bandToBand.exponentialB_V_per_m == Catch::Approx(2.26e9));
    REQUIRE(cfg.bandToBand.sourceIntegration == "transport_node_lumped");
    REQUIRE(cfg.bandToBand.jacobian == "frozen_field");
}

static void requireFiniteNewtonSolution(const NewtonResult& result, Index nodeCount)
{
    REQUIRE(std::isfinite(result.initialResidualNorm));
    REQUIRE(std::isfinite(result.finalResidualNorm));
    for (Index i = 0; i < nodeCount; ++i) {
        const int ii = static_cast<int>(i);
        REQUIRE(std::isfinite(result.solution.psi(ii)));
        REQUIRE(std::isfinite(result.solution.phin(ii)));
        REQUIRE(std::isfinite(result.solution.phip(ii)));
        REQUIRE(std::isfinite(result.solution.n(ii)));
        REQUIRE(std::isfinite(result.solution.p(ii)));
        REQUIRE(result.solution.n(ii) >= 0.0);
        REQUIRE(result.solution.p(ii) >= 0.0);
    }
}

TEST_CASE("NewtonSolver: high doping gradient reverse bias does not diverge", "[newton][stability]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 1.0e21},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 8;
    cfg.reltol = 1.0e-6;
    cfg.abstol = 1.0e-18;
    cfg.lineSearch = true;
    cfg.verbose = false;

    const NewtonResult result = runNewton(
        mesh, matdb, doping, {{"anode", -0.10}, {"cathode", 0.0}}, cfg);

    REQUIRE(result.iters <= cfg.maxIter);
    requireFiniteNewtonSolution(result, mesh.numNodes());
    REQUIRE(result.finalResidualNorm <= result.initialResidualNorm * 1.0e6);
}

TEST_CASE("NewtonSolver: multi-terminal contacted mesh accepts distinct biases", "[newton][contacts]")
{
    DeviceMesh mesh;
    const double L = 1.0e-6;
    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;   n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = L;   n2.y = L;   mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0; n3.y = L;   mesh.addNode(n3);
    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0; c0.node_ids = {0, 1, 2}; mesh.addCell(c0);
    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 0; c1.node_ids = {0, 2, 3}; mesh.addCell(c1);
    Region r0; r0.id = 0; r0.name = "p_region"; r0.material = "Si"; r0.cell_ids = {0, 1}; mesh.addRegion(r0);
    Contact body; body.id = 0; body.name = "body"; body.region_id = 0; body.node_ids = {0}; mesh.addContact(body);
    Contact source; source.id = 1; source.name = "source"; source.region_id = 0; source.node_ids = {1}; mesh.addContact(source);
    Contact gate; gate.id = 2; gate.name = "gate"; gate.region_id = 0; gate.node_ids = {2}; mesh.addContact(gate);
    Contact drain; drain.id = 3; drain.name = "drain"; drain.region_id = 0; drain.node_ids = {3}; mesh.addContact(drain);
    mesh.buildEdges();

    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, std::vector<RegionDopingSpec>{{"p_region", 0.0, 1.0e21}});

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.warmStart = true;

    const NewtonResult result = runNewton(
        mesh,
        matdb,
        doping,
        {{"body", 0.0}, {"source", 0.02}, {"gate", 0.04}, {"drain", 0.06}},
        cfg);

    REQUIRE(result.iters == 0);
    requireFiniteNewtonSolution(result, mesh.numNodes());
    REQUIRE(result.solution.phin(0) == Catch::Approx(0.0));
    REQUIRE(result.solution.phip(0) == Catch::Approx(0.0));
    REQUIRE(result.solution.phin(1) == Catch::Approx(0.02));
    REQUIRE(result.solution.phip(1) == Catch::Approx(0.02));
    REQUIRE(result.solution.phin(2) == Catch::Approx(0.04));
    REQUIRE(result.solution.phip(2) == Catch::Approx(0.04));
    REQUIRE(result.solution.phin(3) == Catch::Approx(0.06));
    REQUIRE(result.solution.phip(3) == Catch::Approx(0.06));
}

TEST_CASE("NewtonSolver: high-bias ohmic contacts keep quasi-Fermi boundary targets when relaxation is disabled",
          "[newton][contacts]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.maxIter = 0;
    cfg.reltol = 0.0;
    cfg.abstol = 0.0;
    cfg.warmStart = true;
    cfg.contactBoundaryMinorityElectronRelaxation = false;

    const Real anodeBias = 0.828125;
    const NewtonResult result = runNewton(
        mesh,
        matdb,
        doping,
        {{"anode", anodeBias}, {"cathode", 0.0}},
        cfg);

    REQUIRE(result.iters == 0);
    requireFiniteNewtonSolution(result, mesh.numNodes());
    for (Index nid : mesh.getContact(1).node_ids) {
        const int ii = static_cast<int>(nid);
        REQUIRE(result.solution.phin(ii) == Catch::Approx(anodeBias).margin(1.0e-10));
        REQUIRE(result.solution.phip(ii) == Catch::Approx(anodeBias).margin(1.0e-10));
    }
    for (Index nid : mesh.getContact(0).node_ids) {
        const int ii = static_cast<int>(nid);
        REQUIRE(result.solution.phin(ii) == Catch::Approx(0.0).margin(1.0e-10));
        REQUIRE(result.solution.phip(ii) == Catch::Approx(0.0).margin(1.0e-10));
    }
}

TEST_CASE("NewtonSolver: carrier row convergence reports unbalanced generation rows", "[newton][carrier_row_convergence]")
{
    CoupledDDCarrierTermDiagnostic row;
    row.nodeId = 7;
    row.electronFlux = 0.0;
    row.electronRecombination = -1.0;
    row.electronImpact = 0.0;
    row.electronResidual = -1.0;

    NewtonCarrierRowConvergenceConfig cfg;
    cfg.mode = "report";
    cfg.epsRow = 1.0e-3;
    cfg.scaleFloor = 1.0e-30;

    const NewtonCarrierRowConvergenceEvaluation evaluation =
        evaluateCarrierRowConvergence({row}, cfg);

    REQUIRE(evaluation.enabled);
    REQUIRE_FALSE(evaluation.enforced);
    REQUIRE_FALSE(evaluation.satisfied);
    REQUIRE(evaluation.violations.size() == 1);
    CHECK(evaluation.violations.front().nodeId == 7);
    CHECK(evaluation.violations.front().carrier == "electron");
    CHECK(evaluation.violations.front().ratio == Catch::Approx(1.0));
}

TEST_CASE("NewtonSolver: carrier row convergence ignores balanced local rows", "[newton][carrier_row_convergence]")
{
    CoupledDDCarrierTermDiagnostic row;
    row.nodeId = 3;
    row.electronFlux = 1.0;
    row.electronRecombination = -1.0;
    row.electronResidual = 0.0;

    NewtonCarrierRowConvergenceConfig cfg;
    cfg.mode = "enforce";
    cfg.epsRow = 1.0e-3;
    cfg.scaleFloor = 1.0e-30;

    const NewtonCarrierRowConvergenceEvaluation evaluation =
        evaluateCarrierRowConvergence({row}, cfg);

    REQUIRE(evaluation.enabled);
    REQUIRE(evaluation.enforced);
    REQUIRE(evaluation.satisfied);
    CHECK(evaluation.violations.empty());
}

TEST_CASE("NewtonSolver: carrier row convergence ignores flux-only rows",
          "[newton][carrier_row_convergence]")
{
    CoupledDDCarrierTermDiagnostic row;
    row.nodeId = 4;
    row.holeFlux = 5.0e-2;
    row.holeFluxAbsSum = 1.0;
    row.holeRecombination = -1.0e-8;
    row.holeImpact = 0.0;
    row.holeResidual = 5.0e-2;

    NewtonCarrierRowConvergenceConfig cfg;
    cfg.mode = "enforce";
    cfg.epsRow = 1.0e-3;
    cfg.scaleFloor = 1.0e-30;
    cfg.minSourceScaleFraction = 1.0e-3;

    const NewtonCarrierRowConvergenceEvaluation evaluation =
        evaluateCarrierRowConvergence({row}, cfg);

    REQUIRE(evaluation.enabled);
    REQUIRE(evaluation.enforced);
    REQUIRE(evaluation.satisfied);
    CHECK(evaluation.violations.empty());
}

TEST_CASE("NewtonSolver: carrier row convergence ignores sources below absolute floor",
          "[newton][carrier_row_convergence]")
{
    CoupledDDCarrierTermDiagnostic row;
    row.nodeId = 5;
    row.electronRecombination = -1.0e-13;
    row.electronResidual = -1.0e-13;

    NewtonCarrierRowConvergenceConfig cfg;
    cfg.mode = "enforce";
    cfg.epsRow = 1.0e-3;
    cfg.scaleFloor = 1.0e-30;
    cfg.minSourceScaleFraction = 0.0;
    cfg.minSourceScale = 1.0e-12;

    const NewtonCarrierRowConvergenceEvaluation ignored =
        evaluateCarrierRowConvergence({row}, cfg);
    REQUIRE(ignored.satisfied);
    CHECK(ignored.violations.empty());

    row.electronRecombination = -1.0e-11;
    row.electronResidual = -1.0e-11;
    const NewtonCarrierRowConvergenceEvaluation enforced =
        evaluateCarrierRowConvergence({row}, cfg);
    REQUIRE_FALSE(enforced.satisfied);
    REQUIRE(enforced.violations.size() == 1);
}

TEST_CASE("NewtonSolver: carrier row convergence qualifies active density and flux rows",
          "[newton][carrier_row_convergence]")
{
    CoupledDDCarrierTermDiagnostic activeDensity;
    activeDensity.nodeId = 10;
    activeDensity.holeDensity_m3 = 1.0e16;
    activeDensity.holeResidual = 2.0e-29;

    CoupledDDCarrierTermDiagnostic activeFlux;
    activeFlux.nodeId = 11;
    activeFlux.holeDensity_m3 = 1.0;
    activeFlux.holeFluxAbsSum = 1.0e-3;
    activeFlux.holeResidual = 2.0e-5;

    CoupledDDCarrierTermDiagnostic inactive;
    inactive.nodeId = 12;
    inactive.holeDensity_m3 = 1.0;
    inactive.holeFluxAbsSum = 1.0e-9;
    inactive.holeResidual = 1.0e-9;

    NewtonCarrierRowConvergenceConfig cfg;
    cfg.mode = "enforce";
    cfg.epsRow = 1.0e-3;
    cfg.minSourceScale = 1.0;
    cfg.minCarrierDensity_m3 = 1.0e16;
    cfg.minFluxScaleFraction = 1.0e-4;

    const auto evaluation = evaluateCarrierRowConvergence(
        {activeDensity, activeFlux, inactive}, cfg);

    REQUIRE_FALSE(evaluation.satisfied);
    REQUIRE(evaluation.violations.size() == 2);
    CHECK(evaluation.qualifiedRowCount == 2);
    CHECK(evaluation.ignoredRowCount == 4);
    CHECK(evaluation.violations[0].densityQualified);
    CHECK_FALSE(evaluation.violations[0].fluxQualified);
    CHECK_FALSE(evaluation.violations[1].densityQualified);
    CHECK(evaluation.violations[1].fluxQualified);
}

TEST_CASE("NewtonSolver: parses carrier row physical qualification controls",
          "[newton][carrier_row_convergence][config]")
{
    const NewtonConfig cfg = newtonConfigFromJson(nlohmann::json{
        {"carrier_row_convergence", {
            {"mode", "report"},
            {"min_source_global_fraction", 1.0e-5},
            {"min_carrier_density_m3", 1.0e16},
            {"min_flux_scale_fraction", 1.0e-6},
            {"min_flux_scale", 1.0e-20},
        }}
    });
    CHECK(cfg.carrierRowConvergence.minSourceGlobalFraction == Catch::Approx(1.0e-5));
    CHECK(cfg.carrierRowConvergence.minCarrierDensity_m3 == Catch::Approx(1.0e16));
    CHECK(cfg.carrierRowConvergence.minFluxScaleFraction == Catch::Approx(1.0e-6));
    CHECK(cfg.carrierRowConvergence.minFluxScale == Catch::Approx(1.0e-20));

    REQUIRE_THROWS_AS(newtonConfigFromJson(nlohmann::json{
        {"carrier_row_convergence", {{"min_carrier_density_m3", -1.0}}}
    }), std::invalid_argument);
}

TEST_CASE("NewtonSolver: carrier row convergence ignores globally negligible sources",
          "[newton][carrier_row_convergence]")
{
    CoupledDDCarrierTermDiagnostic dominant;
    dominant.nodeId = 20;
    dominant.holeRecombination = 1.0;
    dominant.holeResidual = 0.0;

    CoupledDDCarrierTermDiagnostic negligible;
    negligible.nodeId = 21;
    negligible.holeRecombination = 1.0e-9;
    negligible.holeResidual = 1.0e-9;

    NewtonCarrierRowConvergenceConfig cfg;
    cfg.mode = "enforce";
    cfg.epsRow = 1.0e-3;
    cfg.minSourceScaleFraction = 0.0;
    cfg.minSourceGlobalFraction = 1.0e-6;

    const auto evaluation = evaluateCarrierRowConvergence({dominant, negligible}, cfg);
    REQUIRE(evaluation.satisfied);
    CHECK(evaluation.qualifiedRowCount == 1);
    CHECK(evaluation.ignoredRowCount == 3);
}

TEST_CASE("NewtonSolver: global continuity closure compares contact flux with free-node source",
          "[newton][global_continuity_closure]")
{
    std::vector<CoupledDDCarrierTermDiagnostic> rows(3);
    rows[0].nodeId = 0;
    rows[0].electronFlux = -2.0;
    rows[0].holeFlux = -3.0;
    rows[1].nodeId = 1;
    rows[1].electronRecombination = -2.0;
    rows[1].holeRecombination = -3.0;
    rows[2].nodeId = 2;

    NewtonGlobalContinuityClosureConfig cfg;
    cfg.mode = "enforce";
    cfg.tolerance = 1.0e-3;
    cfg.sourceFloor = 1.0e-12;

    const auto balanced =
        evaluateGlobalContinuityClosure(rows, {0, 2}, {0, 2}, cfg);
    REQUIRE(balanced.enabled);
    REQUIRE(balanced.enforced);
    REQUIRE(balanced.satisfied);
    CHECK(balanced.electron.ratio == Catch::Approx(0.0));
    CHECK(balanced.hole.ratio == Catch::Approx(0.0));

    rows[0].holeFlux = -1.0;
    const auto unbalanced =
        evaluateGlobalContinuityClosure(rows, {0, 2}, {0, 2}, cfg);
    REQUIRE_FALSE(unbalanced.satisfied);
    CHECK(unbalanced.electron.ratio == Catch::Approx(0.0));
    CHECK(unbalanced.hole.ratio == Catch::Approx(2.0 / 3.0));
}

TEST_CASE("NewtonSolver: Gummel density recovery jumps a dead carrier row",
          "[newton][carrier_row_recovery]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.recombination = {"srh"};
    cfg.mobility.model = "constant";

    const int n = static_cast<int>(mesh.numNodes());
    DDSolution dead;
    dead.psi = VectorXd::Zero(n);
    dead.phin = VectorXd::Zero(n);
    dead.phip = VectorXd::Zero(n);
    dead.n = VectorXd::Constant(n, 1.0e4);
    dead.p = VectorXd::Constant(n, 1.0e4);
    dead.phin(4) = 8.0;
    dead.n(4) = 1.0e-120;

    NewtonCarrierRowConvergenceViolation violation;
    violation.nodeId = 4;
    violation.carrier = "electron";
    violation.residual = -1.0;
    violation.scale = 1.0;
    violation.ratio = 1.0;

    NewtonCarrierRowRecoveryConfig recovery;
    recovery.mode = "gummel_density";

    const NewtonCarrierRowRecoveryResult recovered =
        recoverCarrierRowsWithGummelDensity(
            mesh, matdb, doping, zeroBias(), cfg, dead, {violation}, recovery);

    REQUIRE(recovered.attempted);
    CHECK(recovered.electronRowsUpdated >= 1);
    CHECK(recovered.solution.n(4) > 1.0e8);
    CHECK(recovered.solution.phin(4) < 1.0);
    CHECK(recovered.maxPsiDelta_V == Catch::Approx(0.0));
}

TEST_CASE("NewtonSolver: Fermi-Dirac Gummel density recovery is available",
          "[newton][carrier_row_recovery][fermi_dirac]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.recombination = {"srh"};
    cfg.mobility.model = "constant";
    cfg.carrierStatistics.model = "fermi_dirac";

    const int n = static_cast<int>(mesh.numNodes());
    DDSolution dead;
    dead.psi = VectorXd::Zero(n);
    dead.phin = VectorXd::Zero(n);
    dead.phip = VectorXd::Zero(n);
    dead.n = VectorXd::Constant(n, 1.0e4);
    dead.p = VectorXd::Constant(n, 1.0e4);
    dead.phin(4) = 8.0;
    dead.n(4) = 1.0e-120;

    NewtonCarrierRowConvergenceViolation violation;
    violation.nodeId = 4;
    violation.carrier = "electron";
    violation.residual = -1.0;
    violation.scale = 1.0;
    violation.ratio = 1.0;

    NewtonCarrierRowRecoveryConfig recovery;
    recovery.mode = "gummel_density";

    const NewtonCarrierRowRecoveryResult recovered =
        recoverCarrierRowsWithGummelDensity(
            mesh, matdb, doping, zeroBias(), cfg, dead, {violation}, recovery);

    REQUIRE(recovered.attempted);
    CHECK(recovered.electronRowsUpdated >= 1);
    CHECK(std::isfinite(recovered.solution.n(4)));
    CHECK(std::isfinite(recovered.solution.phin(4)));
    CHECK(recovered.solution.n(4) > 1.0e8);
    CHECK(recovered.solution.phin(4) < 1.0);
    CHECK(recovered.maxPsiDelta_V == Catch::Approx(0.0));
}

TEST_CASE("NewtonSolver: initial abstol cannot accept an unbalanced carrier row",
          "[newton][carrier_row_convergence][carrier_row_recovery][initial_abstol]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.recombination = {"srh"};
    cfg.mobility.model = "constant";
    cfg.warmStart = true;
    cfg.abstol = 1.0e100;
    cfg.reltol = 0.0;
    cfg.maxIter = 0;
    cfg.verbose = false;
    cfg.carrierRowConvergence.mode = "enforce";
    cfg.carrierRowConvergence.epsRow = 1.0e-3;
    cfg.carrierRowConvergence.minSourceScaleFraction = 0.0;
    cfg.carrierRowConvergence.minSourceScale = 1.0e-30;
    cfg.carrierRowRecovery.mode = "gummel_density";
    cfg.carrierRowRecovery.maxAttempts = 1;
    cfg.carrierRowRecovery.maxCycles = 1;

    const int n = static_cast<int>(mesh.numNodes());
    DDSolution dead;
    dead.psi = VectorXd::Zero(n);
    dead.phin = VectorXd::Zero(n);
    dead.phip = VectorXd::Zero(n);
    dead.phin(4) = 8.0;
    dead.n = VectorXd::Zero(n);
    dead.p = VectorXd::Zero(n);

    const NewtonResult result = runNewton(
        mesh, matdb, doping, zeroBias(), dead, cfg);

    REQUIRE(result.carrierRowRecovery.attempted);
    CHECK(result.carrierRowRecovery.cyclesAttempted >= 1);
    CHECK(result.convergenceReason != "initial_abstol");
}

TEST_CASE("NewtonSolver: Gummel density recovery uses Ohmic contact densities",
          "[newton][carrier_row_recovery][contact_bc]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    NewtonConfig cfg = newtonConfig();
    cfg.recombination = {"srh"};
    cfg.mobility.model = "constant";

    const int n = static_cast<int>(mesh.numNodes());
    DDSolution clean;
    clean.psi = VectorXd::Zero(n);
    clean.phin = VectorXd::Zero(n);
    clean.phip = VectorXd::Zero(n);
    clean.n = VectorXd::Constant(n, 1.0e4);
    clean.p = VectorXd::Constant(n, 1.0e4);
    clean.phin(4) = 8.0;
    clean.n(4) = 1.0e-120;

    DDSolution polluted = clean;
    polluted.n(0) = 1.0e60;
    polluted.p(0) = 1.0e-60;
    polluted.n(1) = 1.0e-60;
    polluted.p(1) = 1.0e60;
    polluted.n(2) = 1.0e-60;
    polluted.p(2) = 1.0e60;
    polluted.n(3) = 1.0e60;
    polluted.p(3) = 1.0e-60;

    NewtonCarrierRowConvergenceViolation violation;
    violation.nodeId = 4;
    violation.carrier = "electron";
    violation.residual = -1.0;
    violation.scale = 1.0;
    violation.ratio = 1.0;

    NewtonCarrierRowRecoveryConfig recovery;
    recovery.mode = "gummel_density";

    const NewtonCarrierRowRecoveryResult fromClean =
        recoverCarrierRowsWithGummelDensity(
            mesh, matdb, doping, zeroBias(), cfg, clean, {violation}, recovery);
    const NewtonCarrierRowRecoveryResult fromPolluted =
        recoverCarrierRowsWithGummelDensity(
            mesh, matdb, doping, zeroBias(), cfg, polluted, {violation}, recovery);

    REQUIRE(fromClean.attempted);
    REQUIRE(fromPolluted.attempted);
    CHECK(fromPolluted.solution.n(4) == Catch::Approx(fromClean.solution.n(4)).epsilon(1.0e-12));
}

TEST_CASE("NewtonConfig unit_scaling default Auger coefficients are TCAD internal",
          "[newton][recombination][auger][scaling]")
{
    const nlohmann::json json = {
        {"recombination", {"srh", "auger"}}
    };

    const NewtonConfig legacy = newtonConfigFromJson(json);
    REQUIRE(legacy.augerCn == Catch::Approx(2.90e-43));
    REQUIRE(legacy.augerCp == Catch::Approx(1.028e-43));

    const NewtonConfig scaled = newtonConfigFromJson(
        json, UnitScalingConfig{UnitScalingMode::UnitScaling});
    REQUIRE(scaled.augerCn == Catch::Approx(2.90e-31));
    REQUIRE(scaled.augerCp == Catch::Approx(1.028e-31));
}
