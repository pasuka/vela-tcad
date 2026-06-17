#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>

#include "vela/material/MaterialDatabase.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/solver/GummelSolver.h"
#include "vela/solver/NewtonSolver.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <vector>

using namespace vela;

static DeviceMesh makePNMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;
    Node n0; n0.id=0; n0.x=0; n0.y=0; mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=L; n1.y=0; mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=L; n2.y=L; mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=0; n3.y=L; mesh.addNode(n3);
    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0; c0.node_ids={0,1,2}; mesh.addCell(c0);
    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=1; c1.node_ids={0,2,3}; mesh.addCell(c1);
    Region r0; r0.id=0; r0.name="n_region"; r0.material="Si"; r0.cell_ids={0}; mesh.addRegion(r0);
    Region r1; r1.id=1; r1.name="p_region"; r1.material="Si"; r1.cell_ids={1}; mesh.addRegion(r1);
    Contact anode; anode.id=0; anode.name="anode"; anode.region_id=1; anode.node_ids={0,3}; mesh.addContact(anode);
    Contact cathode; cathode.id=1; cathode.name="cathode"; cathode.region_id=0; cathode.node_ids={1,2}; mesh.addContact(cathode);
    mesh.buildEdges();
    return mesh;
}

TEST_CASE("Impact ionization none model is zero", "[impact]")
{
    const auto model = makeImpactIonizationModel(impactIonizationModelConfig("none"));
    REQUIRE(model->electronCoefficient(1.0e8) == Catch::Approx(0.0));
    REQUIRE(model->holeCoefficient(1.0e8) == Catch::Approx(0.0));
    REQUIRE(model->generationRate(1.0e8, 1.0e21, 1.0e21) == Catch::Approx(0.0));
}

TEST_CASE("Selberherr impact ionization grows with electric field", "[impact]")
{
    SelberherrImpactIonization model;
    const Real low = model.electronCoefficient(1.0e7);
    const Real high = model.electronCoefficient(5.0e8);
    REQUIRE(low >= 0.0);
    REQUIRE(high > low);
    REQUIRE(model.generationRate(5.0e8, 1.0e20, 2.0e20) > 0.0);
}

TEST_CASE("Van Overstraeten impact ionization matches Sentaurus 2018 silicon defaults",
          "[impact][van_overstraeten]")
{
    const auto model = makeImpactIonizationModel(
        impactIonizationModelConfig("van_overstraeten"));

    const Real lowField = 2.0e7;  // 2e5 V/cm, below E0.
    const Real highField = 5.0e7; // 5e5 V/cm, above E0.

    const Real expectedElectronLow = 7.03e7 * std::exp(-1.231e8 / lowField);
    const Real expectedHoleLow = 1.582e8 * std::exp(-2.036e8 / lowField);
    const Real expectedElectronHigh = 7.03e7 * std::exp(-1.231e8 / highField);
    const Real expectedHoleHigh = 6.71e7 * std::exp(-1.693e8 / highField);

    REQUIRE(model->electronCoefficient(lowField) ==
            Catch::Approx(expectedElectronLow).epsilon(1.0e-12));
    REQUIRE(model->holeCoefficient(lowField) ==
            Catch::Approx(expectedHoleLow).epsilon(1.0e-12));
    REQUIRE(model->electronCoefficient(highField) ==
            Catch::Approx(expectedElectronHigh).epsilon(1.0e-12));
    REQUIRE(model->holeCoefficient(highField) ==
            Catch::Approx(expectedHoleHigh).epsilon(1.0e-12));
}

TEST_CASE("Gummel reverse bias BV regression runs with impact ionization", "[impact][gummel]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);
    const std::unordered_map<std::string, Real> biases = {{"anode", -1.0}, {"cathode", 0.0}};

    GummelConfig cfg;
    cfg.maxIter = 20;
    cfg.reltol = 1.0e-4;
    cfg.abstol = 1.0e12;
    cfg.dampingPsi = 0.3;
    cfg.mobility = mobilityModelConfig("caughey_thomas_field");
    cfg.impactIonization.model = "selberherr";

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));
    REQUIRE(sol.iters >= 1);
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(std::isfinite(sol.psi(i)));
        REQUIRE(std::isfinite(sol.n(i)));
        REQUIRE(std::isfinite(sol.p(i)));
        REQUIRE(sol.n(i) >= 0.0);
        REQUIRE(sol.p(i) >= 0.0);
    }
}


TEST_CASE("Coupled DD residual includes impact-ionization generation", "[impact][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    CoupledDDState state;
    state.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phip = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.psi(1) = 1.0;
    state.psi(2) = 1.0;

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const RecombinationModelConfig recombinationConfig = recombinationModelConfig({"none"});

    CoupledDDAssembler noImpact(
        mesh, matdb, doping, Vt, mobilityConfig, recombinationConfig);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0;
    impactConfig.carrierVelocity = 1.0;
    CoupledDDAssembler withImpact(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityConfig,
        recombinationConfig,
        BandgapNarrowingConfig{},
        impactConfig);

    const VectorXd x = noImpact.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const VectorXd r0 = noImpact.residual(x, bcs);
    const VectorXd r1 = withImpact.residual(x, bcs);

    const int phinOffset = static_cast<int>(mesh.numNodes());
    const int phipOffset = 2 * static_cast<int>(mesh.numNodes());
    bool sawGeneration = false;
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        const Real electronDelta = r1(phinOffset + i) - r0(phinOffset + i);
        const Real holeDelta = r1(phipOffset + i) - r0(phipOffset + i);
        REQUIRE(electronDelta <= 0.0);
        REQUIRE(holeDelta <= 0.0);
        sawGeneration = sawGeneration || electronDelta < 0.0 || holeDelta < 0.0;
    }
    REQUIRE(sawGeneration);
}

TEST_CASE("Quasi-Fermi avalanche driving force ignores built-in electrostatic field",
          "[impact][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    CoupledDDState state;
    state.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phip = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.psi(1) = 1.0;
    state.psi(2) = 1.0;

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const RecombinationModelConfig recombinationConfig = recombinationModelConfig({"none"});
    CoupledDDAssembler noImpact(
        mesh, matdb, doping, Vt, mobilityConfig, recombinationConfig);

    ImpactIonizationModelConfig qfImpact;
    qfImpact.model = "selberherr";
    qfImpact.drivingForce = "quasi_fermi_gradient";
    qfImpact.generation = "current_density";
    qfImpact.electronA = 1.0;
    qfImpact.electronB = 1.0;
    qfImpact.holeA = 1.0;
    qfImpact.holeB = 1.0;
    CoupledDDAssembler withQuasiFermiImpact(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityConfig,
        recombinationConfig,
        BandgapNarrowingConfig{},
        qfImpact);

    const VectorXd x = noImpact.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const VectorXd r0 = noImpact.residual(x, bcs);
    const VectorXd r1 = withQuasiFermiImpact.residual(x, bcs);

    const int phinOffset = static_cast<int>(mesh.numNodes());
    const int phipOffset = 2 * static_cast<int>(mesh.numNodes());
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(r1(phinOffset + i) == Catch::Approx(r0(phinOffset + i)).margin(1.0e-18));
        REQUIRE(r1(phipOffset + i) == Catch::Approx(r0(phipOffset + i)).margin(1.0e-18));
    }
}

TEST_CASE("Coupled DD analytic avalanche Jacobian matches carrier finite differences", "[impact][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    CoupledDDState state;
    state.psi = VectorXd::Zero(static_cast<int>(mesh.numNodes()));
    state.phin = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), 0.01, -0.005);
    state.phip = VectorXd::LinSpaced(static_cast<int>(mesh.numNodes()), -0.008, 0.006);
    state.psi(1) = 0.1;
    state.psi(2) = 0.1;

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0;
    impactConfig.carrierVelocity = 1.0;

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityModelConfig("constant"),
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        impactConfig);

    const VectorXd x = assembler.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const SparseMatrixd analytic = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd finiteDifference = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd denseAnalytic = Eigen::MatrixXd(analytic);
    const Eigen::MatrixXd denseFiniteDifference = Eigen::MatrixXd(finiteDifference);

    const int N = static_cast<int>(mesh.numNodes());
    Real maxAbsDiff = 0.0;
    Real maxAbsRef = 0.0;
    for (int row = 0; row < 3 * N; ++row) {
        for (int col = N; col < 3 * N; ++col) {
            maxAbsDiff = std::max(
                maxAbsDiff,
                std::abs(denseAnalytic(row, col) - denseFiniteDifference(row, col)));
            maxAbsRef = std::max(maxAbsRef, std::abs(denseFiniteDifference(row, col)));
        }
    }

    REQUIRE(maxAbsDiff / std::max<Real>(1.0, maxAbsRef) < 5.0e-5);
}

TEST_CASE("JSON solver config selects impact ionization model", "[impact][json]")
{
    const GummelConfig cfg = gummelConfigFromJson(nlohmann::json{
        {"impact_ionization", {{"model", "selberherr"}, {"electron_A_m_inv", 1.0e6}}}
    });
    REQUIRE(cfg.impactIonization.model == "selberherr");
    REQUIRE(cfg.impactIonization.electronA == Catch::Approx(1.0e6));

    const NewtonConfig stringCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", "selberherr"}
    });
    REQUIRE(stringCfg.impactIonization.model == "selberherr");

    const NewtonConfig objectCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "selberherr"},
            {"electron_A_m_inv", 2.0e6},
            {"electron_B_V_m", 3.0e7},
            {"hole_A_m_inv", 4.0e6},
            {"hole_B_V_m", 5.0e7},
            {"carrier_velocity_m_s", 6.0e4},
        }}
    });
    REQUIRE(objectCfg.impactIonization.model == "selberherr");
    REQUIRE(objectCfg.impactIonization.electronA == Catch::Approx(2.0e6));
    REQUIRE(objectCfg.impactIonization.electronB == Catch::Approx(3.0e7));
    REQUIRE(objectCfg.impactIonization.holeA == Catch::Approx(4.0e6));
    REQUIRE(objectCfg.impactIonization.holeB == Catch::Approx(5.0e7));
    REQUIRE(objectCfg.impactIonization.carrierVelocity == Catch::Approx(6.0e4));

    const NewtonConfig vanOverstraetenCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"electron_a_low_m_inv", 1.0e6},
            {"electron_b_low_V_m", 2.0e7},
            {"hole_a_high_m_inv", 3.0e6},
            {"hole_b_high_V_m", 4.0e7},
            {"switch_field_V_m", 5.0e7},
            {"phonon_energy_eV", 0.063},
            {"temperature_K", 300.0},
        }}
    });
    REQUIRE(vanOverstraetenCfg.impactIonization.model == "van_overstraeten");
    REQUIRE(vanOverstraetenCfg.impactIonization.electronALow == Catch::Approx(1.0e6));
    REQUIRE(vanOverstraetenCfg.impactIonization.electronBLow == Catch::Approx(2.0e7));
    REQUIRE(vanOverstraetenCfg.impactIonization.holeAHigh == Catch::Approx(3.0e6));
    REQUIRE(vanOverstraetenCfg.impactIonization.holeBHigh == Catch::Approx(4.0e7));
    REQUIRE(vanOverstraetenCfg.impactIonization.switchField == Catch::Approx(5.0e7));
    REQUIRE(vanOverstraetenCfg.impactIonization.phononEnergy == Catch::Approx(0.063));
    REQUIRE(vanOverstraetenCfg.impactIonization.temperature_K == Catch::Approx(300.0));

    const NewtonConfig sentaurusCfg = newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "quasi_fermi_gradient"},
            {"generation", "current_density"},
        }}
    });
    REQUIRE(sentaurusCfg.impactIonization.model == "van_overstraeten");
    REQUIRE(sentaurusCfg.impactIonization.drivingForce == "quasi_fermi_gradient");
    REQUIRE(sentaurusCfg.impactIonization.generation == "current_density");

    REQUIRE_THROWS_AS(newtonConfigFromJson(nlohmann::json{
        {"impact_ionization", {
            {"model", "van_overstraeten"},
            {"driving_force", "electrostatic"},
        }}
    }), std::invalid_argument);
}
