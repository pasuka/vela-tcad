#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/GummelSolver.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/physics/RecombinationModel.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using namespace vela;

// ---------------------------------------------------------------------------
// Build a simple PN junction mesh (same geometry as test_poisson.cpp)
//
//   3 -------- 2
//   | p-reg  / |
//   | (T1)  /  |
//   |      /   |
//   |  T0 /    |
//   | n-reg    |
//   0 -------- 1
//
//  Nodes: 0=(0,0), 1=(L,0), 2=(L,L), 3=(0,L)   L = 1 um
//  Cells: T0={0,1,2} n-region,  T1={0,2,3} p-region
//  Contacts:
//    cathode (n): nodes 1, 2   V = 0
//    anode   (p): nodes 0, 3   V = 0
// ---------------------------------------------------------------------------

static DeviceMesh makePNMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;

    Node n0; n0.id=0; n0.x=0;  n0.y=0;  mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=L;  n1.y=0;  mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=L;  n2.y=L;  mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=0;  n3.y=L;  mesh.addNode(n3);

    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0;
    c0.node_ids = {0, 1, 2};  mesh.addCell(c0);

    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=1;
    c1.node_ids = {0, 2, 3};  mesh.addCell(c1);

    Region r0; r0.id=0; r0.name="n_region"; r0.material="Si"; r0.cell_ids={0};
    mesh.addRegion(r0);
    Region r1; r1.id=1; r1.name="p_region"; r1.material="Si"; r1.cell_ids={1};
    mesh.addRegion(r1);

    Contact anode;   anode.id=0;   anode.name="anode";
    anode.region_id=1; anode.node_ids={0,3};
    mesh.addContact(anode);

    Contact cathode; cathode.id=1; cathode.name="cathode";
    cathode.region_id=0; cathode.node_ids={1,2};
    mesh.addContact(cathode);

    mesh.buildEdges();
    return mesh;
}

static DopingModel makePNDoping(const DeviceMesh& mesh)
{
    std::vector<RegionDopingSpec> specs = {
        { "n_region", 1.0e23, 0.0   },
        { "p_region", 0.0,    1.0e23 }
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

static std::vector<Real> readVtkScalar(const std::filesystem::path& path,
                                       const std::string& name,
                                       std::size_t count)
{
    std::ifstream input(path);
    REQUIRE(input.is_open());
    std::string line;
    const std::string marker = "SCALARS " + name + " ";
    while (std::getline(input, line)) {
        if (line.rfind(marker, 0) != 0)
            continue;
        REQUIRE(std::getline(input, line));
        REQUIRE(line == "LOOKUP_TABLE default");
        std::vector<Real> values;
        values.reserve(count);
        Real value = 0.0;
        while (values.size() < count && input >> value)
            values.push_back(value);
        REQUIRE(values.size() == count);
        return values;
    }
    FAIL("VTK scalar not found: " + name);
    return {};
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("GummelSolver: equilibrium (0 V bias) does not crash", "[gummel]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    std::unordered_map<std::string, Real> biases = {
        {"anode",   0.0},
        {"cathode", 0.0}
    };

    GummelConfig cfg;
    cfg.maxIter = 30;
    cfg.reltol  = 1.0e-6;

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));
    REQUIRE(sol.iters >= 1);
}

TEST_CASE("GummelSolver: Fermi-Dirac density block preserves finite equilibrium state",
          "[gummel][fermi_dirac]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    const std::unordered_map<std::string, Real> biases = {
        {"anode", 0.0},
        {"cathode", 0.0}
    };

    GummelConfig cfg;
    cfg.maxIter = 80;
    cfg.reltol = 1.0e-7;
    cfg.dampingPsi = 0.5;
    cfg.mobility.model = "constant";
    cfg.recombination = {"none"};
    cfg.carrierStatistics.model = "fermi_dirac";

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));
    REQUIRE(sol.converged);

    const Material silicon = matdb.getMaterial("Si", cfg.temperature_K);
    REQUIRE(silicon.Nc_m3.has_value());
    REQUIRE(silicon.Nv_m3.has_value());
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(std::isfinite(sol.psi(i)));
        REQUIRE(std::isfinite(sol.phin(i)));
        REQUIRE(std::isfinite(sol.phip(i)));
        REQUIRE(sol.n(i) > 0.0);
        REQUIRE(sol.p(i) > 0.0);
        const Real reconstructedN = electronDensity(
            silicon.ni, *silicon.Nc_m3, sol.psi(i), sol.phin(i),
            constants::Vt_300, cfg.carrierStatistics);
        const Real reconstructedP = holeDensity(
            silicon.ni, *silicon.Nv_m3, sol.psi(i), sol.phip(i),
            constants::Vt_300, cfg.carrierStatistics);
        CHECK(reconstructedN == Catch::Approx(sol.n(i)).epsilon(2.0e-10));
        CHECK(reconstructedP == Catch::Approx(sol.p(i)).epsilon(2.0e-10));
    }
}

TEST_CASE("GummelSolver: n and p are strictly positive", "[gummel]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    std::unordered_map<std::string, Real> biases = {
        {"anode",   0.0},
        {"cathode", 0.0}
    };

    DDSolution sol = runGummel(mesh, matdb, doping, biases);

    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(sol.n(i) > 0.0);
        REQUIRE(sol.p(i) > 0.0);
    }
}

TEST_CASE("GummelSolver: carrier floor config preserves quasi-Fermi consistency", "[gummel]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    std::unordered_map<std::string, Real> biases = {
        {"anode",   0.0},
        {"cathode", 0.0}
    };

    GummelConfig cfg;
    cfg.maxIter = 1;
    cfg.carrierFloor = 1.0;
    DDSolution sol = runGummel(mesh, matdb, doping, biases, cfg);

    const Real ni = matdb.getMaterial("Si").ni;
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        const Real n = electronDensity(ni, sol.psi(i), sol.phin(i), constants::Vt_300);
        const Real p = holeDensity(ni, sol.psi(i), sol.phip(i), constants::Vt_300);
        REQUIRE(sol.n(i) >= cfg.carrierFloor);
        REQUIRE(sol.p(i) >= cfg.carrierFloor);
        REQUIRE(n == Catch::Approx(sol.n(i)).epsilon(1.0e-12));
        REQUIRE(p == Catch::Approx(sol.p(i)).epsilon(1.0e-12));
    }

    const GummelConfig parsed =
        gummelConfigFromJson(nlohmann::json{{"carrier_floor_m3", 2.0}});
    REQUIRE(parsed.carrierFloor == Catch::Approx(2.0));
    REQUIRE_THROWS_AS(
        gummelConfigFromJson(nlohmann::json{{"carrier_floor_m3", -1.0}}),
        std::invalid_argument);
}

TEST_CASE("GummelSolver: no NaN or Inf in any output vector", "[gummel]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    std::unordered_map<std::string, Real> biases = {
        {"anode",   0.0},
        {"cathode", 0.0}
    };

    DDSolution sol = runGummel(mesh, matdb, doping, biases);

    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(std::isfinite(sol.psi (i)));
        REQUIRE(std::isfinite(sol.phin(i)));
        REQUIRE(std::isfinite(sol.phip(i)));
        REQUIRE(std::isfinite(sol.n   (i)));
        REQUIRE(std::isfinite(sol.p   (i)));
    }
}

TEST_CASE("GummelSolver: VTK output is written successfully", "[gummel][vtk]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    std::unordered_map<std::string, Real> biases = {
        {"anode",   0.0},
        {"cathode", 0.0}
    };

    DDSolution sol = runGummel(mesh, matdb, doping, biases);

    const std::string vtkPath =
        (std::filesystem::temp_directory_path() / "test_dd_gummel.vtk").string();

    const MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const RecombinationModelConfig recombinationConfig = recombinationModelConfig({"none"});
    const ImpactIonizationModelConfig impactConfig;
    const BandgapNarrowingConfig bandgapNarrowingConfig;
    REQUIRE_NOTHROW(writeDDSolutionVTK(
        vtkPath,
        mesh,
        matdb,
        doping,
        sol,
        mobilityConfig,
        recombinationConfig,
        impactConfig,
        bandgapNarrowingConfig,
        constants::T0));

    REQUIRE(std::filesystem::exists(vtkPath));
    REQUIRE(std::filesystem::file_size(vtkPath) > 0);

    std::ifstream ifs(vtkPath);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());

    REQUIRE(content.find("Potential")          != std::string::npos);
    REQUIRE(content.find("Electrons")          != std::string::npos);
    REQUIRE(content.find("Holes")              != std::string::npos);
    REQUIRE(content.find("NetDoping")          != std::string::npos);
    REQUIRE(content.find("ElectronQuasiFermi") != std::string::npos);
    REQUIRE(content.find("HoleQuasiFermi")     != std::string::npos);
    REQUIRE(content.find("VECTORS J_n_drift double") != std::string::npos);
    REQUIRE(content.find("VECTORS J_n_diffusion double") != std::string::npos);
    REQUIRE(content.find("VECTORS J_n_total double") != std::string::npos);
    REQUIRE(content.find("VECTORS J_p_drift double") != std::string::npos);
    REQUIRE(content.find("VECTORS J_p_diffusion double") != std::string::npos);
    REQUIRE(content.find("VECTORS J_p_total double") != std::string::npos);
    REQUIRE(content.find("VECTORS TotalCurrentDensityVector double") != std::string::npos);
}

TEST_CASE("GummelSolver: VTK recombination uses Fermi-corrected effective intrinsic density",
          "[gummel][vtk][fermi_dirac][bgn]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const int count = static_cast<int>(mesh.numNodes());
    DDSolution sol;
    sol.psi = VectorXd::Zero(count);
    sol.phin = VectorXd::Zero(count);
    sol.phip = VectorXd::Constant(count, 0.1);
    sol.n = VectorXd::Constant(count, 1.0e22);
    sol.p = VectorXd::Constant(count, 2.0e21);

    BandgapNarrowingConfig bgn = bandgapNarrowingConfig("old_slotboom");
    bgn.fermiStatisticsCorrection = true;
    CarrierStatisticsConfig statistics;
    statistics.model = "fermi_dirac";
    RecombinationModelConfig recombination = recombinationModelConfig({"srh"});
    const std::filesystem::path vtkPath =
        std::filesystem::temp_directory_path() / "test_dd_gummel_fermi_bgn_srh.vtk";
    writeDDSolutionVTK(
        vtkPath.string(), mesh, matdb, doping, sol,
        mobilityModelConfig("constant"), recombination,
        ImpactIonizationModelConfig{}, bgn, constants::T0,
        UnitScalingConfig{}, statistics);

    const auto effectiveNi = detail::buildValidatedEffectiveNodeNi(
        "test", mesh, matdb, doping, bgn, constants::Vt_300);
    const auto Nc = detail::buildNodeDensityOfStates(
        mesh, matdb, constants::T0, true);
    const auto Nv = detail::buildNodeDensityOfStates(
        mesh, matdb, constants::T0, false);
    const RecombinationModel model(recombination);
    const GeneralizedSrhCarrierState state = generalizedSrhCarrierState(
        sol.n(0), sol.p(0), effectiveNi[0], Nc[0], Nv[0],
        sol.phip(0) - sol.phin(0), constants::Vt_300, statistics);
    const Real expected = model.srhRateGeneralizedFromExcessProduct(
        state.excessProduct, sol.n(0), sol.p(0), effectiveNi[0], effectiveNi[0],
        state.electronDegeneracy, state.holeDegeneracy,
        model.srhDopingConcentration(doping.donors(0), doping.acceptors(0))) / 1.0e6;

    const std::vector<Real> exportedNi = readVtkScalar(
        vtkPath, "EffectiveIntrinsicDensity", mesh.numNodes());
    const std::vector<Real> exportedSrh = readVtkScalar(
        vtkPath, "SRHRecombinationCm3PerS", mesh.numNodes());
    REQUIRE(exportedNi[0] == Catch::Approx(effectiveNi[0]).epsilon(1.0e-12));
    REQUIRE(exportedSrh[0] == Catch::Approx(expected).epsilon(1.0e-12));
}

TEST_CASE("GummelSolver: forward bias converges without crash", "[gummel]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    // Small forward bias: anode = 0.3 V, cathode = 0 V
    std::unordered_map<std::string, Real> biases = {
        {"anode",   0.3},
        {"cathode", 0.0}
    };

    GummelConfig cfg;
    cfg.maxIter    = 50;
    cfg.reltol     = 1.0e-5;
    cfg.dampingPsi = 0.5;

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));

    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(std::isfinite(sol.n(i)));
        REQUIRE(std::isfinite(sol.p(i)));
        REQUIRE(sol.n(i) >= 0.0);
        REQUIRE(sol.p(i) >= 0.0);
    }
}

TEST_CASE("GummelSolver: abstol can terminate strongly damped high-doping updates", "[gummel]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    std::vector<RegionDopingSpec> specs = {
        {"n_region", 1.0e24, 0.0},
        {"p_region", 0.0, 1.0e24}
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    std::unordered_map<std::string, Real> biases = {
        {"anode", 0.25},
        {"cathode", 0.0}
    };

    GummelConfig noAbs;
    noAbs.maxIter = 2;
    noAbs.reltol = 0.0;
    noAbs.abstol = 0.0;
    noAbs.dampingPsi = 0.05;

    const DDSolution exhausted = runGummel(mesh, matdb, doping, biases, noAbs);
    REQUIRE_FALSE(exhausted.converged);
    REQUIRE(exhausted.iters == noAbs.maxIter);

    GummelConfig withAbs = noAbs;
    withAbs.abstol = 1.0e40;

    const DDSolution converged = runGummel(mesh, matdb, doping, biases, withAbs);
    REQUIRE(converged.converged);
    REQUIRE(converged.iters == 1);

    const GummelConfig parsed = gummelConfigFromJson(nlohmann::json{{"abstol", 1.0e-12}});
    REQUIRE(parsed.abstol == Catch::Approx(1.0e-12));
}

TEST_CASE("GummelSolver: configured temperature scales ohmic built-in potential", "[gummel][temperature]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const std::unordered_map<std::string, Real> biases = {
        {"anode", 0.0},
        {"cathode", 0.0}
    };

    GummelConfig cfg300;
    cfg300.maxIter = 1;
    cfg300.temperature_K = 300.0;
    const DDSolution sol300 = runGummel(mesh, matdb, doping, biases, cfg300);

    GummelConfig cfg600 = cfg300;
    cfg600.temperature_K = 600.0;
    const DDSolution sol600 = runGummel(mesh, matdb, doping, biases, cfg600);

    const Real builtIn300 = sol300.psi(1) - sol300.psi(0);
    const Real builtIn600 = sol600.psi(1) - sol600.psi(0);
    REQUIRE(builtIn600 > 0.0);
    REQUIRE(builtIn600 < builtIn300);

    const GummelConfig parsed = gummelConfigFromJson(nlohmann::json{{"temperature_K", 325.0}});
    REQUIRE(parsed.temperature_K == Catch::Approx(325.0));
    REQUIRE_THROWS_AS(gummelConfigFromJson(nlohmann::json{{"temperature_K", 0.0}}),
                      std::invalid_argument);
}

TEST_CASE("Gummel solver parses Sentaurus E2 band-to-band parameters",
          "[gummel][btbt][json]")
{
    const GummelConfig cfg = gummelConfigFromJson(nlohmann::json{
        {"band_to_band", {
            {"model", "e2"},
            {"A_m_inv_s_inv_V_inv2", 3.4e23},
            {"B_V_per_m", 2.26e9},
        }},
    });
    REQUIRE(cfg.bandToBand.model == "e2");
    REQUIRE(cfg.bandToBand.prefactorA_SI == Catch::Approx(3.4e23));
    REQUIRE(cfg.bandToBand.exponentialB_V_per_m == Catch::Approx(2.26e9));
    REQUIRE(cfg.bandToBand.jacobian == "frozen_field");
}

static DeviceMesh makeFourTerminalSiliconMesh()
{
    DeviceMesh mesh;
    const double Lx = 1.0e-6;
    const double Ly = 4.0e-7;
    const double xmid = 0.5e-6;

    Node n0; n0.id=0; n0.x=0.0;  n0.y=0.0;  mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=xmid; n1.y=0.0;  mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=Lx;   n2.y=0.0;  mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=0.0;  n3.y=Ly;   mesh.addNode(n3);
    Node n4; n4.id=4; n4.x=xmid; n4.y=Ly;   mesh.addNode(n4);
    Node n5; n5.id=5; n5.x=Lx;   n5.y=Ly;   mesh.addNode(n5);

    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0; c0.node_ids={0,1,4}; mesh.addCell(c0);
    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=0; c1.node_ids={0,4,3}; mesh.addCell(c1);
    Cell c2; c2.id=2; c2.type=CellType::Tri3; c2.region_id=1; c2.node_ids={1,2,5}; mesh.addCell(c2);
    Cell c3; c3.id=3; c3.type=CellType::Tri3; c3.region_id=1; c3.node_ids={1,5,4}; mesh.addCell(c3);

    Region r0; r0.id=0; r0.name="p_channel"; r0.material="Si"; r0.cell_ids={0,1}; mesh.addRegion(r0);
    Region r1; r1.id=1; r1.name="n_drain"; r1.material="Si"; r1.cell_ids={2,3}; mesh.addRegion(r1);

    Contact body; body.id=0; body.name="body"; body.region_id=0; body.node_ids={0,1}; mesh.addContact(body);
    Contact source; source.id=1; source.name="source"; source.region_id=0; source.node_ids={3}; mesh.addContact(source);
    Contact gate; gate.id=2; gate.name="gate"; gate.region_id=0; gate.node_ids={4}; mesh.addContact(gate);
    Contact drain; drain.id=3; drain.name="drain"; drain.region_id=1; drain.node_ids={2,5}; mesh.addContact(drain);

    mesh.buildEdges();
    return mesh;
}

static void requireFiniteDDSolution(const DDSolution& sol, Index nodeCount)
{
    for (Index i = 0; i < nodeCount; ++i) {
        const int ii = static_cast<int>(i);
        REQUIRE(std::isfinite(sol.psi(ii)));
        REQUIRE(std::isfinite(sol.phin(ii)));
        REQUIRE(std::isfinite(sol.phip(ii)));
        REQUIRE(std::isfinite(sol.n(ii)));
        REQUIRE(std::isfinite(sol.p(ii)));
        REQUIRE(sol.n(ii) >= 0.0);
        REQUIRE(sol.p(ii) >= 0.0);
    }
}

TEST_CASE("GummelSolver: high doping gradient reverse bias remains finite", "[gummel][stability]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    std::vector<RegionDopingSpec> specs = {
        {"n_region", 8.0e23, 0.0},
        {"p_region", 0.0, 5.0e21}
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const std::unordered_map<std::string, Real> biases = {
        {"anode", -0.25},
        {"cathode", 0.0}
    };

    GummelConfig cfg;
    cfg.maxIter = 100;
    cfg.reltol = 1.0e-5;
    cfg.abstol = 1.0e8;
    cfg.dampingPsi = 0.25;
    cfg.mobility = mobilityModelConfig("caughey_thomas");

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));
    REQUIRE(sol.iters >= 1);
    REQUIRE(sol.iters <= cfg.maxIter);
    requireFiniteDDSolution(sol, mesh.numNodes());
}

TEST_CASE("GummelSolver: multi-terminal contact biases are imposed", "[gummel][contacts]")
{
    DeviceMesh mesh = makeFourTerminalSiliconMesh();
    MaterialDatabase matdb;
    std::vector<RegionDopingSpec> specs = {
        {"p_channel", 0.0, 1.0e21},
        {"n_drain", 5.0e21, 0.0}
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const std::unordered_map<std::string, Real> biases = {
        {"body", 0.0},
        {"source", 0.0},
        {"gate", 0.05},
        {"drain", 0.10}
    };

    GummelConfig cfg;
    cfg.maxIter = 80;
    cfg.reltol = 1.0e-5;
    cfg.abstol = 1.0e8;
    cfg.dampingPsi = 0.5;

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));
    requireFiniteDDSolution(sol, mesh.numNodes());
    REQUIRE(sol.phin(3) == Catch::Approx(0.0));
    REQUIRE(sol.phip(3) == Catch::Approx(0.0));
    REQUIRE(sol.phin(4) == Catch::Approx(0.05));
    REQUIRE(sol.phip(4) == Catch::Approx(0.05));
    REQUIRE(sol.phin(5) == Catch::Approx(0.10));
    REQUIRE(sol.phip(5) == Catch::Approx(0.10));
}
