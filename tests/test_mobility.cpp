#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScaling.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/solver/GummelSolver.h"

#include <cmath>
#include <unordered_map>
#include <vector>

using namespace vela;

static DeviceMesh makePNMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;

    Node n0; n0.id=0; n0.x=0;  n0.y=0;  mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=L;  n1.y=0;  mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=L;  n2.y=L;  mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=0;  n3.y=L;  mesh.addNode(n3);

    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0; c0.node_ids = {0, 1, 2}; mesh.addCell(c0);
    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=1; c1.node_ids = {0, 2, 3}; mesh.addCell(c1);

    Region r0; r0.id=0; r0.name="n_region"; r0.material="Si"; r0.cell_ids={0}; mesh.addRegion(r0);
    Region r1; r1.id=1; r1.name="p_region"; r1.material="Si"; r1.cell_ids={1}; mesh.addRegion(r1);

    Contact anode; anode.id=0; anode.name="anode"; anode.region_id=1; anode.node_ids={0,3}; mesh.addContact(anode);
    Contact cathode; cathode.id=1; cathode.name="cathode"; cathode.region_id=0; cathode.node_ids={1,2}; mesh.addContact(cathode);

    mesh.buildEdges();
    return mesh;
}

static DopingModel makePNDoping(const DeviceMesh& mesh)
{
    std::vector<RegionDopingSpec> specs = {
        {"n_region", 1.0e23, 0.0},
        {"p_region", 0.0, 1.0e23},
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

TEST_CASE("Caughey-Thomas mobility decreases as doping increases", "[mobility]")
{
    MaterialDatabase matdb;
    const Material& si = matdb.getMaterial("Si");
    DopingDependentMobility mobility;

    const Real lowDopingElectron = mobility.electronMobility(si, 1.0e20, 0.0, 0.0);
    const Real highDopingElectron = mobility.electronMobility(si, 1.0e25, 0.0, 0.0);
    const Real lowDopingHole = mobility.holeMobility(si, 1.0e20, 0.0, 0.0);
    const Real highDopingHole = mobility.holeMobility(si, 1.0e25, 0.0, 0.0);

    REQUIRE(highDopingElectron < lowDopingElectron);
    REQUIRE(highDopingHole < lowDopingHole);
    REQUIRE(lowDopingElectron <= si.mun);
    REQUIRE(lowDopingHole <= si.mup);
}


TEST_CASE("JSON solver config selects mobility and recombination models", "[mobility][json]")
{
    const nlohmann::json json = {
        {"mobility", "caughey_thomas"},
        {"recombination", {"srh", "auger"}},
        {"taun", 2.0e-7},
        {"taup", 3.0e-7},
        {"bandgap_narrowing", {{"model", "slotboom"}, {"coefficient_eV", 0.010}}},
    };

    const GummelConfig cfg = gummelConfigFromJson(json);
    REQUIRE(cfg.mobility.model == "caughey_thomas");
    REQUIRE(cfg.recombination.size() == 2);
    REQUIRE(cfg.recombination[0] == "srh");
    REQUIRE(cfg.recombination[1] == "auger");
    REQUIRE(cfg.taun == Catch::Approx(2.0e-7));
    REQUIRE(cfg.taup == Catch::Approx(3.0e-7));
    REQUIRE(cfg.bandgapNarrowing.model == "slotboom");
    REQUIRE(cfg.bandgapNarrowing.coefficient == Catch::Approx(0.010));
}

TEST_CASE("JSON solver config unit_scaling normalizes mobility and field inputs",
          "[mobility][json][scaling]")
{
    const nlohmann::json json = {
        {"mobility", {
            {"model", "caughey_thomas_field_surface"},
            {"electron_mu_min_m2_V_s", 52.2},
            {"electron_nref_m3", 9.68e16},
            {"hole_mu_min_m2_V_s", 44.9},
            {"hole_nref_m3", 2.23e17},
            {"surface", {
                {"reference_field_V_per_m", 1.5e4},
                {"theta_electron_m_per_V", 1.0e-6},
                {"theta_hole_m_per_V", 2.0e-6}
            }}
        }},
        {"bandgap_narrowing", {
            {"model", "slotboom"},
            {"reference_doping_m3", 1.0e17}
        }},
        {"impact_ionization", {
            {"model", "selberherr"},
            {"electron_A_m_inv", 7.03e5},
            {"electron_B_V_m", 1.231e6},
            {"hole_A_m_inv", 1.582e6},
            {"hole_B_V_m", 2.036e6}
        }}
    };

    const GummelConfig cfg = gummelConfigFromJson(
        json, UnitScalingConfig{UnitScalingMode::UnitScaling});

    REQUIRE(cfg.mobility.electronCT.muMin == Catch::Approx(0.00522));
    REQUIRE(cfg.mobility.electronCT.nRef == Catch::Approx(9.68e22));
    REQUIRE(cfg.mobility.holeCT.muMin == Catch::Approx(0.00449));
    REQUIRE(cfg.mobility.holeCT.nRef == Catch::Approx(2.23e23));
    REQUIRE(cfg.mobility.surface.referenceField == Catch::Approx(1.5e6));
    REQUIRE(cfg.mobility.surface.thetaElectron == Catch::Approx(1.0e-8));
    REQUIRE(cfg.mobility.surface.thetaHole == Catch::Approx(2.0e-8));
    REQUIRE(cfg.bandgapNarrowing.referenceDoping == Catch::Approx(1.0e23));
    REQUIRE(cfg.impactIonization.electronA == Catch::Approx(7.03e7));
    REQUIRE(cfg.impactIonization.electronB == Catch::Approx(1.231e8));
    REQUIRE(cfg.impactIonization.holeA == Catch::Approx(1.582e8));
    REQUIRE(cfg.impactIonization.holeB == Catch::Approx(2.036e8));
}


TEST_CASE("JSON mobility surface interface rejects canonical and alias together",
          "[mobility][json][surface]")
{
    const nlohmann::json json = {
        {"model", "caughey_thomas_surface"},
        {"surface", {
            {"surface_interface", {"channel", "gate_oxide"}},
            {"interface", {"channel", "field_oxide"}}
        }}
    };

    REQUIRE_THROWS_AS(mobilityModelConfigFromJson(json), std::invalid_argument);
}

TEST_CASE("Gummel PN diode runs with configured mobility and recombination", "[mobility][gummel]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);

    GummelConfig cfg;
    cfg.maxIter = 20;
    cfg.reltol = 1.0e-5;
    cfg.dampingPsi = 0.5;
    cfg.mobility = mobilityModelConfig("caughey_thomas");
    cfg.recombination = {"srh", "auger"};

    std::unordered_map<std::string, Real> biases = {
        {"anode", 0.0},
        {"cathode", 0.0},
    };

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));
    REQUIRE(sol.iters >= 1);
    for (int i = 0; i < static_cast<int>(mesh.numNodes()); ++i) {
        REQUIRE(std::isfinite(sol.n(i)));
        REQUIRE(std::isfinite(sol.p(i)));
        REQUIRE(sol.n(i) >= 0.0);
        REQUIRE(sol.p(i) >= 0.0);
    }
}

TEST_CASE("Caughey-Thomas field mobility rolls off toward velocity saturation", "[mobility][field]")
{
    MaterialDatabase matdb;
    const Material& si = matdb.getMaterial("Si");
    MobilityModelConfig config = mobilityModelConfig("caughey_thomas_field");
    DopingDependentMobility mobility(config);

    const Real lowField = mobility.electronMobility(si, 1.0e20, 0.0, 0.0, 0.0);
    const Real highField = mobility.electronMobility(si, 1.0e20, 0.0, 0.0, 1.0e8);

    REQUIRE(highField < lowField);
    REQUIRE(highField * 1.0e8 <= config.electronField.saturationVelocity * 1.01);
}

TEST_CASE("Material temperature path updates intrinsic density and mobility", "[mobility][temperature]")
{
    MaterialDatabase matdb;
    const Material si300 = matdb.getMaterial("Si", 300.0);
    const Material si400 = matdb.getMaterial("Si", 400.0);

    REQUIRE(si300.temperature_K.has_value());
    REQUIRE(*si300.temperature_K == Catch::Approx(300.0));
    REQUIRE(si400.ni > si300.ni);
    REQUIRE(si400.mun < si300.mun);
}

TEST_CASE("surface mobility degradation decreases with vertical field", "[mobility][surface]")
{
    MaterialDatabase matdb;
    const Material& si = matdb.getMaterial("Si");
    MobilityModelConfig config = mobilityModelConfig("caughey_thomas_surface");
    config.surface.thetaElectron = 2.0e-8;
    config.surface.thetaHole = 1.0e-8;
    config.surface.beta = 1.0;
    DopingDependentMobility mobility(config);

    const Real lowField = mobility.electronMobility(si, 1.0e20, 0.0, 0.0, 0.0, 0.0);
    const Real highField = mobility.electronMobility(si, 1.0e20, 0.0, 0.0, 0.0, 1.0e8);
    const Real highFieldHole = mobility.holeMobility(si, 1.0e20, 0.0, 0.0, 0.0, 1.0e8);

    REQUIRE(highField < lowField);
    REQUIRE(highFieldHole < mobility.holeMobility(si, 1.0e20, 0.0, 0.0, 0.0, 0.0));
}

TEST_CASE("surface mobility theta zero preserves baseline", "[mobility][surface]")
{
    MaterialDatabase matdb;
    const Material& si = matdb.getMaterial("Si");
    MobilityModelConfig baselineConfig = mobilityModelConfig("caughey_thomas");
    MobilityModelConfig surfaceConfig = mobilityModelConfig("caughey_thomas_surface");
    surfaceConfig.surface.thetaElectron = 0.0;
    surfaceConfig.surface.thetaHole = 0.0;

    DopingDependentMobility baseline(baselineConfig);
    DopingDependentMobility surface(surfaceConfig);

    REQUIRE(surface.electronMobility(si, 1.0e22, 0.0, 0.0, 0.0, 1.0e9) ==
            Catch::Approx(baseline.electronMobility(si, 1.0e22, 0.0, 0.0, 0.0)));
    REQUIRE(surface.holeMobility(si, 1.0e22, 0.0, 0.0, 0.0, 1.0e9) ==
            Catch::Approx(baseline.holeMobility(si, 1.0e22, 0.0, 0.0, 0.0)));
}

TEST_CASE("surface mobility rejects invalid parameters", "[mobility][surface]")
{
    MaterialDatabase matdb;
    const Material& si = matdb.getMaterial("Si");
    MobilityModelConfig config = mobilityModelConfig("caughey_thomas_surface");
    config.surface.thetaElectron = -1.0e-8;
    DopingDependentMobility mobility(config);

    REQUIRE_THROWS_AS(
        mobility.electronMobility(si, 1.0e20, 0.0, 0.0, 0.0, 1.0e8),
        std::invalid_argument);
}

TEST_CASE("JSON mobility object parses surface settings", "[mobility][json][surface]")
{
    const nlohmann::json json = {
        {"mobility", {
            {"model", "caughey_thomas_field_surface"},
            {"surface", {
                {"theta_electron_m_per_V", 2.0e-8},
                {"theta_hole_m_per_V", 3.0e-8},
                {"beta", 2.0},
                {"reference_field_V_per_m", 1.0e6},
                {"min_factor", 0.1},
                {"surface_region", "p_body"},
                {"surface_interface", {"p_body", "gate_oxide"}}
            }}
        }}
    };

    const GummelConfig cfg = gummelConfigFromJson(json);
    REQUIRE(cfg.mobility.model == "caughey_thomas_field_surface");
    REQUIRE(cfg.mobility.surface.thetaElectron == Catch::Approx(2.0e-8));
    REQUIRE(cfg.mobility.surface.thetaHole == Catch::Approx(3.0e-8));
    REQUIRE(cfg.mobility.surface.beta == Catch::Approx(2.0));
    REQUIRE(cfg.mobility.surface.referenceField == Catch::Approx(1.0e6));
    REQUIRE(cfg.mobility.surface.minFactor == Catch::Approx(0.1));
    REQUIRE(cfg.mobility.surface.surfaceRegion == "p_body");
    REQUIRE(cfg.mobility.surface.surfaceInterface.size() == 2);
}
