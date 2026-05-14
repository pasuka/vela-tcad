#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>

#include "vela/core/PhysicalConstants.h"
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
    REQUIRE(cfg.mobility == "caughey_thomas");
    REQUIRE(cfg.recombination.size() == 2);
    REQUIRE(cfg.recombination[0] == "srh");
    REQUIRE(cfg.recombination[1] == "auger");
    REQUIRE(cfg.taun == Catch::Approx(2.0e-7));
    REQUIRE(cfg.taup == Catch::Approx(3.0e-7));
    REQUIRE(cfg.bandgapNarrowing.model == "slotboom");
    REQUIRE(cfg.bandgapNarrowing.coefficient == Catch::Approx(0.010));
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
    cfg.mobility = "caughey_thomas";
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
