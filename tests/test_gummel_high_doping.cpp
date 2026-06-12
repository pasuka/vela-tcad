#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/GummelSolver.h"

#include <cmath>
#include <unordered_map>
#include <vector>

using namespace vela;

static DeviceMesh makeHighDopingPNMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;

    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;   n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = L;   n2.y = L;   mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0; n3.y = L;   mesh.addNode(n3);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0; c0.node_ids = {0, 1, 2}; mesh.addCell(c0);
    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 1; c1.node_ids = {0, 2, 3}; mesh.addCell(c1);

    Region r0; r0.id = 0; r0.name = "n_region"; r0.material = "Si"; r0.cell_ids = {0}; mesh.addRegion(r0);
    Region r1; r1.id = 1; r1.name = "p_region"; r1.material = "Si"; r1.cell_ids = {1}; mesh.addRegion(r1);

    Contact anode; anode.id = 0; anode.name = "anode"; anode.region_id = 1; anode.node_ids = {0, 3}; mesh.addContact(anode);
    Contact cathode; cathode.id = 1; cathode.name = "cathode"; cathode.region_id = 0; cathode.node_ids = {1, 2}; mesh.addContact(cathode);

    mesh.buildEdges();
    return mesh;
}

static DopingModel makeOneE24PNDoping(const DeviceMesh& mesh)
{
    std::vector<RegionDopingSpec> specs = {
        {"n_region", 1.0e24, 0.0},
        {"p_region", 0.0, 1.0e24},
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

static void requireFinitePositiveCarriers(const DDSolution& sol, Index nodeCount)
{
    for (Index i = 0; i < nodeCount; ++i) {
        const int ii = static_cast<int>(i);
        REQUIRE(std::isfinite(sol.psi(ii)));
        REQUIRE(std::isfinite(sol.phin(ii)));
        REQUIRE(std::isfinite(sol.phip(ii)));
        REQUIRE(std::isfinite(sol.n(ii)));
        REQUIRE(std::isfinite(sol.p(ii)));
        REQUIRE(sol.n(ii) > 0.0);
        REQUIRE(sol.p(ii) > 0.0);
    }
}

TEST_CASE("Gummel high-doping PN diode remains stable at 1e24 m^-3", "[gummel][high-doping]")
{
    DeviceMesh mesh = makeHighDopingPNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeOneE24PNDoping(mesh);

    const std::unordered_map<std::string, Real> biases = {
        {"anode", 0.0},
        {"cathode", 0.0},
    };

    GummelConfig cfg;
    cfg.maxIter = 80;
    cfg.reltol = 1.0e-6;
    cfg.abstol = 1.0e-8;
    cfg.dampingPsi = 0.25;

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));

    REQUIRE(sol.iters >= 1);
    REQUIRE(sol.iters <= cfg.maxIter);
    REQUIRE(sol.converged);
    requireFinitePositiveCarriers(sol, mesh.numNodes());

    const Material& si = matdb.getMaterial("Si");
    const Real Vt = constants::kb * cfg.temperature_K / constants::q;
    const Real ni = si.ni;
    const Real minority = ni * ni / 1.0e24;
    const Real builtIn = Vt * std::log(1.0e24 / minority);

    REQUIRE(sol.n(1) == Catch::Approx(1.0e24).epsilon(1.0e-12));
    REQUIRE(sol.p(1) == Catch::Approx(minority).epsilon(1.0e-12));
    REQUIRE(sol.p(3) == Catch::Approx(1.0e24).epsilon(1.0e-12));
    REQUIRE(sol.n(3) == Catch::Approx(minority).epsilon(1.0e-12));
    REQUIRE(sol.psi(1) - sol.psi(3) == Catch::Approx(builtIn).epsilon(1.0e-12));
}

TEST_CASE("Gummel high-doping forward bias does not produce invalid state", "[gummel][high-doping]")
{
    DeviceMesh mesh = makeHighDopingPNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeOneE24PNDoping(mesh);

    const std::unordered_map<std::string, Real> biases = {
        {"anode", 0.2},
        {"cathode", 0.0},
    };

    GummelConfig cfg;
    cfg.maxIter = 100;
    cfg.reltol = 1.0e-5;
    cfg.abstol = 1.0e-8;
    cfg.dampingPsi = 0.2;
    cfg.mobility = mobilityModelConfig("caughey_thomas");
    cfg.recombination = {"srh", "auger"};

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));

    REQUIRE(sol.iters >= 1);
    REQUIRE(sol.iters <= cfg.maxIter);
    REQUIRE(sol.converged);
    requireFinitePositiveCarriers(sol, mesh.numNodes());

    REQUIRE(sol.phin(0) == Catch::Approx(0.2));
    REQUIRE(sol.phip(0) == Catch::Approx(0.2));
    REQUIRE(sol.phin(1) == Catch::Approx(0.0));
    REQUIRE(sol.phip(1) == Catch::Approx(0.0));
}

TEST_CASE("Gummel high-doping contacts use Slotboom bandgap narrowing", "[gummel][high-doping][bgn]")
{
    DeviceMesh mesh = makeHighDopingPNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeOneE24PNDoping(mesh);

    const std::unordered_map<std::string, Real> biases = {
        {"anode", 0.0},
        {"cathode", 0.0},
    };

    GummelConfig cfg;
    cfg.maxIter = 80;
    cfg.reltol = 1.0e-6;
    cfg.abstol = 1.0e-8;
    cfg.dampingPsi = 0.25;
    cfg.bandgapNarrowing.model = "slotboom";

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));

    REQUIRE(sol.converged);
    requireFinitePositiveCarriers(sol, mesh.numNodes());

    const Material& si = matdb.getMaterial("Si");
    const Real Vt = constants::kb * cfg.temperature_K / constants::q;
    const SlotboomBandgapNarrowing bgn(cfg.bandgapNarrowing);
    const Real niEff = effectiveIntrinsicDensity(
        si.ni, Vt, bgn.deltaEg(1.0e24, 0.0, 0.0));
    const Real minority = niEff * niEff / 1.0e24;
    const Real builtIn = Vt * std::log(1.0e24 / minority);

    REQUIRE(niEff > si.ni);
    REQUIRE(sol.p(1) == Catch::Approx(minority).epsilon(1.0e-12));
    REQUIRE(sol.n(3) == Catch::Approx(minority).epsilon(1.0e-12));
    REQUIRE(sol.psi(1) - sol.psi(3) == Catch::Approx(builtIn).epsilon(1.0e-12));
}

TEST_CASE("Gummel Slotboom BGN uses total impurity at compensated contacts",
          "[gummel][high-doping][bgn][doping]")
{
    DeviceMesh mesh = makeHighDopingPNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e24, 1.0e24},
        {"p_region", 1.0e24, 1.0e24},
    });

    const std::unordered_map<std::string, Real> biases = {
        {"anode", 0.0},
        {"cathode", 0.0},
    };

    GummelConfig cfg;
    cfg.maxIter = 20;
    cfg.reltol = 1.0e-8;
    cfg.abstol = 1.0e-10;
    cfg.bandgapNarrowing.model = "slotboom";

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));

    const Material& si = matdb.getMaterial("Si");
    const Real Vt = constants::kb * cfg.temperature_K / constants::q;
    const SlotboomBandgapNarrowing bgn(cfg.bandgapNarrowing);
    const Real niEff = effectiveIntrinsicDensity(
        si.ni, Vt, bgn.deltaEg(2.0e24, 0.0, 0.0));

    REQUIRE(sol.converged);
    REQUIRE(doping.netDoping(1) == Catch::Approx(0.0));
    REQUIRE(doping.totalImpurity(1) == Catch::Approx(2.0e24));
    REQUIRE(niEff > si.ni);
    REQUIRE(sol.n(1) == Catch::Approx(niEff).epsilon(1.0e-12));
    REQUIRE(sol.p(1) == Catch::Approx(niEff).epsilon(1.0e-12));
    REQUIRE(sol.n(3) == Catch::Approx(niEff).epsilon(1.0e-12));
    REQUIRE(sol.p(3) == Catch::Approx(niEff).epsilon(1.0e-12));
}

TEST_CASE("Gummel high-doping asymmetric reverse bias does not diverge", "[gummel][high-doping][stability]")
{
    DeviceMesh mesh = makeHighDopingPNMesh();
    MaterialDatabase matdb;
    std::vector<RegionDopingSpec> specs = {
        {"n_region", 1.0e24, 0.0},
        {"p_region", 0.0, 2.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const std::unordered_map<std::string, Real> biases = {
        {"anode", -0.20},
        {"cathode", 0.0},
    };

    GummelConfig cfg;
    cfg.maxIter = 120;
    cfg.reltol = 1.0e-5;
    cfg.abstol = 1.0e8;
    cfg.dampingPsi = 0.2;
    cfg.mobility = mobilityModelConfig("caughey_thomas");
    cfg.bandgapNarrowing.model = "slotboom";

    DDSolution sol;
    REQUIRE_NOTHROW(sol = runGummel(mesh, matdb, doping, biases, cfg));
    REQUIRE(sol.iters >= 1);
    REQUIRE(sol.iters <= cfg.maxIter);
    requireFinitePositiveCarriers(sol, mesh.numNodes());
    REQUIRE(sol.phin(0) == Catch::Approx(-0.20));
    REQUIRE(sol.phip(0) == Catch::Approx(-0.20));
    REQUIRE(sol.phin(1) == Catch::Approx(0.0));
    REQUIRE(sol.phip(1) == Catch::Approx(0.0));
}
