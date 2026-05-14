#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>

#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/solver/GummelSolver.h"

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
    cfg.mobility = "caughey_thomas_field";
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

TEST_CASE("JSON solver config selects impact ionization model", "[impact][json]")
{
    const GummelConfig cfg = gummelConfigFromJson(nlohmann::json{
        {"impact_ionization", {{"model", "selberherr"}, {"electron_A_m_inv", 1.0e6}}}
    });
    REQUIRE(cfg.impactIonization.model == "selberherr");
    REQUIRE(cfg.impactIonization.electronA == Catch::Approx(1.0e6));
}
