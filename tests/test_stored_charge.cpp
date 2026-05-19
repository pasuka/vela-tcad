#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include "vela/mesh/DeviceMesh.h"
#include "vela/post/StoredCharge.h"

TEST_CASE("StoredCharge sums finite non-negative mobile charge", "[post]") {
    vela::DeviceMesh mesh;
    mesh.addNode(vela::Node{0, 0.0, 0.0, 0.0});
    mesh.addNode(vela::Node{1, 1.0, 0.0, 0.0});
    mesh.addNode(vela::Node{2, 0.0, 1.0, 0.0});
    mesh.addCell(vela::Cell{0, vela::CellType::Tri3, 0, {0, 1, 2}});
    mesh.addRegion(vela::Region{0, "bulk", "Si", {0}});
    mesh.buildEdges();

    vela::DDSolution sol;
    sol.n = vela::VectorXd::Constant(3, 1.0e20);
    sol.p = vela::VectorXd::Constant(3, 2.0e20);

    vela::StoredCharge sc(mesh);
    vela::StoredChargeConfig cfg;
    cfg.regions = {"bulk"};
    cfg.perMeter = true;

    const auto out = sc.compute(sol, cfg);
    REQUIRE(out.perMeter);
    REQUIRE(out.charge >= 0.0);
    REQUIRE(std::isfinite(out.charge));
}
