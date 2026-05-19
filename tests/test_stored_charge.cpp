#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "vela/core/PhysicalConstants.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/post/StoredCharge.h"

namespace {

vela::DeviceMesh makeTwoRegionMesh()
{
    vela::DeviceMesh mesh;
    mesh.addNode(vela::Node{0, 0.0, 0.0, 0.0});
    mesh.addNode(vela::Node{1, 1.0, 0.0, 0.0});
    mesh.addNode(vela::Node{2, 1.0, 1.0, 0.0});
    mesh.addNode(vela::Node{3, 0.0, 1.0, 0.0});
    mesh.addCell(vela::Cell{0, vela::CellType::Tri3, 0, {0, 1, 2}});
    mesh.addCell(vela::Cell{1, vela::CellType::Tri3, 1, {0, 2, 3}});
    mesh.addRegion(vela::Region{0, "left", "Si", {0}});
    mesh.addRegion(vela::Region{1, "right", "Si", {1}});
    mesh.buildEdges();
    return mesh;
}

} // namespace

TEST_CASE("StoredCharge computes positive proxy for positive carriers", "[post][stored_charge]")
{
    const vela::DeviceMesh mesh = makeTwoRegionMesh();
    vela::DDSolution sol;
    sol.n = vela::VectorXd::Constant(4, 1.0e20);
    sol.p = vela::VectorXd::Constant(4, 2.0e20);

    vela::StoredCharge sc(mesh);
    vela::StoredChargeConfig cfg;
    cfg.regions = {"left", "right"};
    const auto out = sc.compute(sol, cfg);

    REQUIRE(out.perMeter);
    REQUIRE(std::isfinite(out.charge));
    REQUIRE(out.charge > 0.0);
}

TEST_CASE("StoredCharge returns zero when carriers are zero", "[post][stored_charge]")
{
    const vela::DeviceMesh mesh = makeTwoRegionMesh();
    vela::DDSolution sol;
    sol.n = vela::VectorXd::Zero(4);
    sol.p = vela::VectorXd::Zero(4);

    vela::StoredCharge sc(mesh);
    vela::StoredChargeConfig cfg;
    cfg.regions = {"left", "right"};
    const auto out = sc.compute(sol, cfg);

    REQUIRE(out.charge == 0.0);
}

TEST_CASE("StoredCharge region selection filters contribution", "[post][stored_charge]")
{
    const vela::DeviceMesh mesh = makeTwoRegionMesh();
    vela::DDSolution sol;
    sol.n = vela::VectorXd::Constant(4, 1.0e20);
    sol.p = vela::VectorXd::Constant(4, 1.0e20);

    vela::StoredCharge sc(mesh);
    vela::StoredChargeConfig allCfg;
    allCfg.regions = {"left", "right"};
    const auto all = sc.compute(sol, allCfg);

    vela::StoredChargeConfig oneCfg;
    oneCfg.regions = {"left"};
    const auto one = sc.compute(sol, oneCfg);

    REQUIRE(one.charge > 0.0);
    REQUIRE(all.charge > one.charge);
}
