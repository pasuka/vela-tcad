#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/preprocess/GeometryRegionBuilder.h"

#include <algorithm>
#include <cmath>
#include <string>

using namespace Catch::Matchers;
using namespace vela;

namespace {

DeviceIr2D buildTwoRegionIr()
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Left", "Silicon", {0}});
    ir.regions.push_back(RegionIr{"R.Right", "Silicon", {1}});

    GeometryPrimitiveIr left;
    left.kind = GeometryPrimitiveKind::Rectangle;
    left.name = "rect_0";
    left.region = "R.Left";
    left.material = "Silicon";
    left.points = {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}};
    ir.geometry.push_back(left);

    GeometryPrimitiveIr right;
    right.kind = GeometryPrimitiveKind::Rectangle;
    right.name = "rect_1";
    right.region = "R.Right";
    right.material = "Silicon";
    right.points = {Point2D{1.0, 0.0}, Point2D{2.0, 1.0}};
    ir.geometry.push_back(right);

    ContactIr contact;
    contact.name = "anode";
    contact.ownerRegion = "R.Left";
    contact.boundaryRefs = {"rect_0#edge0"};
    ir.contacts.push_back(contact);
    return ir;
}

} // namespace

TEST_CASE("geometry region builder tags regions and boundaries deterministically", "[preprocess][geometry]")
{
    const DeviceIr2D ir = buildTwoRegionIr();
    GeometryRegionBuilder builder;
    const GeometryRegionTopology topology = builder.build(ir);

    REQUIRE(topology.regions.size() == 2);
    CHECK(topology.regions[0].id == 0);
    CHECK(topology.regions[0].name == "R.Left");
    CHECK(topology.regions[1].id == 1);
    CHECK(topology.regions[1].name == "R.Right");

    REQUIRE(topology.boundaries.size() == 8);
    std::size_t internalCount = 0;
    std::size_t externalCount = 0;
    std::size_t contactCount = 0;
    for (const auto& edge : topology.boundaries) {
        if (edge.isExternal) {
            ++externalCount;
        } else {
            ++internalCount;
            CHECK(!edge.adjacentRegion.empty());
        }
        if (!edge.contactName.empty()) {
            ++contactCount;
            CHECK(edge.contactName == "anode");
        }
    }
    CHECK(internalCount == 2);
    CHECK(externalCount == 6);
    CHECK(contactCount == 1);
}

TEST_CASE("geometry region builder rejects unknown contact boundary refs", "[preprocess][geometry]")
{
    DeviceIr2D ir = buildTwoRegionIr();
    ir.contacts.front().boundaryRefs = {"missing#edge0"};

    GeometryRegionBuilder builder;
    REQUIRE_THROWS_WITH(
        builder.build(ir),
        ContainsSubstring("references unknown boundary"));
}

TEST_CASE("geometry region builder recognizes reversed shared edges", "[preprocess][geometry]")
{
    const GeometryRegionTopology topology = GeometryRegionBuilder{}.build(buildTwoRegionIr());
    std::size_t internalCount = 0;
    for (const auto& boundary : topology.boundaries) {
        if (!boundary.isExternal) {
            ++internalCount;
            CHECK(!boundary.adjacentRegion.empty());
        }
    }
    CHECK(internalCount == 2);
}

TEST_CASE("geometry region builder splits long edge against two short edges", "[preprocess][geometry]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Left", "Silicon", {0}});
    ir.regions.push_back(RegionIr{"R.Right", "Silicon", {1}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "left", "R.Left", "Silicon",
        {Point2D{0.0, 0.0}, Point2D{1.0, 2.0}}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Polygon, "right", "R.Right", "Silicon",
        {Point2D{1.0, 0.0}, Point2D{2.0, 0.0}, Point2D{2.0, 2.0},
         Point2D{1.0, 2.0}, Point2D{1.0, 1.0}}});

    const GeometryRegionTopology topology = GeometryRegionBuilder{}.build(ir);
    std::size_t leftInterfaceParts = 0;
    for (const auto& boundary : topology.boundaries) {
        if (boundary.regionId == 0 && !boundary.isExternal) {
            ++leftInterfaceParts;
            CHECK(boundary.adjacentRegion == "R.Right");
        }
    }
    CHECK(leftInterfaceParts == 2);
}

TEST_CASE("geometry region builder classifies only positive partial overlap as internal", "[preprocess][geometry]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Left", "Silicon", {0}});
    ir.regions.push_back(RegionIr{"R.Right", "Silicon", {1}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "left", "R.Left", "Silicon",
        {Point2D{0.0, 0.0}, Point2D{1.0, 2.0}}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "right", "R.Right", "Silicon",
        {Point2D{1.0, 0.5}, Point2D{2.0, 1.5}}});

    const GeometryRegionTopology topology = GeometryRegionBuilder{}.build(ir);
    std::size_t leftExternalParts = 0;
    std::size_t leftInternalParts = 0;
    for (const auto& boundary : topology.boundaries) {
        if (boundary.regionId != 0 ||
            std::abs(boundary.startUm.x_um - 1.0) > 1e-12 ||
            std::abs(boundary.endUm.x_um - 1.0) > 1e-12) {
            continue;
        }
        boundary.isExternal ? ++leftExternalParts : ++leftInternalParts;
    }
    CHECK(leftExternalParts == 2);
    CHECK(leftInternalParts == 1);
}

TEST_CASE("geometry region builder does not classify T point contact as interface", "[preprocess][geometry]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.A", "Silicon", {0}});
    ir.regions.push_back(RegionIr{"R.B", "Silicon", {1}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "a", "R.A", "Silicon",
        {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "b", "R.B", "Silicon",
        {Point2D{1.0, 1.0}, Point2D{2.0, 2.0}}});

    const GeometryRegionTopology topology = GeometryRegionBuilder{}.build(ir);
    CHECK(std::all_of(topology.boundaries.begin(), topology.boundaries.end(),
                      [](const auto& boundary) { return boundary.isExternal; }));
}

TEST_CASE("geometry region builder keeps short separated parallel edges external", "[preprocess][geometry]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.A", "Silicon", {0}});
    ir.regions.push_back(RegionIr{"R.B", "Silicon", {1}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "a", "R.A", "Silicon",
        {Point2D{0.0, 0.0}, Point2D{1.0e-6, 1.0e-6}}});
    ir.geometry.push_back(GeometryPrimitiveIr{
        GeometryPrimitiveKind::Rectangle, "b", "R.B", "Silicon",
        {Point2D{0.0, 1.0e-4}, Point2D{1.0e-6, 1.01e-4}}});

    const GeometryRegionTopology topology = GeometryRegionBuilder{}.build(ir);
    CHECK(std::all_of(topology.boundaries.begin(), topology.boundaries.end(),
                      [](const auto& boundary) { return boundary.isExternal; }));
}

TEST_CASE("geometry region builder rejects non-manifold shared subsegment", "[preprocess][geometry]")
{
    DeviceIr2D ir;
    for (int i = 0; i < 3; ++i) {
        ir.regions.push_back(RegionIr{"R." + std::to_string(i), "Silicon", {static_cast<std::size_t>(i)}});
        ir.geometry.push_back(GeometryPrimitiveIr{
            GeometryPrimitiveKind::Rectangle, "r" + std::to_string(i), "R." + std::to_string(i),
            "Silicon", {Point2D{static_cast<Real>(i == 0 ? 0 : 1), 0.0},
                         Point2D{static_cast<Real>(i == 0 ? 1 : 2), 1.0}}});
    }
    REQUIRE_THROWS_WITH(
        GeometryRegionBuilder{}.build(ir),
        ContainsSubstring("non-manifold"));
}

TEST_CASE("geometry region builder rejects contact covering internal interface", "[preprocess][geometry]")
{
    DeviceIr2D ir = buildTwoRegionIr();
    ir.contacts.clear();
    ir.contacts.push_back(ContactIr{"bad", "R.Left", {"rect_0#edge1"}});
    REQUIRE_THROWS_WITH(
        GeometryRegionBuilder{}.build(ir),
        ContainsSubstring("covers an internal interface"));
}
