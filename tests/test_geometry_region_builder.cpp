#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/preprocess/GeometryRegionBuilder.h"

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

