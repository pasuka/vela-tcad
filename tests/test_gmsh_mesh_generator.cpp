#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/mesh/GmshMeshGenerator.h"

using namespace Catch::Matchers;
using namespace vela;

TEST_CASE("gmsh generator currently reports not wired backend", "[mesh][gmsh][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {}});

    GmshMeshGenerator generator;
    REQUIRE_THROWS_WITH(
        generator.generate(ir),
        ContainsSubstring("backend not wired yet"));
}

