#include <catch2/catch_test_macros.hpp>

#include "vela/io/SdeScriptReader.h"
#include "vela/mesh/GmshMeshGenerator.h"

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>

using namespace vela;

namespace {

void checkSample(const std::filesystem::path& root, const std::string& sampleName)
{
    SdeScriptReader reader;
    GmshMeshGenerator generator;

    const auto sampleDir = root / "samples" / "sde2d" / sampleName;
    const auto scriptPath = sampleDir / "device.sde";
    const auto expectedPath = sampleDir / "expected.json";

    const DeviceIr2D ir = reader.parseFile(scriptPath.string());
    const MeshBundle2D mesh = generator.generate(ir);

    std::ifstream in(expectedPath);
    REQUIRE(in.is_open());
    nlohmann::json expected;
    in >> expected;

    if (expected.contains("nodes")) {
        CHECK(mesh.mesh.numNodes() == expected.at("nodes").get<Index>());
    }
    if (expected.contains("nodes_min")) {
        CHECK(mesh.mesh.numNodes() >= expected.at("nodes_min").get<Index>());
    }
    if (expected.contains("cells")) {
        CHECK(mesh.mesh.numCells() == expected.at("cells").get<Index>());
    }
    if (expected.contains("cells_min")) {
        CHECK(mesh.mesh.numCells() >= expected.at("cells_min").get<Index>());
    }
    CHECK(mesh.mesh.numRegions() == expected.at("regions").get<Index>());
    CHECK(mesh.mesh.numContacts() == expected.at("contacts").get<Index>());
}

} // namespace

TEST_CASE("sde2d sample corpus baseline", "[sde][samples][preprocess]")
{
    const std::filesystem::path repoRoot = VELA_SOURCE_DIR;
    checkSample(repoRoot, "minimal");
    checkSample(repoRoot, "doping");
    checkSample(repoRoot, "regions_contacts");
    checkSample(repoRoot, "refine");
}

