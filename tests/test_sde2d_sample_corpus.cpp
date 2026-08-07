#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/io/SdeScriptReader.h"
#include "vela/mesh/GmshMeshGenerator.h"
#include "vela/physics/DopingProfileEvaluator.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
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
    const DopingModel doping = DopingProfileEvaluator::evaluate(mesh.mesh, ir);

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

    if (expected.contains("doping_checks")) {
        for (const auto& check : expected.at("doping_checks")) {
            const Real x = check.at("x_um").get<Real>();
            const Real y = check.at("y_um").get<Real>();
            const std::string kind = check.at("kind").get<std::string>();
            const std::string field = check.at("field").get<std::string>();
            const Real expectedValue = check.at("expected").get<Real>();
            const Real tolerance = check.at("relative_tolerance").get<Real>();
            Real bestDistance = std::numeric_limits<Real>::max();
            Index bestNodeId = 0;
            bool found = false;
            for (Index nodeId = 0; nodeId < mesh.mesh.numNodes(); ++nodeId) {
                const auto& node = mesh.mesh.getNode(nodeId);
                const Real distance = std::abs(node.x - x) + std::abs(node.y - y);
                if (distance < bestDistance) {
                    bestDistance = distance;
                    bestNodeId = nodeId;
                    found = true;
                }
            }
            CHECK(found);
            const Real actual = field == "donors_cm3" ? doping.donors(bestNodeId) : doping.acceptors(bestNodeId);
            if (kind == "peak") {
                CHECK(actual == Catch::Approx(expectedValue).epsilon(tolerance));
            } else {
                CHECK(actual == Catch::Approx(expectedValue).margin(std::max<Real>(1.0, expectedValue * tolerance)));
            }
        }
    }
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
