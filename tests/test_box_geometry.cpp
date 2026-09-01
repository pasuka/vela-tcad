#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/mesh/BoxGeometryBuilder.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/simulation/ConfigParsing.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <streambuf>

using namespace vela;

namespace {

DeviceMesh makeSingleEquilateralTriangle()
{
    DeviceMesh mesh;
    const Real h = std::sqrt(3.0) / 2.0;

    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0; n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.5; n2.y = h; mesh.addNode(n2);

    Cell c; c.id = 0; c.type = CellType::Tri3; c.region_id = 0;
    c.node_ids = {0, 1, 2};
    mesh.addCell(c);

    Region r; r.id = 0; r.name = "body"; r.material = "Si"; r.cell_ids = {0};
    mesh.addRegion(r);

    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeUnitSquareTwoTriangles()
{
    DeviceMesh mesh;

    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0; n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 1.0; n2.y = 1.0; mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0; n3.y = 1.0; mesh.addNode(n3);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0;
    c0.node_ids = {0, 1, 2};
    mesh.addCell(c0);

    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 0;
    c1.node_ids = {0, 2, 3};
    mesh.addCell(c1);

    Region r; r.id = 0; r.name = "body"; r.material = "Si"; r.cell_ids = {0, 1};
    mesh.addRegion(r);

    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeObtuseTriangle(bool reverseWinding = false)
{
    DeviceMesh mesh;

    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 2.0; n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.2; n2.y = 0.1; mesh.addNode(n2);

    Cell c; c.id = 0; c.type = CellType::Tri3; c.region_id = 0;
    c.node_ids = reverseWinding ? std::vector<Index>{2, 1, 0}
                                : std::vector<Index>{0, 1, 2};
    mesh.addCell(c);

    Region r; r.id = 0; r.name = "body"; r.material = "Si"; r.cell_ids = {0};
    mesh.addRegion(r);

    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeRefinementTransitionPatch()
{
    DeviceMesh mesh;

    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 2.0; n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.0; n2.y = 1.0; mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 2.0; n3.y = 1.0; mesh.addNode(n3);
    Node n4; n4.id = 4; n4.x = 1.0; n4.y = 0.0; mesh.addNode(n4);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0;
    c0.node_ids = {0, 4, 2};
    mesh.addCell(c0);

    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 0;
    c1.node_ids = {4, 1, 3};
    mesh.addCell(c1);

    Cell c2; c2.id = 2; c2.type = CellType::Tri3; c2.region_id = 0;
    c2.node_ids = {4, 3, 2};
    mesh.addCell(c2);

    Region r; r.id = 0; r.name = "body"; r.material = "Si"; r.cell_ids = {0, 1, 2};
    mesh.addRegion(r);

    mesh.buildEdges();
    return mesh;
}

struct CerrRedirectGuard {
    explicit CerrRedirectGuard(std::ostream& stream, std::streambuf* replacement)
        : stream_(stream), original_(stream.rdbuf(replacement))
    {}

    ~CerrRedirectGuard()
    {
        stream_.rdbuf(original_);
    }

    std::ostream& stream_;
    std::streambuf* original_;
};

} // namespace

TEST_CASE("BoxGeometryBuilder: equilateral triangle area is correct", "[box_geometry]")
{
    DeviceMesh mesh = makeSingleEquilateralTriangle();
    const Real expectedArea = std::sqrt(3.0) / 4.0;

    REQUIRE(BoxGeometryBuilder::triangleArea(mesh.getNode(0), mesh.getNode(1), mesh.getNode(2)) ==
            Catch::Approx(expectedArea));

    Real volumeSum = 0.0;
    for (const Node& node : mesh.nodes())
        volumeSum += node.volume;
    REQUIRE(volumeSum == Catch::Approx(expectedArea));
}

TEST_CASE("BoxGeometryBuilder: square node volumes sum to total area", "[box_geometry]")
{
    DeviceMesh mesh = makeUnitSquareTwoTriangles();

    Real volumeSum = 0.0;
    for (const Node& node : mesh.nodes())
        volumeSum += node.volume;

    REQUIRE(volumeSum == Catch::Approx(1.0));
}

TEST_CASE("BoxGeometryBuilder: default node volumes remain barycentric", "[box_geometry]")
{
    DeviceMesh mesh = makeRefinementTransitionPatch();

    REQUIRE(mesh.getNode(0).volume == Catch::Approx(1.0 / 6.0));
    REQUIRE(mesh.getNode(1).volume == Catch::Approx(1.0 / 6.0));
    REQUIRE(mesh.getNode(2).volume == Catch::Approx(0.5));
    REQUIRE(mesh.getNode(3).volume == Catch::Approx(0.5));
    REQUIRE(mesh.getNode(4).volume == Catch::Approx(2.0 / 3.0));
}

TEST_CASE("BoxGeometryBuilder: mixed Voronoi node volumes are opt-in and conservative",
          "[box_geometry]")
{
    DeviceMesh mesh = makeRefinementTransitionPatch();

    BoxGeometryBuilder::Options options;
    options.nodeVolumePolicy = BoxGeometryBuilder::NodeVolumePolicy::MixedVoronoi;
    mesh.buildBoxGeometry(options);

    REQUIRE(mesh.getNode(0).volume == Catch::Approx(0.25));
    REQUIRE(mesh.getNode(1).volume == Catch::Approx(0.25));
    REQUIRE(mesh.getNode(2).volume == Catch::Approx(0.375));
    REQUIRE(mesh.getNode(3).volume == Catch::Approx(0.375));
    REQUIRE(mesh.getNode(4).volume == Catch::Approx(0.75));

    Real volumeSum = 0.0;
    for (const Node& node : mesh.nodes())
        volumeSum += node.volume;
    REQUIRE(volumeSum == Catch::Approx(2.0));
}

TEST_CASE("BoxGeometryBuilder: mixed Voronoi obtuse shares are positive and winding invariant",
          "[box_geometry][mixed_voronoi]")
{
    for (const bool reverseWinding : {false, true}) {
        DeviceMesh mesh = makeObtuseTriangle(reverseWinding);
        BoxGeometryBuilder::Options options;
        options.nodeVolumePolicy = BoxGeometryBuilder::NodeVolumePolicy::MixedVoronoi;
        mesh.buildBoxGeometry(options);

        REQUIRE(mesh.getNode(0).volume == Catch::Approx(0.025));
        REQUIRE(mesh.getNode(1).volume == Catch::Approx(0.025));
        REQUIRE(mesh.getNode(2).volume == Catch::Approx(0.050));

        Real volumeSum = 0.0;
        for (const Node& node : mesh.nodes()) {
            REQUIRE(node.volume > 0.0);
            volumeSum += node.volume;
        }
        REQUIRE(volumeSum == Catch::Approx(0.1));
    }
}

TEST_CASE("ConfigParsing: omitted node volume policy preserves explicit barycentric legacy behavior",
          "[box_geometry][config]")
{
    const auto omitted = parseBoxGeometryOptions(nlohmann::json::object());
    const auto emptyObject = parseBoxGeometryOptions(
        nlohmann::json{{"mesh_geometry", nlohmann::json::object()}});
    const auto explicitBarycentric = parseBoxGeometryOptions(
        nlohmann::json{{"mesh_geometry", {{"node_volume_policy", "barycentric"}}}});
    const auto explicitMixed = parseBoxGeometryOptions(
        nlohmann::json{{"mesh_geometry", {{"node_volume_policy", "mixed_voronoi"}}}});
    const auto qualifiedMixed = parseBoxGeometryOptions(
        nlohmann::json{{"mesh_geometry",
                        {{"node_volume_policy", "mixed_voronoi"},
                         {"require_non_obtuse", true}}}});
    const auto noNegativeFallback = parseBoxGeometryOptions(
        nlohmann::json{{"mesh_geometry",
                        {{"fallback_negative_cotangent", false}}}});

    REQUIRE(omitted.nodeVolumePolicy == BoxGeometryBuilder::NodeVolumePolicy::Barycentric);
    REQUIRE(emptyObject.nodeVolumePolicy == BoxGeometryBuilder::NodeVolumePolicy::Barycentric);
    REQUIRE(explicitBarycentric.nodeVolumePolicy ==
            BoxGeometryBuilder::NodeVolumePolicy::Barycentric);
    REQUIRE(explicitMixed.nodeVolumePolicy ==
            BoxGeometryBuilder::NodeVolumePolicy::MixedVoronoi);
    REQUIRE_FALSE(omitted.requireNonObtuse);
    REQUIRE_FALSE(explicitMixed.requireNonObtuse);
    REQUIRE(qualifiedMixed.requireNonObtuse);
    REQUIRE(omitted.fallbackNegativeCotangent);
    REQUIRE_FALSE(noNegativeFallback.fallbackNegativeCotangent);
    REQUIRE_THROWS_WITH(
        parseBoxGeometryOptions(
            nlohmann::json{{"mesh_geometry", {{"node_volume_policy", "unknown"}}}}),
        "ConfigParsing: mesh_geometry.node_volume_policy must be 'barycentric' or "
        "'mixed_voronoi'.");
    REQUIRE_THROWS_WITH(
        parseBoxGeometryOptions(
            nlohmann::json{{"mesh_geometry", {{"require_non_obtuse", "yes"}}}}),
        "ConfigParsing: mesh_geometry.require_non_obtuse must be boolean.");
    REQUIRE_THROWS_WITH(
        parseBoxGeometryOptions(nlohmann::json{
            {"mesh_geometry", {{"fallback_negative_cotangent", "yes"}}}}),
        "ConfigParsing: mesh_geometry.fallback_negative_cotangent must be boolean.");
}

TEST_CASE("External AverageBox carrier profile is explicit and leaves Poisson couples unchanged",
          "[box_geometry][carrier_transport_profile]")
{
    DeviceMesh mesh = makeSingleEquilateralTriangle();
    const Real original = mesh.getEdge(0).couple;
    const auto path = std::filesystem::temp_directory_path() /
        "vela_external_averagebox_transport_profile.csv";
    {
        std::ofstream output(path);
        REQUIRE(output.is_open());
        output << "node0,node1,couple_m\n";
        output << mesh.getEdge(0).n0 << ',' << mesh.getEdge(0).n1
               << ",2e-7\n";
    }
    const nlohmann::json cfg{{"mesh_geometry", {
        {"carrier_transport_couple_profile", "templates_ldmos_external_averagebox"},
        {"external_averagebox_couples_file", path.string()},
        {"external_averagebox_expected_edges", 1},
    }}};
    const CarrierTransportCoupleProfileReport report =
        applyCarrierTransportCoupleProfile(
            mesh, cfg, path.parent_path(), UnitScalingConfig{UnitScalingMode::UnitScaling});
    std::filesystem::remove(path);

    REQUIRE(report.profile == "templates_ldmos_external_averagebox");
    REQUIRE(report.records == 1);
    REQUIRE(mesh.getEdge(0).couple == Catch::Approx(original));
    REQUIRE(mesh.getEdge(0).transport_couple == Catch::Approx(0.2));
    REQUIRE(detail::computeEdgeCouplings(mesh).at(0) == Catch::Approx(original));
    REQUIRE(detail::computeTransportEdgeCouplings(mesh).at(0) == Catch::Approx(0.2));
    REQUIRE(mesh.getEdge(1).transport_couple < 0.0);
}

TEST_CASE("External AverageBox fields cannot be silently ignored by the default profile",
          "[box_geometry][carrier_transport_profile]")
{
    DeviceMesh mesh = makeSingleEquilateralTriangle();
    REQUIRE_THROWS_WITH(
        applyCarrierTransportCoupleProfile(
            mesh,
            nlohmann::json{{"mesh_geometry", {
                {"external_averagebox_couples_file", "unused.csv"},
            }}},
            std::filesystem::current_path()),
        Catch::Matchers::ContainsSubstring(
            "external AverageBox fields require"));
}

TEST_CASE("Templates LDMOS region AverageBox changes Poisson but preserves transport override",
          "[box_geometry][poisson_profile]")
{
    DeviceMesh mesh = makeSingleEquilateralTriangle();
    const auto transportPath = std::filesystem::temp_directory_path() /
        "vela_external_averagebox_transport_for_poisson.csv";
    const auto poissonPath = std::filesystem::temp_directory_path() /
        "vela_external_averagebox_poisson_profile.csv";
    {
        std::ofstream output(transportPath);
        REQUIRE(output.is_open());
        output << "node0,node1,couple_m\n";
        output << mesh.getEdge(0).n0 << ',' << mesh.getEdge(0).n1
               << ",2e-7\n";
    }
    {
        std::ofstream output(poissonPath);
        REQUIRE(output.is_open());
        output << "node0,node1,couple_m\n";
        output << mesh.getEdge(0).n0 << ',' << mesh.getEdge(0).n1
               << ",3e-7\n";
    }
    const nlohmann::json cfg{{"mesh_geometry", {
        {"node_volume_policy", "barycentric"},
        {"carrier_transport_couple_profile", "templates_ldmos_external_averagebox"},
        {"external_averagebox_couples_file", transportPath.string()},
        {"external_averagebox_expected_edges", 1},
        {"poisson_couple_profile", "templates_ldmos_region_averagebox"},
        {"external_averagebox_poisson_couples_file", poissonPath.string()},
        {"external_averagebox_poisson_expected_edges", 1},
    }}};
    const UnitScalingConfig scaling{UnitScalingMode::UnitScaling};
    applyCarrierTransportCoupleProfile(
        mesh, cfg, transportPath.parent_path(), scaling);
    const PoissonCoupleProfileReport report = applyPoissonCoupleProfile(
        mesh, cfg, poissonPath.parent_path(), scaling);
    std::filesystem::remove(transportPath);
    std::filesystem::remove(poissonPath);

    REQUIRE(report.profile == "templates_ldmos_region_averagebox");
    REQUIRE(report.records == 1);
    REQUIRE(mesh.getEdge(0).couple == Catch::Approx(0.3));
    REQUIRE(mesh.getEdge(0).transport_couple == Catch::Approx(0.2));
    REQUIRE(detail::computeEdgeCouplings(mesh).at(0) == Catch::Approx(0.3));
    REQUIRE(detail::computeTransportEdgeCouplings(mesh).at(0) == Catch::Approx(0.2));
}

TEST_CASE("Templates LDMOS region AverageBox requires qualified transport profile",
          "[box_geometry][poisson_profile]")
{
    DeviceMesh mesh = makeSingleEquilateralTriangle();
    REQUIRE_THROWS_WITH(
        applyPoissonCoupleProfile(
            mesh,
            nlohmann::json{{"mesh_geometry", {
                {"poisson_couple_profile", "templates_ldmos_region_averagebox"},
                {"external_averagebox_poisson_couples_file", "unused.csv"},
                {"external_averagebox_poisson_expected_edges", 1},
            }}},
            std::filesystem::current_path()),
        Catch::Matchers::ContainsSubstring(
            "requires the qualified templates_ldmos_external_averagebox"));
}

TEST_CASE("Templates LDMOS AverageBox profile rejects non-barycentric source geometry",
          "[box_geometry][carrier_transport_profile]")
{
    DeviceMesh mesh = makeSingleEquilateralTriangle();
    REQUIRE_THROWS_WITH(
        applyCarrierTransportCoupleProfile(
            mesh,
            nlohmann::json{{"mesh_geometry", {
                {"node_volume_policy", "mixed_voronoi"},
                {"carrier_transport_couple_profile",
                 "templates_ldmos_external_averagebox"},
                {"external_averagebox_couples_file", "unused.csv"},
                {"external_averagebox_expected_edges", 1},
            }}},
            std::filesystem::current_path()),
        Catch::Matchers::ContainsSubstring(
            "requires mesh_geometry.node_volume_policy='barycentric'"));
}

TEST_CASE("BoxGeometryBuilder: non-obtuse qualification accepts eligible meshes and rejects obtuse cells",
          "[box_geometry][mixed_voronoi][qualification]")
{
    BoxGeometryBuilder::Options options;
    options.nodeVolumePolicy = BoxGeometryBuilder::NodeVolumePolicy::MixedVoronoi;
    options.requireNonObtuse = true;

    DeviceMesh acute = makeSingleEquilateralTriangle();
    REQUIRE_NOTHROW(acute.buildBoxGeometry(options));

    DeviceMesh right = makeUnitSquareTwoTriangles();
    REQUIRE_NOTHROW(right.buildBoxGeometry(options));

    DeviceMesh obtuse = makeObtuseTriangle();
    REQUIRE_THROWS_WITH(
        obtuse.buildBoxGeometry(options),
        Catch::Matchers::ContainsSubstring("mesh_geometry.require_non_obtuse rejected cell 0"));
}

TEST_CASE("BoxGeometryBuilder: all edge couplings are non-negative", "[box_geometry]")
{
    DeviceMesh square = makeUnitSquareTwoTriangles();
    for (const Edge& edge : square.edges())
        REQUIRE(edge.couple >= 0.0);

    DeviceMesh obtuse = makeObtuseTriangle();
    for (const Edge& edge : obtuse.edges())
        REQUIRE(edge.couple >= 0.0);
}

TEST_CASE("BoxGeometryBuilder: warnings are opt-in", "[box_geometry]")
{
    REQUIRE_FALSE(BoxGeometryBuilder::Options{}.warnOnNegativeCotangent);

    std::ostringstream capturedErr;
    CerrRedirectGuard guard(std::cerr, capturedErr.rdbuf());

    DeviceMesh mesh;
    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 2.0; n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.2; n2.y = 0.1; mesh.addNode(n2);
    Cell c; c.id = 0; c.type = CellType::Tri3; c.region_id = 0;
    c.node_ids = {0, 1, 2};
    mesh.addCell(c);
    Region r; r.id = 0; r.name = "body"; r.material = "Si"; r.cell_ids = {0};
    mesh.addRegion(r);
    mesh.buildEdges();

    REQUIRE(mesh.numEdges() == 3);
    REQUIRE(mesh.getEdge(0).couple >= 0.0);
    REQUIRE(capturedErr.str().empty());
}

TEST_CASE("BoxGeometryBuilder: obtuse triangle reports fallback and keeps couplings non-negative", "[box_geometry]")
{
    DeviceMesh mesh = makeObtuseTriangle();
    const GeometryBuildReport& report = mesh.lastGeometryBuildReport();

    REQUIRE(report.totalCells == 1);
    REQUIRE(report.degenerateCells == 0);
    REQUIRE(report.negativeCotangentCount == 1);
    REQUIRE(report.fallbackCount == 1);
    REQUIRE(report.minAngleDegrees > 0.0);
    REQUIRE(report.maxAngleDegrees > 90.0);
    REQUIRE(report.minEdgeLength > 0.0);

    for (const Edge& edge : mesh.edges())
        REQUIRE(edge.couple >= 0.0);
}
