#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/io/NeutralMeshWriter.h"
#include "vela/mesh/GmshMeshGenerator.h"
#include "vela/physics/DopingProfileEvaluator.h"
#include "vela/simulation/PoissonSimulation.h"

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

using namespace vela;

TEST_CASE("poisson simulation can run from neutral preprocess outputs", "[poisson][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Si", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Si";
    primitive.material = "Si";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}};
    ir.geometry.push_back(primitive);
    ir.contacts.push_back(ContactIr{"anode", "R.Si", {"rect_0#edge0"}});
    ir.contacts.push_back(ContactIr{"cathode", "R.Si", {"rect_0#edge2"}});
    ir.dopingProfiles.push_back(DopingProfileIr{
        "ConstN",
        DopingProfileKind::Constant,
        "R.Si",
        0,
        1.0e17,
        0.0});

    GmshMeshGenerator generator;
    const MeshBundle2D meshBundle = generator.generate(ir);
    const DopingModel doping = DopingProfileEvaluator::evaluate(meshBundle.mesh, ir);

    const std::filesystem::path root = std::filesystem::temp_directory_path() / "vela_preprocess_solver_test";
    const std::filesystem::path neutralDir = root / "neutral";
    const std::filesystem::path outputVtk = root / "result.vtk";
    const std::filesystem::path configPath = root / "sim.json";
    std::error_code ec;
    std::filesystem::remove_all(root, ec);
    std::filesystem::create_directories(root);

    NeutralMeshWriter::write(meshBundle, neutralDir, &doping);

    nlohmann::json contacts = nlohmann::json::array();
    contacts.push_back({{"name", "anode"}, {"bias", 0.0}});
    contacts.push_back({{"name", "cathode"}, {"bias", 0.0}});

    nlohmann::json cfg = {
        {"neutral_mesh_dir", neutralDir.string()},
        {"scaling", {{"mode", "unit_scaling"}}},
        {"output_vtk", outputVtk.string()},
        {"contacts", contacts},
        {"doping", nlohmann::json::array()},
    };
    std::ofstream out(configPath);
    out << cfg.dump(2);
    out.close();

    PoissonSimulation sim;
    const PoissonResult result = sim.runWithResult(configPath.string());
    CHECK(result.mesh.numNodes() == meshBundle.mesh.numNodes());
    CHECK(result.potential.size() == static_cast<int>(meshBundle.mesh.numNodes()));

    std::filesystem::remove_all(root, ec);
}

TEST_CASE("poisson simulation rejects neutral mesh without unit_scaling", "[poisson][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Si", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Si";
    primitive.material = "Si";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}};
    ir.geometry.push_back(primitive);

    GmshMeshGenerator generator;
    const MeshBundle2D meshBundle = generator.generate(ir);
    const DopingModel doping = DopingProfileEvaluator::evaluate(meshBundle.mesh, ir);

    const std::filesystem::path root = std::filesystem::temp_directory_path() / "vela_preprocess_solver_test_fail";
    const std::filesystem::path neutralDir = root / "neutral";
    const std::filesystem::path configPath = root / "sim.json";
    std::error_code ec;
    std::filesystem::remove_all(root, ec);
    std::filesystem::create_directories(root);

    NeutralMeshWriter::write(meshBundle, neutralDir, &doping);

    nlohmann::json cfg = {
        {"neutral_mesh_dir", neutralDir.string()},
        {"output_vtk", (root / "result.vtk").string()},
        {"contacts", nlohmann::json::array()},
        {"doping", nlohmann::json::array()},
    };
    std::ofstream out(configPath);
    out << cfg.dump(2);
    out.close();

    PoissonSimulation sim;
    REQUIRE_THROWS_WITH(sim.runWithResult(configPath.string()), Catch::Matchers::ContainsSubstring("unit_scaling"));

    std::filesystem::remove_all(root, ec);
}

TEST_CASE("poisson simulation rejects both mesh_file and neutral_mesh_dir", "[poisson][preprocess]")
{
    const std::filesystem::path root = std::filesystem::temp_directory_path() / "vela_preprocess_solver_dual_mesh_test";
    const std::filesystem::path configPath = root / "sim.json";
    std::error_code ec;
    std::filesystem::remove_all(root, ec);
    std::filesystem::create_directories(root);

    nlohmann::json cfg = {
        {"mesh_file", "mesh.json"},
        {"neutral_mesh_dir", "neutral"},
        {"output_vtk", (root / "result.vtk").string()},
        {"contacts", nlohmann::json::array()},
        {"doping", nlohmann::json::array()},
    };
    std::ofstream out(configPath);
    out << cfg.dump(2);
    out.close();

    PoissonSimulation sim;
    REQUIRE_THROWS_WITH(sim.runWithResult(configPath.string()), Catch::Matchers::ContainsSubstring("exactly one of 'mesh_file' or 'neutral_mesh_dir'"));

    std::filesystem::remove_all(root, ec);
}

TEST_CASE("poisson simulation rejects missing mesh source", "[poisson][preprocess]")
{
    const std::filesystem::path root = std::filesystem::temp_directory_path() / "vela_preprocess_solver_missing_mesh_test";
    const std::filesystem::path configPath = root / "sim.json";
    std::error_code ec;
    std::filesystem::remove_all(root, ec);
    std::filesystem::create_directories(root);

    nlohmann::json cfg = {
        {"output_vtk", (root / "result.vtk").string()},
        {"contacts", nlohmann::json::array()},
        {"doping", nlohmann::json::array()},
    };
    std::ofstream out(configPath);
    out << cfg.dump(2);
    out.close();

    PoissonSimulation sim;
    REQUIRE_THROWS_WITH(sim.runWithResult(configPath.string()), Catch::Matchers::ContainsSubstring("exactly one of 'mesh_file' or 'neutral_mesh_dir'"));

    std::filesystem::remove_all(root, ec);
}

TEST_CASE("poisson simulation prefers explicit node_doping_file over neutral doping csv", "[poisson][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Si", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Si";
    primitive.material = "Si";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}};
    ir.geometry.push_back(primitive);
    ir.contacts.push_back(ContactIr{"anode", "R.Si", {"rect_0#edge0"}});
    ir.contacts.push_back(ContactIr{"cathode", "R.Si", {"rect_0#edge2"}});

    GmshMeshGenerator generator;
    const MeshBundle2D meshBundle = generator.generate(ir);
    DopingModel neutralDoping(meshBundle.mesh.numNodes());
    for (Index i = 0; i < neutralDoping.numNodes(); ++i) {
        neutralDoping.setNodeDoping(i, 1.0e17, 0.0);
    }

    const std::filesystem::path root = std::filesystem::temp_directory_path() / "vela_preprocess_solver_doping_priority_test";
    const std::filesystem::path neutralDir = root / "neutral";
    const std::filesystem::path outputVtk = root / "result.vtk";
    const std::filesystem::path configPath = root / "sim.json";
    const std::filesystem::path explicitDopingPath = root / "explicit_doping.csv";
    std::error_code ec;
    std::filesystem::remove_all(root, ec);
    std::filesystem::create_directories(root);

    NeutralMeshWriter::write(meshBundle, neutralDir, &neutralDoping);
    std::ofstream explicitDoping(explicitDopingPath);
    explicitDoping << "node_id,donors_cm3,acceptors_cm3\n";
    for (Index i = 0; i < meshBundle.mesh.numNodes(); ++i) {
        explicitDoping << i << ",0,1e17\n";
    }
    explicitDoping.close();

    nlohmann::json contacts = nlohmann::json::array();
    contacts.push_back({{"name", "anode"}, {"bias", 0.0}});
    contacts.push_back({{"name", "cathode"}, {"bias", 0.0}});

    nlohmann::json cfg = {
        {"neutral_mesh_dir", neutralDir.string()},
        {"node_doping_file", explicitDopingPath.string()},
        {"scaling", {{"mode", "unit_scaling"}}},
        {"output_vtk", outputVtk.string()},
        {"contacts", contacts},
        {"doping", nlohmann::json::array()},
    };
    std::ofstream out(configPath);
    out << cfg.dump(2);
    out.close();

    PoissonSimulation sim;
    const PoissonResult result = sim.runWithResult(configPath.string());
    CHECK(result.netDoping[0] < 0.0);

    std::filesystem::remove_all(root, ec);
}
