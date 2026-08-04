#include <catch2/catch_test_macros.hpp>

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
