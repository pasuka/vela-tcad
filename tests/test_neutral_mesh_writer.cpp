#include <catch2/catch_test_macros.hpp>

#include "vela/io/NeutralMeshWriter.h"
#include "vela/mesh/GmshMeshGenerator.h"
#include "vela/physics/DopingProfileEvaluator.h"

#include <filesystem>
#include <fstream>

using namespace vela;

TEST_CASE("neutral mesh writer exports expected files", "[io][neutral][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {0}});
    GeometryPrimitiveIr primitive;
    primitive.kind = GeometryPrimitiveKind::Rectangle;
    primitive.name = "rect_0";
    primitive.region = "R.Si";
    primitive.material = "Silicon";
    primitive.points = {Point2D{0.0, 0.0}, Point2D{1.0, 1.0}};
    ir.geometry.push_back(primitive);

    DopingProfileIr profile;
    profile.name = "ConstN";
    profile.kind = DopingProfileKind::Constant;
    profile.targetRegion = "R.Si";
    profile.donors_cm3 = 1.0e17;
    ir.dopingProfiles.push_back(profile);

    GmshMeshGenerator generator;
    const MeshBundle2D mesh = generator.generate(ir);
    const DopingModel doping = DopingProfileEvaluator::evaluate(mesh.mesh, ir);

    const std::filesystem::path outDir =
        std::filesystem::temp_directory_path() / "vela_neutral_writer_test";
    std::error_code ec;
    std::filesystem::remove_all(outDir, ec);

    NeutralMeshWriter::write(mesh, outDir, &doping);

    CHECK(std::filesystem::exists(outDir / "nodes.csv"));
    CHECK(std::filesystem::exists(outDir / "elements.csv"));
    CHECK(std::filesystem::exists(outDir / "regions.csv"));
    CHECK(std::filesystem::exists(outDir / "boundaries.csv"));
    CHECK(std::filesystem::exists(outDir / "contacts.csv"));
    CHECK(std::filesystem::exists(outDir / "doping.csv"));

    std::ifstream nodesFile(outDir / "nodes.csv");
    std::string header;
    std::getline(nodesFile, header);
    CHECK(header == "id,x_um,y_um");

    std::filesystem::remove_all(outDir, ec);
}

