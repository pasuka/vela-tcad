#include <catch2/catch_test_macros.hpp>

#include "vela/io/NeutralMeshReader.h"
#include "vela/io/NeutralMeshWriter.h"
#include "vela/mesh/GmshMeshGenerator.h"
#include "vela/physics/DopingProfileEvaluator.h"
#include "vela/core/UnitScaling.h"

#include <filesystem>

using namespace vela;

TEST_CASE("neutral mesh reader round-trip from writer output", "[io][neutral][preprocess]")
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
    ir.contacts.push_back(ContactIr{"anode", "R.Si", {"rect_0#edge0"}});

    GmshMeshGenerator generator;
    const MeshBundle2D generated = generator.generate(ir);
    const DopingModel doping = DopingProfileEvaluator::evaluate(generated.mesh, ir);

    const auto outDir = std::filesystem::temp_directory_path() / "vela_neutral_reader_test";
    std::error_code ec;
    std::filesystem::remove_all(outDir, ec);
    NeutralMeshWriter::write(generated, outDir, &doping);

    NeutralMeshReader reader;
    const UnitScalingConfig scaling{};
    const DeviceMesh mesh = reader.readDirectory(outDir, scaling);
    const DopingModel loadedDoping = reader.readDopingCsv(outDir / "doping.csv", mesh.numNodes(), scaling);

    CHECK(mesh.numNodes() == generated.mesh.numNodes());
    CHECK(mesh.numCells() == generated.mesh.numCells());
    CHECK(mesh.numRegions() == generated.mesh.numRegions());
    CHECK(mesh.numContacts() == generated.mesh.numContacts());
    CHECK(loadedDoping.numNodes() == mesh.numNodes());

    std::filesystem::remove_all(outDir, ec);
}
