#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/io/NeutralMeshReader.h"
#include "vela/io/NeutralMeshWriter.h"
#include "vela/mesh/GmshMeshGenerator.h"
#include "vela/physics/DopingProfileEvaluator.h"
#include "vela/core/UnitScaling.h"

#include <filesystem>
#include <fstream>

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
    CHECK(loadedDoping.donors(0) == doping.donors(0));

    std::filesystem::remove_all(outDir, ec);
}

TEST_CASE("neutral mesh reader rejects invalid doping rows", "[io][neutral][preprocess]")
{
    const std::filesystem::path outDir = std::filesystem::temp_directory_path() / "vela_neutral_reader_invalid_test";
    std::error_code ec;
    std::filesystem::remove_all(outDir, ec);
    std::filesystem::create_directories(outDir);

    std::ofstream(outDir / "nodes.csv") << "id,x_um,y_um\n0,0,0\n1,1,0\n2,1,1\n3,0,1\n";
    std::ofstream(outDir / "elements.csv") << "id,node0,node1,node2,region_id\n0,0,1,2,0\n1,0,2,3,0\n";
    std::ofstream(outDir / "regions.csv") << "region_id,region_name,material\n0,R.Si,Silicon\n";
    std::ofstream(outDir / "contacts.csv") << "id,name,region_id,node_ids\n";
    std::ofstream(outDir / "doping.csv") << "node_id,donors_cm3,acceptors_cm3\n0,1e17,0\n0,1e17,0\n";

    NeutralMeshReader reader;
    const UnitScalingConfig scaling{};
    REQUIRE_THROWS_WITH(reader.readDopingCsv(outDir / "doping.csv", 4, scaling), Catch::Matchers::ContainsSubstring("duplicate"));

    std::filesystem::remove_all(outDir, ec);
}

TEST_CASE("neutral mesh reader rejects missing duplicate malformed and invalid doping ids", "[io][neutral][preprocess]")
{
    const std::filesystem::path outDir = std::filesystem::temp_directory_path() / "vela_neutral_reader_doping_contract_test";
    std::error_code ec;
    std::filesystem::remove_all(outDir, ec);
    std::filesystem::create_directories(outDir);

    auto writeBase = [&]() {
        std::ofstream(outDir / "nodes.csv") << "id,x_um,y_um\n0,0,0\n1,1,0\n2,1,1\n3,0,1\n";
        std::ofstream(outDir / "elements.csv") << "id,node0,node1,node2,region_id\n0,0,1,2,0\n1,0,2,3,0\n";
        std::ofstream(outDir / "regions.csv") << "region_id,region_name,material\n0,R.Si,Silicon\n";
        std::ofstream(outDir / "contacts.csv") << "id,name,region_id,node_ids\n";
    };

    NeutralMeshReader reader;
    const UnitScalingConfig scaling{};

    writeBase();
    std::ofstream(outDir / "doping.csv") << "node_id,donors_cm3,acceptors_cm3\n0,0,0\n1,0,0\n3,0,0\n";
    REQUIRE_THROWS_WITH(reader.readDopingCsv(outDir / "doping.csv", 4, scaling), Catch::Matchers::ContainsSubstring("missing node id"));

    writeBase();
    std::ofstream(outDir / "doping.csv") << "node_id,donors_cm3,acceptors_cm3\n0,0,0\n1,0,0\n1,0,0\n2,0,0\n3,0,0\n";
    REQUIRE_THROWS_WITH(reader.readDopingCsv(outDir / "doping.csv", 4, scaling), Catch::Matchers::ContainsSubstring("duplicate"));

    writeBase();
    std::ofstream(outDir / "doping.csv") << "node_id,donors_cm3,acceptors_cm3\n0,0,0\n1abc,0,0\n2,0,0\n3,0,0\n";
    REQUIRE_THROWS_WITH(reader.readDopingCsv(outDir / "doping.csv", 4, scaling), Catch::Matchers::ContainsSubstring("invalid"));

    writeBase();
    std::ofstream(outDir / "doping.csv") << "node_id,donors_cm3,acceptors_cm3\n0,0,0\n4,0,0\n2,0,0\n3,0,0\n";
    REQUIRE_THROWS_WITH(reader.readDopingCsv(outDir / "doping.csv", 4, scaling), Catch::Matchers::ContainsSubstring("invalid"));

    writeBase();
    std::ofstream(outDir / "doping.csv") << "node_id,donors_cm3,acceptors_cm3\n0,0,0\n1,nan,0\n2,0,0\n3,0,0\n";
    REQUIRE_THROWS_WITH(reader.readDopingCsv(outDir / "doping.csv", 4, scaling), Catch::Matchers::ContainsSubstring("invalid"));

    writeBase();
    std::ofstream(outDir / "doping.csv") << "node_id,donors_cm3,acceptors_cm3\n0,0,0\n1,-1,0\n2,0,0\n3,0,0\n";
    REQUIRE_THROWS_WITH(reader.readDopingCsv(outDir / "doping.csv", 4, scaling), Catch::Matchers::ContainsSubstring("invalid"));

    writeBase();
    std::ofstream(outDir / "doping.csv") << "node_id,donors_cm3,acceptors_cm3\n0,0,0\n1,0\n2,0,0\n3,0,0\n";
    REQUIRE_THROWS_WITH(reader.readDopingCsv(outDir / "doping.csv", 4, scaling), Catch::Matchers::ContainsSubstring("malformed"));

    std::filesystem::remove_all(outDir, ec);
}

TEST_CASE("neutral mesh reader rejects duplicate and missing ids", "[io][neutral][preprocess]")
{
    const std::filesystem::path outDir = std::filesystem::temp_directory_path() / "vela_neutral_reader_ids_test";
    std::error_code ec;
    std::filesystem::remove_all(outDir, ec);
    std::filesystem::create_directories(outDir);

    std::ofstream(outDir / "nodes.csv") << "id,x_um,y_um\n0,0,0\n1,1,0\n2,1,1\n3,0,1\n";
    std::ofstream(outDir / "elements.csv") << "id,node0,node1,node2,region_id\n0,0,1,2,0\n1,0,2,3,0\n";
    std::ofstream(outDir / "regions.csv") << "region_id,region_name,material\n0,R.Si,Silicon\n";
    std::ofstream(outDir / "contacts.csv") << "id,name,region_id,node_ids\n0,anode,0,0;1\n2,cathode,0,2;3\n";
    std::ofstream(outDir / "doping.csv") << "node_id,donors_cm3,acceptors_cm3\n0,1e17,0\n1,1e17,0\n2,1e17,0\n3,1e17,0\n";

    NeutralMeshReader reader;
    const UnitScalingConfig scaling{};
    REQUIRE_THROWS_WITH(reader.readDirectory(outDir, scaling), Catch::Matchers::ContainsSubstring("out-of-order"));

    std::filesystem::remove_all(outDir, ec);
}

TEST_CASE("neutral mesh reader rejects out-of-order nodes regions cells and contacts", "[io][neutral][preprocess]")
{
    const std::filesystem::path outDir = std::filesystem::temp_directory_path() / "vela_neutral_reader_order_test";
    std::error_code ec;
    std::filesystem::remove_all(outDir, ec);
    std::filesystem::create_directories(outDir);

    NeutralMeshReader reader;
    const UnitScalingConfig scaling{};

    std::ofstream(outDir / "nodes.csv") << "id,x_um,y_um\n1,1,0\n0,0,0\n2,1,1\n3,0,1\n";
    std::ofstream(outDir / "elements.csv") << "id,node0,node1,node2,region_id\n0,0,1,2,0\n1,0,2,3,0\n";
    std::ofstream(outDir / "regions.csv") << "region_id,region_name,material\n0,R.Si,Silicon\n";
    std::ofstream(outDir / "contacts.csv") << "id,name,region_id,node_ids\n0,anode,0,0;1\n";
    REQUIRE_THROWS_WITH(reader.readDirectory(outDir, scaling), Catch::Matchers::ContainsSubstring("out-of-order"));

    std::ofstream(outDir / "nodes.csv") << "id,x_um,y_um\n0,0,0\n1,1,0\n2,1,1\n3,0,1\n";
    std::ofstream(outDir / "regions.csv") << "region_id,region_name,material\n1,R.Si,Silicon\n0,R.Ox,Oxide\n";
    REQUIRE_THROWS_WITH(reader.readDirectory(outDir, scaling), Catch::Matchers::ContainsSubstring("out-of-order"));

    std::ofstream(outDir / "regions.csv") << "region_id,region_name,material\n0,R.Si,Silicon\n";
    std::ofstream(outDir / "elements.csv") << "id,node0,node1,node2,region_id\n1,0,1,2,0\n0,0,2,3,0\n";
    REQUIRE_THROWS_WITH(reader.readDirectory(outDir, scaling), Catch::Matchers::ContainsSubstring("out-of-order"));

    std::ofstream(outDir / "elements.csv") << "id,node0,node1,node2,region_id\n0,0,1,2,0\n1,0,2,3,0\n";
    std::ofstream(outDir / "contacts.csv") << "id,name,region_id,node_ids\n1,anode,0,0;1\n0,cathode,0,2;3\n";
    REQUIRE_THROWS_WITH(reader.readDirectory(outDir, scaling), Catch::Matchers::ContainsSubstring("out-of-order"));

    std::filesystem::remove_all(outDir, ec);
}
