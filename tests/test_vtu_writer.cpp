#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/io/VtuWriter.h"
#include "vela/mesh/DeviceMesh.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using namespace vela;

namespace {

DeviceMesh makeTwoTriangleMesh()
{
    DeviceMesh mesh;
    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0; n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 1.0; n2.y = 1.0; mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0; n3.y = 1.0; mesh.addNode(n3);
    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0; c0.node_ids = {0, 1, 2}; mesh.addCell(c0);
    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 1; c1.node_ids = {0, 2, 3}; mesh.addCell(c1);
    Region r0; r0.id = 0; r0.name = "n_region"; r0.material = "Si"; r0.cell_ids = {0}; mesh.addRegion(r0);
    Region r1; r1.id = 1; r1.name = "p_region"; r1.material = "Si"; r1.cell_ids = {1}; mesh.addRegion(r1);
    mesh.buildEdges();
    return mesh;
}

/// Minimal independent reader for the subset of VTK XML the writer emits.
struct ParsedVtu {
    int numberOfPoints = 0;
    int numberOfCells = 0;
    std::vector<double> points;
    std::vector<int> connectivity;
    std::vector<int> offsets;
    std::vector<int> types;
    struct Array {
        std::string type;
        std::string format;
        int components = 1;
        std::string body;
    };
    std::map<std::string, Array> pointData;
    std::map<std::string, Array> cellData;
    std::map<std::string, Array> fieldData;
    std::vector<std::string> pointDataOrder;
};

std::string sectionBetween(const std::string& text, const std::string& open, const std::string& close)
{
    const std::size_t start = text.find(open);
    REQUIRE(start != std::string::npos);
    const std::size_t end = text.find(close, start);
    REQUIRE(end != std::string::npos);
    return text.substr(start, end - start);
}

std::map<std::string, ParsedVtu::Array> parseArrays(const std::string& section, std::vector<std::string>* order = nullptr)
{
    static const std::regex arrayRe(
        R"(<DataArray\s+([^>]*)>([\s\S]*?)</DataArray>)");
    static const std::regex attrRe(R"((\w+)=\"([^\"]*)\")");
    std::map<std::string, ParsedVtu::Array> out;
    for (auto it = std::sregex_iterator(section.begin(), section.end(), arrayRe); it != std::sregex_iterator(); ++it) {
        ParsedVtu::Array array;
        std::string name;
        const std::string attrs = (*it)[1].str();
        for (auto a = std::sregex_iterator(attrs.begin(), attrs.end(), attrRe); a != std::sregex_iterator(); ++a) {
            const std::string key = (*a)[1].str();
            const std::string value = (*a)[2].str();
            if (key == "Name") name = value;
            else if (key == "type") array.type = value;
            else if (key == "format") array.format = value;
            else if (key == "NumberOfComponents") array.components = std::stoi(value);
        }
        array.body = (*it)[2].str();
        REQUIRE_FALSE(name.empty());
        out[name] = array;
        if (order != nullptr)
            order->push_back(name);
    }
    return out;
}

template <typename T>
std::vector<T> asciiValues(const std::string& body)
{
    std::vector<T> values;
    std::istringstream stream(body);
    std::string token;
    while (stream >> token) {
        if constexpr (std::is_integral_v<T>)
            values.push_back(static_cast<T>(std::stoll(token)));
        else
            values.push_back(static_cast<T>(std::stod(token)));
    }
    return values;
}

std::vector<unsigned char> base64Decode(const std::string& text)
{
    auto value = [](char c) -> int {
        if (c >= 'A' && c <= 'Z') return c - 'A';
        if (c >= 'a' && c <= 'z') return c - 'a' + 26;
        if (c >= '0' && c <= '9') return c - '0' + 52;
        if (c == '+') return 62;
        if (c == '/') return 63;
        return -1;
    };
    std::vector<unsigned char> out;
    int buffer = 0;
    int bits = 0;
    for (const char c : text) {
        const int v = value(c);
        if (v < 0)
            continue;
        buffer = (buffer << 6) | v;
        bits += 6;
        if (bits >= 8) {
            bits -= 8;
            out.push_back(static_cast<unsigned char>((buffer >> bits) & 0xFF));
        }
    }
    return out;
}

std::vector<double> binaryDoubles(const std::string& body)
{
    const std::vector<unsigned char> bytes = base64Decode(body);
    REQUIRE(bytes.size() >= 4);
    std::uint32_t count = 0;
    std::memcpy(&count, bytes.data(), 4);
    REQUIRE(bytes.size() == 4 + count);
    std::vector<double> values(count / sizeof(double));
    std::memcpy(values.data(), bytes.data() + 4, count);
    return values;
}

std::vector<int> binaryInts(const std::string& body)
{
    const std::vector<unsigned char> bytes = base64Decode(body);
    REQUIRE(bytes.size() >= 4);
    std::uint32_t count = 0;
    std::memcpy(&count, bytes.data(), 4);
    REQUIRE(bytes.size() == 4 + count);
    std::vector<std::int32_t> raw(count / sizeof(std::int32_t));
    std::memcpy(raw.data(), bytes.data() + 4, count);
    return std::vector<int>(raw.begin(), raw.end());
}

std::vector<int> intArray(const ParsedVtu::Array& array)
{
    return array.format == "ascii" ? asciiValues<int>(array.body) : binaryInts(array.body);
}

std::string stringFieldValue(const ParsedVtu::Array& array)
{
    std::string text;
    for (const int code : asciiValues<int>(array.body)) {
        if (code == 0)
            break;
        text += static_cast<char>(code);
    }
    return text;
}

ParsedVtu parse(const std::string& text)
{
    ParsedVtu parsed;
    REQUIRE(text.rfind("<?xml version=\"1.0\"?>", 0) == 0);
    REQUIRE(text.find("<VTKFile type=\"UnstructuredGrid\"") != std::string::npos);
    const std::regex pieceRe(R"(<Piece NumberOfPoints=\"(\d+)\" NumberOfCells=\"(\d+)\">)");
    std::smatch piece;
    REQUIRE(std::regex_search(text, piece, pieceRe));
    parsed.numberOfPoints = std::stoi(piece[1]);
    parsed.numberOfCells = std::stoi(piece[2]);
    const auto points = parseArrays(sectionBetween(text, "<Points>", "</Points>"));
    REQUIRE(points.count("Points") == 1);
    parsed.points = points.at("Points").format == "ascii"
        ? asciiValues<double>(points.at("Points").body)
        : binaryDoubles(points.at("Points").body);
    const auto cells = parseArrays(sectionBetween(text, "<Cells>", "</Cells>"));
    parsed.connectivity = intArray(cells.at("connectivity"));
    parsed.offsets = intArray(cells.at("offsets"));
    parsed.types = intArray(cells.at("types"));
    parsed.pointData = parseArrays(sectionBetween(text, "<PointData>", "</PointData>"), &parsed.pointDataOrder);
    parsed.cellData = parseArrays(sectionBetween(text, "<CellData>", "</CellData>"));
    parsed.fieldData = parseArrays(sectionBetween(text, "<FieldData>", "</FieldData>"));
    return parsed;
}

} // namespace

TEST_CASE("VtuWriter writes an XML unstructured grid that an independent reader recovers",
          "[vtu_writer]")
{
    const DeviceMesh mesh = makeTwoTriangleMesh();
    VtuWriter writer(mesh);
    writer.addPointScalar("eDensity", "cm^-3", {1.0e10, 2.0e10, 3.0e10, 4.0e10});
    writer.addPointVector("ElectricField", "V/cm", {
        Point3{1.0, 0.0, 0.0}, Point3{0.0, 2.0, 0.0}, Point3{3.0, 4.0, 0.0}, Point3{-1.0, -2.0, 0.0}});
    writer.addCellScalar("CellMarker", "1", {7.0, 8.0});
    writer.addFieldString("VelaVersion", "0.1.0-test");
    writer.addFieldValue("Bias_V", 0.25);
    writer.addFieldValue("goal_reached", 0.0);

    const ParsedVtu parsed = parse(writer.render());
    REQUIRE(parsed.numberOfPoints == 4);
    REQUIRE(parsed.numberOfCells == 2);
    REQUIRE(parsed.points == std::vector<double>{0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0});
    REQUIRE(parsed.connectivity == std::vector<int>{0, 1, 2, 0, 2, 3});
    REQUIRE(parsed.offsets == std::vector<int>{3, 6});
    REQUIRE(parsed.types == std::vector<int>{5, 5});

    REQUIRE(parsed.pointDataOrder == std::vector<std::string>{"eDensity", "ElectricField"});
    const auto& density = parsed.pointData.at("eDensity");
    REQUIRE(density.type == "Float64");
    REQUIRE(density.format == "ascii");
    REQUIRE(density.components == 1);
    REQUIRE(asciiValues<double>(density.body) == std::vector<double>{1.0e10, 2.0e10, 3.0e10, 4.0e10});
    const auto& field = parsed.pointData.at("ElectricField");
    REQUIRE(field.components == 3);
    REQUIRE(asciiValues<double>(field.body) == std::vector<double>{1, 0, 0, 0, 2, 0, 3, 4, 0, -1, -2, 0});

    REQUIRE(parsed.cellData.at("region_id").type == "Int32");
    REQUIRE(asciiValues<int>(parsed.cellData.at("region_id").body) == std::vector<int>{0, 1});
    REQUIRE(asciiValues<double>(parsed.cellData.at("CellMarker").body) == std::vector<double>{7.0, 8.0});

    REQUIRE(stringFieldValue(parsed.fieldData.at("VelaVersion")) == "0.1.0-test");
    REQUIRE(asciiValues<double>(parsed.fieldData.at("Bias_V").body) == std::vector<double>{0.25});
    REQUIRE(asciiValues<double>(parsed.fieldData.at("goal_reached").body) == std::vector<double>{0.0});
    REQUIRE(stringFieldValue(parsed.fieldData.at(VtuWriter::unitArrayName("eDensity"))) == "cm^-3");
    REQUIRE(stringFieldValue(parsed.fieldData.at(VtuWriter::unitArrayName("ElectricField"))) == "V/cm");
    REQUIRE(stringFieldValue(parsed.fieldData.at(VtuWriter::unitArrayName("CellMarker"))) == "1");
}

TEST_CASE("VtuWriter binary encoding round-trips Float64 arrays through base64", "[vtu_writer]")
{
    const DeviceMesh mesh = makeTwoTriangleMesh();
    VtuWriter writer(mesh, VtuWriter::Encoding::Binary);
    const std::vector<Real> values{-1.5, 2.25, 1.0e-30, 6.02214076e23};
    writer.addPointScalar("ElectrostaticPotential", "V", values);
    const ParsedVtu parsed = parse(writer.render());
    const auto& array = parsed.pointData.at("ElectrostaticPotential");
    REQUIRE(array.format == "binary");
    const std::vector<double> decoded = binaryDoubles(array.body);
    REQUIRE(decoded.size() == values.size());
    for (std::size_t i = 0; i < values.size(); ++i)
        REQUIRE(decoded[i] == values[i]);
    REQUIRE(parsed.points.size() == 12);
    REQUIRE(parsed.points[3] == 1.0);
}

TEST_CASE("VtuWriter rejects arrays whose length does not match the mesh", "[vtu_writer]")
{
    const DeviceMesh mesh = makeTwoTriangleMesh();
    VtuWriter writer(mesh);
    REQUIRE_THROWS_AS(writer.addPointScalar("eDensity", "cm^-3", {1.0, 2.0}), std::invalid_argument);
    REQUIRE_THROWS_AS(writer.addCellScalar("x", "1", {1.0}), std::invalid_argument);
    REQUIRE_THROWS_AS(writer.addPointVector("ElectricField", "V/cm", {Point3::Zero()}), std::invalid_argument);
}

TEST_CASE("TDR field catalogue lists unique names with manual units", "[vtu_writer][tdr_catalogue]")
{
    const auto& catalogue = tdrFieldCatalogue();
    REQUIRE(catalogue.size() >= 40);
    std::vector<std::string> names;
    for (const TdrFieldSpec& spec : catalogue) {
        REQUIRE_FALSE(spec.tdrName.empty());
        REQUIRE_FALSE(spec.velaSource.empty());
        REQUIRE_FALSE(spec.unit.empty());
        REQUIRE((spec.components == 1 || spec.components == 3));
        names.push_back(spec.tdrName);
    }
    std::sort(names.begin(), names.end());
    REQUIRE(std::adjacent_find(names.begin(), names.end()) == names.end());

    const std::map<std::string, std::pair<std::string, int>> expected = {
        {"ElectrostaticPotential", {"V", 1}},
        {"eQuasiFermiPotential", {"V", 1}},
        {"eDensity", {"cm^-3", 1}},
        {"DopingConcentration", {"cm^-3", 1}},
        {"DonorConcentration", {"cm^-3", 1}},
        {"SpaceCharge", {"cm^-3", 1}},
        {"ElectricField", {"V/cm", 3}},
        {"eCurrentDensity", {"A/cm^2", 3}},
        {"TotalCurrentDensity", {"A/cm^2", 3}},
        {"eMobility", {"cm^2/(V*s)", 1}},
        {"eVelocity", {"cm/s", 1}},
        {"srhRecombination", {"cm^-3*s^-1", 1}},
        {"AvalancheGeneration", {"cm^-3*s^-1", 1}},
        {"eAlphaAvalanche", {"cm^-1", 1}},
        {"eIonIntegral", {"1", 1}},
        {"BandGap", {"eV", 1}},
        {"ConductionBandEnergy", {"eV", 1}},
        {"eEparallel", {"V/cm", 1}},
        {"eGradQuasiFermi", {"V/cm", 3}},
    };
    for (const auto& [name, unitAndComponents] : expected) {
        const TdrFieldSpec* spec = tdrFieldSpec(name);
        REQUIRE(spec != nullptr);
        REQUIRE(spec->unit == unitAndComponents.first);
        REQUIRE(spec->components == unitAndComponents.second);
    }
    REQUIRE(tdrFieldSpec("eTemperature") == nullptr);
    REQUIRE(tdrFieldSpec("LatticeTemperature") == nullptr);
    REQUIRE(tdrFieldSpec("eQuantumPotential") == nullptr);
}

TEST_CASE("TDR field selection resolves Plot requests, vector suffixes and unavailable names",
          "[vtu_writer][tdr_catalogue]")
{
    TdrFieldSelection selection({"eDensity", "ElectricField/Vector", "eCurrentDensity",
                                 "TotalCurrentDensity", "TotalCurrentDensity/Vector",
                                 "eTemperature", "eDensity/Vector"});
    REQUIRE_FALSE(selection.selectsAll());
    REQUIRE(selection.wantsScalar("eDensity"));
    REQUIRE(selection.wantsVector("eDensity"));  // requested (and reported unavailable below)
    REQUIRE(selection.wantsVector("ElectricField"));
    REQUIRE_FALSE(selection.wantsScalar("ElectricField"));
    REQUIRE(selection.wantsScalar("eCurrentDensity"));
    REQUIRE_FALSE(selection.wantsVector("eCurrentDensity"));
    REQUIRE(selection.wantsScalar("TotalCurrentDensity"));
    REQUIRE(selection.wantsVector("TotalCurrentDensity"));
    REQUIRE_FALSE(selection.wants("hDensity"));
    // eTemperature is not provided; eDensity/Vector asks for a vector form of a scalar.
    REQUIRE(selection.unavailable() == std::vector<std::string>{"eTemperature", "eDensity/Vector"});

    selection.markProvided("eDensity", false);
    selection.markProvided("ElectricField", true);
    selection.markProvided("eCurrentDensity", false);
    REQUIRE(selection.undelivered() == std::vector<std::string>{"TotalCurrentDensity", "TotalCurrentDensity/Vector"});

    const TdrFieldSelection everything;
    REQUIRE(everything.selectsAll());
    REQUIRE(everything.wantsScalar("eDensity"));
    REQUIRE(everything.wantsVector("ElectricField"));
    REQUIRE_FALSE(everything.wantsScalar("ElectricField"));
    REQUIRE(everything.unavailable().empty());
}
