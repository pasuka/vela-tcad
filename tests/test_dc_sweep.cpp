#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/simulation/DCSweep.h"
#include "vela/simulation/DCSweepStepControl.h"
#include "vela/post/TerminalCharge.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

using namespace vela;

namespace {

constexpr int kMaxDirCreationAttempts = 8;

struct ScopedDirectoryCleanup {
    std::filesystem::path path;

    ~ScopedDirectoryCleanup()
    {
        std::error_code ec;
        std::filesystem::remove_all(path, ec);
        if (ec)
            std::cerr << "Failed to remove test directory '" << path.string()
                      << "': " << ec.message() << '\n';
    }
};

std::filesystem::path makeUniqueSweepDir()
{
    const auto base = std::filesystem::temp_directory_path();
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto tid = std::hash<std::thread::id>{}(std::this_thread::get_id());
    std::mt19937_64 rng(static_cast<std::mt19937_64::result_type>(stamp ^ tid));
    std::uniform_int_distribution<unsigned long long> dist;

    for (int attempt = 0; attempt < kMaxDirCreationAttempts; ++attempt) {
        const auto dir = base /
            ("vela_dc_sweep_test_" + std::to_string(stamp) + "_" + std::to_string(dist(rng)));
        if (!std::filesystem::exists(dir))
            return dir;
    }

    throw std::runtime_error("Failed to create a unique temp directory for DCSweep test.");
}

std::filesystem::path writePNMesh(const std::filesystem::path& dir)
{
    nlohmann::json mesh = {
        {"nodes", {
            {{"id", 0}, {"x", 0.0e-6}, {"y", 0.0e-6}},
            {{"id", 1}, {"x", 1.0e-6}, {"y", 0.0e-6}},
            {{"id", 2}, {"x", 1.0e-6}, {"y", 1.0e-6}},
            {{"id", 3}, {"x", 0.0e-6}, {"y", 1.0e-6}}
        }},
        {"triangles", {
            {{"id", 0}, {"region_id", 0}, {"node_ids", {0, 1, 2}}},
            {{"id", 1}, {"region_id", 1}, {"node_ids", {0, 2, 3}}}
        }},
        {"regions", {
            {{"id", 0}, {"name", "n_region"}, {"material", "Si"}, {"cell_ids", {0}}},
            {{"id", 1}, {"name", "p_region"}, {"material", "Si"}, {"cell_ids", {1}}}
        }},
        {"contacts", {
            {{"id", 0}, {"name", "anode"}, {"region_id", 1}, {"node_ids", {0, 3}}},
            {{"id", 1}, {"name", "cathode"}, {"region_id", 0}, {"node_ids", {1, 2}}}
        }}
    };

    const auto meshPath = dir / "pn_mesh.json";
    std::ofstream(meshPath) << mesh.dump(2);
    return meshPath;
}

std::filesystem::path writePNMeshMicrometers(const std::filesystem::path& dir)
{
    nlohmann::json mesh = {
        {"nodes", {
            {{"id", 0}, {"x", 0.0}, {"y", 0.0}},
            {{"id", 1}, {"x", 1.0}, {"y", 0.0}},
            {{"id", 2}, {"x", 1.0}, {"y", 1.0}},
            {{"id", 3}, {"x", 0.0}, {"y", 1.0}}
        }},
        {"triangles", {
            {{"id", 0}, {"region_id", 0}, {"node_ids", {0, 1, 2}}},
            {{"id", 1}, {"region_id", 1}, {"node_ids", {0, 2, 3}}}
        }},
        {"regions", {
            {{"id", 0}, {"name", "n_region"}, {"material", "Si"}, {"cell_ids", {0}}},
            {{"id", 1}, {"name", "p_region"}, {"material", "Si"}, {"cell_ids", {1}}}
        }},
        {"contacts", {
            {{"id", 0}, {"name", "anode"}, {"region_id", 1}, {"node_ids", {0, 3}}},
            {{"id", 1}, {"name", "cathode"}, {"region_id", 0}, {"node_ids", {1, 2}}}
        }}
    };

    const auto meshPath = dir / "pn_mesh_um.json";
    std::ofstream(meshPath) << mesh.dump(2);
    return meshPath;
}

nlohmann::json baseSweepConfig(const std::filesystem::path& dir,
                               const std::filesystem::path& meshPath,
                               const std::filesystem::path& csvPath)
{
    return {
        {"mesh_file", meshPath.string()},
        {"output_csv", csvPath.string()},
        {"doping", {
            {{"region", "n_region"}, {"donors", 1.0e23}, {"acceptors", 0.0}},
            {{"region", "p_region"}, {"donors", 0.0}, {"acceptors", 1.0e23}}
        }},
        {"contacts", {
            {{"name", "anode"}, {"bias", 0.0}},
            {{"name", "cathode"}, {"bias", 0.0}}
        }},
        {"solver", {
            {"max_iter", 80},
            {"reltol", 1.0e-5},
            {"damping_psi", 0.5}
        }},
        {"sweep", {
            {"contact", "anode"},
            {"start", 0.0},
            {"stop", 0.5},
            {"step", 0.25},
            {"current_contact", "anode"},
            {"write_vtk", true},
            {"vtk_prefix", (dir / "pn_sweep").string()}
        }}
    };
}

std::filesystem::path writeSweepConfig(const std::filesystem::path& dir,
                                       const std::filesystem::path& meshPath,
                                       const std::filesystem::path& csvPath,
                                       const nlohmann::json& sweepOverrides = {},
                                       const nlohmann::json& solverOverrides = {})
{
    nlohmann::json cfg = baseSweepConfig(dir, meshPath, csvPath);
    for (auto it = sweepOverrides.begin(); it != sweepOverrides.end(); ++it)
        cfg["sweep"][it.key()] = it.value();
    for (auto it = solverOverrides.begin(); it != solverOverrides.end(); ++it)
        cfg["solver"][it.key()] = it.value();

    const auto cfgPath = dir / "pn_sweep.json";
    std::ofstream(cfgPath) << cfg.dump(2);
    return cfgPath;
}

std::filesystem::path writeUnitScalingSweepConfig(
    const std::filesystem::path& dir,
    const std::filesystem::path& meshPath,
    const std::filesystem::path& csvPath,
    const nlohmann::json& sweepOverrides = {},
    const nlohmann::json& solverOverrides = {})
{
    nlohmann::json cfg = baseSweepConfig(dir, meshPath, csvPath);
    cfg["scaling"] = {{"mode", "unit_scaling"}};
    cfg["doping"] = {
        {{"region", "n_region"}, {"donors", 1.0e17}, {"acceptors", 0.0}},
        {{"region", "p_region"}, {"donors", 0.0}, {"acceptors", 1.0e17}}
    };
    for (auto it = sweepOverrides.begin(); it != sweepOverrides.end(); ++it)
        cfg["sweep"][it.key()] = it.value();
    for (auto it = solverOverrides.begin(); it != solverOverrides.end(); ++it)
        cfg["solver"][it.key()] = it.value();

    const auto cfgPath = dir / "pn_sweep_unit_scaling.json";
    std::ofstream(cfgPath) << cfg.dump(2);
    return cfgPath;
}

std::vector<std::vector<std::string>> readCsvRows(const std::filesystem::path& csvPath)
{
    std::ifstream input(csvPath);
    std::vector<std::vector<std::string>> rows;
    std::string line;
    while (std::getline(input, line)) {
        std::vector<std::string> columns;
        std::stringstream ss(line);
        std::string column;
        while (std::getline(ss, column, ','))
            columns.push_back(column);
        rows.push_back(columns);
    }
    return rows;
}

void writeNodeDopingCsv(const std::filesystem::path& csvPath,
                        const std::vector<std::tuple<Index, Real, Real>>& rows)
{
    std::ofstream output(csvPath);
    output << "node_id,donors_cm3,acceptors_cm3\n";
    for (const auto& [nodeId, donors, acceptors] : rows)
        output << nodeId << ',' << donors << ',' << acceptors << '\n';
}

void convertMeshToMicrometersInPlace(const std::filesystem::path& meshPath)
{
    std::ifstream input(meshPath);
    nlohmann::json mesh;
    input >> mesh;
    for (auto& node : mesh["nodes"]) {
        node["x"] = node.at("x").get<Real>() * 1.0e6;
        node["y"] = node.at("y").get<Real>() * 1.0e6;
    }
    std::ofstream(meshPath) << mesh.dump(2);
}

void convertDopingToCm3InPlace(nlohmann::json& cfg)
{
    for (auto& region : cfg["doping"]) {
        region["donors"] = region.at("donors").get<Real>() / 1.0e6;
        region["acceptors"] = region.at("acceptors").get<Real>() / 1.0e6;
        if (region.contains("fixed_charge_m3"))
            region["fixed_charge_m3"] = region.at("fixed_charge_m3").get<Real>() / 1.0e6;
    }
}

std::size_t csvColumnIndex(const std::vector<std::string>& header,
                           const std::string& column)
{
    const auto it = std::find(header.begin(), header.end(), column);
    REQUIRE(it != header.end());
    return static_cast<std::size_t>(std::distance(header.begin(), it));
}

Real csvReal(const std::vector<std::string>& row, std::size_t column)
{
    REQUIRE(column < row.size());
    return std::stod(row.at(column));
}



DeviceMesh makeTwoRegionUnitSquareMesh()
{
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0.0, 0.0, 0.0});
    mesh.addNode(Node{1, 1.0, 0.0, 0.0});
    mesh.addNode(Node{2, 1.0, 1.0, 0.0});
    mesh.addNode(Node{3, 0.0, 1.0, 0.0});
    mesh.addCell(Cell{0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.addCell(Cell{1, CellType::Tri3, 1, {0, 2, 3}});
    mesh.addRegion(Region{0, "right", "Si", {0}});
    mesh.addRegion(Region{1, "left", "Si", {1}});
    mesh.addContact(Contact{0, "left_contact", 1, {0, 3}});
    mesh.addContact(Contact{1, "right_contact", 0, {1, 2}});
    mesh.buildEdges();
    return mesh;
}

DDSolution uniformCarrierSolution(Index numNodes, Real electrons, Real holes)
{
    DDSolution solution;
    solution.psi = VectorXd::Zero(static_cast<int>(numNodes));
    solution.phin = VectorXd::Zero(static_cast<int>(numNodes));
    solution.phip = VectorXd::Zero(static_cast<int>(numNodes));
    solution.n = VectorXd::Constant(static_cast<int>(numNodes), electrons);
    solution.p = VectorXd::Constant(static_cast<int>(numNodes), holes);
    solution.converged = true;
    return solution;
}

Real runMosExampleDrainCurrentAtGate(const std::string& exampleName, Real gateBias, Real drainBias)
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    const std::filesystem::path src = std::filesystem::path(VELA_SOURCE_DIR) / "examples" / exampleName;
    std::filesystem::copy(src, dir, std::filesystem::copy_options::recursive);
    std::filesystem::create_directories(dir / "outputs");

    const auto cfgPath = dir / "simulation_iv.json";
    std::ifstream input(cfgPath);
    nlohmann::json cfg;
    input >> cfg;
    bool foundGateContact = false;
    for (auto& contact : cfg["contacts"]) {
        if (contact.at("name").get<std::string>() == "gate") {
            contact["bias"] = gateBias;
            foundGateContact = true;
        }
    }
    REQUIRE(foundGateContact);
    cfg["output_csv"] = "outputs/mos_idvd_test.csv";
    cfg["sweep"]["start"] = drainBias;
    cfg["sweep"]["stop"] = drainBias;
    cfg["sweep"]["step"] = (drainBias >= 0.0) ? 0.05 : -0.05;
    cfg["sweep"]["write_vtk"] = false;
    std::ofstream(cfgPath) << cfg.dump(2);

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 1);
    const DCSweepPoint& point = result.points.front();
    REQUIRE(point.converged);
    REQUIRE(std::isfinite(point.totalCurrent));
    REQUIRE(point.iterations > 0);
    return point.totalCurrent;
}

} // namespace

TEST_CASE("DCSweep: PN diode forward sweep writes CSV and finite monotonic IV data", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "iv.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath);

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    const std::vector<DCSweepPoint>& points = result.points;

    REQUIRE(points.size() == 3);
    REQUIRE(result.mesh.numNodes() == 4);
    REQUIRE(result.mesh.lastGeometryBuildReport().totalCells == 2);
    REQUIRE(std::filesystem::exists(csvPath));
    REQUIRE(std::filesystem::file_size(csvPath) > 0);

    for (const DCSweepPoint& point : points) {
        REQUIRE(point.converged);
        REQUIRE(std::isfinite(point.electronCurrent));
        REQUIRE(std::isfinite(point.holeCurrent));
        REQUIRE(std::isfinite(point.totalCurrent));
    }

    REQUIRE(points[0].attemptedStep == Catch::Approx(0.0));
    REQUIRE(points[1].attemptedStep == Catch::Approx(0.25));
    REQUIRE(points[1].acceptedStep == Catch::Approx(0.25));
    REQUIRE(points[1].retryCount == 0);
    REQUIRE(std::abs(points.back().totalCurrent) >= std::abs(points.front().totalCurrent));
    REQUIRE(std::filesystem::exists(dir / "pn_sweep_0000_0V.vtk"));

    const auto rows = readCsvRows(csvPath);
    REQUIRE(rows.front() == std::vector<std::string>{"mode", "bias_contact", "bias_V",
                                                     "current_contact", "current_electron", "current_hole",
                                                     "current_total", "converged", "iterations",
                                                     "step_diagnostics", "validation_diagnostics"});
}

TEST_CASE("DCSweep: unit_scaling CSV appends per-micron currents and V-per-cm field",
          "[dc_sweep][scaling]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "bv_unit_scaling.csv";
    const auto cfgPath = writeUnitScalingSweepConfig(dir, meshPath, csvPath, {
        {"mode", "bv_reverse"},
        {"start", 0.0},
        {"stop", 0.25},
        {"step", 0.25},
        {"write_vtk", false},
        {"breakdown", {
            {"max_electric_field_V_per_m", 1.0e12},
            {"current_jump_ratio", 1.0e12},
            {"non_convergence", true}
        }}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 2);
    REQUIRE(result.points.back().converged);

    const auto rows = readCsvRows(csvPath);
    REQUIRE(rows.size() == 3);
    const auto& header = rows.front();
    const std::size_t currentTotal = csvColumnIndex(header, "current_total");
    const std::size_t currentElectron = csvColumnIndex(header, "current_electron");
    const std::size_t currentHole = csvColumnIndex(header, "current_hole");
    const std::size_t maxField = csvColumnIndex(header, "max_electric_field_V_per_m");
    const std::size_t currentTotalUm = csvColumnIndex(header, "current_total_A_per_um");
    const std::size_t currentElectronUm = csvColumnIndex(header, "current_electron_A_per_um");
    const std::size_t currentHoleUm = csvColumnIndex(header, "current_hole_A_per_um");
    const std::size_t maxFieldCm = csvColumnIndex(header, "max_electric_field_V_per_cm");

    for (std::size_t r = 1; r < rows.size(); ++r) {
        const auto& row = rows.at(r);
        REQUIRE(csvReal(row, currentTotalUm) ==
                Catch::Approx(csvReal(row, currentTotal) / 1.0e6).epsilon(1.0e-12));
        REQUIRE(csvReal(row, currentElectronUm) ==
                Catch::Approx(csvReal(row, currentElectron) / 1.0e6).epsilon(1.0e-12));
        REQUIRE(csvReal(row, currentHoleUm) ==
                Catch::Approx(csvReal(row, currentHole) / 1.0e6).epsilon(1.0e-12));
        REQUIRE(csvReal(row, maxFieldCm) ==
                Catch::Approx(csvReal(row, maxField) / 100.0).epsilon(1.0e-12));
    }
}

TEST_CASE("DCSweep: PN forward IV unit_scaling remains physically equivalent to legacy SI",
          "[dc_sweep][scaling][dd_gummel]")
{
    const auto legacyDir = makeUniqueSweepDir();
    const auto scaledDir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanupLegacy{legacyDir};
    const ScopedDirectoryCleanup cleanupScaled{scaledDir};
    std::filesystem::create_directories(legacyDir);
    std::filesystem::create_directories(scaledDir);

    const auto legacyMesh = writePNMesh(legacyDir);
    const auto scaledMesh = writePNMeshMicrometers(scaledDir);
    const auto legacyCsv = legacyDir / "iv_legacy.csv";
    const auto scaledCsv = scaledDir / "iv_unit_scaling.csv";
    const auto legacyCfg = writeSweepConfig(legacyDir, legacyMesh, legacyCsv, {
        {"start", 0.0}, {"stop", 0.5}, {"step", 0.25}, {"write_vtk", false}
    });
    const auto scaledCfg = writeUnitScalingSweepConfig(scaledDir, scaledMesh, scaledCsv, {
        {"start", 0.0}, {"stop", 0.5}, {"step", 0.25}, {"write_vtk", false}
    });

    DCSweep sweep;
    const DCSweepResult legacy = sweep.runWithResult(legacyCfg.string());
    const DCSweepResult scaled = sweep.runWithResult(scaledCfg.string());

    REQUIRE(legacy.points.size() == scaled.points.size());
    REQUIRE_FALSE(legacy.points.empty());
    for (std::size_t i = 0; i < legacy.points.size(); ++i) {
        REQUIRE(legacy.points[i].converged);
        REQUIRE(scaled.points[i].converged);
        REQUIRE(std::abs(scaled.points[i].totalCurrent)
                == Catch::Approx(std::abs(legacy.points[i].totalCurrent)).epsilon(5.0e-2));
    }

    const Real legacyEnd = std::abs(legacy.points.back().totalCurrent);
    const Real scaledEnd = std::abs(scaled.points.back().totalCurrent);
    REQUIRE(scaledEnd == Catch::Approx(legacyEnd).epsilon(5.0e-2));
}

TEST_CASE("DCSweep reads node_doping_file before region averages", "[dc_sweep][doping]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "node_doping_iv.csv";
    const auto cfgPath = writeUnitScalingSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.1},
        {"write_vtk", false}
    });

    std::ifstream input(cfgPath);
    nlohmann::json cfg;
    input >> cfg;
    cfg["node_doping_file"] = "doping.csv";
    cfg["doping"] = {
        {{"region", "n_region"}, {"donors", 0.0}, {"acceptors", 0.0}},
        {{"region", "p_region"}, {"donors", 0.0}, {"acceptors", 0.0}}
    };
    std::ofstream(cfgPath) << cfg.dump(2);

    writeNodeDopingCsv(dir / "doping.csv", {
        {0, 0.0, 1.0e17},
        {1, 1.0e17, 0.0},
        {2, 1.0e17, 0.0},
        {3, 0.0, 1.0e17},
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);
    REQUIRE(std::isfinite(result.points.front().electronCurrent));
    REQUIRE(std::isfinite(result.points.front().holeCurrent));
    REQUIRE(std::isfinite(result.points.front().totalCurrent));

    writeNodeDopingCsv(dir / "doping.csv", {
        {0, 0.0, 1.0e17},
        {7, 1.0e17, 0.0},
    });
    REQUIRE_THROWS_WITH(
        sweep.runWithResult(cfgPath.string()),
        Catch::Matchers::ContainsSubstring(
            "DCSweep: node_doping_file references missing node id 7"));

    {
        std::ofstream malformed(dir / "doping.csv");
        malformed << "node_id,donors_cm3,acceptors_cm3\n";
        malformed << "1abc,1.0e17,0.0\n";
    }
    REQUIRE_THROWS_WITH(
        sweep.runWithResult(cfgPath.string()),
        Catch::Matchers::ContainsSubstring(
            "DCSweep: node_doping_file has invalid node id '1abc'"));

    writeNodeDopingCsv(dir / "doping.csv", {
        {0, 0.0, 1.0e17},
        {1, 1.0e17, 0.0},
        {2, 1.0e17, 0.0},
    });
    REQUIRE_THROWS_WITH(
        sweep.runWithResult(cfgPath.string()),
        Catch::Matchers::ContainsSubstring(
            "DCSweep: node_doping_file missing row for node id 3"));

    writeNodeDopingCsv(dir / "doping.csv", {
        {0, 0.0, 1.0e17},
        {1, 1.0e17, 0.0},
        {1, 1.0e17, 0.0},
        {2, 1.0e17, 0.0},
        {3, 0.0, 1.0e17},
    });
    REQUIRE_THROWS_WITH(
        sweep.runWithResult(cfgPath.string()),
        Catch::Matchers::ContainsSubstring(
            "DCSweep: node_doping_file has duplicate row for node id 1"));

    {
        std::ofstream malformed(dir / "doping.csv");
        malformed << "node_id,donors_cm3,acceptors_cm3\n";
        malformed << "0,0.0,1.0e17\n";
        malformed << "1,1.0e17abc,0.0\n";
        malformed << "2,1.0e17,0.0\n";
        malformed << "3,0.0,1.0e17\n";
    }
    REQUIRE_THROWS_WITH(
        sweep.runWithResult(cfgPath.string()),
        Catch::Matchers::ContainsSubstring(
            "DCSweep: node_doping_file has invalid donors_cm3 '1.0e17abc' for node id 1"));

    {
        std::ofstream malformed(dir / "doping.csv");
        malformed << "node_id,donors_cm3,acceptors_cm3\n";
        malformed << "0,0.0,1.0e17\n";
        malformed << "1,nan,0.0\n";
        malformed << "2,1.0e17,0.0\n";
        malformed << "3,0.0,1.0e17\n";
    }
    REQUIRE_THROWS_WITH(
        sweep.runWithResult(cfgPath.string()),
        Catch::Matchers::ContainsSubstring(
            "DCSweep: node_doping_file has non-finite donors_cm3 'nan' for node id 1"));

    {
        std::ofstream malformed(dir / "doping.csv");
        malformed << "node_id,donors_cm3,acceptors_cm3\n";
        malformed << "0,0.0,1.0e17\n";
        malformed << "\"1\",1.0e17,0.0\n";
        malformed << "2,1.0e17,0.0\n";
        malformed << "3,0.0,1.0e17\n";
    }
    REQUIRE_THROWS_WITH(
        sweep.runWithResult(cfgPath.string()),
        Catch::Matchers::ContainsSubstring(
            "DCSweep: node_doping_file does not support quoted fields"));
}

TEST_CASE("DCSweep: NMOS and PMOS unit_scaling low-bias smoke sweeps converge",
          "[dc_sweep][scaling][dd_gummel]")
{
    const std::vector<std::string> devices = {"nmos2d_dd", "pmos2d_dd"};

    for (const std::string& device : devices) {
        INFO(device);
        const auto dir = makeUniqueSweepDir();
        const ScopedDirectoryCleanup cleanup{dir};
        const std::filesystem::path src = std::filesystem::path(VELA_SOURCE_DIR) / "examples" / device;
        std::filesystem::copy(src, dir, std::filesystem::copy_options::recursive);
        std::filesystem::create_directories(dir / "outputs");

        const auto meshPath = dir / "mesh.json";
        const auto cfgPath = dir / "simulation_iv.json";
        convertMeshToMicrometersInPlace(meshPath);

        std::ifstream cfgIn(cfgPath);
        nlohmann::json cfg;
        cfgIn >> cfg;
        cfg["scaling"] = { {"mode", "unit_scaling"} };
        convertDopingToCm3InPlace(cfg);
        cfg["sweep"]["write_vtk"] = false;
        cfg["output_csv"] = (dir / "outputs" / (device + "_unit_scaling_iv.csv")).string();
        std::ofstream(cfgPath) << cfg.dump(2);

        DCSweep sweep;
        const DCSweepResult result = sweep.runWithResult(cfgPath.string());
        REQUIRE_FALSE(result.points.empty());
        for (const DCSweepPoint& point : result.points) {
            REQUIRE(point.converged);
            REQUIRE(std::isfinite(point.totalCurrent));
        }
    }
}

TEST_CASE("DCSweep: unit_scaling CV CSV appends per-micron charge and capacitance",
          "[dc_sweep][scaling]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "cv_unit_scaling.csv";
    const auto cfgPath = writeUnitScalingSweepConfig(dir, meshPath, csvPath, {
        {"mode", "cv_quasistatic"},
        {"start", 0.0},
        {"stop", 0.25},
        {"step", 0.25},
        {"write_vtk", false},
        {"terminal_charge", {
            {"contact", "anode"},
            {"regions", {"p_region"}},
            {"per_meter", true}
        }}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 2);
    REQUIRE(result.points.back().converged);

    const auto rows = readCsvRows(csvPath);
    REQUIRE(rows.size() == 3);
    const auto& header = rows.front();
    const std::size_t charge = csvColumnIndex(header, "charge_C_per_m");
    const std::size_t capacitance = csvColumnIndex(header, "capacitance_F_per_m");
    const std::size_t chargeUm = csvColumnIndex(header, "charge_C_per_um");
    const std::size_t capacitanceUm = csvColumnIndex(header, "capacitance_F_per_um");

    for (std::size_t r = 1; r < rows.size(); ++r) {
        const auto& row = rows.at(r);
        REQUIRE(csvReal(row, chargeUm) ==
                Catch::Approx(csvReal(row, charge) / 1.0e6).epsilon(1.0e-12));
        REQUIRE(csvReal(row, capacitanceUm) ==
                Catch::Approx(csvReal(row, capacitance) / 1.0e6).epsilon(1.0e-12));
    }
}


TEST_CASE("TerminalCharge: region selections use region-local cell volume", "[terminal_charge]")
{
    DeviceMesh mesh = makeTwoRegionUnitSquareMesh();
    DopingModel doping(mesh.numNodes());
    const DDSolution solution = uniformCarrierSolution(mesh.numNodes(), 0.0, 1.0);

    TerminalChargeConfig config;
    config.regions = {"left"};
    config.includeIonizedDopants = false;

    const TerminalChargeResult result = TerminalCharge::compute(mesh, doping, solution, config);

    REQUIRE(result.charge / constants::q == Catch::Approx(0.5));
}

TEST_CASE("TerminalCharge: unknown region selections are rejected", "[terminal_charge]")
{
    DeviceMesh mesh = makeTwoRegionUnitSquareMesh();
    DopingModel doping(mesh.numNodes());
    const DDSolution solution = uniformCarrierSolution(mesh.numNodes(), 0.0, 1.0);

    TerminalChargeConfig config;
    config.regions = {"missing"};

    REQUIRE_THROWS_AS(TerminalCharge::compute(mesh, doping, solution, config),
                      std::invalid_argument);
}


TEST_CASE("DCSweep: curve output schemas distinguish IV, CV, and BV modes", "[dc_sweep][curve]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);

    SECTION("CV quasistatic adds terminal charge and capacitance columns")
    {
        const auto csvPath = dir / "cv.csv";
        const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
            {"mode", "cv_quasistatic"},
            {"start", 0.0},
            {"stop", 0.25},
            {"step", 0.25},
            {"write_vtk", false},
            {"terminal_charge", {
                {"contact", "anode"},
                {"regions", {"p_region"}},
                {"per_meter", true}
            }}
        });

        DCSweep sweep;
        const DCSweepResult result = sweep.runWithResult(cfgPath.string());
        REQUIRE(result.points.size() == 2);
        REQUIRE(result.points[0].converged);
        REQUIRE(result.points[1].converged);
        REQUIRE(std::isfinite(result.points[1].terminalCharge));
        REQUIRE(std::isfinite(result.points[1].capacitance));

        const auto rows = readCsvRows(csvPath);
        REQUIRE(rows.front() == std::vector<std::string>{"mode", "bias_contact", "bias_V",
                                                         "current_contact", "current_electron", "current_hole",
                                                         "current_total", "converged", "iterations",
                                                         "step_diagnostics", "validation_diagnostics", "charge_C_per_m",
                                                         "capacitance_F_per_m"});
        REQUIRE(rows.at(1).at(0) == "cv_quasistatic");
    }

    SECTION("CV quasistatic adds multi-terminal charge and capacitance columns")
    {
        const auto csvPath = dir / "cv_multi.csv";
        const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
            {"mode", "cv_quasistatic"},
            {"start", 0.0},
            {"stop", 0.25},
            {"step", 0.25},
            {"write_vtk", false},
            {"terminal_charges", {
                {
                    {"name", "gate"},
                    {"contact", "anode"},
                    {"regions", {"p_region"}},
                    {"per_meter", true}
                },
                {
                    {"name", "source"},
                    {"contact", "cathode"},
                    {"regions", {"n_region"}},
                    {"per_meter", true}
                },
                {
                    {"name", "substrate"},
                    {"contact", "cathode"},
                    {"regions", {"n_region"}},
                    {"per_meter", true}
                }
            }}
        });

        DCSweep sweep;
        const DCSweepResult result = sweep.runWithResult(cfgPath.string());
        REQUIRE(result.points.size() == 2);
        REQUIRE(result.points[1].terminalChargeValues.size() == 3);
        REQUIRE(result.points[1].terminalCapacitanceValues.size() == 3);
        REQUIRE(std::isfinite(result.points[1].terminalChargeValues[0].second));
        REQUIRE(std::isfinite(result.points[1].terminalCapacitanceValues[0].second));
        REQUIRE(std::isfinite(result.points[1].extraFields[0].second));

        const auto rows = readCsvRows(csvPath);
        REQUIRE(rows.front() == std::vector<std::string>{"mode", "bias_contact", "bias_V",
                                                         "current_contact", "current_electron", "current_hole",
                                                         "current_total", "converged", "iterations",
                                                         "step_diagnostics", "validation_diagnostics", "charge_C_per_m",
                                                         "capacitance_F_per_m", "charge_gate_C_per_m",
                                                         "capacitance_Canode_gate_F_per_m", "charge_source_C_per_m",
                                                         "capacitance_Canode_source_F_per_m", "charge_substrate_C_per_m",
                                                         "capacitance_Canode_substrate_F_per_m"});
    }

    SECTION("CV quasistatic rejects an empty terminal_charges array")
    {
        const auto csvPath = dir / "cv_empty_multi.csv";
        const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
            {"mode", "cv_quasistatic"},
            {"start", 0.0},
            {"stop", 0.25},
            {"step", 0.25},
            {"write_vtk", false},
            {"terminal_charges", nlohmann::json::array()}
        });

        DCSweep sweep;
        REQUIRE_THROWS_WITH(sweep.runWithResult(cfgPath.string()),
                            Catch::Matchers::ContainsSubstring("sweep.terminal_charges must not be empty"));
    }

    SECTION("BV reverse adds breakdown diagnostic columns")
    {
        const auto csvPath = dir / "bv.csv";
        const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
            {"mode", "bv_reverse"},
            {"start", 0.0},
            {"stop", 0.25},
            {"step", 0.25},
            {"write_vtk", false},
            {"breakdown", {
                {"max_electric_field_V_per_m", 1.0},
                {"current_jump_ratio", 1.0e12},
                {"non_convergence", true}
            }}
        });

        DCSweep sweep;
        const DCSweepResult result = sweep.runWithResult(cfgPath.string());
        REQUIRE(result.points.size() == 2);
        REQUIRE(result.points.back().converged);
        REQUIRE(result.points.back().breakdownDetected);
        REQUIRE(result.points.back().breakdownCriterion == "max_electric_field");

        const auto rows = readCsvRows(csvPath);
        REQUIRE(rows.front() == std::vector<std::string>{"mode", "bias_contact", "bias_V",
                                                         "current_contact", "current_electron", "current_hole",
                                                         "current_total", "converged", "iterations",
                                                         "step_diagnostics", "validation_diagnostics", "max_electric_field_V_per_m",
                                                         "current_jump_ratio", "breakdown_detected",
                                                         "breakdown_voltage", "criterion", "last_stable_bias",
                                                         "failed_bias", "failure_reason"});
        REQUIRE(rows.at(1).at(0) == "bv_reverse");
    }
}




TEST_CASE("DCSweep: LDMOS BV diagnostic deck writes complete schema", "[dc_sweep][ldmos]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    const std::filesystem::path src = std::filesystem::path(VELA_SOURCE_DIR) / "examples" / "ldmos2d";
    std::filesystem::copy(src, dir, std::filesystem::copy_options::recursive);
    std::filesystem::create_directories(dir / "outputs");

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult((dir / "simulation_bv.json").string());

    REQUIRE(result.points.size() == 3);
    Real previousMaxField = -1.0;
    for (const DCSweepPoint& point : result.points) {
        REQUIRE(point.converged);
        REQUIRE(std::isfinite(point.maxElectricField));
        REQUIRE(point.maxElectricField >= 0.0);
        REQUIRE(point.maxElectricField + 1.0e-9 >= previousMaxField);
        previousMaxField = point.maxElectricField;
    }

    const auto rows = readCsvRows(dir / "outputs" / "ldmos2d_bv.csv");
    REQUIRE(rows.size() == 4);
    REQUIRE(rows.front() == std::vector<std::string>{"mode", "bias_contact", "bias_V",
                                                     "current_contact", "current_electron", "current_hole",
                                                     "current_total", "converged", "iterations",
                                                     "step_diagnostics", "validation_diagnostics", "max_electric_field_V_per_m",
                                                     "current_jump_ratio", "breakdown_detected",
                                                     "breakdown_voltage", "criterion", "last_stable_bias",
                                                     "failed_bias", "failure_reason"});
    REQUIRE(rows.at(1).at(0) == "bv_reverse");
    REQUIRE(rows.at(1).at(1) == "drain");
    REQUIRE(rows.at(1).at(3) == "drain");
}


TEST_CASE("DCSweep: BV reverse start failure records failed diagnostic row", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "bv_nonconvergence.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"mode", "bv_reverse"},
        {"start", 0.0},
        {"stop", -0.5},
        {"step", -0.5},
        {"min_step", 0.2},
        {"max_step", 0.5},
        {"shrink_factor", 0.5},
        {"growth_factor", 1.0},
        {"max_retries", 3},
        {"stop_on_failure", true},
        {"write_vtk", false},
        {"breakdown", {
            {"max_electric_field_V_per_m", 0.0},
            {"current_jump_ratio", 0.0},
            {"non_convergence", true}
        }}
    }, {
        {"max_iter", 0},
        {"reltol", 1.0e-30}
    });

    DCSweep sweep;
    const std::vector<DCSweepPoint> points = sweep.run(cfgPath.string());

    REQUIRE(points.size() == 1);
    const DCSweepPoint& point = points.back();
    REQUIRE_FALSE(point.converged);
    REQUIRE(point.failed);
    REQUIRE_FALSE(point.breakdownDetected);
    REQUIRE(point.breakdownCriterion.empty());
    REQUIRE(point.failedBias == Catch::Approx(0.0));
    REQUIRE(point.lastStableBias == Catch::Approx(0.0));
    REQUIRE(point.failureReason == "non_convergence");

    const auto rows = readCsvRows(csvPath);
    REQUIRE(rows.front() == std::vector<std::string>{"mode", "bias_contact", "bias_V",
                                                     "current_contact", "current_electron", "current_hole",
                                                     "current_total", "converged", "iterations",
                                                     "step_diagnostics", "validation_diagnostics", "max_electric_field_V_per_m",
                                                     "current_jump_ratio", "breakdown_detected",
                                                     "breakdown_voltage", "criterion", "last_stable_bias",
                                                     "failed_bias", "failure_reason"});
    REQUIRE(rows.at(1).at(15).empty());
    REQUIRE(rows.at(1).at(18) == "non_convergence");
}

TEST_CASE("DCSweep: PN diode reverse sweep reaches descending targets", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "reverse.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.5},
        {"stop", 0.0},
        {"step", -0.25},
        {"write_vtk", false}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    const std::vector<DCSweepPoint>& points = result.points;

    REQUIRE(points.size() == 3);
    REQUIRE(result.mesh.numNodes() == 4);
    REQUIRE(result.mesh.lastGeometryBuildReport().totalCells == 2);
    REQUIRE(points[0].voltage == Catch::Approx(0.5));
    REQUIRE(points[1].voltage == Catch::Approx(0.25));
    REQUIRE(points[2].voltage == Catch::Approx(0.0));
    REQUIRE(points[1].attemptedStep == Catch::Approx(-0.25));
    REQUIRE(points[1].acceptedStep == Catch::Approx(-0.25));
}


TEST_CASE("DCSweep step control: invalid direct-call config fails fast", "[dc_sweep]")
{
    const auto attempt = [](Real, Real, int) { return true; };
    const auto record = [](const detail::DCSweepStepControlEvent&) {};

    SECTION("default config has a zero step")
    {
        REQUIRE_THROWS_AS(detail::runDCSweepStepControl({}, attempt, record),
                          std::invalid_argument);
    }

    SECTION("zero maxStep is rejected before attempting a solve")
    {
        detail::DCSweepStepControlConfig cfg;
        cfg.start = 0.0;
        cfg.stop = 1.0;
        cfg.step = 0.25;
        cfg.minStep = 0.125;
        cfg.maxStep = 0.0;
        cfg.growthFactor = 1.0;
        cfg.shrinkFactor = 0.5;
        cfg.maxRetries = 1;

        bool attempted = false;
        REQUIRE_THROWS_AS(
            detail::runDCSweepStepControl(
                cfg,
                [&](Real, Real, int) {
                    attempted = true;
                    return true;
                },
                record),
            std::invalid_argument);
        REQUIRE_FALSE(attempted);
    }

    SECTION("step direction must move toward stop")
    {
        detail::DCSweepStepControlConfig cfg;
        cfg.start = 1.0;
        cfg.stop = 0.0;
        cfg.step = 0.25;
        cfg.minStep = 0.125;
        cfg.maxStep = 0.25;
        cfg.growthFactor = 1.0;
        cfg.shrinkFactor = 0.5;
        cfg.maxRetries = 1;

        REQUIRE_THROWS_AS(detail::runDCSweepStepControl(cfg, attempt, record),
                          std::invalid_argument);
    }
}

TEST_CASE("DCSweep step control: failure after growth shrinks and retries", "[dc_sweep]")
{
    detail::DCSweepStepControlConfig cfg;
    cfg.start = 0.0;
    cfg.stop = 0.5;
    cfg.step = 0.5;
    cfg.minStep = 0.0625;
    cfg.maxStep = 0.5;
    cfg.growthFactor = 2.0;
    cfg.shrinkFactor = 0.5;
    cfg.maxRetries = 4;
    cfg.stopOnFailure = true;

    std::vector<detail::DCSweepStepControlEvent> events;
    std::vector<Real> attempts;

    detail::runDCSweepStepControl(
        cfg,
        [&](Real voltage, Real, int) {
            attempts.push_back(voltage);
            return attempts.size() == 2 || attempts.size() == 4 || attempts.size() == 5;
        },
        [&](const detail::DCSweepStepControlEvent& event) {
            events.push_back(event);
        });

    REQUIRE(attempts.size() == 5);
    REQUIRE(attempts[0] == Catch::Approx(0.5));
    REQUIRE(attempts[1] == Catch::Approx(0.25));
    REQUIRE(attempts[2] == Catch::Approx(0.5));
    REQUIRE(attempts[3] == Catch::Approx(0.375));
    REQUIRE(attempts[4] == Catch::Approx(0.5));

    REQUIRE(events.size() == 3);
    REQUIRE(events[0].converged);
    REQUIRE(events[0].voltage == Catch::Approx(0.25));
    REQUIRE(events[0].attemptedStep == Catch::Approx(0.25));
    REQUIRE(events[0].acceptedStep == Catch::Approx(0.25));
    REQUIRE(events[0].retryCount == 1);

    REQUIRE(events[1].converged);
    REQUIRE(events[1].voltage == Catch::Approx(0.375));
    REQUIRE(events[1].attemptedStep == Catch::Approx(0.125));
    REQUIRE(events[1].acceptedStep == Catch::Approx(0.125));
    REQUIRE(events[1].retryCount == 1);

    REQUIRE(events[2].converged);
    REQUIRE(events[2].voltage == Catch::Approx(0.5));
    REQUIRE(events[2].attemptedStep == Catch::Approx(0.125));
    REQUIRE(events[2].acceptedStep == Catch::Approx(0.125));
    REQUIRE(events[2].retryCount == 0);
}

TEST_CASE("DCSweep step control: minStep boundary records aborting failed attempt", "[dc_sweep]")
{
    detail::DCSweepStepControlConfig cfg;
    cfg.start = 0.0;
    cfg.stop = 0.5;
    cfg.step = 0.5;
    cfg.minStep = 0.2;
    cfg.maxStep = 0.5;
    cfg.growthFactor = 1.0;
    cfg.shrinkFactor = 0.5;
    cfg.maxRetries = 5;
    cfg.stopOnFailure = true;

    std::vector<detail::DCSweepStepControlEvent> events;
    std::vector<Real> attempts;

    detail::runDCSweepStepControl(
        cfg,
        [&](Real voltage, Real, int) {
            attempts.push_back(voltage);
            return false;
        },
        [&](const detail::DCSweepStepControlEvent& event) {
            events.push_back(event);
        });

    REQUIRE(attempts.size() == 2);
    REQUIRE(attempts[0] == Catch::Approx(0.5));
    REQUIRE(attempts[1] == Catch::Approx(0.25));

    REQUIRE(events.size() == 1);
    REQUIRE_FALSE(events[0].converged);
    REQUIRE(events[0].failureReason == "min_step_exhausted");
    REQUIRE(events[0].voltage == Catch::Approx(0.25));
    REQUIRE(events[0].attemptedStep == Catch::Approx(0.25));
    REQUIRE(events[0].acceptedStep == Catch::Approx(0.0));
    REQUIRE(events[0].retryCount == 1);
}

TEST_CASE("DCSweep: failed solve records retry diagnostics", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "retry_failure.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"stop", 0.5},
        {"step", 0.5},
        {"min_step", 0.0625},
        {"max_step", 0.5},
        {"shrink_factor", 0.5},
        {"growth_factor", 1.0},
        {"max_retries", 3},
        {"stop_on_failure", true},
        {"write_vtk", false}
    }, {
        {"max_iter", 0},
        {"reltol", 1.0e-30}
    });

    DCSweep sweep;
    const std::vector<DCSweepPoint> points = sweep.run(cfgPath.string());

    REQUIRE(points.size() >= 1);
    REQUIRE_FALSE(points.back().converged);
    REQUIRE(points.back().retryCount <= 3);
    REQUIRE(points.back().attemptedStep == Catch::Approx(0.0));
    REQUIRE(points.back().acceptedStep == Catch::Approx(0.0));
}

TEST_CASE("DCSweep: final stop point is reached exactly without overshoot", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "exact_stop.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"stop", 0.55},
        {"step", 0.2},
        {"min_step", 0.05},
        {"max_step", 0.2},
        {"growth_factor", 2.0},
        {"shrink_factor", 0.5},
        {"write_vtk", false}
    });

    DCSweep sweep;
    const std::vector<DCSweepPoint> points = sweep.run(cfgPath.string());

    REQUIRE(points.size() == 4);
    REQUIRE(points[1].voltage == Catch::Approx(0.2));
    REQUIRE(points[2].voltage == Catch::Approx(0.4));
    REQUIRE(points[3].voltage == Catch::Approx(0.55));
    REQUIRE(points[3].attemptedStep == Catch::Approx(0.15));
    REQUIRE(points[3].acceptedStep == Catch::Approx(0.15));
}


TEST_CASE("DCSweep: Gummel method ignores Newton-only solver fields", "[dc_sweep][gummel]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "gummel_with_newton_fields.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false}
    }, {
        {"method", "gummel"},
        {"jacobian", "ignored_by_gummel"}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);
    REQUIRE(std::filesystem::exists(csvPath));
}

TEST_CASE("DCSweep: invalid solver type message mentions both solver keys", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "invalid_solver.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false}
    }, {
        {"type", "gummle"}
    });

    DCSweep sweep;
    try {
        (void)sweep.runWithResult(cfgPath.string());
        FAIL("Expected invalid solver type to throw");
    } catch (const std::invalid_argument& ex) {
        const std::string message = ex.what();
        REQUIRE(message.find("solver.method/type") != std::string::npos);
        REQUIRE(message.find("gummel") != std::string::npos);
        REQUIRE(message.find("newton") != std::string::npos);
        REQUIRE(message.find("gummel_newton") != std::string::npos);
    }
}

TEST_CASE("DCSweep: explicit Newton solver method is reachable from config", "[dc_sweep][newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "newton_start.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false}
    }, {
        {"method", "newton"},
        {"max_iter", 10},
        {"reltol", 1.0e-8},
        {"abstol", 1.0e-18},
        {"damping_factor", 1.0},
        {"line_search", true},
        {"verbose", false}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    const std::vector<DCSweepPoint>& points = result.points;

    REQUIRE(points.size() == 1);
    REQUIRE(points.front().converged);
    REQUIRE(points.front().voltage == Catch::Approx(0.0));
    REQUIRE(points.front().attemptedStep == Catch::Approx(0.0));
    REQUIRE(std::filesystem::exists(csvPath));
}

TEST_CASE("DCSweep: hybrid Gummel-Newton method is reachable from config",
          "[dc_sweep][gummel_newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "gummel_newton_start.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 12},
        {"reltol", 1.0e-8},
        {"abstol", 1.0e-18},
        {"damping_psi", 0.35},
        {"damping_factor", 1.0},
        {"line_search", true},
        {"verbose", false}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    const DCSweepPoint& point = result.points.front();
    REQUIRE(point.converged);
    REQUIRE(point.solverMethod == "gummel_newton");
    REQUIRE(point.gummelIterations > 0);
    REQUIRE(point.newtonIterations >= 0);
    REQUIRE(point.handoffStage == "newton");
    REQUIRE(std::filesystem::exists(csvPath));
}

TEST_CASE("DCSweep: hybrid path uses Gummel iterations before Newton handoff",
          "[dc_sweep][gummel_newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "gummel_newton_forward.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.2},
        {"step", 0.2},
        {"write_vtk", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 20},
        {"reltol", 1.0e-8},
        {"abstol", 1.0e-18},
        {"damping_psi", 0.35},
        {"damping_factor", 1.0},
        {"line_search", true},
        {"warm_start", true},
        {"verbose", false}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 2);
    for (const DCSweepPoint& point : result.points) {
        REQUIRE(point.converged);
        REQUIRE(point.solverMethod == "gummel_newton");
        REQUIRE(point.gummelIterations > 0);
        REQUIRE(point.handoffStage == "newton");
    }
}

TEST_CASE("DCSweep: hybrid validates Gummel initializer before Newton handoff",
          "[dc_sweep][gummel_newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "gummel_newton_gummel_validation.csv";

    nlohmann::json cfg = baseSweepConfig(dir, meshPath, csvPath);
    cfg["sweep"]["start"] = 0.0;
    cfg["sweep"]["stop"] = 0.0;
    cfg["sweep"]["step"] = 0.25;
    cfg["sweep"]["write_vtk"] = false;
    cfg["solver"] = {
        {"method", "gummel_newton"},
        {"max_iter", 12},
        {"reltol", 1.0e-8},
        {"abstol", 1.0e-18},
        {"damping_psi", 0.35},
        {"line_search", true},
        {"verbose", false}
    };
    cfg["validation"] = {
        {"enforce_minimum_carrier_density", true},
        {"minimum_carrier_density", 1.0e40}
    };
    const auto cfgPath = dir / "gummel_newton_gummel_validation.json";
    std::ofstream(cfgPath) << cfg.dump(2);

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    const DCSweepPoint& point = result.points.front();
    REQUIRE_FALSE(point.converged);
    REQUIRE(point.solverMethod == "gummel_newton");
    REQUIRE(point.gummelIterations > 0);
    REQUIRE(point.newtonIterations == 0);
    REQUIRE(point.handoffStage == "gummel_validation_failed");
    REQUIRE(point.failureReason == "gummel_validation_failed");
}


TEST_CASE("DCSweep: NMOS and PMOS DD examples increase drain current with stronger gate drive", "[dc_sweep][mos]")
{
    const Real nmosOff = std::abs(runMosExampleDrainCurrentAtGate("nmos2d_dd", 0.0, 0.1));
    const Real nmosOn = std::abs(runMosExampleDrainCurrentAtGate("nmos2d_dd", 0.05, 0.1));
    REQUIRE(nmosOn > nmosOff);

    const Real pmosOff = std::abs(runMosExampleDrainCurrentAtGate("pmos2d_dd", 0.0, -0.1));
    const Real pmosOn = std::abs(runMosExampleDrainCurrentAtGate("pmos2d_dd", -0.05, -0.1));
    REQUIRE(pmosOn > pmosOff);
}
