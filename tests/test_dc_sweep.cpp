#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/io/DDSolutionCsv.h"
#include "vela/simulation/DCSweep.h"
#include "vela/simulation/DCSweepPredictor.h"
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

std::filesystem::path writePNMeshWithInterior(const std::filesystem::path& dir)
{
    nlohmann::json mesh = {
        {"nodes", {
            {{"id", 0}, {"x", 0.0e-6}, {"y", 0.0e-6}},
            {{"id", 1}, {"x", 1.0e-6}, {"y", 0.0e-6}},
            {{"id", 2}, {"x", 1.0e-6}, {"y", 1.0e-6}},
            {{"id", 3}, {"x", 0.0e-6}, {"y", 1.0e-6}},
            {{"id", 4}, {"x", 0.5e-6}, {"y", 0.5e-6}}
        }},
        {"triangles", {
            {{"id", 0}, {"region_id", 0}, {"node_ids", {0, 1, 4}}},
            {{"id", 1}, {"region_id", 0}, {"node_ids", {1, 2, 4}}},
            {{"id", 2}, {"region_id", 1}, {"node_ids", {2, 3, 4}}},
            {{"id", 3}, {"region_id", 1}, {"node_ids", {3, 0, 4}}}
        }},
        {"regions", {
            {{"id", 0}, {"name", "n_region"}, {"material", "Si"}, {"cell_ids", {0, 1}}},
            {{"id", 1}, {"name", "p_region"}, {"material", "Si"}, {"cell_ids", {2, 3}}}
        }},
        {"contacts", {
            {{"id", 0}, {"name", "anode"}, {"region_id", 1}, {"node_ids", {0, 3}}},
            {{"id", 1}, {"name", "cathode"}, {"region_id", 0}, {"node_ids", {1, 2}}}
        }}
    };

    const auto meshPath = dir / "pn_mesh_with_interior.json";
    std::ofstream(meshPath) << mesh.dump(2);
    return meshPath;
}

std::filesystem::path writeRefinementTransitionMesh(const std::filesystem::path& dir)
{
    nlohmann::json mesh = {
        {"nodes", {
            {{"id", 0}, {"x", 0.0e-6}, {"y", 0.0e-6}},
            {{"id", 1}, {"x", 2.0e-6}, {"y", 0.0e-6}},
            {{"id", 2}, {"x", 0.0e-6}, {"y", 1.0e-6}},
            {{"id", 3}, {"x", 2.0e-6}, {"y", 1.0e-6}},
            {{"id", 4}, {"x", 1.0e-6}, {"y", 0.0e-6}}
        }},
        {"triangles", {
            {{"id", 0}, {"region_id", 0}, {"node_ids", {0, 4, 2}}},
            {{"id", 1}, {"region_id", 1}, {"node_ids", {4, 1, 3}}},
            {{"id", 2}, {"region_id", 1}, {"node_ids", {4, 3, 2}}}
        }},
        {"regions", {
            {{"id", 0}, {"name", "n_region"}, {"material", "Si"}, {"cell_ids", {0}}},
            {{"id", 1}, {"name", "p_region"}, {"material", "Si"}, {"cell_ids", {1, 2}}}
        }},
        {"contacts", {
            {{"id", 0}, {"name", "anode"}, {"region_id", 0}, {"node_ids", {0, 2}}},
            {{"id", 1}, {"name", "cathode"}, {"region_id", 1}, {"node_ids", {1, 3}}}
        }}
    };

    const auto meshPath = dir / "refinement_transition_mesh.json";
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

std::string readTextFile(const std::filesystem::path& path)
{
    std::ifstream input(path);
    REQUIRE(input.is_open());
    std::ostringstream ss;
    ss << input.rdbuf();
    return ss.str();
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

bool hasExtraField(const DCSweepPoint& point, const std::string& name, Real* outValue = nullptr)
{
    const auto it = std::find_if(point.extraFields.begin(), point.extraFields.end(),
        [&](const auto& entry) { return entry.first == name; });
    if (it == point.extraFields.end())
        return false;
    if (outValue != nullptr)
        *outValue = it->second;
    return true;
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
                                                     "current_contact", "current_electron", "current_electron_drift",
                                                     "current_electron_diffusion", "current_hole", "current_hole_drift",
                                                     "current_hole_diffusion", "current_total", "converged", "iterations",
                                                     "solver_method", "gummel_iterations", "newton_iterations",
                                                     "handoff_stage", "step_diagnostics", "validation_diagnostics",
                                                     "failure_reason", "newton_failure_class",
                                                     "newton_failure_diagnostics_json"});
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
    const std::size_t currentElectronDrift = csvColumnIndex(header, "current_electron_drift");
    const std::size_t currentElectronDiffusion = csvColumnIndex(header, "current_electron_diffusion");
    const std::size_t currentHole = csvColumnIndex(header, "current_hole");
    const std::size_t currentHoleDrift = csvColumnIndex(header, "current_hole_drift");
    const std::size_t currentHoleDiffusion = csvColumnIndex(header, "current_hole_diffusion");
    const std::size_t maxField = csvColumnIndex(header, "max_electric_field_V_per_m");
    const std::size_t currentTotalUm = csvColumnIndex(header, "current_total_A_per_um");
    const std::size_t currentElectronUm = csvColumnIndex(header, "current_electron_A_per_um");
    const std::size_t currentElectronDriftUm = csvColumnIndex(header, "current_electron_drift_A_per_um");
    const std::size_t currentElectronDiffusionUm = csvColumnIndex(header, "current_electron_diffusion_A_per_um");
    const std::size_t currentHoleUm = csvColumnIndex(header, "current_hole_A_per_um");
    const std::size_t currentHoleDriftUm = csvColumnIndex(header, "current_hole_drift_A_per_um");
    const std::size_t currentHoleDiffusionUm = csvColumnIndex(header, "current_hole_diffusion_A_per_um");
    const std::size_t maxFieldCm = csvColumnIndex(header, "max_electric_field_V_per_cm");

    for (std::size_t r = 1; r < rows.size(); ++r) {
        const auto& row = rows.at(r);
        REQUIRE(csvReal(row, currentTotalUm) ==
                Catch::Approx(csvReal(row, currentTotal) / 1.0e6).epsilon(1.0e-12));
        REQUIRE(csvReal(row, currentElectronUm) ==
                Catch::Approx(csvReal(row, currentElectron) / 1.0e6).epsilon(1.0e-12));
        REQUIRE(csvReal(row, currentElectronDriftUm) ==
            Catch::Approx(csvReal(row, currentElectronDrift) / 1.0e6).epsilon(1.0e-12));
        REQUIRE(csvReal(row, currentElectronDiffusionUm) ==
            Catch::Approx(csvReal(row, currentElectronDiffusion) / 1.0e6).epsilon(1.0e-12));
        REQUIRE(csvReal(row, currentHoleUm) ==
                Catch::Approx(csvReal(row, currentHole) / 1.0e6).epsilon(1.0e-12));
        REQUIRE(csvReal(row, currentHoleDriftUm) ==
            Catch::Approx(csvReal(row, currentHoleDrift) / 1.0e6).epsilon(1.0e-12));
        REQUIRE(csvReal(row, currentHoleDiffusionUm) ==
            Catch::Approx(csvReal(row, currentHoleDiffusion) / 1.0e6).epsilon(1.0e-12));
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

TEST_CASE("DCSweep: mesh_geometry node_volume_policy selects mixed Voronoi volumes",
          "[dc_sweep][mesh_geometry]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "mixed_volume_iv.csv";
    nlohmann::json cfg = baseSweepConfig(dir, meshPath, csvPath);
    cfg["mesh_geometry"] = {{"node_volume_policy", "mixed_voronoi"}};
    cfg["sweep"]["start"] = 0.0;
    cfg["sweep"]["stop"] = 0.0;
    cfg["sweep"]["step"] = 0.1;
    cfg["sweep"]["write_vtk"] = false;

    const auto cfgPath = dir / "mixed_volume_sweep.json";
    std::ofstream(cfgPath) << cfg.dump(2);

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);
    REQUIRE(result.mesh.getNode(0).volume == Catch::Approx(0.25e-12));
    REQUIRE(result.mesh.getNode(1).volume == Catch::Approx(0.25e-12));
    REQUIRE(result.mesh.getNode(2).volume == Catch::Approx(0.25e-12));
    REQUIRE(result.mesh.getNode(3).volume == Catch::Approx(0.25e-12));
}

TEST_CASE("DCSweep: high-doping node-level PN diode converges with hybrid handoff",
          "[dc_sweep][gummel_newton][doping]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "node_doping_hybrid_iv.csv";
    const auto cfgPath = writeUnitScalingSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.1},
        {"write_vtk", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 60},
        {"reltol", 1.0e-7},
        {"abstol", 1.0e-18},
        {"damping_psi", 0.2},
        {"line_search", true},
        {"warm_start", true},
        {"verbose", false},
        {"handoff", {{"fallback", "none"}}}
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

    REQUIRE(result.points.size() >= 1);
    REQUIRE(result.points.front().converged);
    REQUIRE(result.points.front().solverMethod == "gummel_newton");
    REQUIRE(result.points.front().handoffStage == "newton");
    REQUIRE(std::isfinite(result.points.front().totalCurrent));
}

TEST_CASE("DCSweep: recombination diagnostics are opt-in for hybrid handoff",
          "[dc_sweep][gummel_newton][diagnostics]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "node_doping_hybrid_diag.csv";
    const auto cfgPath = writeUnitScalingSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.1},
        {"write_vtk", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 60},
        {"reltol", 1.0e-7},
        {"abstol", 1.0e-18},
        {"damping_psi", 0.2},
        {"line_search", true},
        {"warm_start", true},
        {"verbose", false},
        {"handoff", {{"fallback", "none"}}}
    });

    std::ifstream input(cfgPath);
    nlohmann::json cfg;
    input >> cfg;
    cfg["node_doping_file"] = "doping.csv";
    cfg["doping"] = {
        {{"region", "n_region"}, {"donors", 0.0}, {"acceptors", 0.0}},
        {{"region", "p_region"}, {"donors", 0.0}, {"acceptors", 0.0}}
    };
    cfg["solver"]["diagnostics"] = false;
    std::ofstream(cfgPath) << cfg.dump(2);

    writeNodeDopingCsv(dir / "doping.csv", {
        {0, 0.0, 1.0e17},
        {1, 1.0e17, 0.0},
        {2, 1.0e17, 0.0},
        {3, 0.0, 1.0e17},
    });

    DCSweep sweep;
    const DCSweepResult noDiag = sweep.runWithResult(cfgPath.string());
    REQUIRE(noDiag.points.size() == 1);
    REQUIRE(noDiag.points.front().converged);
    REQUIRE_FALSE(hasExtraField(noDiag.points.front(), "recombination_max_abs_rate_m3_per_s"));

    const auto rowsNoDiag = readCsvRows(csvPath);
    REQUIRE(std::find(rowsNoDiag.front().begin(), rowsNoDiag.front().end(),
                      "recombination_max_abs_rate_m3_per_s") == rowsNoDiag.front().end());

    cfg["solver"]["diagnostics"] = true;
    std::ofstream(cfgPath) << cfg.dump(2);

    const DCSweepResult withDiag = sweep.runWithResult(cfgPath.string());
    REQUIRE(withDiag.points.size() == 1);
    REQUIRE(withDiag.points.front().converged);

    Real maxAbsRate = 0.0;
    Real meanAbsRate = 0.0;
    Real maxNpOverNi2 = 0.0;
    REQUIRE(hasExtraField(withDiag.points.front(), "recombination_max_abs_rate_m3_per_s", &maxAbsRate));
    REQUIRE(hasExtraField(withDiag.points.front(), "recombination_mean_abs_rate_m3_per_s", &meanAbsRate));
    REQUIRE(hasExtraField(withDiag.points.front(), "carrier_product_max_np_over_ni2", &maxNpOverNi2));
    REQUIRE(std::isfinite(maxAbsRate));
    REQUIRE(std::isfinite(meanAbsRate));
    REQUIRE(std::isfinite(maxNpOverNi2));

    const auto rowsWithDiag = readCsvRows(csvPath);
    const std::size_t maxAbsRateColumn =
        csvColumnIndex(rowsWithDiag.front(), "recombination_max_abs_rate_m3_per_s");
    const std::size_t meanAbsRateColumn =
        csvColumnIndex(rowsWithDiag.front(), "recombination_mean_abs_rate_m3_per_s");
    const std::size_t maxNpOverNi2Column =
        csvColumnIndex(rowsWithDiag.front(), "carrier_product_max_np_over_ni2");

    REQUIRE(std::isfinite(csvReal(rowsWithDiag.at(1), maxAbsRateColumn)));
    REQUIRE(std::isfinite(csvReal(rowsWithDiag.at(1), meanAbsRateColumn)));
    REQUIRE(std::isfinite(csvReal(rowsWithDiag.at(1), maxNpOverNi2Column)));
}

TEST_CASE("DCSweep: transport diagnostics append mobility and current-driver columns",
          "[dc_sweep][diagnostics][transport]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "iv_transport.csv";
    const auto cfgPath = writeUnitScalingSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.05},
        {"write_vtk", false},
        {"diagnostics", {
            {"transport", {
                {"enabled", true},
            }},
        }},
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);
    const auto rows = readCsvRows(csvPath);
    REQUIRE(rows.size() == 2);
    const auto& header = rows.front();
    csvColumnIndex(header, "mean_electron_mobility_m2_V_s");
    csvColumnIndex(header, "mean_hole_mobility_m2_V_s");
    csvColumnIndex(header, "min_electron_mobility_m2_V_s");
    csvColumnIndex(header, "min_hole_mobility_m2_V_s");
    csvColumnIndex(header, "max_electric_field_V_per_cm");
    csvColumnIndex(header, "mean_electron_qf_gradient_V_per_cm");
    csvColumnIndex(header, "mean_hole_qf_gradient_V_per_cm");
    const std::size_t meanElectronDriveCol =
        csvColumnIndex(header, "mean_electron_high_field_drive_V_per_cm");
    const std::size_t meanHoleDriveCol =
        csvColumnIndex(header, "mean_hole_high_field_drive_V_per_cm");
    const std::size_t minElectronLimiterCol =
        csvColumnIndex(header, "min_electron_mobility_limiter");
    const std::size_t minHoleLimiterCol =
        csvColumnIndex(header, "min_hole_mobility_limiter");
    const std::size_t meanElectronLimiterCol =
        csvColumnIndex(header, "mean_electron_mobility_limiter");
    const std::size_t meanHoleLimiterCol =
        csvColumnIndex(header, "mean_hole_mobility_limiter");
    const auto& data = rows.at(1);
    REQUIRE(std::isfinite(csvReal(data, meanElectronDriveCol)));
    REQUIRE(std::isfinite(csvReal(data, meanHoleDriveCol)));
    REQUIRE(std::isfinite(csvReal(data, minElectronLimiterCol)));
    REQUIRE(std::isfinite(csvReal(data, minHoleLimiterCol)));
    REQUIRE(std::isfinite(csvReal(data, meanElectronLimiterCol)));
    REQUIRE(std::isfinite(csvReal(data, meanHoleLimiterCol)));
}

TEST_CASE("DCSweep: VTK transport diagnostics include mobility decomposition fields",
          "[dc_sweep][diagnostics][vtk]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "iv_transport_vtk.csv";
    const auto vtkPrefix = dir / "iv_transport_vtk";
    const auto cfgPath = writeUnitScalingSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.05},
        {"write_vtk", true},
        {"vtk_prefix", vtkPrefix.string()},
        {"solver", {
            {"mobility", {
                {"model", "masetti_field"},
                {"high_field_driving_force", "quasi_fermi_gradient"},
            }},
        }},
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);
    REQUIRE_FALSE(result.points.front().outputVtk.empty());

    std::ifstream ifs(result.points.front().outputVtk);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());

    REQUIRE(content.find("ElectronLowFieldMobility") != std::string::npos);
    REQUIRE(content.find("HoleLowFieldMobility") != std::string::npos);
    REQUIRE(content.find("ElectronHighFieldDrive") != std::string::npos);
    REQUIRE(content.find("HoleHighFieldDrive") != std::string::npos);
    REQUIRE(content.find("ElectronMobilityLimiter") != std::string::npos);
    REQUIRE(content.find("HoleMobilityLimiter") != std::string::npos);
}

TEST_CASE("DCSweep: contact-edge diagnostics are opt-in and write per-edge rows",
          "[dc_sweep][diagnostics][contact_edge]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "iv_unit_scaling.csv";
    const auto edgeDiagPath = dir / "iv_contact_edges.csv";
    const auto cfgPath = writeUnitScalingSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.1},
        {"write_vtk", false},
        {"diagnostics", {
            {"contact_edge", {
                {"enabled", true},
                {"csv_file", edgeDiagPath.string()}
            }}
        }}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);

    REQUIRE(std::filesystem::exists(edgeDiagPath));
    const auto rows = readCsvRows(edgeDiagPath);
    REQUIRE(rows.size() >= 1);

    const auto& header = rows.front();
    const std::size_t pointIndexCol = csvColumnIndex(header, "point_index");
    const std::size_t electronBranchCol = csvColumnIndex(header, "electron_branch");
    const std::size_t holeBranchCol = csvColumnIndex(header, "hole_branch");
    const std::size_t edgeCurrentCol = csvColumnIndex(header, "current_total");
    const std::size_t edgeCurrentUmCol = csvColumnIndex(header, "current_total_A_per_um");

    if (rows.size() > 1) {
        for (std::size_t i = 1; i < rows.size(); ++i) {
            const auto& row = rows.at(i);
            REQUIRE(csvReal(row, pointIndexCol) == Catch::Approx(0.0));
            const std::string electronBranch = row.at(electronBranchCol);
            const std::string holeBranch = row.at(holeBranchCol);
            REQUIRE((electronBranch == "density" || electronBranch == "quasi_fermi"));
            REQUIRE((holeBranch == "density" || holeBranch == "quasi_fermi"));
            REQUIRE(std::isfinite(csvReal(row, edgeCurrentCol)));
            REQUIRE(csvReal(row, edgeCurrentUmCol) ==
                    Catch::Approx(csvReal(row, edgeCurrentCol) / 1.0e6).epsilon(1.0e-12));
        }
    }
}

TEST_CASE("DCSweep: SG avalanche edge diagnostics write assembled source rows",
          "[dc_sweep][diagnostics][sg_avalanche_edges]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "bv_sg_source.csv";
    const auto edgeSourcePath = dir / "bv_sg_source_edges.csv";
    const auto cfgPath = writeUnitScalingSweepConfig(dir, meshPath, csvPath, {
        {"mode", "bv_reverse"},
        {"start", 0.0},
        {"stop", 0.0},
        {"step", -0.05},
        {"write_vtk", false},
        {"diagnostics", {
            {"sg_avalanche_edges", {
                {"enabled", true},
                {"csv_file", edgeSourcePath.string()}
            }}
        }}
    }, {
        {"method", "gummel_newton"},
        {"handoff", {
            {"gummel_max_iter", 0},
            {"newton_max_iter", 80},
            {"require_gummel_convergence", false}
        }},
        {"impact_ionization", {
            {"model", "selberherr"},
            {"driving_force", "electric_field"},
            {"generation", "current_density"},
            {"current_approximation", "density_gradient"},
            {"electron_A_m_inv", 1.0},
            {"electron_B_V_m", 1.0e-30},
            {"hole_A_m_inv", 1.0},
            {"hole_B_V_m", 1.0e-30}
        }}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);

    REQUIRE(std::filesystem::exists(edgeSourcePath));
    const auto rows = readCsvRows(edgeSourcePath);
    REQUIRE(rows.size() > 1);

    const auto& header = rows.front();
    const std::size_t pointIndexCol = csvColumnIndex(header, "point_index");
    const std::size_t edgeSourceCol = csvColumnIndex(header, "edge_source_integral");
    const std::size_t node0SourceCol = csvColumnIndex(header, "node0_source_integral");
    const std::size_t node1SourceCol = csvColumnIndex(header, "node1_source_integral");
    (void)csvColumnIndex(header, "edge_area_proxy_m2");
    (void)csvColumnIndex(header, "electric_field_V_per_m");
    (void)csvColumnIndex(header, "electron_impact_field_V_per_m");
    (void)csvColumnIndex(header, "hole_impact_field_V_per_m");
    (void)csvColumnIndex(header, "electron_alpha_m_inv");
    (void)csvColumnIndex(header, "hole_alpha_m_inv");
    (void)csvColumnIndex(header, "electron_flux_proxy");
    (void)csvColumnIndex(header, "hole_flux_proxy");
    (void)csvColumnIndex(header, "edge_class");

    for (std::size_t i = 1; i < rows.size(); ++i) {
        const auto& row = rows.at(i);
        REQUIRE(csvReal(row, pointIndexCol) == Catch::Approx(0.0));
        const Real edgeSource = csvReal(row, edgeSourceCol);
        REQUIRE(edgeSource >= 0.0);
        REQUIRE(csvReal(row, node0SourceCol) == Catch::Approx(0.5 * edgeSource));
        REQUIRE(csvReal(row, node1SourceCol) == Catch::Approx(0.5 * edgeSource));
    }
}

TEST_CASE("DCSweep: continuity-balance diagnostics write contact-adjacent residual rows",
          "[dc_sweep][diagnostics][continuity_balance]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "iv_unit_scaling.csv";
    const auto balancePath = dir / "iv_continuity_balance.csv";
    const auto cfgPath = writeUnitScalingSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.1},
        {"write_vtk", false},
        {"diagnostics", {
            {"continuity_balance", {
                {"enabled", true},
                {"contacts", {"anode", "cathode"}},
                {"csv_file", balancePath.string()}
            }}
        }}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);

    REQUIRE(std::filesystem::exists(balancePath));
    const auto rows = readCsvRows(balancePath);
    REQUIRE(rows.size() > 1);

    const auto& header = rows.front();
    const std::size_t pointIndexCol = csvColumnIndex(header, "point_index");
    const std::size_t contactCol = csvColumnIndex(header, "contact");
    const std::size_t carrierCol = csvColumnIndex(header, "carrier");
    const std::size_t residualCol = csvColumnIndex(header, "continuity_residual");
    const std::size_t contactFluxCol = csvColumnIndex(header, "contact_edge_flux");
    const std::size_t neighborFluxCol = csvColumnIndex(header, "neighbor_edge_flux");
    const std::size_t recombinationCol = csvColumnIndex(header, "recombination_term");
    (void)csvColumnIndex(header, "contact_node");
    (void)csvColumnIndex(header, "interior_node");
    (void)csvColumnIndex(header, "interior_volume_m2");
    (void)csvColumnIndex(header, "qf_contact_V");
    (void)csvColumnIndex(header, "qf_interior_V");

    bool sawElectron = false;
    bool sawHole = false;
    for (std::size_t i = 1; i < rows.size(); ++i) {
        const auto& row = rows.at(i);
        REQUIRE(csvReal(row, pointIndexCol) == Catch::Approx(0.0));
        REQUIRE((row.at(contactCol) == "anode" || row.at(contactCol) == "cathode"));
        sawElectron = sawElectron || row.at(carrierCol) == "electron";
        sawHole = sawHole || row.at(carrierCol) == "hole";
        REQUIRE(std::isfinite(csvReal(row, residualCol)));
        REQUIRE(std::isfinite(csvReal(row, contactFluxCol)));
        REQUIRE(std::isfinite(csvReal(row, neighborFluxCol)));
        REQUIRE(std::isfinite(csvReal(row, recombinationCol)));
    }
    REQUIRE(sawElectron);
    REQUIRE(sawHole);
}

TEST_CASE("DCSweep: Newton history diagnostic writes accepted iteration block residuals",
          "[dc_sweep][diagnostics][newton_history]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshWithInterior(dir);
    const auto csvPath = dir / "newton_history.csv";
    const auto historyPath = dir / "newton_history_iterations.csv";
    const auto initialStatePath = dir / "newton_history_initial_state.csv";
    {
        std::ofstream state(initialStatePath);
        state << "node_id,psi,phin,phip,electrons_m3,holes_m3\n";
        for (int node = 0; node < 5; ++node)
            state << node << ",0,0,0,1e10,1e10\n";
    }
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.05},
        {"stop", 0.05},
        {"step", 0.25},
        {"write_vtk", false},
        {"initial_state_file", initialStatePath.string()},
        {"diagnostics", {
            {"newton_history", {
                {"enabled", true},
                {"csv_file", historyPath.string()}
            }}
        }}
    }, {
        {"method", "newton"},
        {"line_search", true},
        {"warm_start", true},
        {"verbose", false},
        {"max_iter", 80}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);
    REQUIRE(std::filesystem::exists(historyPath));

    const auto rows = readCsvRows(historyPath);
    REQUIRE(rows.size() > 1);
    const auto& header = rows.front();
    const std::size_t pointIndexCol = csvColumnIndex(header, "point_index");
    const std::size_t biasCol = csvColumnIndex(header, "bias_V");
    const std::size_t iterCol = csvColumnIndex(header, "iteration");
    const std::size_t residualCol = csvColumnIndex(header, "residual_norm");
    const std::size_t relativeCol = csvColumnIndex(header, "relative_residual_norm");
    const std::size_t rawStepCol = csvColumnIndex(header, "raw_step_norm");
    const std::size_t appliedStepCol = csvColumnIndex(header, "applied_step_norm");
    const std::size_t dampingCol = csvColumnIndex(header, "damping_factor");
    const std::size_t attemptsCol = csvColumnIndex(header, "line_search_attempts");
    const std::size_t psiBlockCol = csvColumnIndex(header, "block_psi");
    const std::size_t phinBlockCol = csvColumnIndex(header, "block_phin");
    const std::size_t phipBlockCol = csvColumnIndex(header, "block_phip");
    const std::size_t combinedBlockCol = csvColumnIndex(header, "block_combined");

    REQUIRE(rows.at(1).at(pointIndexCol) == "0");
    REQUIRE(std::stod(rows.at(1).at(biasCol)) == Catch::Approx(0.05));
    REQUIRE(std::stoi(rows.at(1).at(iterCol)) >= 1);
    REQUIRE(std::stod(rows.at(1).at(residualCol)) > 0.0);
    REQUIRE(std::stod(rows.at(1).at(relativeCol)) >= 0.0);
    REQUIRE(std::stod(rows.at(1).at(rawStepCol)) >=
            std::stod(rows.at(1).at(appliedStepCol)));
    REQUIRE(std::stod(rows.at(1).at(dampingCol)) > 0.0);
    REQUIRE(std::stoi(rows.at(1).at(attemptsCol)) >= 1);
    REQUIRE(std::stod(rows.at(1).at(psiBlockCol)) >= 0.0);
    REQUIRE(std::stod(rows.at(1).at(phinBlockCol)) >= 0.0);
    REQUIRE(std::stod(rows.at(1).at(phipBlockCol)) >= 0.0);
    REQUIRE(std::stod(rows.at(1).at(combinedBlockCol)) >=
            std::stod(rows.at(1).at(psiBlockCol)));
}

TEST_CASE("DCSweep: continuation predictor config is validated",
          "[dc_sweep][continuation][predictor]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "continuation_predictor.csv";

    auto writeConfigWithContinuation = [&](const nlohmann::json& continuation) {
        const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
            {"start", 0.0},
            {"stop", 0.0},
            {"step", 0.25},
            {"write_vtk", false},
            {"continuation", continuation}
        });
        return cfgPath;
    };

    DCSweep sweep;

    SECTION("valid predictor modes and branch acceptance parse")
    {
        for (const std::string mode : {"none", "constant", "linear", "secant"}) {
            INFO(mode);
            const auto cfgPath = writeConfigWithContinuation({
                {"predictor", {
                    {"mode", mode},
                    {"fields", {"psi", "phin", "phip"}},
                    {"max_extrapolation_ratio", 2.0}
                }},
                {"branch_acceptance", {
                    {"terminal_current_consistency", true},
                    {"min_terminal_current_ratio", 1.0e-6}
                }}
            });
            const DCSweepResult result = sweep.runWithResult(cfgPath.string());
            REQUIRE(result.points.size() == 1);
        }
    }

    SECTION("invalid predictor mode is rejected")
    {
        const auto cfgPath = writeConfigWithContinuation({
            {"predictor", {
                {"mode", "quadratic"}
            }}
        });
        REQUIRE_THROWS_WITH(
            sweep.runWithResult(cfgPath.string()),
            Catch::Matchers::ContainsSubstring(
                "DCSweep: sweep.continuation.predictor.mode must be"));
    }

    SECTION("invalid predictor field is rejected")
    {
        const auto cfgPath = writeConfigWithContinuation({
            {"predictor", {
                {"mode", "linear"},
                {"fields", {"psi", "electrons"}}
            }}
        });
        REQUIRE_THROWS_WITH(
            sweep.runWithResult(cfgPath.string()),
            Catch::Matchers::ContainsSubstring(
                "DCSweep: sweep.continuation.predictor.fields entries must be"));
    }

    SECTION("invalid predictor ratio is rejected")
    {
        const auto cfgPath = writeConfigWithContinuation({
            {"predictor", {
                {"mode", "linear"},
                {"max_extrapolation_ratio", 0.5}
            }}
        });
        REQUIRE_THROWS_WITH(
            sweep.runWithResult(cfgPath.string()),
            Catch::Matchers::ContainsSubstring(
                "DCSweep: sweep.continuation.predictor.max_extrapolation_ratio"));
    }

    SECTION("invalid terminal current threshold is rejected")
    {
        const auto cfgPath = writeConfigWithContinuation({
            {"branch_acceptance", {
                {"min_terminal_current_ratio", -1.0}
            }}
        });
        REQUIRE_THROWS_WITH(
            sweep.runWithResult(cfgPath.string()),
            Catch::Matchers::ContainsSubstring(
                "DCSweep: sweep.continuation.branch_acceptance.min_terminal_current_ratio"));
    }

    SECTION("invalid psi-phin jump threshold is rejected")
    {
        const auto cfgPath = writeConfigWithContinuation({
            {"branch_acceptance", {
                {"psi_phin_jump", true},
                {"max_psi_phin_jump_V", -1.0}
            }}
        });
        REQUIRE_THROWS_WITH(
            sweep.runWithResult(cfgPath.string()),
            Catch::Matchers::ContainsSubstring(
                "DCSweep: sweep.continuation.branch_acceptance.max_psi_phin_jump_V"));
    }

    SECTION("invalid p95 electron density jump threshold is rejected")
    {
        const auto cfgPath = writeConfigWithContinuation({
            {"branch_acceptance", {
                {"carrier_density_jump", true},
                {"max_electron_density_jump_dex", 100.0},
                {"max_electron_density_jump_p95_abs_dex", -0.1}
            }}
        });
        REQUIRE_THROWS_WITH(
            sweep.runWithResult(cfgPath.string()),
            Catch::Matchers::ContainsSubstring(
                "DCSweep: sweep.continuation.branch_acceptance."
                "max_electron_density_jump_p95_abs_dex"));
    }
}

TEST_CASE("DCSweep: terminal balance diagnostics reuse one solution for two contacts",
          "[dc_sweep][diagnostics][terminal_balance][contact_edge]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "iv_unit_scaling.csv";
    const auto balancePath = dir / "iv_terminal_balance.csv";
    const auto edgeDiagPath = dir / "iv_contact_edges.csv";
    const auto cfgPath = writeUnitScalingSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.1},
        {"write_vtk", false},
        {"diagnostics", {
            {"terminal_balance", {
                {"enabled", true},
                {"contacts", {"anode", "cathode"}},
                {"csv_file", balancePath.string()}
            }},
            {"contact_edge", {
                {"enabled", true},
                {"contacts", {"anode", "cathode"}},
                {"csv_file", edgeDiagPath.string()}
            }}
        }}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);

    REQUIRE(std::filesystem::exists(balancePath));
    const auto balanceRows = readCsvRows(balancePath);
    REQUIRE(balanceRows.size() == 3);
    const auto& balanceHeader = balanceRows.front();
    const std::size_t balanceContactCol = csvColumnIndex(balanceHeader, "contact");
    const std::size_t electronCol = csvColumnIndex(balanceHeader, "current_electron");
    const std::size_t holeCol = csvColumnIndex(balanceHeader, "current_hole");
    const std::size_t minusCol = csvColumnIndex(balanceHeader, "electron_minus_hole");
    const std::size_t plusCol = csvColumnIndex(balanceHeader, "electron_plus_hole");
    const std::size_t minusUmCol = csvColumnIndex(balanceHeader, "electron_minus_hole_A_per_um");
    const std::size_t plusUmCol = csvColumnIndex(balanceHeader, "electron_plus_hole_A_per_um");
    Real minusUmPairSum = 0.0;

    for (std::size_t i = 1; i < balanceRows.size(); ++i) {
        const auto& row = balanceRows.at(i);
        const Real electron = csvReal(row, electronCol);
        const Real hole = csvReal(row, holeCol);
        REQUIRE(csvReal(row, minusCol) == Catch::Approx(electron - hole).epsilon(1.0e-12));
        REQUIRE(csvReal(row, plusCol) == Catch::Approx(electron + hole).epsilon(1.0e-12));
        REQUIRE(csvReal(row, minusUmCol) == Catch::Approx(csvReal(row, minusCol) / 1.0e6).epsilon(1.0e-12));
        REQUIRE(csvReal(row, plusUmCol) == Catch::Approx(csvReal(row, plusCol) / 1.0e6).epsilon(1.0e-12));
        minusUmPairSum += csvReal(row, minusUmCol);
    }
    REQUIRE(std::abs(minusUmPairSum) <= 1.0e-24);

    REQUIRE(std::filesystem::exists(edgeDiagPath));
    const auto edgeRows = readCsvRows(edgeDiagPath);
    REQUIRE(edgeRows.size() > 2);
    const auto& edgeHeader = edgeRows.front();
    const std::size_t edgeContactCol = csvColumnIndex(edgeHeader, "current_contact");
    const std::size_t edgeTotalCol = csvColumnIndex(edgeHeader, "current_total");
    (void)csvColumnIndex(edgeHeader, "psi0");
    (void)csvColumnIndex(edgeHeader, "phin0");
    (void)csvColumnIndex(edgeHeader, "phip0");
    (void)csvColumnIndex(edgeHeader, "n0");
    (void)csvColumnIndex(edgeHeader, "p0");
    (void)csvColumnIndex(edgeHeader, "ni0");
    (void)csvColumnIndex(edgeHeader, "mun");
    (void)csvColumnIndex(edgeHeader, "electron_continuity_flux");
    (void)csvColumnIndex(edgeHeader, "hole_continuity_flux");

    for (const std::string contact : {"anode", "cathode"}) {
        Real edgeSum = 0.0;
        int edgeCount = 0;
        for (std::size_t i = 1; i < edgeRows.size(); ++i) {
            if (edgeRows.at(i).at(edgeContactCol) == contact) {
                edgeSum += csvReal(edgeRows.at(i), edgeTotalCol);
                ++edgeCount;
            }
        }
        REQUIRE(edgeCount > 0);

        bool foundTerminal = false;
        for (std::size_t i = 1; i < balanceRows.size(); ++i) {
            if (balanceRows.at(i).at(balanceContactCol) == contact) {
                REQUIRE(edgeSum == Catch::Approx(csvReal(balanceRows.at(i), minusCol)).epsilon(1.0e-12));
                foundTerminal = true;
            }
        }
        REQUIRE(foundTerminal);
    }
}

TEST_CASE("DCSweep: contact current QF floor reporting uses initial edge drops only when enabled",
          "[dc_sweep][diagnostics][contact_current_qf_floor]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshWithInterior(dir);
    const auto defaultCsvPath = dir / "default.csv";
    const auto floorCsvPath = dir / "floor.csv";
    const auto defaultEdgesPath = dir / "default_edges.csv";
    const auto floorEdgesPath = dir / "floor_edges.csv";
    const auto initialStatePath = dir / "initial_state.csv";
    {
        std::ofstream state(initialStatePath);
        state << "node_id,psi,phin,phip,electrons_m3,holes_m3\n";
        state << "0,0,0,-1e-6,1e10,1e23\n";
        state << "1,0,0,0,1e23,1e10\n";
        state << "2,0,0,0,1e23,1e10\n";
        state << "3,0,0,-1e-6,1e10,1e23\n";
        state << "4,0,0,0,1e12,1e12\n";
    }

    const nlohmann::json commonSweep = {
        {"start", 1.0e-6},
        {"stop", 1.0e-6},
        {"step", 1.0e-6},
        {"write_vtk", false},
        {"initial_state_file", initialStatePath.string()},
    };
    nlohmann::json defaultSweep = commonSweep;
    defaultSweep["diagnostics"] = {
        {"contact_edge", {
            {"enabled", true},
            {"contacts", {"anode"}},
            {"csv_file", defaultEdgesPath.string()}
        }}
    };
    nlohmann::json floorSweep = commonSweep;
    floorSweep["csv_file"] = floorCsvPath.string();
    floorSweep["diagnostics"] = {
        {"contact_edge", {
            {"enabled", true},
            {"contacts", {"anode"}},
            {"csv_file", floorEdgesPath.string()}
        }},
        {"contact_current_qf_floor", {
            {"enabled", true},
            {"contacts", {"anode"}}
        }}
    };
    const nlohmann::json solverOverrides = {
        {"method", "newton"},
        {"warm_start", true},
        {"line_search", true},
        {"reltol", 1.0e-4},
        {"max_iter", 80}
    };

    DCSweep sweep;
    const auto defaultCfg = writeSweepConfig(
        dir, meshPath, defaultCsvPath, defaultSweep, solverOverrides);
    const DCSweepResult defaultResult = sweep.runWithResult(defaultCfg.string());
    const auto floorCfg = writeSweepConfig(
        dir, meshPath, floorCsvPath, floorSweep, solverOverrides);
    const DCSweepResult floorResult = sweep.runWithResult(floorCfg.string());
    REQUIRE(defaultResult.points.size() == 1);
    REQUIRE(floorResult.points.size() == 1);
    REQUIRE(defaultResult.points.front().converged);
    REQUIRE(floorResult.points.front().converged);

    const auto defaultEdgeRows = readCsvRows(defaultEdgesPath);
    const auto floorEdgeRows = readCsvRows(floorEdgesPath);
    REQUIRE(defaultEdgeRows.size() == floorEdgeRows.size());
    const auto& defaultEdgeHeader = defaultEdgeRows.front();
    const auto& floorEdgeHeader = floorEdgeRows.front();
    const std::size_t defaultOverrideCol =
        csvColumnIndex(defaultEdgeHeader, "hole_qf_drop_override_applied");
    const std::size_t floorOverrideCol =
        csvColumnIndex(floorEdgeHeader, "hole_qf_drop_override_applied");
    const std::size_t phip0Col = csvColumnIndex(floorEdgeHeader, "phip0");
    const std::size_t phip1Col = csvColumnIndex(floorEdgeHeader, "phip1");
    const std::size_t holeCurrentCol = csvColumnIndex(floorEdgeHeader, "current_hole");

    bool sawOverride = false;
    Real floorEdgeHoleCurrent = 0.0;
    for (std::size_t i = 1; i < floorEdgeRows.size(); ++i) {
        REQUIRE(defaultEdgeRows.at(i).at(defaultOverrideCol) == "0");
        floorEdgeHoleCurrent += csvReal(floorEdgeRows.at(i), holeCurrentCol);
        if (floorEdgeRows.at(i).at(floorOverrideCol) == "1") {
            sawOverride = true;
            REQUIRE(std::abs(csvReal(floorEdgeRows.at(i), phip1Col) -
                             csvReal(floorEdgeRows.at(i), phip0Col)) > 0.0);
        }
    }
    REQUIRE(sawOverride);
    REQUIRE(floorResult.points.front().holeCurrent ==
            Catch::Approx(floorEdgeHoleCurrent).margin(1.0e-18));
}

TEST_CASE("DCSweep: contact current QF floor reporting ignores continuation states",
          "[dc_sweep][diagnostics][contact_current_qf_floor]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshWithInterior(dir);
    const auto csvPath = dir / "sweep.csv";
    const auto edgePath = dir / "contact_edges.csv";

    const nlohmann::json sweepOverrides = {
        {"start", 0.0},
        {"stop", -0.05},
        {"step", -0.05},
        {"write_vtk", false},
        {"diagnostics", {
            {"contact_edge", {
                {"enabled", true},
                {"contacts", {"anode"}},
                {"csv_file", edgePath.string()}
            }},
            {"contact_current_qf_floor", {
                {"enabled", true},
                {"contacts", {"anode"}}
            }}
        }}
    };
    const nlohmann::json solverOverrides = {
        {"method", "newton"},
        {"warm_start", true},
        {"line_search", true},
        {"reltol", 1.0e-4},
        {"max_iter", 80}
    };

    const auto cfgPath = writeSweepConfig(
        dir, meshPath, csvPath, sweepOverrides, solverOverrides);
    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 2);
    REQUIRE(result.points.at(0).converged);
    REQUIRE(result.points.at(1).converged);

    const auto rows = readCsvRows(edgePath);
    const std::size_t overrideCol =
        csvColumnIndex(rows.front(), "hole_qf_drop_override_applied");
    for (std::size_t i = 1; i < rows.size(); ++i)
        REQUIRE(rows.at(i).at(overrideCol) == "0");
}

TEST_CASE("DCSweep: contact current reporting policy preserves initial endpoint QF drops",
          "[dc_sweep][contact_current_reporting]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMeshWithInterior(dir);
    const auto csvPath = dir / "reporting_policy.csv";
    const auto edgePath = dir / "reporting_policy_edges.csv";
    const auto initialStatePath = dir / "initial_state.csv";
    {
        std::ofstream state(initialStatePath);
        state << "node_id,psi,phin,phip,electrons_m3,holes_m3\n";
        state << "0,0,0,-1e-6,1e10,1e23\n";
        state << "1,0,0,0,1e23,1e10\n";
        state << "2,0,0,0,1e23,1e10\n";
        state << "3,0,0,-1e-6,1e10,1e23\n";
        state << "4,0,0,0,1e12,1e12\n";
    }

    const nlohmann::json sweepOverrides = {
        {"start", 1.0e-6},
        {"stop", 1.0e-6},
        {"step", 1.0e-6},
        {"write_vtk", false},
        {"initial_state_file", initialStatePath.string()},
        {"diagnostics", {
            {"contact_edge", {
                {"enabled", true},
                {"contacts", {"anode"}},
                {"csv_file", edgePath.string()}
            }}
        }},
        {"contact_current_reporting", {
            {"endpoint_qf_floor", {
                {"enabled", true},
                {"contacts", {"anode"}}
            }}
        }},
        {"continuation", {
            {"predictor", {
                {"mode", "constant"},
                {"fields", {"psi", "phin", "phip"}}
            }}
        }}
    };
    const nlohmann::json solverOverrides = {
        {"method", "newton"},
        {"warm_start", true},
        {"line_search", true},
        {"reltol", 1.0e-4},
        {"max_iter", 80}
    };

    const auto cfgPath = writeSweepConfig(
        dir, meshPath, csvPath, sweepOverrides, solverOverrides);
    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);

    const auto rows = readCsvRows(edgePath);
    const std::size_t overrideCol =
        csvColumnIndex(rows.front(), "hole_qf_drop_override_applied");
    const std::size_t phip0Col = csvColumnIndex(rows.front(), "phip0");
    const std::size_t phip1Col = csvColumnIndex(rows.front(), "phip1");
    const std::size_t holeCurrentCol =
        csvColumnIndex(rows.front(), "current_hole");
    bool sawOverride = false;
    Real edgeHoleCurrent = 0.0;
    for (std::size_t i = 1; i < rows.size(); ++i) {
        edgeHoleCurrent += csvReal(rows.at(i), holeCurrentCol);
        if (rows.at(i).at(overrideCol) == "1") {
            sawOverride = true;
            REQUIRE(csvReal(rows.at(i), phip1Col) - csvReal(rows.at(i), phip0Col) ==
                    Catch::Approx(1.0e-6).margin(1.0e-15));
        }
    }
    REQUIRE(sawOverride);
    REQUIRE(result.points.front().holeCurrent ==
            Catch::Approx(edgeHoleCurrent).margin(1.0e-18));
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
                                                         "current_contact", "current_electron", "current_electron_drift",
                                                         "current_electron_diffusion", "current_hole", "current_hole_drift",
                                                         "current_hole_diffusion", "current_total", "converged", "iterations",
                                                         "solver_method", "gummel_iterations", "newton_iterations",
                                                         "handoff_stage",
                                                         "step_diagnostics", "validation_diagnostics",
                                                         "failure_reason", "newton_failure_class",
                                                         "newton_failure_diagnostics_json", "charge_C_per_m",
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
                                                         "current_contact", "current_electron", "current_electron_drift",
                                                         "current_electron_diffusion", "current_hole", "current_hole_drift",
                                                         "current_hole_diffusion", "current_total", "converged", "iterations",
                                                         "solver_method", "gummel_iterations", "newton_iterations",
                                                         "handoff_stage",
                                                         "step_diagnostics", "validation_diagnostics",
                                                         "failure_reason", "newton_failure_class",
                                                         "newton_failure_diagnostics_json", "charge_C_per_m",
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
                                                         "current_contact", "current_electron", "current_electron_drift",
                                                         "current_electron_diffusion", "current_hole", "current_hole_drift",
                                                         "current_hole_diffusion", "current_total", "converged", "iterations",
                                                         "solver_method", "gummel_iterations", "newton_iterations",
                                                         "handoff_stage",
                                                         "step_diagnostics", "validation_diagnostics",
                                                         "failure_reason", "newton_failure_class",
                                                         "newton_failure_diagnostics_json", "max_electric_field_V_per_m",
                                                         "current_jump_ratio", "breakdown_detected",
                                                         "breakdown_voltage", "criterion", "last_stable_bias",
                                                         "failed_bias", "breakdown_failure_reason"});
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
                                                     "current_contact", "current_electron", "current_electron_drift",
                                                     "current_electron_diffusion", "current_hole", "current_hole_drift",
                                                     "current_hole_diffusion", "current_total", "converged", "iterations",
                                                     "solver_method", "gummel_iterations", "newton_iterations",
                                                     "handoff_stage",
                                                     "step_diagnostics", "validation_diagnostics",
                                                     "failure_reason", "newton_failure_class",
                                                     "newton_failure_diagnostics_json", "max_electric_field_V_per_m",
                                                     "current_jump_ratio", "breakdown_detected",
                                                     "breakdown_voltage", "criterion", "last_stable_bias",
                                                     "failed_bias", "breakdown_failure_reason"});
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
                                                     "current_contact", "current_electron", "current_electron_drift",
                                                     "current_electron_diffusion", "current_hole", "current_hole_drift",
                                                     "current_hole_diffusion", "current_total", "converged", "iterations",
                                                     "solver_method", "gummel_iterations", "newton_iterations",
                                                     "handoff_stage",
                                                     "step_diagnostics", "validation_diagnostics",
                                                     "failure_reason", "newton_failure_class",
                                                     "newton_failure_diagnostics_json", "max_electric_field_V_per_m",
                                                     "current_jump_ratio", "breakdown_detected",
                                                     "breakdown_voltage", "criterion", "last_stable_bias",
                                                     "failed_bias", "breakdown_failure_reason"});
    const std::size_t criterionColumn = csvColumnIndex(rows.front(), "criterion");
    const std::size_t failureReasonColumn = csvColumnIndex(rows.front(), "breakdown_failure_reason");
    REQUIRE(rows.at(1).at(criterionColumn).empty());
    REQUIRE(rows.at(1).at(failureReasonColumn) == "non_convergence");
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

TEST_CASE("DCSweep: explicit bias_points solve only requested biases", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "bias_points.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.5},
        {"step", 0.25},
        {"bias_points", {0.0, 0.125, 0.4}},
        {"write_vtk", false}
    });

    DCSweep sweep;
    const std::vector<DCSweepPoint> points = sweep.run(cfgPath.string());

    REQUIRE(points.size() == 3);
    REQUIRE(points[0].voltage == Catch::Approx(0.0));
    REQUIRE(points[1].voltage == Catch::Approx(0.125));
    REQUIRE(points[2].voltage == Catch::Approx(0.4));
    REQUIRE(points[1].attemptedStep == Catch::Approx(0.125));
    REQUIRE(points[1].acceptedStep == Catch::Approx(0.125));
    REQUIRE(points[2].attemptedStep == Catch::Approx(0.025));
    REQUIRE(points[2].acceptedStep == Catch::Approx(0.025));
}

TEST_CASE("DCSweep: write_state_file stores latest converged restart state", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "state_writer.csv";
    const auto statePath = dir / "latest_state.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.25},
        {"step", 0.25},
        {"write_vtk", false},
        {"write_state_file", statePath.string()}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 2);
    REQUIRE(std::filesystem::exists(statePath));
    const auto rows = readCsvRows(statePath);
    REQUIRE(rows.size() == result.mesh.numNodes() + 1);
    REQUIRE(rows.front() == std::vector<std::string>{
        "node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"});
    for (std::size_t row = 1; row < rows.size(); ++row) {
        REQUIRE(rows[row].size() == 6);
        REQUIRE(std::stoul(rows[row][0]) == row - 1);
        for (std::size_t column = 1; column < rows[row].size(); ++column)
            REQUIRE(std::isfinite(std::stod(rows[row][column])));
    }
}

TEST_CASE("DDSolution CSV shared IO roundtrips restart state", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto path = dir / "state.csv";

    DDSolution solution;
    solution.psi = VectorXd::LinSpaced(3, -0.1, 0.1);
    solution.phin = VectorXd::LinSpaced(3, 0.2, 0.4);
    solution.phip = VectorXd::LinSpaced(3, -0.4, -0.2);
    solution.n = VectorXd::Constant(3, 1.0e16);
    solution.p = VectorXd::Constant(3, 2.0e16);

    writeDDSolutionStateCsv(path, solution);
    const DDSolution loaded = readDDSolutionStateCsv(path, 3);

    REQUIRE((loaded.psi - solution.psi).norm() == Catch::Approx(0.0));
    REQUIRE((loaded.phin - solution.phin).norm() == Catch::Approx(0.0));
    REQUIRE((loaded.phip - solution.phip).norm() == Catch::Approx(0.0));
    REQUIRE((loaded.n - solution.n).norm() == Catch::Approx(0.0));
    REQUIRE((loaded.p - solution.p).norm() == Catch::Approx(0.0));
}

TEST_CASE("DCSweep: write_state_every_point_prefix stores accepted states", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "point_states.csv";
    const auto prefix = dir / "states" / "bv_state";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", -0.1},
        {"step", -0.05},
        {"write_vtk", false},
        {"write_state_every_point_prefix", prefix.string()}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 3);
    REQUIRE(std::filesystem::exists(dir / "states" / "bv_state_bias_0p000000.csv"));
    REQUIRE(std::filesystem::exists(dir / "states" / "bv_state_bias_m0p050000.csv"));
    REQUIRE(std::filesystem::exists(dir / "states" / "bv_state_bias_m0p100000.csv"));
}

TEST_CASE("DCSweep: initial_state_file validates restart node coverage", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "bad_restart.csv";
    const auto statePath = dir / "bad_state.csv";
    {
        std::ofstream state(statePath);
        state << "node_id,psi,phin,phip,electrons_m3,holes_m3\n";
        state << "0,0,0,0,1e10,1e10\n";
    }
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false},
        {"initial_state_file", statePath.string()}
    });

    DCSweep sweep;
    REQUIRE_THROWS_WITH(
        sweep.run(cfgPath.string()),
        Catch::Matchers::ContainsSubstring("DCSweep: initial_state_file missing row for node id 1"));
}

TEST_CASE("DCSweep predictor: extrapolates selected coupled variables",
          "[dc_sweep][continuation][predictor]")
{
    DDSolution previous;
    previous.psi = VectorXd::LinSpaced(3, 1.0, 3.0);
    previous.phin = VectorXd::LinSpaced(3, 10.0, 12.0);
    previous.phip = VectorXd::LinSpaced(3, 20.0, 22.0);
    previous.n = VectorXd::Constant(3, 100.0);
    previous.p = VectorXd::Constant(3, 200.0);
    previous.iters = 4;
    previous.converged = true;

    DDSolution current;
    current.psi = VectorXd::LinSpaced(3, 2.0, 4.0);
    current.phin = VectorXd::LinSpaced(3, 12.0, 14.0);
    current.phip = VectorXd::LinSpaced(3, 23.0, 25.0);
    current.n = VectorXd::Constant(3, 300.0);
    current.p = VectorXd::Constant(3, 400.0);
    current.iters = 5;
    current.converged = true;

    SECTION("none and constant return current state")
    {
        for (const std::string mode : {"none", "constant"}) {
            SweepPredictorConfig config;
            config.mode = mode;
            const DDSolution predicted = detail::predictDCSweepInitialState(
                config, &previous, current, -12.65, -12.70, -12.75);

            REQUIRE(predicted.psi.isApprox(current.psi));
            REQUIRE(predicted.phin.isApprox(current.phin));
            REQUIRE(predicted.phip.isApprox(current.phip));
            REQUIRE(predicted.n.isApprox(current.n));
            REQUIRE(predicted.p.isApprox(current.p));
        }
    }

    SECTION("linear extrapolates selected fields and leaves carriers from current")
    {
        SweepPredictorConfig config;
        config.mode = "linear";
        config.fields = {"psi", "phin"};
        config.maxExtrapolationRatio = 2.0;

        const DDSolution predicted = detail::predictDCSweepInitialState(
            config, &previous, current, -12.65, -12.70, -12.75);

        REQUIRE(predicted.psi.isApprox(current.psi + (current.psi - previous.psi)));
        REQUIRE(predicted.phin.isApprox(current.phin + (current.phin - previous.phin)));
        REQUIRE(predicted.phip.isApprox(current.phip));
        REQUIRE(predicted.n.isApprox(current.n));
        REQUIRE(predicted.p.isApprox(current.p));
    }

    SECTION("secant currently uses the same bounded extrapolation")
    {
        SweepPredictorConfig config;
        config.mode = "secant";
        config.fields = {"phip"};
        config.maxExtrapolationRatio = 2.0;

        const DDSolution predicted = detail::predictDCSweepInitialState(
            config, &previous, current, -12.65, -12.70, -12.75);

        REQUIRE(predicted.psi.isApprox(current.psi));
        REQUIRE(predicted.phin.isApprox(current.phin));
        REQUIRE(predicted.phip.isApprox(current.phip + (current.phip - previous.phip)));
    }

    SECTION("linear extrapolation ratio is clamped")
    {
        SweepPredictorConfig config;
        config.mode = "linear";
        config.fields = {"psi"};
        config.maxExtrapolationRatio = 1.5;

        const DDSolution predicted = detail::predictDCSweepInitialState(
            config, &previous, current, -12.65, -12.70, -12.85);

        REQUIRE(predicted.psi.isApprox(current.psi + 1.5 * (current.psi - previous.psi)));
    }

    SECTION("linear predictor is disabled for shrunken retry attempts")
    {
        SweepPredictorConfig config;
        config.mode = "linear";
        config.fields = {"psi", "phin", "phip"};
        config.maxExtrapolationRatio = 2.0;

        const int retryCount = 1;
        const DDSolution predicted = detail::predictDCSweepInitialState(
            config, &previous, current, -12.65, -12.70, -12.75, retryCount);

        REQUIRE(predicted.psi.isApprox(current.psi));
        REQUIRE(predicted.phin.isApprox(current.phin));
        REQUIRE(predicted.phip.isApprox(current.phip));
        REQUIRE(predicted.n.isApprox(current.n));
        REQUIRE(predicted.p.isApprox(current.p));
    }
}

TEST_CASE("DCSweep: continuation predictor writes branch diagnostics",
          "[dc_sweep][continuation][predictor]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "continuation_predictor_diagnostics.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.5},
        {"step", 0.25},
        {"write_vtk", false},
        {"continuation", {
            {"predictor", {
                {"mode", "linear"},
                {"fields", {"psi", "phin", "phip"}},
                {"max_extrapolation_ratio", 2.0}
            }},
            {"branch_acceptance", {
                {"terminal_current_consistency", true},
                {"min_terminal_current_ratio", 0.0}
            }}
        }}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 3);
    REQUIRE(result.points.at(2).converged);

    const auto rows = readCsvRows(csvPath);
    REQUIRE(rows.size() == 4);
    const auto& header = rows.front();
    const std::size_t predictorModeCol = csvColumnIndex(header, "predictor_mode");
    const std::size_t predictedStateCol = csvColumnIndex(header, "predicted_initial_state");
    const std::size_t branchStatusCol = csvColumnIndex(header, "branch_acceptance_status");
    const std::size_t branchReasonCol = csvColumnIndex(header, "branch_acceptance_reason");
    const std::size_t ratioCol = csvColumnIndex(header, "terminal_current_consistency_ratio");

    REQUIRE(rows.at(1).at(predictorModeCol) == "linear");
    REQUIRE(rows.at(1).at(predictedStateCol) == "0");
    REQUIRE(rows.at(2).at(predictedStateCol) == "1");
    REQUIRE(rows.at(3).at(predictedStateCol) == "1");
    REQUIRE(rows.at(3).at(branchStatusCol) == "accepted");
    REQUIRE(rows.at(3).at(branchReasonCol).empty());
    REQUIRE(std::isfinite(std::stod(rows.at(3).at(ratioCol))));
}

TEST_CASE("DCSweep branch acceptance: measures psi-phin exponent jumps",
          "[dc_sweep][continuation][branch_acceptance]")
{
    DDSolution previous;
    previous.psi = VectorXd::Zero(3);
    previous.phin = VectorXd::Zero(3);
    previous.phip = VectorXd::Zero(3);
    previous.n = VectorXd::Constant(3, 1.0e10);
    previous.p = VectorXd::Constant(3, 1.0e10);

    DDSolution current = previous;
    current.psi(0) = 0.01;
    current.phin(0) = 0.01;
    current.psi(1) = 0.18;
    current.phin(1) = 0.02;
    current.psi(2) = -0.04;
    current.phin(2) = -0.01;

    REQUIRE(detail::maxPsiPhinJump(previous, current) == Catch::Approx(0.16));
}

TEST_CASE("DCSweep branch acceptance: measures electron density jump statistics",
          "[dc_sweep][continuation][branch_acceptance]")
{
    DDSolution previous;
    previous.psi = VectorXd::Zero(4);
    previous.phin = VectorXd::Zero(4);
    previous.phip = VectorXd::Zero(4);
    previous.n = VectorXd::Constant(4, 1.0e10);
    previous.p = VectorXd::Constant(4, 1.0e10);

    DDSolution current = previous;
    current.n(0) = 1.0e10;
    current.n(1) = 1.0e11;
    current.n(2) = 1.0e12;
    current.n(3) = 1.0e9;

    const auto stats = detail::electronDensityJumpStats(previous, current);

    REQUIRE(stats.medianDex == Catch::Approx(0.5));
    REQUIRE(stats.maxAbsDex == Catch::Approx(2.0));
    REQUIRE(stats.maxSignedDex == Catch::Approx(2.0));
    REQUIRE(stats.maxNode == 2);
}

TEST_CASE("DCSweep branch acceptance: rejects invalid electron density jump threshold",
          "[dc_sweep][continuation][branch_acceptance]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "bad_carrier_branch_guard.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.25},
        {"step", 0.25},
        {"write_vtk", false},
        {"continuation", {
            {"branch_acceptance", {
                {"carrier_density_jump", true},
                {"max_electron_density_jump_dex", -0.1}
            }}
        }}
    });

    DCSweep sweep;
    REQUIRE_THROWS_WITH(
        sweep.run(cfgPath.string()),
        Catch::Matchers::ContainsSubstring(
            "DCSweep: sweep.continuation.branch_acceptance.max_electron_density_jump_dex"));
}

TEST_CASE("DCSweep: psi-phin branch guard writes jump diagnostics",
          "[dc_sweep][continuation][branch_acceptance]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "psi_phin_branch_guard.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.25},
        {"step", 0.25},
        {"min_step", 0.125},
        {"max_step", 0.25},
        {"max_retries", 0},
        {"stop_on_failure", true},
        {"write_vtk", false},
        {"continuation", {
            {"branch_acceptance", {
                {"psi_phin_jump", true},
                {"max_psi_phin_jump_V", 1.0}
            }}
        }}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 2);
    REQUIRE(result.points.front().converged);
    REQUIRE(result.points.back().converged);
    REQUIRE(result.points.back().branchAcceptanceStatus == "accepted");
    REQUIRE(result.points.back().branchAcceptanceReason.empty());
    REQUIRE(result.points.back().psiPhinMaxJump_V >= 0.0);

    const auto rows = readCsvRows(csvPath);
    REQUIRE(rows.size() == 3);
    const auto& header = rows.front();
    const std::size_t statusCol = csvColumnIndex(header, "branch_acceptance_status");
    const std::size_t reasonCol = csvColumnIndex(header, "branch_acceptance_reason");
    const std::size_t jumpCol = csvColumnIndex(header, "psi_phin_max_jump_V");

    REQUIRE(rows.at(2).at(statusCol) == "accepted");
    REQUIRE(rows.at(2).at(reasonCol).empty());
    REQUIRE(std::stod(rows.at(2).at(jumpCol)) >= 0.0);
}

TEST_CASE("DCSweep: carrier density branch guard writes jump diagnostics",
          "[dc_sweep][continuation][branch_acceptance]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "carrier_density_branch_guard.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.25},
        {"step", 0.25},
        {"min_step", 0.125},
        {"max_step", 0.25},
        {"max_retries", 0},
        {"stop_on_failure", true},
        {"write_vtk", false},
        {"continuation", {
            {"branch_acceptance", {
                {"carrier_density_jump", true},
                {"max_electron_density_jump_dex", 100.0}
            }}
        }}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 2);
    REQUIRE(result.points.front().converged);
    REQUIRE(result.points.back().converged);
    REQUIRE(result.points.back().branchAcceptanceStatus == "accepted");
    REQUIRE(result.points.back().branchAcceptanceReason.empty());
    REQUIRE(std::isfinite(result.points.back().electronDensityJumpMedianDex));
    REQUIRE(result.points.back().electronDensityJumpP95AbsDex >= 0.0);
    REQUIRE(result.points.back().electronDensityJumpMaxAbsDex >= 0.0);
    REQUIRE(result.points.back().electronDensityJumpMaxNode >= 0);

    const auto rows = readCsvRows(csvPath);
    REQUIRE(rows.size() == 3);
    const auto& header = rows.front();
    const std::size_t statusCol = csvColumnIndex(header, "branch_acceptance_status");
    const std::size_t reasonCol = csvColumnIndex(header, "branch_acceptance_reason");
    const std::size_t medianCol =
        csvColumnIndex(header, "electron_density_jump_median_dex");
    const std::size_t p95Col =
        csvColumnIndex(header, "electron_density_jump_p95_abs_dex");
    const std::size_t maxCol =
        csvColumnIndex(header, "electron_density_jump_max_abs_dex");
    const std::size_t nodeCol =
        csvColumnIndex(header, "electron_density_jump_max_node");

    REQUIRE(rows.at(2).at(statusCol) == "accepted");
    REQUIRE(rows.at(2).at(reasonCol).empty());
    REQUIRE(std::isfinite(std::stod(rows.at(2).at(medianCol))));
    REQUIRE(std::stod(rows.at(2).at(p95Col)) >= 0.0);
    REQUIRE(std::stod(rows.at(2).at(maxCol)) >= 0.0);
    REQUIRE(std::stoi(rows.at(2).at(nodeCol)) >= 0);
}

TEST_CASE("DCSweep branch acceptance: classifies p95 electron density jumps",
          "[dc_sweep][continuation][branch_acceptance]")
{
    SweepBranchAcceptanceConfig cfg;
    cfg.carrierDensityJump = true;
    cfg.maxElectronDensityJumpDex = 100.0;
    cfg.maxElectronDensityJumpP95AbsDex = 0.15;

    detail::ElectronDensityJumpStats stats;
    stats.medianDex = 0.01;
    stats.p95AbsDex = 0.20;
    stats.maxAbsDex = 0.30;
    stats.maxNode = 7;

    REQUIRE(detail::electronDensityJumpAcceptanceFailure(cfg, stats) ==
            "electron_density_p95_jump_exceeded");
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

TEST_CASE("DCSweep: CSV records hybrid solver handoff provenance",
          "[dc_sweep][gummel_newton][csv]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "handoff_columns.csv";
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
        {"line_search", true},
        {"warm_start", true},
        {"verbose", false},
        {"handoff", {{"fallback", "none"}}}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().handoffStage == "newton");

    const std::string csv = readTextFile(csvPath);
    REQUIRE(csv.find("solver_method,gummel_iterations,newton_iterations,handoff_stage") !=
            std::string::npos);
    REQUIRE(csv.find("gummel_newton") != std::string::npos);
    REQUIRE(csv.find(",newton,") != std::string::npos);
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

TEST_CASE("DCSweep: hybrid fallback can accept converged Gummel when Newton fails",
          "[dc_sweep][gummel_newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "gummel_newton_fallback.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 20},
        {"handoff", {
            {"fallback", "gummel_on_newton_failure"},
            {"newton_max_iter", 0}
        }},
        {"verbose", false}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    const DCSweepPoint& point = result.points.front();
    REQUIRE(point.converged);
    REQUIRE(point.handoffStage == "gummel_fallback");
    REQUIRE(point.gummelIterations > 0);
    REQUIRE(point.newtonIterations == 0);
}

TEST_CASE("DCSweep: hybrid handoff has separate Gummel and Newton iteration budgets",
          "[dc_sweep][gummel_newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "hybrid_budget.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 0},
        {"reltol", 1.0e-8},
        {"abstol", 1.0e-18},
        {"damping_psi", 0.35},
        {"line_search", true},
        {"warm_start", true},
        {"verbose", false},
        {"handoff", {
            {"fallback", "none"},
            {"gummel_max_iter", 20},
            {"newton_max_iter", 12}
        }}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);
    REQUIRE(result.points.front().gummelIterations > 0);
    REQUIRE(result.points.front().newtonIterations <= 12);
    REQUIRE(result.points.front().handoffStage == "newton");
}

TEST_CASE("DCSweep: hybrid strict policy rejects Newton failure",
          "[dc_sweep][gummel_newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "gummel_newton_strict.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false},
        {"stop_on_failure", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 20},
        {"handoff", {
            {"fallback", "none"},
            {"newton_max_iter", 0}
        }},
        {"verbose", false}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    const DCSweepPoint& point = result.points.front();
    REQUIRE_FALSE(point.converged);
    REQUIRE(point.handoffStage == "newton_failed");
    REQUIRE(point.failureReason == "newton_non_convergence");
    REQUIRE(point.newtonFailureClass.empty());
    REQUIRE(point.failureDiagnosticsJson.empty());

    const auto rows = readCsvRows(csvPath);
    const std::size_t failureReasonColumn = csvColumnIndex(rows.front(), "failure_reason");
    const std::size_t newtonFailureColumn = csvColumnIndex(rows.front(), "newton_failure_class");
    REQUIRE(rows.at(1).at(failureReasonColumn) == "newton_non_convergence");
    REQUIRE(rows.at(1).at(newtonFailureColumn).empty());
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
