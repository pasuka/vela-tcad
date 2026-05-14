#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/simulation/DCSweep.h"
#include "vela/simulation/DCSweepStepControl.h"

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
    REQUIRE(rows.front() == std::vector<std::string>{"voltage", "electron_current", "hole_current",
                                                     "total_current", "converged", "iterations",
                                                     "attempted_step", "accepted_step", "retry_count"});
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
