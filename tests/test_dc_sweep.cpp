#include <catch2/catch_test_macros.hpp>

#include "vela/simulation/DCSweep.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <random>
#include <thread>

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

std::filesystem::path writeSweepConfig(const std::filesystem::path& dir,
                                       const std::filesystem::path& meshPath,
                                       const std::filesystem::path& csvPath)
{
    nlohmann::json cfg = {
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

    const auto cfgPath = dir / "pn_sweep.json";
    std::ofstream(cfgPath) << cfg.dump(2);
    return cfgPath;
}

} // namespace

TEST_CASE("DCSweep: PN diode sweep writes CSV and finite monotonic IV data", "[dc_sweep]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "iv.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath);

    DCSweep sweep;
    const std::vector<DCSweepPoint> points = sweep.run(cfgPath.string());

    REQUIRE(points.size() == 3);
    REQUIRE(std::filesystem::exists(csvPath));
    REQUIRE(std::filesystem::file_size(csvPath) > 0);

    for (const DCSweepPoint& point : points) {
        REQUIRE(point.converged);
        REQUIRE(std::isfinite(point.electronCurrent));
        REQUIRE(std::isfinite(point.holeCurrent));
        REQUIRE(std::isfinite(point.totalCurrent));
    }

    REQUIRE(std::abs(points.back().totalCurrent) >= std::abs(points.front().totalCurrent));
    REQUIRE(std::filesystem::exists(dir / "pn_sweep_0000_0V.vtk"));
}
