#include "vela/simulation/DCSweep.h"
#include "vela/simulation/PoissonSimulation.h"
#include "vela/io/MeshReader.h"

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>

namespace {

void usage(const char* argv0)
{
    std::cerr << "Usage: " << argv0 << " --config <simulation.json> [--mesh-report]\n";
}

std::filesystem::path configDirectory(const std::string& configFile)
{
    const std::filesystem::path path(configFile);
    const std::filesystem::path parent = path.parent_path();
    return parent.empty() ? std::filesystem::current_path() : parent;
}

std::string resolvePath(const std::filesystem::path& baseDir, const std::string& path)
{
    std::filesystem::path resolved(path);
    if (resolved.is_relative())
        resolved = baseDir / resolved;
    return resolved.string();
}

nlohmann::json meshReportJson(const vela::GeometryBuildReport& report)
{
    return {
        {"total_cells", report.totalCells},
        {"degenerate_cells", report.degenerateCells},
        {"negative_cotangent_count", report.negativeCotangentCount},
        {"fallback_count", report.fallbackCount},
        {"min_angle_degrees", report.minAngleDegrees},
        {"max_angle_degrees", report.maxAngleDegrees},
        {"min_edge_length", report.minEdgeLength},
    };
}

} // namespace

int main(int argc, char** argv)
{
    std::string configFile;
    bool includeMeshReport = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--config" && i + 1 < argc) {
            configFile = argv[++i];
        } else if (arg == "--mesh-report") {
            includeMeshReport = true;
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return 0;
        } else {
            usage(argv[0]);
            return 2;
        }
    }

    if (configFile.empty()) {
        usage(argv[0]);
        return 2;
    }

    try {
        std::ifstream ifs(configFile);
        if (!ifs.is_open()) {
            std::cerr << "Cannot open config file: " << configFile << '\n';
            return 1;
        }

        nlohmann::json cfg;
        ifs >> cfg;
        const std::string type = cfg.value("simulation_type", cfg.contains("sweep") ? "dc_sweep" : "poisson");

        nlohmann::json status;
        status["config"] = configFile;
        status["simulation_type"] = type;
        status["converged"] = true;

        if (includeMeshReport) {
            vela::JsonMeshReader reader;
            const std::filesystem::path cfgDir = configDirectory(configFile);
            vela::DeviceMesh mesh = reader.read(resolvePath(cfgDir, cfg.at("mesh_file").get<std::string>()));
            status["mesh_report"] = meshReportJson(mesh.lastGeometryBuildReport());
        }

        if (type == "dc_sweep") {
            vela::DCSweep sweep;
            const auto points = sweep.run(configFile);
            bool allConverged = !points.empty();
            for (const auto& point : points)
                allConverged = allConverged && point.converged;
            status["converged"] = allConverged;
            status["points"] = points.size();
        } else if (type == "poisson") {
            vela::PoissonSimulation sim;
            const auto psi = sim.run(configFile);
            status["nodes"] = psi.size();
        } else {
            std::cerr << "Unknown simulation_type: " << type << '\n';
            return 2;
        }

        std::cout << status.dump() << '\n';
        return status.value("converged", false) ? 0 : 1;
    } catch (const std::exception& ex) {
        std::cerr << "vela_example_runner failed: " << ex.what() << '\n';
        return 1;
    }
}
