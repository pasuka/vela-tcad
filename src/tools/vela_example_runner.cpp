#include "vela/simulation/DCSweep.h"
#include "vela/simulation/PoissonSimulation.h"
#include <exception>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>

namespace {

void usage(const char* argv0)
{
    std::cerr << "Usage: " << argv0 << " --config <simulation.json> [--mesh-report]\n";
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

        if (type == "dc_sweep") {
            vela::DCSweep sweep;
            const auto result = sweep.runWithResult(configFile);
            bool allConverged = !result.points.empty();
            for (const auto& point : result.points)
                allConverged = allConverged && point.converged;
            status["converged"] = allConverged;
            status["points"] = result.points.size();
            if (includeMeshReport)
                status["mesh_report"] = meshReportJson(result.mesh.lastGeometryBuildReport());
        } else if (type == "poisson") {
            vela::PoissonSimulation sim;
            const auto result = sim.runWithResult(configFile);
            status["nodes"] = result.potential.size();
            if (includeMeshReport)
                status["mesh_report"] = meshReportJson(result.mesh.lastGeometryBuildReport());
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
