#include "vela/io/MeshReader.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/NewtonSolver.h"
#include "vela/simulation/ConfigParsing.h"
#include "vela/simulation/DCSweep.h"
#include "vela/simulation/PoissonSimulation.h"
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

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

std::unordered_map<std::string, vela::Real> contactBiasesFromJson(const nlohmann::json& cfg)
{
    std::unordered_map<std::string, vela::Real> biases;
    for (const auto& contact : cfg.at("contacts")) {
        biases[contact.at("name").get<std::string>()] =
            contact.at("bias").get<vela::Real>();
    }
    return biases;
}

struct NewtonCliResult {
    vela::DeviceMesh mesh;
    vela::NewtonResult result;
};

NewtonCliResult runNewtonConfig(const std::string& configFile, const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    const vela::UnitScalingConfig scaling = vela::parseUnitScalingConfig(cfg);

    vela::JsonMeshReader reader;
    vela::DeviceMesh mesh = reader.read(
        resolvePath(cfgDir, cfg.at("mesh_file").get<std::string>()),
        scaling);

    vela::MaterialDatabase matdb;
    if (cfg.contains("materials_file"))
        matdb.loadJson(resolvePath(cfgDir, cfg.at("materials_file").get<std::string>()), scaling);

    vela::DopingModel doping =
        vela::DopingModel::fromMeshAndRegions(mesh, vela::parseDopingSpecs(cfg, scaling));
    const auto biases = contactBiasesFromJson(cfg);
    vela::NewtonConfig newton = cfg.contains("solver")
        ? vela::newtonConfigFromJson(cfg.at("solver"), scaling)
        : vela::NewtonConfig{};

    vela::NewtonResult result = vela::runNewton(mesh, matdb, doping, biases, newton);

    if (cfg.contains("output_vtk")) {
        vela::writeDDSolutionVTK(
            resolvePath(cfgDir, cfg.at("output_vtk").get<std::string>()),
            mesh,
            doping,
            result.solution);
    }

    return NewtonCliResult{std::move(mesh), std::move(result)};
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
        } else if (type == "newton") {
            const auto result = runNewtonConfig(configFile, cfg);
            status["converged"] = result.result.converged;
            status["nodes"] = result.mesh.numNodes();
            status["iterations"] = result.result.iters;
            status["initial_residual"] = result.result.initialResidualNorm;
            status["final_residual"] = result.result.finalResidualNorm;
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
