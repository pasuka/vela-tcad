#include "vela/io/MeshReader.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/NewtonSolver.h"
#include "vela/simulation/ConfigParsing.h"
#include "vela/simulation/DCSweep.h"
#include "vela/simulation/PoissonSimulation.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <sstream>
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

struct NewtonProblem {
    vela::DeviceMesh mesh;
    vela::MaterialDatabase matdb;
    vela::DopingModel doping;
    std::unordered_map<std::string, vela::Real> biases;
    vela::NewtonConfig newton;
};

vela::DopingModel readNodeDopingCsv(const std::filesystem::path& path,
                                    vela::Index nodeCount,
                                    vela::UnitScalingConfig scaling);

NewtonProblem loadNewtonProblem(const std::string& configFile, const nlohmann::json& cfg)
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

    vela::DopingModel doping = cfg.contains("node_doping_file")
        ? readNodeDopingCsv(
            resolvePath(cfgDir, cfg.at("node_doping_file").get<std::string>()),
            mesh.numNodes(),
            scaling)
        : vela::DopingModel::fromMeshAndRegions(mesh, vela::parseDopingSpecs(cfg, scaling));
    const auto biases = contactBiasesFromJson(cfg);
    vela::NewtonConfig newton = cfg.contains("solver")
        ? vela::newtonConfigFromJson(cfg.at("solver"), scaling)
        : vela::NewtonConfig{};

    return NewtonProblem{
        std::move(mesh),
        std::move(matdb),
        std::move(doping),
        std::move(biases),
        std::move(newton)};
}

NewtonCliResult runNewtonConfig(const std::string& configFile, const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);

    vela::NewtonResult result = vela::runNewton(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);

    if (cfg.contains("output_vtk")) {
        vela::writeDDSolutionVTK(
            resolvePath(cfgDir, cfg.at("output_vtk").get<std::string>()),
            problem.mesh,
            problem.doping,
            result.solution);
    }

    return NewtonCliResult{std::move(problem.mesh), std::move(result)};
}

std::vector<std::string> splitCsvLine(const std::string& line)
{
    std::vector<std::string> cells;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ','))
        cells.push_back(cell);
    return cells;
}

std::vector<vela::Real> readNodeScalarCsv(const std::filesystem::path& path,
                                          vela::Index nodeCount)
{
    std::ifstream input(path);
    if (!input.is_open())
        throw std::runtime_error("Cannot open scalar field CSV: " + path.string());
    std::string headerLine;
    if (!std::getline(input, headerLine))
        throw std::runtime_error("Empty scalar field CSV: " + path.string());
    const std::vector<std::string> header = splitCsvLine(headerLine);
    std::size_t nodeCol = header.size();
    std::size_t valueCol = header.size();
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i] == "node_id")
            nodeCol = i;
        if (header[i] == "component0")
            valueCol = i;
    }
    if (nodeCol == header.size() || valueCol == header.size())
        throw std::runtime_error(
            "Scalar field CSV must contain node_id and component0 columns: " + path.string());

    std::vector<vela::Real> values(static_cast<std::size_t>(nodeCount), 0.0);
    std::vector<bool> seen(static_cast<std::size_t>(nodeCount), false);
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> cells = splitCsvLine(line);
        if (cells.size() <= std::max(nodeCol, valueCol))
            throw std::runtime_error("Malformed scalar field CSV row: " + path.string());
        const auto node = static_cast<vela::Index>(std::stoll(cells[nodeCol]));
        if (node >= nodeCount)
            throw std::runtime_error("Scalar field CSV node_id out of range: " + path.string());
        values[static_cast<std::size_t>(node)] = std::stod(cells[valueCol]);
        seen[static_cast<std::size_t>(node)] = true;
    }
    for (vela::Index node = 0; node < nodeCount; ++node) {
        if (!seen[static_cast<std::size_t>(node)])
            throw std::runtime_error("Scalar field CSV is missing a node row: " + path.string());
    }
    return values;
}

vela::DopingModel readNodeDopingCsv(const std::filesystem::path& path,
                                    vela::Index nodeCount,
                                    vela::UnitScalingConfig scaling)
{
    std::ifstream input(path);
    if (!input.is_open())
        throw std::runtime_error("Cannot open node_doping_file: " + path.string());
    std::string headerLine;
    if (!std::getline(input, headerLine))
        throw std::runtime_error("node_doping_file is empty: " + path.string());
    const std::vector<std::string> header = splitCsvLine(headerLine);
    std::size_t nodeCol = header.size();
    std::size_t donorsCol = header.size();
    std::size_t acceptorsCol = header.size();
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i] == "node_id")
            nodeCol = i;
        if (header[i] == "donors_cm3")
            donorsCol = i;
        if (header[i] == "acceptors_cm3")
            acceptorsCol = i;
    }
    if (nodeCol == header.size() || donorsCol == header.size() || acceptorsCol == header.size()) {
        throw std::runtime_error(
            "node_doping_file must contain node_id, donors_cm3, and acceptors_cm3 columns: "
            + path.string());
    }

    vela::DopingModel model(nodeCount);
    std::vector<bool> seen(static_cast<std::size_t>(nodeCount), false);
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> cells = splitCsvLine(line);
        if (cells.size() <= std::max({nodeCol, donorsCol, acceptorsCol}))
            throw std::runtime_error("Malformed node_doping_file row: " + path.string());
        const auto node = static_cast<vela::Index>(std::stoll(cells[nodeCol]));
        if (node >= nodeCount)
            throw std::runtime_error("node_doping_file node_id out of range: " + path.string());
        if (seen[static_cast<std::size_t>(node)])
            throw std::runtime_error("node_doping_file has duplicate node_id: " + path.string());
        model.setNodeDoping(
            node,
            scaling.concentrationToSI(std::stod(cells[donorsCol])),
            scaling.concentrationToSI(std::stod(cells[acceptorsCol])));
        seen[static_cast<std::size_t>(node)] = true;
    }
    for (vela::Index node = 0; node < nodeCount; ++node) {
        if (!seen[static_cast<std::size_t>(node)])
            throw std::runtime_error("node_doping_file is missing a node row: " + path.string());
    }
    return model;
}

vela::DDSolution readExternalState(const std::filesystem::path& cfgDir,
                                   const nlohmann::json& cfg,
                                   vela::Index nodeCount)
{
    const std::filesystem::path fieldsDir =
        resolvePath(cfgDir, cfg.at("state_fields_dir").get<std::string>());
    const auto read = [&](const char* field) {
        return readNodeScalarCsv(fieldsDir / (std::string(field) + "_region0.csv"), nodeCount);
    };
    const std::vector<vela::Real> psi = read("ElectrostaticPotential");
    const std::vector<vela::Real> phin = read("eQuasiFermiPotential");
    const std::vector<vela::Real> phip = read("hQuasiFermiPotential");

    const int n = static_cast<int>(nodeCount);
    vela::DDSolution state;
    state.psi.resize(n);
    state.phin.resize(n);
    state.phip.resize(n);
    state.n = vela::VectorXd::Zero(n);
    state.p = vela::VectorXd::Zero(n);
    for (int i = 0; i < n; ++i) {
        state.psi(i) = psi[static_cast<std::size_t>(i)];
        state.phin(i) = phin[static_cast<std::size_t>(i)];
        state.phip(i) = phip[static_cast<std::size_t>(i)];
    }
    return state;
}

void writeResidualProbeCsv(const std::filesystem::path& path,
                           const vela::DeviceMesh& mesh,
                           const vela::DopingModel& doping,
                           const vela::DDSolution& state,
                           const vela::NewtonResidualEvaluation& residual)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write residual probe CSV: " + path.string());
    out << "node_id,x,y,psi,phin,phip,psi_residual,phin_residual,phip_residual,"
        << "abs_psi_residual,abs_phin_residual,abs_phip_residual,"
        << "donors_m3,acceptors_m3,net_doping_m3,ni_eff_m3\n";
    const int n = static_cast<int>(mesh.numNodes());
    for (int i = 0; i < n; ++i) {
        const auto nodeId = static_cast<vela::Index>(i);
        const vela::Node& node = mesh.getNode(nodeId);
        const vela::Real rPsi = residual.raw(i);
        const vela::Real rPhin = residual.raw(n + i);
        const vela::Real rPhip = residual.raw(2 * n + i);
        out << nodeId << ','
            << node.x << ','
            << node.y << ','
            << state.psi(i) << ','
            << state.phin(i) << ','
            << state.phip(i) << ','
            << rPsi << ','
            << rPhin << ','
            << rPhip << ','
            << std::abs(rPsi) << ','
            << std::abs(rPhin) << ','
            << std::abs(rPhip) << ','
            << doping.donors(nodeId) << ','
            << doping.acceptors(nodeId) << ','
            << doping.netDoping(nodeId) << ','
            << residual.intrinsicDensity[static_cast<std::size_t>(nodeId)]
            << '\n';
    }
}

nlohmann::json runNewtonResidualProbe(const std::string& configFile, const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state = readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const vela::NewtonResidualEvaluation residual = solver.evaluateResidual(state);
    writeResidualProbeCsv(
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>()),
        problem.mesh,
        problem.doping,
        state,
        residual);

    return {
        {"nodes", problem.mesh.numNodes()},
        {"scaled_state", residual.scaledState},
        {"potential_scale", residual.potentialScale},
        {"block_residuals", {
            {"psi", residual.blockNorms.psi},
            {"phin", residual.blockNorms.phin},
            {"phip", residual.blockNorms.phip},
            {"combined", residual.blockNorms.combined},
        }},
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
        } else if (type == "newton") {
            const auto result = runNewtonConfig(configFile, cfg);
            status["converged"] = result.result.converged;
            status["nodes"] = result.mesh.numNodes();
            status["iterations"] = result.result.iters;
            status["initial_residual"] = result.result.initialResidualNorm;
            status["final_residual"] = result.result.finalResidualNorm;
            if (includeMeshReport)
                status["mesh_report"] = meshReportJson(result.mesh.lastGeometryBuildReport());
        } else if (type == "newton_residual_probe") {
            status.update(runNewtonResidualProbe(configFile, cfg));
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
