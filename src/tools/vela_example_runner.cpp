#include "vela/io/CsvUtils.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/io/DDSolutionCsv.h"
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
#include <iomanip>
#include <iostream>
#include <memory>
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
    mesh.buildBoxGeometry(vela::parseBoxGeometryOptions(cfg));

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

std::vector<vela::Real> readNodeScalarCsv(const std::filesystem::path& path,
                                          vela::Index nodeCount)
{
    std::ifstream input(path);
    if (!input.is_open())
        throw std::runtime_error("Cannot open scalar field CSV: " + path.string());
    std::string headerLine;
    if (!std::getline(input, headerLine))
        throw std::runtime_error("Empty scalar field CSV: " + path.string());
    const std::vector<std::string> header = vela::splitCsvLine(headerLine);
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
        const std::vector<std::string> cells = vela::splitCsvLine(line);
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
    const std::vector<std::string> header = vela::splitCsvLine(headerLine);
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
        const std::vector<std::string> cells = vela::splitCsvLine(line);
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

void writeNewtonStepProbeCsv(const std::filesystem::path& path,
                             const vela::DeviceMesh& mesh,
                             const vela::DopingModel& doping,
                             const vela::DDSolution& state,
                             const vela::NewtonStepEvaluation& step)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write Newton step probe CSV: " + path.string());
    out << "node_id,x,y,psi,phin,phip,delta_psi_V,delta_phin_V,delta_phip_V,"
        << "delta_psi_minus_phin_V,delta_phip_minus_psi_V,"
        << "trial_psi,trial_phin,trial_phip,trial_electron_density_m3,trial_hole_density_m3,"
        << "psi_residual,phin_residual,phip_residual,"
        << "trial_psi_residual,trial_phin_residual,trial_phip_residual,"
        << "donors_m3,acceptors_m3,net_doping_m3,ni_eff_m3\n";
    const int n = static_cast<int>(mesh.numNodes());
    for (int i = 0; i < n; ++i) {
        const auto nodeId = static_cast<vela::Index>(i);
        const vela::Node& node = mesh.getNode(nodeId);
        out << nodeId << ','
            << node.x << ','
            << node.y << ','
            << state.psi(i) << ','
            << state.phin(i) << ','
            << state.phip(i) << ','
            << step.deltaPsi(i) << ','
            << step.deltaPhin(i) << ','
            << step.deltaPhip(i) << ','
            << (step.deltaPsi(i) - step.deltaPhin(i)) << ','
            << (step.deltaPhip(i) - step.deltaPsi(i)) << ','
            << step.trialSolution.psi(i) << ','
            << step.trialSolution.phin(i) << ','
            << step.trialSolution.phip(i) << ','
            << step.trialSolution.n(i) << ','
            << step.trialSolution.p(i) << ','
            << step.residual.raw(i) << ','
            << step.residual.raw(n + i) << ','
            << step.residual.raw(2 * n + i) << ','
            << step.trialResidual.raw(i) << ','
            << step.trialResidual.raw(n + i) << ','
            << step.trialResidual.raw(2 * n + i) << ','
            << doping.donors(nodeId) << ','
            << doping.acceptors(nodeId) << ','
            << doping.netDoping(nodeId) << ','
            << step.residual.intrinsicDensity[static_cast<std::size_t>(nodeId)]
            << '\n';
    }
}

nlohmann::json runNewtonStepProbe(const std::string& configFile, const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state = readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const vela::NewtonStepEvaluation step = solver.evaluateStep(state);
    writeNewtonStepProbeCsv(
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>()),
        problem.mesh,
        problem.doping,
        state,
        step);

    return {
        {"nodes", problem.mesh.numNodes()},
        {"scaled_state", step.residual.scaledState},
        {"potential_scale", step.residual.potentialScale},
        {"raw_step_norm", step.rawStepNorm},
        {"step_norm", step.stepNorm},
        {"block_residuals", {
            {"psi", step.residual.blockNorms.psi},
            {"phin", step.residual.blockNorms.phin},
            {"phip", step.residual.blockNorms.phip},
            {"combined", step.residual.blockNorms.combined},
        }},
        {"trial_block_residuals", {
            {"psi", step.trialResidual.blockNorms.psi},
            {"phin", step.trialResidual.blockNorms.phin},
            {"phip", step.trialResidual.blockNorms.phip},
            {"combined", step.trialResidual.blockNorms.combined},
        }},
    };
}

std::vector<bool> contactNodeMask(const vela::DeviceMesh& mesh)
{
    std::vector<bool> mask(static_cast<std::size_t>(mesh.numNodes()), false);
    for (const auto& contact : mesh.contacts()) {
        for (const vela::Index nodeId : contact.node_ids) {
            if (nodeId < mask.size())
                mask[static_cast<std::size_t>(nodeId)] = true;
        }
    }
    return mask;
}

bool coordinateInRange(const nlohmann::json& direction,
                       const std::string& axis,
                       vela::Real value)
{
    const std::string minKey = axis + "_min";
    const std::string maxKey = axis + "_max";
    const std::string minUmKey = axis + "_min_um";
    const std::string maxUmKey = axis + "_max_um";
    if (direction.contains(minKey) && value < direction.at(minKey).get<vela::Real>())
        return false;
    if (direction.contains(maxKey) && value > direction.at(maxKey).get<vela::Real>())
        return false;
    if (direction.contains(minUmKey) && value < 1.0e-6 * direction.at(minUmKey).get<vela::Real>())
        return false;
    if (direction.contains(maxUmKey) && value > 1.0e-6 * direction.at(maxUmKey).get<vela::Real>())
        return false;
    return true;
}

struct JvpProbeDirection {
    std::string name;
    std::string mode;
    vela::Real amplitude_V = 0.0;
    vela::DDSolution perturbation;
    int selectedNodes = 0;
};

JvpProbeDirection makeJvpProbeDirection(const vela::DeviceMesh& mesh,
                                        const nlohmann::json& direction)
{
    JvpProbeDirection probe;
    probe.name = direction.at("name").get<std::string>();
    probe.mode = direction.at("mode").get<std::string>();
    probe.amplitude_V = direction.value("amplitude_V", 1.0e-6);
    if (!std::isfinite(probe.amplitude_V) || probe.amplitude_V == 0.0)
        throw std::invalid_argument("newton_jvp_probe direction amplitude_V must be finite and non-zero.");

    const int n = static_cast<int>(mesh.numNodes());
    probe.perturbation.psi = vela::VectorXd::Zero(n);
    probe.perturbation.phin = vela::VectorXd::Zero(n);
    probe.perturbation.phip = vela::VectorXd::Zero(n);

    const bool excludeContacts = direction.value("exclude_contacts", true);
    const std::vector<bool> contactNodes = contactNodeMask(mesh);
    for (int i = 0; i < n; ++i) {
        const vela::Node& node = mesh.getNode(static_cast<vela::Index>(i));
        if (excludeContacts && contactNodes[static_cast<std::size_t>(i)])
            continue;
        if (!coordinateInRange(direction, "x", node.x) ||
            !coordinateInRange(direction, "y", node.y)) {
            continue;
        }

        if (probe.mode == "psi") {
            probe.perturbation.psi(i) = probe.amplitude_V;
        } else if (probe.mode == "phin") {
            probe.perturbation.phin(i) = probe.amplitude_V;
        } else if (probe.mode == "phip") {
            probe.perturbation.phip(i) = probe.amplitude_V;
        } else if (probe.mode == "psi_minus_phin") {
            probe.perturbation.psi(i) = 0.5 * probe.amplitude_V;
            probe.perturbation.phin(i) = -0.5 * probe.amplitude_V;
        } else if (probe.mode == "phip_minus_psi") {
            probe.perturbation.phip(i) = 0.5 * probe.amplitude_V;
            probe.perturbation.psi(i) = -0.5 * probe.amplitude_V;
        } else {
            throw std::invalid_argument(
                "newton_jvp_probe direction mode must be one of psi, phin, phip, "
                "psi_minus_phin, phip_minus_psi.");
        }
        ++probe.selectedNodes;
    }
    if (probe.selectedNodes == 0)
        throw std::runtime_error("newton_jvp_probe direction selected no nodes: " + probe.name);
    return probe;
}

vela::Real vectorBlockNorm(const vela::VectorXd& values, int offset, int count)
{
    return values.segment(offset, count).norm();
}

void writeNewtonJvpProbeCsv(const std::filesystem::path& path,
                            const std::vector<nlohmann::json>& rows)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write Newton JVP probe CSV: " + path.string());
    out << "direction,mode,amplitude_V,selected_nodes,perturbation_norm,"
        << "analytic_norm,finite_difference_norm,absolute_error,relative_error,"
        << "analytic_psi_norm,analytic_phin_norm,analytic_phip_norm,"
        << "finite_difference_psi_norm,finite_difference_phin_norm,finite_difference_phip_norm,"
        << "psi_relative_error,phin_relative_error,phip_relative_error\n";
    for (const auto& row : rows) {
        out << row.at("direction").get<std::string>() << ','
            << row.at("mode").get<std::string>() << ','
            << row.at("amplitude_V").get<vela::Real>() << ','
            << row.at("selected_nodes").get<int>() << ','
            << row.at("perturbation_norm").get<vela::Real>() << ','
            << row.at("analytic_norm").get<vela::Real>() << ','
            << row.at("finite_difference_norm").get<vela::Real>() << ','
            << row.at("absolute_error").get<vela::Real>() << ','
            << row.at("relative_error").get<vela::Real>() << ','
            << row.at("analytic_psi_norm").get<vela::Real>() << ','
            << row.at("analytic_phin_norm").get<vela::Real>() << ','
            << row.at("analytic_phip_norm").get<vela::Real>() << ','
            << row.at("finite_difference_psi_norm").get<vela::Real>() << ','
            << row.at("finite_difference_phin_norm").get<vela::Real>() << ','
            << row.at("finite_difference_phip_norm").get<vela::Real>() << ','
            << row.at("psi_relative_error").get<vela::Real>() << ','
            << row.at("phin_relative_error").get<vela::Real>() << ','
            << row.at("phip_relative_error").get<vela::Real>()
            << '\n';
    }
}

nlohmann::json runNewtonJvpProbe(const std::string& configFile, const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state = readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);

    if (!cfg.contains("directions") || !cfg.at("directions").is_array())
        throw std::invalid_argument("newton_jvp_probe requires a directions array.");

    const int n = static_cast<int>(problem.mesh.numNodes());
    std::vector<nlohmann::json> rows;
    vela::Real maxRelativeError = 0.0;
    for (const auto& directionConfig : cfg.at("directions")) {
        const JvpProbeDirection direction = makeJvpProbeDirection(
            problem.mesh, directionConfig);
        const vela::NewtonDirectionalDerivativeEvaluation jvp =
            solver.evaluateDirectionalDerivative(state, direction.perturbation);
        const vela::VectorXd error = jvp.analyticJv - jvp.finiteDifferenceJv;
        const vela::Real psiFd = vectorBlockNorm(jvp.finiteDifferenceJv, 0, n);
        const vela::Real phinFd = vectorBlockNorm(jvp.finiteDifferenceJv, n, n);
        const vela::Real phipFd = vectorBlockNorm(jvp.finiteDifferenceJv, 2 * n, n);
        const vela::Real psiError = vectorBlockNorm(error, 0, n);
        const vela::Real phinError = vectorBlockNorm(error, n, n);
        const vela::Real phipError = vectorBlockNorm(error, 2 * n, n);
        maxRelativeError = std::max(maxRelativeError, jvp.relativeError);
        rows.push_back({
            {"direction", direction.name},
            {"mode", direction.mode},
            {"amplitude_V", direction.amplitude_V},
            {"selected_nodes", direction.selectedNodes},
            {"perturbation_norm", jvp.perturbationNorm},
            {"analytic_norm", jvp.analyticNorm},
            {"finite_difference_norm", jvp.finiteDifferenceNorm},
            {"absolute_error", jvp.absoluteError},
            {"relative_error", jvp.relativeError},
            {"analytic_psi_norm", vectorBlockNorm(jvp.analyticJv, 0, n)},
            {"analytic_phin_norm", vectorBlockNorm(jvp.analyticJv, n, n)},
            {"analytic_phip_norm", vectorBlockNorm(jvp.analyticJv, 2 * n, n)},
            {"finite_difference_psi_norm", psiFd},
            {"finite_difference_phin_norm", phinFd},
            {"finite_difference_phip_norm", phipFd},
            {"psi_relative_error", psiError / std::max<vela::Real>(1.0, psiFd)},
            {"phin_relative_error", phinError / std::max<vela::Real>(1.0, phinFd)},
            {"phip_relative_error", phipError / std::max<vela::Real>(1.0, phipFd)},
        });
    }

    writeNewtonJvpProbeCsv(
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>()),
        rows);
    return {
        {"nodes", problem.mesh.numNodes()},
        {"direction_count", rows.size()},
        {"max_relative_error", maxRelativeError},
        {"directions", rows},
    };
}

void writeNewtonBlockStepProbeCsv(const std::filesystem::path& path,
                                  const vela::DeviceMesh& mesh,
                                  const vela::DopingModel& doping,
                                  const vela::DDSolution& state,
                                  const std::vector<vela::NewtonBlockStepEvaluation>& steps)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write Newton block-step probe CSV: " + path.string());
    out << "mode,node_id,x,y,psi,phin,phip,delta_psi_V,delta_phin_V,delta_phip_V,"
        << "delta_psi_minus_phin_V,delta_phip_minus_psi_V,"
        << "trial_psi,trial_phin,trial_phip,trial_electron_density_m3,trial_hole_density_m3,"
        << "psi_residual,phin_residual,phip_residual,"
        << "trial_psi_residual,trial_phin_residual,trial_phip_residual,"
        << "donors_m3,acceptors_m3,net_doping_m3,ni_eff_m3\n";
    const int n = static_cast<int>(mesh.numNodes());
    for (const auto& step : steps) {
        for (int i = 0; i < n; ++i) {
            const auto nodeId = static_cast<vela::Index>(i);
            const vela::Node& node = mesh.getNode(nodeId);
            out << step.mode << ','
                << nodeId << ','
                << node.x << ','
                << node.y << ','
                << state.psi(i) << ','
                << state.phin(i) << ','
                << state.phip(i) << ','
                << step.deltaPsi(i) << ','
                << step.deltaPhin(i) << ','
                << step.deltaPhip(i) << ','
                << (step.deltaPsi(i) - step.deltaPhin(i)) << ','
                << (step.deltaPhip(i) - step.deltaPsi(i)) << ','
                << step.trialSolution.psi(i) << ','
                << step.trialSolution.phin(i) << ','
                << step.trialSolution.phip(i) << ','
                << step.trialSolution.n(i) << ','
                << step.trialSolution.p(i) << ','
                << step.residual.raw(i) << ','
                << step.residual.raw(n + i) << ','
                << step.residual.raw(2 * n + i) << ','
                << step.trialResidual.raw(i) << ','
                << step.trialResidual.raw(n + i) << ','
                << step.trialResidual.raw(2 * n + i) << ','
                << doping.donors(nodeId) << ','
                << doping.acceptors(nodeId) << ','
                << doping.netDoping(nodeId) << ','
                << step.residual.intrinsicDensity[static_cast<std::size_t>(nodeId)]
                << '\n';
        }
    }
}

nlohmann::json blockStepSummaryJson(const vela::NewtonBlockStepEvaluation& step)
{
    return {
        {"mode", step.mode},
        {"raw_step_norm", step.rawStepNorm},
        {"step_norm", step.stepNorm},
        {"block_residuals", {
            {"psi", step.residual.blockNorms.psi},
            {"phin", step.residual.blockNorms.phin},
            {"phip", step.residual.blockNorms.phip},
            {"combined", step.residual.blockNorms.combined},
        }},
        {"trial_block_residuals", {
            {"psi", step.trialResidual.blockNorms.psi},
            {"phin", step.trialResidual.blockNorms.phin},
            {"phip", step.trialResidual.blockNorms.phip},
            {"combined", step.trialResidual.blockNorms.combined},
        }},
    };
}

nlohmann::json runNewtonBlockStepProbe(const std::string& configFile, const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state = readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);

    std::vector<std::string> modes = {"poisson_only", "carrier_only"};
    if (cfg.contains("block_modes")) {
        modes.clear();
        for (const auto& value : cfg.at("block_modes"))
            modes.push_back(value.get<std::string>());
    }
    if (modes.empty())
        throw std::invalid_argument("newton_block_step_probe requires at least one block mode.");

    std::vector<vela::NewtonBlockStepEvaluation> steps;
    std::vector<nlohmann::json> summaries;
    for (const std::string& mode : modes) {
        steps.push_back(solver.evaluateBlockStep(state, mode));
        summaries.push_back(blockStepSummaryJson(steps.back()));
    }

    writeNewtonBlockStepProbeCsv(
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>()),
        problem.mesh,
        problem.doping,
        state,
        steps);

    return {
        {"nodes", problem.mesh.numNodes()},
        {"block_step_count", steps.size()},
        {"block_steps", summaries},
    };
}

void writeNewtonRegularizedCarrierStepProbeCsv(
    const std::filesystem::path& path,
    const vela::DeviceMesh& mesh,
    const vela::DopingModel& doping,
    const vela::DDSolution& state,
    const std::vector<vela::NewtonRegularizedCarrierStepEvaluation>& steps)
{
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error(
            "Cannot write Newton regularized carrier-step probe CSV: " + path.string());
    }
    out << "regularization_scale,node_id,x,y,psi,phin,phip,"
        << "delta_psi_V,delta_phin_V,delta_phip_V,"
        << "delta_psi_minus_phin_V,delta_phip_minus_psi_V,"
        << "trial_psi,trial_phin,trial_phip,trial_electron_density_m3,"
        << "trial_hole_density_m3,"
        << "psi_residual,phin_residual,phip_residual,"
        << "trial_psi_residual,trial_phin_residual,trial_phip_residual,"
        << "donors_m3,acceptors_m3,net_doping_m3,ni_eff_m3,"
        << "raw_step_norm,step_norm,regularization_diagonal_norm\n";
    const int n = static_cast<int>(mesh.numNodes());
    for (const auto& step : steps) {
        for (int i = 0; i < n; ++i) {
            const auto nodeId = static_cast<vela::Index>(i);
            const vela::Node& node = mesh.getNode(nodeId);
            out << step.regularizationScale << ','
                << nodeId << ','
                << node.x << ','
                << node.y << ','
                << state.psi(i) << ','
                << state.phin(i) << ','
                << state.phip(i) << ','
                << step.deltaPsi(i) << ','
                << step.deltaPhin(i) << ','
                << step.deltaPhip(i) << ','
                << (step.deltaPsi(i) - step.deltaPhin(i)) << ','
                << (step.deltaPhip(i) - step.deltaPsi(i)) << ','
                << step.trialSolution.psi(i) << ','
                << step.trialSolution.phin(i) << ','
                << step.trialSolution.phip(i) << ','
                << step.trialSolution.n(i) << ','
                << step.trialSolution.p(i) << ','
                << step.residual.raw(i) << ','
                << step.residual.raw(n + i) << ','
                << step.residual.raw(2 * n + i) << ','
                << step.trialResidual.raw(i) << ','
                << step.trialResidual.raw(n + i) << ','
                << step.trialResidual.raw(2 * n + i) << ','
                << doping.donors(nodeId) << ','
                << doping.acceptors(nodeId) << ','
                << doping.netDoping(nodeId) << ','
                << step.residual.intrinsicDensity[static_cast<std::size_t>(nodeId)] << ','
                << step.rawStepNorm << ','
                << step.stepNorm << ','
                << step.regularizationDiagonalNorm
                << '\n';
        }
    }
}

nlohmann::json regularizedCarrierStepSummaryJson(
    const vela::NewtonRegularizedCarrierStepEvaluation& step)
{
    return {
        {"regularization_scale", step.regularizationScale},
        {"raw_step_norm", step.rawStepNorm},
        {"step_norm", step.stepNorm},
        {"regularization_diagonal_norm", step.regularizationDiagonalNorm},
        {"block_residuals", {
            {"psi", step.residual.blockNorms.psi},
            {"phin", step.residual.blockNorms.phin},
            {"phip", step.residual.blockNorms.phip},
            {"combined", step.residual.blockNorms.combined},
        }},
        {"trial_block_residuals", {
            {"psi", step.trialResidual.blockNorms.psi},
            {"phin", step.trialResidual.blockNorms.phin},
            {"phip", step.trialResidual.blockNorms.phip},
            {"combined", step.trialResidual.blockNorms.combined},
        }},
    };
}

nlohmann::json runNewtonRegularizedCarrierStepProbe(
    const std::string& configFile,
    const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state = readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);

    if (!cfg.contains("regularization_scales"))
        throw std::invalid_argument(
            "newton_regularized_carrier_step_probe requires regularization_scales.");

    std::vector<vela::Real> scales;
    for (const auto& value : cfg.at("regularization_scales"))
        scales.push_back(value.get<vela::Real>());
    if (scales.empty()) {
        throw std::invalid_argument(
            "newton_regularized_carrier_step_probe requires at least one scale.");
    }

    std::vector<vela::NewtonRegularizedCarrierStepEvaluation> steps;
    std::vector<nlohmann::json> summaries;
    for (const vela::Real scale : scales) {
        steps.push_back(solver.evaluateRegularizedCarrierStep(state, scale));
        summaries.push_back(regularizedCarrierStepSummaryJson(steps.back()));
    }

    writeNewtonRegularizedCarrierStepProbeCsv(
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>()),
        problem.mesh,
        problem.doping,
        state,
        steps);

    return {
        {"nodes", problem.mesh.numNodes()},
        {"regularized_step_count", steps.size()},
        {"regularized_steps", summaries},
    };
}

void writeNewtonCarrierRowProbeCsv(
    const std::filesystem::path& path,
    const vela::DeviceMesh& mesh,
    const vela::DopingModel& doping,
    const vela::NewtonCarrierRowDiagnosticsEvaluation& diagnostics)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write Newton carrier-row probe CSV: " + path.string());
    out << "node_id,x,y,electron_residual,hole_residual,"
        << "electron_diagonal,hole_diagonal,"
        << "electron_row_abs_sum,hole_row_abs_sum,"
        << "electron_offdiag_abs_sum,hole_offdiag_abs_sum,"
        << "electron_row_l2_norm,hole_row_l2_norm,"
        << "raw_delta_phin_V,raw_delta_phip_V,"
        << "capped_delta_phin_V,capped_delta_phip_V,"
        << "donors_m3,acceptors_m3,net_doping_m3,ni_eff_m3\n";
    for (const auto& row : diagnostics.rows) {
        const vela::Node& node = mesh.getNode(row.nodeId);
        out << row.nodeId << ','
            << node.x << ','
            << node.y << ','
            << row.electronResidual << ','
            << row.holeResidual << ','
            << row.electronDiagonal << ','
            << row.holeDiagonal << ','
            << row.electronRowAbsSum << ','
            << row.holeRowAbsSum << ','
            << row.electronOffdiagAbsSum << ','
            << row.holeOffdiagAbsSum << ','
            << row.electronRowL2Norm << ','
            << row.holeRowL2Norm << ','
            << row.rawDeltaPhin_V << ','
            << row.rawDeltaPhip_V << ','
            << row.cappedDeltaPhin_V << ','
            << row.cappedDeltaPhip_V << ','
            << doping.donors(row.nodeId) << ','
            << doping.acceptors(row.nodeId) << ','
            << doping.netDoping(row.nodeId) << ','
            << diagnostics.residual.intrinsicDensity[static_cast<std::size_t>(row.nodeId)]
            << '\n';
    }
}

nlohmann::json runNewtonCarrierRowProbe(const std::string& configFile, const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state = readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const vela::NewtonCarrierRowDiagnosticsEvaluation diagnostics =
        solver.evaluateCarrierRowDiagnostics(state);
    writeNewtonCarrierRowProbeCsv(
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>()),
        problem.mesh,
        problem.doping,
        diagnostics);
    return {
        {"nodes", problem.mesh.numNodes()},
        {"row_count", diagnostics.rows.size()},
        {"scaled_state", diagnostics.residual.scaledState},
        {"potential_scale", diagnostics.potentialScale},
        {"raw_carrier_step_norm", diagnostics.rawCarrierStepNorm},
        {"capped_carrier_step_norm", diagnostics.cappedCarrierStepNorm},
        {"block_residuals", {
            {"psi", diagnostics.residual.blockNorms.psi},
            {"phin", diagnostics.residual.blockNorms.phin},
            {"phip", diagnostics.residual.blockNorms.phip},
            {"combined", diagnostics.residual.blockNorms.combined},
        }},
    };
}

void writeNewtonCarrierTermProbeCsv(
    const std::filesystem::path& path,
    const vela::DeviceMesh& mesh,
    const vela::DopingModel& doping,
    const vela::NewtonCarrierTermDiagnosticsEvaluation& diagnostics,
    vela::Real electronImpactScale,
    vela::Real holeImpactScale)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write Newton carrier-term probe CSV: " + path.string());
    out << "node_id,x,y,"
        << "electron_flux,electron_recombination,electron_impact,electron_gauge,electron_boundary,"
        << "electron_term_sum,electron_residual,electron_adjusted_impact,"
        << "electron_adjusted_term_sum,electron_adjusted_residual,"
        << "hole_flux,hole_recombination,hole_impact,hole_gauge,hole_boundary,"
        << "hole_term_sum,hole_residual,hole_adjusted_impact,"
        << "hole_adjusted_term_sum,hole_adjusted_residual,"
        << "impact_electron_source,impact_hole_source,impact_combined_source,"
        << "donors_m3,acceptors_m3,net_doping_m3,ni_eff_m3\n";
    for (const auto& row : diagnostics.rows) {
        const vela::Node& node = mesh.getNode(row.nodeId);
        const vela::Real electronSum = row.electronFlux
            + row.electronRecombination
            + row.electronImpact
            + row.electronGauge
            + row.electronBoundary;
        const vela::Real holeSum = row.holeFlux
            + row.holeRecombination
            + row.holeImpact
            + row.holeGauge
            + row.holeBoundary;
        const vela::Real electronAdjustedImpact = row.electronImpact * electronImpactScale;
        const vela::Real holeAdjustedImpact = row.holeImpact * holeImpactScale;
        const vela::Real electronAdjustedSum = row.electronFlux
            + row.electronRecombination
            + electronAdjustedImpact
            + row.electronGauge
            + row.electronBoundary;
        const vela::Real holeAdjustedSum = row.holeFlux
            + row.holeRecombination
            + holeAdjustedImpact
            + row.holeGauge
            + row.holeBoundary;
        out << row.nodeId << ','
            << node.x << ','
            << node.y << ','
            << row.electronFlux << ','
            << row.electronRecombination << ','
            << row.electronImpact << ','
            << row.electronGauge << ','
            << row.electronBoundary << ','
            << electronSum << ','
            << row.electronResidual << ','
            << electronAdjustedImpact << ','
            << electronAdjustedSum << ','
            << electronAdjustedSum << ','
            << row.holeFlux << ','
            << row.holeRecombination << ','
            << row.holeImpact << ','
            << row.holeGauge << ','
            << row.holeBoundary << ','
            << holeSum << ','
            << row.holeResidual << ','
            << holeAdjustedImpact << ','
            << holeAdjustedSum << ','
            << holeAdjustedSum << ','
            << row.impactElectronSource << ','
            << row.impactHoleSource << ','
            << row.impactCombinedSource << ','
            << doping.donors(row.nodeId) << ','
            << doping.acceptors(row.nodeId) << ','
            << doping.netDoping(row.nodeId) << ','
            << diagnostics.residual.intrinsicDensity[static_cast<std::size_t>(row.nodeId)]
            << '\n';
    }
}

struct CarrierTermProbeOptions {
    vela::Real electronImpactScale = 1.0;
    vela::Real holeImpactScale = 1.0;
};

CarrierTermProbeOptions carrierTermProbeOptionsFromJson(const nlohmann::json& cfg)
{
    CarrierTermProbeOptions options;
    if (!cfg.contains("carrier_term_probe"))
        return options;
    const auto& probe = cfg.at("carrier_term_probe");
    if (!probe.is_object())
        throw std::runtime_error("carrier_term_probe must be an object.");
    options.electronImpactScale = probe.value("electron_impact_scale", 1.0);
    options.holeImpactScale = probe.value("hole_impact_scale", 1.0);
    if (!std::isfinite(options.electronImpactScale) || !std::isfinite(options.holeImpactScale))
        throw std::runtime_error("carrier_term_probe impact scales must be finite.");
    return options;
}

vela::NewtonBlockResidualInfo adjustedCarrierTermBlocks(
    const vela::NewtonCarrierTermDiagnosticsEvaluation& diagnostics,
    vela::Real electronImpactScale,
    vela::Real holeImpactScale)
{
    vela::Real phinSq = 0.0;
    vela::Real phipSq = 0.0;
    for (const auto& row : diagnostics.rows) {
        const vela::Real electronResidual = row.electronFlux
            + row.electronRecombination
            + row.electronImpact * electronImpactScale
            + row.electronGauge
            + row.electronBoundary;
        const vela::Real holeResidual = row.holeFlux
            + row.holeRecombination
            + row.holeImpact * holeImpactScale
            + row.holeGauge
            + row.holeBoundary;
        phinSq += electronResidual * electronResidual;
        phipSq += holeResidual * holeResidual;
    }
    const vela::Real psi = diagnostics.residual.blockNorms.psi;
    const vela::Real phin = std::sqrt(phinSq);
    const vela::Real phip = std::sqrt(phipSq);
    return {psi, phin, phip, std::sqrt(psi * psi + phin * phin + phip * phip)};
}

nlohmann::json runNewtonCarrierTermProbe(const std::string& configFile, const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state = readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const vela::NewtonCarrierTermDiagnosticsEvaluation diagnostics =
        solver.evaluateCarrierTermDiagnostics(state);
    const CarrierTermProbeOptions options = carrierTermProbeOptionsFromJson(cfg);
    writeNewtonCarrierTermProbeCsv(
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>()),
        problem.mesh,
        problem.doping,
        diagnostics,
        options.electronImpactScale,
        options.holeImpactScale);
    const vela::NewtonBlockResidualInfo adjustedBlocks = adjustedCarrierTermBlocks(
        diagnostics,
        options.electronImpactScale,
        options.holeImpactScale);
    return {
        {"nodes", problem.mesh.numNodes()},
        {"row_count", diagnostics.rows.size()},
        {"scaled_state", diagnostics.residual.scaledState},
        {"potential_scale", diagnostics.residual.potentialScale},
        {"carrier_term_probe", {
            {"electron_impact_scale", options.electronImpactScale},
            {"hole_impact_scale", options.holeImpactScale},
        }},
        {"block_residuals", {
            {"psi", diagnostics.residual.blockNorms.psi},
            {"phin", diagnostics.residual.blockNorms.phin},
            {"phip", diagnostics.residual.blockNorms.phip},
            {"combined", diagnostics.residual.blockNorms.combined},
        }},
        {"adjusted_block_residuals", {
            {"psi", adjustedBlocks.psi},
            {"phin", adjustedBlocks.phin},
            {"phip", adjustedBlocks.phip},
            {"combined", adjustedBlocks.combined},
        }},
    };
}

void writeSgEdgeFluxProbeCsv(
    const std::filesystem::path& path,
    const std::vector<vela::CoupledDDEdgeFluxDiagnostic>& edges)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write SG edge flux probe CSV: " + path.string());
    out << "edge_id,node0,node1,x0,y0,x1,y1,length_m,couple_m,"
        << "net_doping_avg_m3,ni0_m3,ni1_m3,psi0_V,psi1_V,phin0_V,phin1_V,"
        << "phip0_V,phip1_V,electric_field_V_m,electron_mobility_m2_V_s,"
        << "hole_mobility_m2_V_s,electron_flux,hole_flux\n";
    out << std::setprecision(17);
    for (const auto& edge : edges) {
        out << edge.edgeId << ','
            << edge.node0 << ','
            << edge.node1 << ','
            << edge.x0 << ','
            << edge.y0 << ','
            << edge.x1 << ','
            << edge.y1 << ','
            << edge.length_m << ','
            << edge.couple_m << ','
            << edge.netDopingAvg_m3 << ','
            << edge.ni0_m3 << ','
            << edge.ni1_m3 << ','
            << edge.psi0_V << ','
            << edge.psi1_V << ','
            << edge.phin0_V << ','
            << edge.phin1_V << ','
            << edge.phip0_V << ','
            << edge.phip1_V << ','
            << edge.electricField_V_m << ','
            << edge.electronMobility_m2_V_s << ','
            << edge.holeMobility_m2_V_s << ','
            << edge.electronFlux << ','
            << edge.holeFlux << '\n';
    }
}

nlohmann::json runSgEdgeFluxProbe(const std::string& configFile, const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state = readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const std::vector<vela::CoupledDDEdgeFluxDiagnostic> edges =
        solver.evaluateSgEdgeFluxDiagnostics(state);
    writeSgEdgeFluxProbeCsv(
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>()),
        edges);
    return {
        {"nodes", problem.mesh.numNodes()},
        {"edge_count", edges.size()},
    };
}

void writeEdgeMobilityProbeCsv(const std::filesystem::path& path,
                               const vela::DeviceMesh& mesh,
                               const vela::DopingModel& doping,
                               const vela::DDSolution& state,
                               const vela::MaterialDatabase& matdb,
                               const vela::NewtonConfig& newton)
{
    const auto edgeCells = vela::detail::buildEdgeCellMap(mesh);
    const auto cellMaterials = vela::detail::buildCellMaterials(
        mesh, matdb, newton.temperature_K);
    const std::unique_ptr<vela::MobilityModel> mobility =
        vela::makeMobilityModel(newton.mobility);

    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write edge mobility probe CSV: " + path.string());
    out << "edge_id,node0,node1,x0,y0,x1,y1,length_m,couple_m,"
        << "net_doping_avg_m3,electric_field_V_m,electron_qf_field_V_m,"
        << "hole_qf_field_V_m,electron_mobility_field_V_m,"
        << "hole_mobility_field_V_m,electron_low_field_mobility_m2_V_s,"
        << "hole_low_field_mobility_m2_V_s,electron_final_mobility_m2_V_s,"
        << "hole_final_mobility_m2_V_s,electron_mobility_limiter,"
        << "hole_mobility_limiter,adjacent_cell_count\n";

    for (vela::Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const vela::Edge& edge = mesh.getEdge(edgeId);
        const vela::Node& n0 = mesh.getNode(edge.n0);
        const vela::Node& n1 = mesh.getNode(edge.n1);
        const vela::Real length = edge.length;
        if (length <= 0.0)
            continue;
        const int i0 = static_cast<int>(edge.n0);
        const int i1 = static_cast<int>(edge.n1);
        const vela::Real electricField = std::abs(state.psi(i1) - state.psi(i0)) / length;
        const vela::Real electronQfField = std::abs(state.phin(i1) - state.phin(i0)) / length;
        const vela::Real holeQfField = std::abs(state.phip(i1) - state.phip(i0)) / length;
        const vela::Real electronMobilityField =
            newton.mobility.highFieldDrivingForce == "quasi_fermi_gradient"
            ? electronQfField
            : electricField;
        const vela::Real holeMobilityField =
            newton.mobility.highFieldDrivingForce == "quasi_fermi_gradient"
            ? holeQfField
            : electricField;
        const vela::Real electronLowField = vela::detail::edgeMobility(
            edgeCells, mesh, doping, *mobility, cellMaterials, edgeId,
            vela::CarrierType::Electron, 0.0, &newton.mobility, nullptr);
        const vela::Real holeLowField = vela::detail::edgeMobility(
            edgeCells, mesh, doping, *mobility, cellMaterials, edgeId,
            vela::CarrierType::Hole, 0.0, &newton.mobility, nullptr);
        const vela::Real electronFinal = vela::detail::edgeMobility(
            edgeCells, mesh, doping, *mobility, cellMaterials, edgeId,
            vela::CarrierType::Electron, electronMobilityField, &newton.mobility, &state.psi);
        const vela::Real holeFinal = vela::detail::edgeMobility(
            edgeCells, mesh, doping, *mobility, cellMaterials, edgeId,
            vela::CarrierType::Hole, holeMobilityField, &newton.mobility, &state.psi);
        const vela::Real electronLimiter =
            electronLowField > 0.0 ? electronFinal / electronLowField : 0.0;
        const vela::Real holeLimiter =
            holeLowField > 0.0 ? holeFinal / holeLowField : 0.0;
        const vela::Real netDoping = 0.5 * (
            doping.netDoping(edge.n0) + doping.netDoping(edge.n1));

        out << edgeId << ','
            << edge.n0 << ','
            << edge.n1 << ','
            << n0.x << ','
            << n0.y << ','
            << n1.x << ','
            << n1.y << ','
            << length << ','
            << edge.couple << ','
            << netDoping << ','
            << electricField << ','
            << electronQfField << ','
            << holeQfField << ','
            << electronMobilityField << ','
            << holeMobilityField << ','
            << electronLowField << ','
            << holeLowField << ','
            << electronFinal << ','
            << holeFinal << ','
            << electronLimiter << ','
            << holeLimiter << ','
            << edgeCells[edgeId].size()
            << '\n';
    }
}

nlohmann::json runEdgeMobilityProbe(const std::string& configFile, const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state = readExternalState(cfgDir, cfg, problem.mesh.numNodes());

    writeEdgeMobilityProbeCsv(
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>()),
        problem.mesh,
        problem.doping,
        state,
        problem.matdb,
        problem.newton);

    return {
        {"nodes", problem.mesh.numNodes()},
        {"edge_count", problem.mesh.numEdges()},
        {"temperature_K", problem.newton.temperature_K},
        {"mobility_model", problem.newton.mobility.model},
        {"high_field_driving_force", problem.newton.mobility.highFieldDrivingForce},
    };
}

nlohmann::json runNewtonJacobianBlockProbe(const std::string& configFile,
                                           const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const std::filesystem::path statePath =
        resolvePath(cfgDir, cfg.at("state_file").get<std::string>());
    const vela::DDSolution state =
        vela::readDDSolutionStateCsv(statePath, problem.mesh.numNodes());
    const vela::Real fdStep = cfg.value("finite_difference_step", 1.0e-7);

    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const auto rows = solver.evaluateJacobianBlockAudit(state, fdStep);

    const std::filesystem::path outputPath =
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>());
    if (!outputPath.parent_path().empty())
        std::filesystem::create_directories(outputPath.parent_path());
    std::ofstream out(outputPath);
    if (!out.is_open())
        throw std::runtime_error("Cannot write jacobian block probe CSV: " + outputPath.string());
    out << "block,analytic_norm,fd_norm,diff_norm,rel_diff\n";
    out << std::setprecision(17);
    for (const auto& row : rows) {
        out << row.block << ','
            << row.analyticNorm << ','
            << row.fdNorm << ','
            << row.diffNorm << ','
            << row.relDiff << '\n';
    }

    return {
        {"nodes", problem.mesh.numNodes()},
        {"blocks", rows.size()},
        {"output_csv", outputPath.string()},
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
        } else if (type == "newton_step_probe") {
            status.update(runNewtonStepProbe(configFile, cfg));
        } else if (type == "newton_jvp_probe") {
            status.update(runNewtonJvpProbe(configFile, cfg));
        } else if (type == "newton_block_step_probe") {
            status.update(runNewtonBlockStepProbe(configFile, cfg));
        } else if (type == "newton_regularized_carrier_step_probe") {
            status.update(runNewtonRegularizedCarrierStepProbe(configFile, cfg));
        } else if (type == "newton_carrier_row_probe") {
            status.update(runNewtonCarrierRowProbe(configFile, cfg));
        } else if (type == "newton_carrier_term_probe") {
            status.update(runNewtonCarrierTermProbe(configFile, cfg));
        } else if (type == "sg_edge_flux_probe") {
            status.update(runSgEdgeFluxProbe(configFile, cfg));
        } else if (type == "edge_mobility_probe") {
            status.update(runEdgeMobilityProbe(configFile, cfg));
        } else if (type == "newton_jacobian_block_probe") {
            status.update(runNewtonJacobianBlockProbe(configFile, cfg));
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
