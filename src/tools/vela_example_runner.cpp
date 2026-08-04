#include "vela/io/CsvUtils.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/io/DDSolutionCsv.h"
#include "vela/io/MeshReader.h"
#include "vela/io/NeutralMeshReader.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/post/ContactCurrent.h"
#include "vela/core/RuntimeLog.h"
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
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

void usage(const char* argv0)
{
    std::cerr << "Usage: " << argv0
              << " --config <simulation.json> [--mesh-report] [--log <auto|off|path>] [--log-profile <minimal|default|debug>]\n";
}

struct RuntimeLogOverrideGuard {
    explicit RuntimeLogOverrideGuard(const vela::RuntimeLogCliOverrides& overrides)
    {
        vela::setRuntimeLogCliOverrides(overrides);
    }
    ~RuntimeLogOverrideGuard()
    {
        vela::clearRuntimeLogCliOverrides();
    }
};

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

nlohmann::json releaseBVConfigAuditJson(const vela::ReleaseBVConfigAuditMetadata& metadata)
{
    return {
        {"model", metadata.model},
        {"coupling_mode", metadata.couplingMode},
        {"driving_force", metadata.drivingForce},
        {"parameter_set", metadata.parameterSet},
        {"A_scale", metadata.aScale},
        {"B_scale", metadata.bScale},
        {"switchField_V_per_cm", metadata.switchField_V_per_cm},
        {"cutoff_minField_V_per_cm", metadata.minimumField_V_per_cm},
        {"smoothing", metadata.smoothing},
        {"RefDens_electron_cm_minus3", metadata.electronRefDens_cm3},
        {"RefDens_hole_cm_minus3", metadata.holeRefDens_cm3},
        {"source_mapping_mode", metadata.sourceMappingMode},
        {"current_magnitude_mode", metadata.currentMagnitudeMode},
        {"lambda_ava", metadata.lambdaAva},
        {"depth_2D_um", metadata.depth2D_um},
        {"current_normalization", metadata.currentNormalization},
        {"qG_normalization", metadata.qGNormalization},
        {"audit_csv", metadata.auditCsvFile},
        {"audit_summary", metadata.auditSummaryFile},
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
nlohmann::json carrierRowConvergenceJson(
    const vela::NewtonCarrierRowConvergenceEvaluation& evaluation)
{
    nlohmann::json violations = nlohmann::json::array();
    for (const auto& row : evaluation.violations) {
        violations.push_back({
            {"node_id", row.nodeId},
            {"carrier", row.carrier},
            {"residual", row.residual},
            {"scale", row.scale},
            {"ratio", row.ratio},
            {"flux", row.flux},
            {"srh", row.recombination},
            {"impact", row.impact},
        });
    }
    return {
        {"enabled", evaluation.enabled},
        {"enforced", evaluation.enforced},
        {"satisfied", evaluation.satisfied},
        {"eps_row", evaluation.epsRow},
        {"violation_count", evaluation.violations.size()},
        {"max_ratio", evaluation.maxRatio},
        {"max_ratio_node", evaluation.maxRatioNode},
        {"max_ratio_carrier", evaluation.maxRatioCarrier},
        {"violations", std::move(violations)},
    };
}
nlohmann::json carrierRowRecoveryJson(
    const vela::NewtonCarrierRowRecoveryResult& recovery)
{
    return {
        {"attempted", recovery.attempted},
        {"mode", recovery.mode},
        {"electron_rows_updated", recovery.electronRowsUpdated},
        {"hole_rows_updated", recovery.holeRowsUpdated},
        {"density_passes", recovery.densityPasses},
        {"cycles_attempted", recovery.cyclesAttempted},
        {"density_converged", recovery.densityConverged},
        {"max_density_relative_change", recovery.maxDensityRelativeChange},
        {"max_psi_delta_V", recovery.maxPsiDelta_V},
        {"max_carrier_density_ratio", recovery.maxCarrierDensityRatio},
    };
}

nlohmann::json blockResidualsJson(const vela::NewtonBlockResidualInfo& blocks)
{
    return {
        {"psi", blocks.psi},
        {"phin", blocks.phin},
        {"phip", blocks.phip},
        {"combined", blocks.combined},
    };
}

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

vela::DDSolution readExternalState(const std::filesystem::path& cfgDir,
                                   const nlohmann::json& cfg,
                                   vela::Index nodeCount);

NewtonProblem loadNewtonProblem(const std::string& configFile, const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    const vela::UnitScalingConfig scaling = vela::parseUnitScalingConfig(cfg);

    const bool useNeutralMesh = cfg.contains("neutral_mesh_dir");
    if (useNeutralMesh && !scaling.isUnitScaling()) {
        throw std::runtime_error(
            "neutral_mesh_dir requires scaling.mode = 'unit_scaling' because neutral mesh uses "
            "x_um/y_um coordinates and cm^-3 doping values.");
    }
    vela::DeviceMesh mesh;
    if (useNeutralMesh) {
        vela::NeutralMeshReader reader;
        mesh = reader.readDirectory(resolvePath(cfgDir, cfg.at("neutral_mesh_dir").get<std::string>()),
                                    scaling);
    } else {
        vela::JsonMeshReader reader;
        mesh = reader.read(
            resolvePath(cfgDir, cfg.at("mesh_file").get<std::string>()),
            scaling);
    }
    mesh.buildBoxGeometry(vela::parseBoxGeometryOptions(cfg));

    vela::MaterialDatabase matdb(scaling);
    if (cfg.contains("materials_file"))
        matdb.loadJson(resolvePath(cfgDir, cfg.at("materials_file").get<std::string>()), scaling);

    vela::DopingModel doping(mesh.numNodes());
    if (cfg.contains("node_doping_file")) {
        doping = readNodeDopingCsv(
            resolvePath(cfgDir, cfg.at("node_doping_file").get<std::string>()),
            mesh.numNodes(),
            scaling);
    } else if (useNeutralMesh &&
               std::filesystem::exists(
                   std::filesystem::path(resolvePath(cfgDir, cfg.at("neutral_mesh_dir").get<std::string>())) /
                   "doping.csv")) {
        vela::NeutralMeshReader reader;
        doping = reader.readDopingCsv(
            std::filesystem::path(resolvePath(cfgDir, cfg.at("neutral_mesh_dir").get<std::string>())) /
                "doping.csv",
            mesh.numNodes(),
            scaling);
    } else {
        doping = vela::DopingModel::fromMeshAndRegions(mesh, vela::parseDopingSpecs(cfg, scaling));
    }
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
            result.solution,
            problem.newton.inputScaling);
    }

    return NewtonCliResult{std::move(problem.mesh), std::move(result)};
}

nlohmann::json runNewtonSolveFromState(const std::string& configFile,
                                       const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution initial =
        readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    vela::NewtonResult result = solver.solve(initial);

    if (cfg.contains("output_state_file")) {
        vela::writeDDSolutionStateCsv(
            resolvePath(cfgDir, cfg.at("output_state_file").get<std::string>()),
            result.solution,
            problem.newton.inputScaling);
    }
    if (cfg.contains("output_vtk")) {
        vela::writeDDSolutionVTK(
            resolvePath(cfgDir, cfg.at("output_vtk").get<std::string>()),
            problem.mesh,
            problem.doping,
            result.solution,
            problem.newton.inputScaling);
    }

    // On-state terminal current extraction (no re-solve). Computing ContactCurrent
    // directly on the converged solution avoids the state relaxation that a
    // dc_sweep restart introduces, so the reported current belongs to the same
    // state whose node densities are inspected by the branch gate. The scaling and
    // model setup mirror DCSweep so the value is directly comparable to a dc_sweep
    // terminal current.
    vela::DDScalingSpec ddScaling;
    if (problem.newton.inputScaling.isUnitScaling()) {
        const vela::UnitScalingSystem sc = vela::UnitScalingSystem::fromInputs(
            problem.newton.temperature_K,
            (11.7 * vela::constants::eps0),
            vela::UnitScalingSystem::autoInputsFrom(
                problem.mesh, problem.doping, problem.matdb, 1e10),
            vela::UnitScalingReferenceConfig{},
            problem.newton.inputScaling.unitSystem());
        ddScaling.enabled = true;
        ddScaling.V0 = sc.V0();
        ddScaling.C0 = sc.C0();
        ddScaling.mu0 = sc.mu0();
        ddScaling.D0 = sc.D0();
        ddScaling.L0 = sc.L0();
        ddScaling.permittivityReference_F_per_m = (11.7 * vela::constants::eps0);
        ddScaling.unitSystem = problem.newton.inputScaling.unitSystem();
        ddScaling.chargeAreaFactor = problem.newton.inputScaling.unitSystem().chargeAreaFactor();
        ddScaling.chargeLineFactor = problem.newton.inputScaling.unitSystem().chargeLineFactor();
        ddScaling.fieldFromCoordinateDeltaFactor = problem.newton.inputScaling.unitSystem().fieldFromCoordinateDeltaFactor();
        ddScaling.currentDensityLineIntegralFactor =
            problem.newton.inputScaling.unitSystem().currentDensityAM2PerInternal() *
            problem.newton.inputScaling.unitSystem().lengthMPerInternal();
    }
    const vela::ContactCurrent contactCurrent(
        problem.mesh,
        problem.matdb,
        problem.doping,
        problem.newton.mobility,
        problem.newton.temperature_K,
        ddScaling,
        problem.newton.bandgapNarrowing);

    constexpr vela::Real perMeterToPerMicron = 1.0e-6;
    nlohmann::json contactCurrentsJson = nlohmann::json::object();
    vela::Real drivenCurrentPerMicron = 0.0;
    vela::Real drivenAbsBias = -1.0;
    for (const auto& [name, bias] : problem.biases) {
        const vela::ContactCurrentResult cc =
            contactCurrent.compute(result.solution, name);
        const vela::Real currentPerMicron = cc.totalCurrent * perMeterToPerMicron;
        contactCurrentsJson[name] = currentPerMicron;
        if (std::abs(bias) > drivenAbsBias) {
            drivenAbsBias = std::abs(bias);
            drivenCurrentPerMicron = currentPerMicron;
        }
    }

    return {
        {"nodes", problem.mesh.numNodes()},
        {"converged", result.converged},
        {"iterations", result.iters},
        {"initial_residual", result.initialResidualNorm},
        {"final_residual", result.finalResidualNorm},
        {"convergence_reason", result.convergenceReason},
        {"failure_reason", result.failureDiagnostics.failureReason},
        {"final_block_residuals", blockResidualsJson(result.finalBlockNorms)},
        {"carrier_row_convergence", carrierRowConvergenceJson(result.finalCarrierRowConvergence)},
        {"carrier_row_recovery", carrierRowRecoveryJson(result.carrierRowRecovery)},
        {"contact_currents_A_per_um", contactCurrentsJson},
        {"current_total_A_per_um", drivenCurrentPerMicron},
    };
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
        if (seen[static_cast<std::size_t>(node)])
            throw std::runtime_error("Scalar field CSV has duplicate node_id: " + path.string());
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
            scaling.concentrationToInternal(std::stod(cells[donorsCol])),
            scaling.concentrationToInternal(std::stod(cells[acceptorsCol])));
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

vela::DDSolution readFeedbackReplacementState(
    const std::filesystem::path& cfgDir,
    const nlohmann::json& cfg,
    vela::Index nodeCount,
    const vela::UnitScalingConfig& scaling)
{
    const std::filesystem::path fieldsDir =
        resolvePath(cfgDir, cfg.at("feedback_state_fields_dir").get<std::string>());
    const auto read = [&](const char* field) {
        return readNodeScalarCsv(fieldsDir / (std::string(field) + "_region0.csv"), nodeCount);
    };
    const std::vector<vela::Real> phin = read("eQuasiFermiPotential");
    const std::vector<vela::Real> phip = read("hQuasiFermiPotential");
    const std::vector<vela::Real> electrons_m3 = read("eDensity_m3");
    const std::vector<vela::Real> holes_m3 = read("hDensity_m3");
    const auto& units = scaling.unitSystem();

    const int n = static_cast<int>(nodeCount);
    vela::DDSolution state;
    state.psi = vela::VectorXd::Zero(n);
    state.phin.resize(n);
    state.phip.resize(n);
    state.n.resize(n);
    state.p.resize(n);
    for (int i = 0; i < n; ++i) {
        state.phin(i) = phin[static_cast<std::size_t>(i)];
        state.phip(i) = phip[static_cast<std::size_t>(i)];
        state.n(i) = units.m3ToInternalConcentration(
            electrons_m3[static_cast<std::size_t>(i)]);
        state.p(i) = units.m3ToInternalConcentration(
            holes_m3[static_cast<std::size_t>(i)]);
    }
    return state;
}

void writeResidualProbeCsv(const std::filesystem::path& path,
                           const vela::DeviceMesh& mesh,
                           const vela::DopingModel& doping,
                           const vela::DDSolution& state,
                           const vela::NewtonResidualEvaluation& residual,
                           const vela::UnitScalingConfig& scaling)
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
        const vela::PhysicalUnitSystem& units = scaling.unitSystem();
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
            << units.internalConcentrationToM3(doping.donors(nodeId)) << ','
            << units.internalConcentrationToM3(doping.acceptors(nodeId)) << ','
            << units.internalConcentrationToM3(doping.netDoping(nodeId)) << ','
            << units.internalConcentrationToM3(
                   residual.intrinsicDensity[static_cast<std::size_t>(nodeId)])
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
        residual,
        problem.newton.inputScaling);

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
                             const vela::NewtonStepEvaluation& step,
                             const vela::UnitScalingConfig& scaling)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write Newton step probe CSV: " + path.string());
    out << std::setprecision(17);
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
        const vela::PhysicalUnitSystem& units = scaling.unitSystem();
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
            << units.internalConcentrationToM3(step.trialSolution.n(i)) << ','
            << units.internalConcentrationToM3(step.trialSolution.p(i)) << ','
            << step.residual.raw(i) << ','
            << step.residual.raw(n + i) << ','
            << step.residual.raw(2 * n + i) << ','
            << step.trialResidual.raw(i) << ','
            << step.trialResidual.raw(n + i) << ','
            << step.trialResidual.raw(2 * n + i) << ','
            << units.internalConcentrationToM3(doping.donors(nodeId)) << ','
            << units.internalConcentrationToM3(doping.acceptors(nodeId)) << ','
            << units.internalConcentrationToM3(doping.netDoping(nodeId)) << ','
            << units.internalConcentrationToM3(
                   step.residual.intrinsicDensity[static_cast<std::size_t>(nodeId)])
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
        step,
        problem.newton.inputScaling);

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

void writeNewtonFeedbackSubstitutionProbeCsv(
    const std::filesystem::path& path,
    const vela::DeviceMesh& mesh,
    const vela::DDSolution& state,
    const vela::DDSolution& replacement,
    const std::vector<vela::NewtonFeedbackSubstitutionEvaluation>& evaluations,
    const vela::UnitScalingConfig& scaling)
{
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error(
            "Cannot write Newton feedback-substitution probe CSV: " + path.string());
    }
    out << std::setprecision(17);
    out << "variant,node_id,x,y,is_contact,"
        << "baseline_phin_V,baseline_phip_V,replacement_phin_V,replacement_phip_V,"
        << "replacement_electron_density_m3,replacement_hole_density_m3,"
        << "psi_residual,electron_residual,hole_residual,"
        << "desired_psi_residual,desired_electron_residual,desired_hole_residual,"
        << "delta_psi_V,delta_phin_V,delta_phip_V,"
        << "trial_phin_V,trial_phip_V,"
        << "carrier_only_delta_phin_V,carrier_only_delta_phip_V,"
        << "carrier_only_trial_phin_V,carrier_only_trial_phip_V,"
        << "production_trial_psi_residual,production_trial_electron_residual,"
        << "production_trial_hole_residual,"
        << "electron_flux,electron_recombination,electron_impact,electron_gauge,"
        << "electron_boundary,electron_term_sum,electron_closure_error,"
        << "hole_flux,hole_recombination,hole_impact,hole_gauge,hole_boundary,"
        << "hole_term_sum,hole_closure_error\n";
    const int n = static_cast<int>(mesh.numNodes());
    const std::vector<bool> contacts = contactNodeMask(mesh);
    const auto& units = scaling.unitSystem();
    for (const auto& evaluation : evaluations) {
        if (evaluation.carrierTerms.size() != static_cast<std::size_t>(n))
            throw std::runtime_error("feedback-substitution carrier-term support mismatch.");
        for (int i = 0; i < n; ++i) {
            const auto nodeId = static_cast<vela::Index>(i);
            const vela::Node& node = mesh.getNode(nodeId);
            const auto& term = evaluation.carrierTerms[static_cast<std::size_t>(i)];
            const vela::Real electronTermSum =
                term.electronFlux + term.electronRecombination +
                term.electronImpact + term.electronGauge + term.electronBoundary;
            const vela::Real holeTermSum =
                term.holeFlux + term.holeRecombination +
                term.holeImpact + term.holeGauge + term.holeBoundary;
            const vela::Real electronResidual = evaluation.residual.raw(n + i);
            const vela::Real holeResidual = evaluation.residual.raw(2 * n + i);
            out << evaluation.variant << ','
                << nodeId << ','
                << node.x << ','
                << node.y << ','
                << static_cast<int>(contacts[static_cast<std::size_t>(i)]) << ','
                << state.phin(i) << ','
                << state.phip(i) << ','
                << replacement.phin(i) << ','
                << replacement.phip(i) << ','
                << units.internalConcentrationToM3(replacement.n(i)) << ','
                << units.internalConcentrationToM3(replacement.p(i)) << ','
                << evaluation.residual.raw(i) << ','
                << electronResidual << ','
                << holeResidual << ','
                << evaluation.desiredResidual(i) << ','
                << evaluation.desiredResidual(n + i) << ','
                << evaluation.desiredResidual(2 * n + i) << ','
                << evaluation.deltaPsi(i) << ','
                << evaluation.deltaPhin(i) << ','
                << evaluation.deltaPhip(i) << ','
                << state.phin(i) + evaluation.deltaPhin(i) << ','
                << state.phip(i) + evaluation.deltaPhip(i) << ','
                << evaluation.carrierOnlyDeltaPhin(i) << ','
                << evaluation.carrierOnlyDeltaPhip(i) << ','
                << state.phin(i) + evaluation.carrierOnlyDeltaPhin(i) << ','
                << state.phip(i) + evaluation.carrierOnlyDeltaPhip(i) << ','
                << evaluation.productionTrialResidual.raw(i) << ','
                << evaluation.productionTrialResidual.raw(n + i) << ','
                << evaluation.productionTrialResidual.raw(2 * n + i) << ','
                << term.electronFlux << ','
                << term.electronRecombination << ','
                << term.electronImpact << ','
                << term.electronGauge << ','
                << term.electronBoundary << ','
                << electronTermSum << ','
                << electronTermSum - electronResidual << ','
                << term.holeFlux << ','
                << term.holeRecombination << ','
                << term.holeImpact << ','
                << term.holeGauge << ','
                << term.holeBoundary << ','
                << holeTermSum << ','
                << holeTermSum - holeResidual
                << '\n';
        }
    }
}

nlohmann::json runNewtonFeedbackSubstitutionProbe(
    const std::string& configFile,
    const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state =
        readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::DDSolution replacement = readFeedbackReplacementState(
        cfgDir,
        cfg,
        problem.mesh.numNodes(),
        problem.newton.inputScaling);
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const auto evaluations =
        solver.evaluateFeedbackSubstitutions(state, replacement);
    writeNewtonFeedbackSubstitutionProbeCsv(
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>()),
        problem.mesh,
        state,
        replacement,
        evaluations,
        problem.newton.inputScaling);

    nlohmann::json variants = nlohmann::json::array();
    for (const auto& evaluation : evaluations) {
        variants.push_back({
            {"variant", evaluation.variant},
            {"replaces_density", evaluation.replacesDensity},
            {"replaces_quasi_fermi", evaluation.replacesQuasiFermi},
            {"raw_step_norm", evaluation.rawStepNorm},
            {"step_norm", evaluation.stepNorm},
            {"carrier_only_raw_step_norm", evaluation.carrierOnlyRawStepNorm},
            {"carrier_only_step_norm", evaluation.carrierOnlyStepNorm},
            {"block_residuals", blockResidualsJson(evaluation.residual.blockNorms)},
            {"production_trial_block_residuals",
             blockResidualsJson(evaluation.productionTrialResidual.blockNorms)},
        });
    }
    return {
        {"nodes", problem.mesh.numNodes()},
        {"variants", variants},
        {"contract", {
            {"baseline_state", "converged_vela_avalanche_on_exact_state"},
            {"replacement_state", "sentaurus_avalanche_on_exact_state"},
            {"jacobian", "single_production_baseline_jacobian"},
            {"boundary_rows", "baseline_preserved"},
            {"substitution", "frozen_residual_operator_inputs"},
            {"production_defaults_changed", false},
        }},
    };
}

void writeNewtonPoissonQfpCrossBlockCsv(
    const std::filesystem::path& path,
    const vela::DeviceMesh& mesh,
    const vela::DDSolution& state,
    const vela::DDSolution& replacement,
    const vela::NewtonPoissonQfpCrossBlockEvaluation& evaluation)
{
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error(
            "Cannot write Poisson-QFP cross-block probe CSV: " + path.string());
    }
    out << std::setprecision(17);
    out << "node_id,x,y,is_contact,baseline_phin_V,baseline_phip_V,"
        << "replacement_phin_V,replacement_phip_V,target_delta_phin_V,"
        << "target_delta_phip_V,psi_residual,electron_residual,hole_residual,"
        << "independent_delta_psi_V,independent_delta_phin_V,"
        << "independent_delta_phip_V,no_psi_qfp_delta_psi_V,"
        << "no_psi_qfp_delta_phin_V,no_psi_qfp_delta_phip_V,"
        << "no_qfp_psi_delta_psi_V,no_qfp_psi_delta_phin_V,"
        << "no_qfp_psi_delta_phip_V,schur_delta_psi_V,schur_delta_phin_V,"
        << "schur_delta_phip_V,full_raw_delta_psi_V,full_raw_delta_phin_V,"
        << "full_raw_delta_phip_V,full_capped_delta_psi_V,"
        << "full_capped_delta_phin_V,full_capped_delta_phip_V,"
        << "psi_qfp_product,qfp_psi_electron_product,qfp_psi_hole_product,"
        << "full_minus_independent_delta_phin_V,"
        << "full_minus_independent_delta_phip_V,"
        << "qfp_fd_direction_phin_V,qfp_fd_direction_phip_V,"
        << "psi_fd_direction_V,B_directional_analytic,"
        << "B_directional_finite_difference,C_directional_electron_analytic,"
        << "C_directional_electron_finite_difference,"
        << "C_directional_hole_analytic,"
        << "C_directional_hole_finite_difference,"
        << "leave_out_transport_boundary_delta_phin_V,"
        << "leave_out_transport_boundary_delta_phip_V,"
        << "leave_out_srh_auger_delta_phin_V,"
        << "leave_out_srh_auger_delta_phip_V,"
        << "leave_out_sg_avalanche_delta_phin_V,"
        << "leave_out_sg_avalanche_delta_phip_V,"
        << "only_transport_boundary_delta_phin_V,"
        << "only_transport_boundary_delta_phip_V,"
        << "only_srh_auger_delta_phin_V,"
        << "only_srh_auger_delta_phip_V,"
        << "only_sg_avalanche_delta_phin_V,"
        << "only_sg_avalanche_delta_phip_V\n";
    const int n = static_cast<int>(mesh.numNodes());
    const std::vector<bool> contacts = contactNodeMask(mesh);
    const auto& transport = evaluation.loopComponents.at(0);
    const auto& recombination = evaluation.loopComponents.at(1);
    const auto& avalanche = evaluation.loopComponents.at(2);
    for (int i = 0; i < n; ++i) {
        const auto nodeId = static_cast<vela::Index>(i);
        const vela::Node& node = mesh.getNode(nodeId);
        out << nodeId << ','
            << node.x << ','
            << node.y << ','
            << static_cast<int>(contacts[static_cast<std::size_t>(i)]) << ','
            << state.phin(i) << ','
            << state.phip(i) << ','
            << replacement.phin(i) << ','
            << replacement.phip(i) << ','
            << evaluation.targetDeltaPhin(i) << ','
            << evaluation.targetDeltaPhip(i) << ','
            << evaluation.residual.raw(i) << ','
            << evaluation.residual.raw(n + i) << ','
            << evaluation.residual.raw(2 * n + i) << ','
            << evaluation.independentDeltaPsi(i) << ','
            << evaluation.independentDeltaPhin(i) << ','
            << evaluation.independentDeltaPhip(i) << ','
            << evaluation.noPsiQfpDeltaPsi(i) << ','
            << evaluation.noPsiQfpDeltaPhin(i) << ','
            << evaluation.noPsiQfpDeltaPhip(i) << ','
            << evaluation.noQfpPsiDeltaPsi(i) << ','
            << evaluation.noQfpPsiDeltaPhin(i) << ','
            << evaluation.noQfpPsiDeltaPhip(i) << ','
            << evaluation.schurDeltaPsi(i) << ','
            << evaluation.schurDeltaPhin(i) << ','
            << evaluation.schurDeltaPhip(i) << ','
            << evaluation.fullRawDeltaPsi(i) << ','
            << evaluation.fullRawDeltaPhin(i) << ','
            << evaluation.fullRawDeltaPhip(i) << ','
            << evaluation.fullCappedDeltaPsi(i) << ','
            << evaluation.fullCappedDeltaPhin(i) << ','
            << evaluation.fullCappedDeltaPhip(i) << ','
            << evaluation.psiQfpProduct(i) << ','
            << evaluation.qfpPsiProduct(i) << ','
            << evaluation.qfpPsiProduct(n + i) << ','
            << evaluation.fullRawDeltaPhin(i)
                - evaluation.independentDeltaPhin(i) << ','
            << evaluation.fullRawDeltaPhip(i)
                - evaluation.independentDeltaPhip(i)
            << ',' << evaluation.qfpFiniteDifferenceDirectionPhin(i)
            << ',' << evaluation.qfpFiniteDifferenceDirectionPhip(i)
            << ',' << evaluation.psiFiniteDifferenceDirection(i)
            << ',' << evaluation.analyticPsiQfpDirectionalDerivative(i)
            << ',' << evaluation.finiteDifferencePsiQfpDirectionalDerivative(i)
            << ',' << evaluation.analyticQfpPsiDirectionalDerivative(i)
            << ',' << evaluation.finiteDifferenceQfpPsiDirectionalDerivative(i)
            << ',' << evaluation.analyticQfpPsiDirectionalDerivative(n + i)
            << ',' << evaluation.finiteDifferenceQfpPsiDirectionalDerivative(n + i)
            << ',' << transport.leaveOutDeltaPhin(i)
            << ',' << transport.leaveOutDeltaPhip(i)
            << ',' << recombination.leaveOutDeltaPhin(i)
            << ',' << recombination.leaveOutDeltaPhip(i)
            << ',' << avalanche.leaveOutDeltaPhin(i)
            << ',' << avalanche.leaveOutDeltaPhip(i)
            << ',' << transport.onlyDeltaPhin(i)
            << ',' << transport.onlyDeltaPhip(i)
            << ',' << recombination.onlyDeltaPhin(i)
            << ',' << recombination.onlyDeltaPhip(i)
            << ',' << avalanche.onlyDeltaPhin(i)
            << ',' << avalanche.onlyDeltaPhip(i)
            << '\n';
    }
}

void writeJacobianBlockEntries(
    std::ofstream& out,
    const std::string& name,
    const vela::SparseMatrixd& block,
    int rowOffset,
    int columnOffset)
{
    for (int column = 0; column < block.outerSize(); ++column) {
        for (vela::SparseMatrixd::InnerIterator it(block, column); it; ++it) {
            out << name << ','
                << it.row() << ','
                << it.col() << ','
                << rowOffset + it.row() << ','
                << columnOffset + it.col() << ','
                << it.value() << '\n';
        }
    }
}

void writeNewtonPoissonQfpJacobianBlocksCsv(
    const std::filesystem::path& path,
    const vela::NewtonPoissonQfpCrossBlockEvaluation& evaluation,
    int nodeCount)
{
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error(
            "Cannot write Poisson-QFP Jacobian block CSV: " + path.string());
    }
    out << std::setprecision(17);
    out << "block,row_local,col_local,row_global,col_global,value\n";
    writeJacobianBlockEntries(
        out, "J_psi_psi", evaluation.jacobianPsiPsi, 0, 0);
    writeJacobianBlockEntries(
        out, "J_psi_qfp", evaluation.jacobianPsiQfp, 0, nodeCount);
    writeJacobianBlockEntries(
        out, "J_qfp_psi", evaluation.jacobianQfpPsi, nodeCount, 0);
    writeJacobianBlockEntries(
        out, "J_qfp_qfp", evaluation.jacobianQfpQfp, nodeCount, nodeCount);
}

nlohmann::json matrixConditionJson(
    const vela::NewtonMatrixConditionEstimate& condition)
{
    return {
        {"rows", condition.rows},
        {"columns", condition.columns},
        {"numerical_rank", condition.numericalRank},
        {"largest_singular_value", condition.largestSingularValue},
        {"smallest_resolved_singular_value",
         condition.smallestResolvedSingularValue},
        {"resolved_condition_number", condition.resolvedConditionNumber},
    };
}

void writeNewtonPoissonQfpSchurLoopCsv(
    const std::filesystem::path& path,
    const vela::DeviceMesh& mesh,
    const vela::NewtonPoissonQfpCrossBlockEvaluation& evaluation)
{
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error(
            "Cannot write Poisson-QFP Schur-loop CSV: " + path.string());
    }
    out << std::setprecision(17);
    out << "matrix,component,row_carrier,row_node,row_x,row_y,"
        << "col_variable,col_carrier,col_node,col_x,col_y,value,sign\n";
    const int n = static_cast<int>(mesh.numNodes());
    const auto writeEntries =
        [&](const std::string& matrixName,
            const std::string& componentName,
            const vela::SparseMatrixd& matrix,
            bool qfpColumns) {
            for (int column = 0; column < matrix.outerSize(); ++column) {
                for (vela::SparseMatrixd::InnerIterator it(
                         matrix, column);
                     it;
                     ++it) {
                    const int rowNode = it.row() % n;
                    const int colNode = it.col() % n;
                    const vela::Node& row =
                        mesh.getNode(static_cast<vela::Index>(rowNode));
                    const vela::Node& columnNode =
                        mesh.getNode(static_cast<vela::Index>(colNode));
                    const std::string rowCarrier =
                        it.row() < n ? "electron" : "hole";
                    const std::string colCarrier = !qfpColumns
                        ? ""
                        : (it.col() < n ? "electron" : "hole");
                    const double value = it.value();
                    out << matrixName << ','
                        << componentName << ','
                        << rowCarrier << ','
                        << rowNode << ','
                        << row.x << ','
                        << row.y << ','
                        << (qfpColumns ? "qfp" : "psi") << ','
                        << colCarrier << ','
                        << colNode << ','
                        << columnNode.x << ','
                        << columnNode.y << ','
                        << value << ','
                        << (value > 0.0 ? "positive"
                                        : (value < 0.0 ? "negative" : "zero"))
                        << '\n';
                }
            }
        };
    writeEntries(
        "C_Ainv_B", "all", evaluation.effectiveSchurLoop, true);
    for (const auto& component : evaluation.loopComponents) {
        writeEntries(
            "C_component",
            component.name,
            component.jacobianQfpPsi,
            false);
        writeEntries(
            "C_Ainv_B_component",
            component.name,
            component.effectiveLoop,
            true);
    }
}

nlohmann::json runNewtonPoissonQfpCrossBlockProbe(
    const std::string& configFile,
    const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state =
        readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::DDSolution replacement = readFeedbackReplacementState(
        cfgDir,
        cfg,
        problem.mesh.numNodes(),
        problem.newton.inputScaling);
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const auto evaluation =
        solver.evaluatePoissonQfpCrossBlockDecomposition(state, replacement);
    const std::filesystem::path outputPath =
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>());
    const std::filesystem::path blockPath = cfg.contains("jacobian_blocks_csv")
        ? std::filesystem::path(
            resolvePath(cfgDir, cfg.at("jacobian_blocks_csv").get<std::string>()))
        : outputPath.parent_path() /
            (outputPath.stem().string() + "_jacobian_blocks.csv");
    const std::filesystem::path schurLoopPath = cfg.contains("schur_loop_csv")
        ? std::filesystem::path(
            resolvePath(cfgDir, cfg.at("schur_loop_csv").get<std::string>()))
        : outputPath.parent_path() /
            (outputPath.stem().string() + "_schur_loop.csv");
    writeNewtonPoissonQfpCrossBlockCsv(
        outputPath, problem.mesh, state, replacement, evaluation);
    writeNewtonPoissonQfpJacobianBlocksCsv(
        blockPath,
        evaluation,
        static_cast<int>(problem.mesh.numNodes()));
    writeNewtonPoissonQfpSchurLoopCsv(
        schurLoopPath, problem.mesh, evaluation);

    nlohmann::json componentNorms = nlohmann::json::object();
    for (const auto& component : evaluation.loopComponents) {
        componentNorms[component.name] = {
            {"J_qfp_psi_norm", component.jacobianQfpPsi.norm()},
            {"C_Ainv_B_norm", component.effectiveLoop.norm()},
        };
    }

    return {
        {"nodes", problem.mesh.numNodes()},
        {"output_csv", outputPath.string()},
        {"jacobian_blocks_csv", blockPath.string()},
        {"schur_loop_csv", schurLoopPath.string()},
        {"jacobian_block_norms", {
            {"J_psi_psi", evaluation.jacobianPsiPsiNorm},
            {"J_psi_qfp", evaluation.jacobianPsiQfpNorm},
            {"J_qfp_psi", evaluation.jacobianQfpPsiNorm},
            {"J_qfp_qfp", evaluation.jacobianQfpQfpNorm},
        }},
        {"full_linear_closure_norm", evaluation.fullLinearClosureNorm},
        {"schur_closure_norm", evaluation.schurClosureNorm},
        {"schur_relative_closure", evaluation.schurRelativeClosure},
        {"loop_component_closure_norm",
         evaluation.loopComponentClosureNorm},
        {"loop_component_norms", componentNorms},
        {"condition_estimates", {
            {"J_psi_psi",
             matrixConditionJson(evaluation.jacobianPsiPsiCondition)},
            {"J_psi_psi_l2_equilibrated",
             matrixConditionJson(
                 evaluation.jacobianPsiPsiEquilibratedCondition)},
            {"J_qfp_qfp",
             matrixConditionJson(evaluation.jacobianQfpQfpCondition)},
            {"J_qfp_qfp_l2_equilibrated",
             matrixConditionJson(
                 evaluation.jacobianQfpQfpEquilibratedCondition)},
            {"schur",
             matrixConditionJson(evaluation.schurCondition)},
            {"schur_l2_equilibrated",
             matrixConditionJson(evaluation.schurEquilibratedCondition)},
            {"C_Ainv_B",
             matrixConditionJson(
                 evaluation.effectiveSchurLoopCondition)},
        }},
        {"directional_derivative_check", {
            {"relative_step",
             evaluation.finiteDifferenceRelativeStep},
            {"J_psi_qfp_relative_error",
             evaluation.psiQfpDirectionalDerivativeRelativeError},
            {"J_psi_qfp_analytic_norm",
             evaluation.analyticPsiQfpDirectionalDerivative.norm()},
            {"J_psi_qfp_finite_difference_norm",
             evaluation.finiteDifferencePsiQfpDirectionalDerivative.norm()},
            {"J_qfp_psi_relative_error",
             evaluation.qfpPsiDirectionalDerivativeRelativeError},
            {"J_qfp_psi_analytic_norm",
             evaluation.analyticQfpPsiDirectionalDerivative.norm()},
            {"J_qfp_psi_finite_difference_norm",
             evaluation.finiteDifferenceQfpPsiDirectionalDerivative.norm()},
            {"residual", "production_baseline"},
            {"scheme", "double_symmetric"},
        }},
        {"contract", {
            {"residual", "qfp_only_frozen_substitution"},
            {"jacobian", "single_production_baseline_jacobian"},
            {"partition", "psi_vs_electron_and_hole_qfp"},
            {"effective_loop", "C_A_inverse_B"},
            {"model_components", {
                "transport_boundary",
                "srh_auger",
                "sg_avalanche",
            }},
            {"counterfactuals", {
                "independent_blocks",
                "remove_J_psi_qfp",
                "remove_J_qfp_psi",
                "full_schur",
            }},
            {"boundary_rows", "baseline_preserved"},
            {"production_defaults_changed", false},
        }},
    };
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
    out << std::setprecision(17);
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

void writeNewtonCarrierBlockColumnsCsv(
    const std::filesystem::path& path,
    const vela::DeviceMesh& mesh,
    const vela::NewtonCarrierBlockDecompositionEvaluation& diagnostics)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write carrier-block columns CSV: " + path.string());
    out << std::setprecision(17);
    out << "reduced_column,carrier,node_id,x,y,diagonal,column_l2_norm,"
        << "electron_row_l2_norm,hole_row_l2_norm,diagonal_fraction,"
        << "cross_carrier_row_fraction,continuity_row_weight,residual,"
        << "full_delta_qfp_V\n";
    for (const auto& column : diagnostics.columns) {
        const vela::Node& node = mesh.getNode(column.nodeId);
        out << column.reducedColumn << ',' << column.carrier << ','
            << column.nodeId << ',' << node.x << ',' << node.y << ','
            << column.diagonal << ',' << column.columnL2Norm << ','
            << column.electronRowL2Norm << ',' << column.holeRowL2Norm << ','
            << column.diagonalFraction << ',' << column.crossCarrierRowFraction << ','
            << column.continuityRowWeight << ',' << column.residual << ','
            << column.fullDeltaQfp_V << '\n';
    }
}

void writeNewtonCarrierBlockSingularModesCsv(
    const std::filesystem::path& path,
    const vela::NewtonCarrierBlockDecompositionEvaluation& diagnostics)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write carrier-block singular modes CSV: " + path.string());
    out << std::setprecision(17);
    out << "mode_index,singular_value,relative_singular_value,rhs_projection,"
        << "rhs_energy_fraction,step_amplitude,step_energy_fraction,"
        << "transport_jacobian_projection,recombination_jacobian_projection,"
        << "avalanche_diagonal_jacobian_projection,"
        << "avalanche_cross_jacobian_projection,jacobian_projection_closure,"
        << "transport_rhs_projection,recombination_rhs_projection,"
        << "avalanche_rhs_projection,rhs_projection_closure,"
        << "no_cross_carrier_step_amplitude,"
        << "no_recombination_step_amplitude,no_avalanche_step_amplitude,"
        << "transport_only_step_amplitude,"
        << "right_electron_fraction,left_electron_fraction,"
        << "top_right_carrier,top_right_node,top_right_value,"
        << "top_left_carrier,top_left_node,top_left_value\n";
    for (const auto& mode : diagnostics.singularModes) {
        out << mode.modeIndex << ',' << mode.singularValue << ','
            << mode.relativeSingularValue << ',' << mode.rhsProjection << ','
            << mode.rhsEnergyFraction << ',' << mode.stepAmplitude << ','
            << mode.stepEnergyFraction << ','
            << mode.transportJacobianProjection << ','
            << mode.recombinationJacobianProjection << ','
            << mode.avalancheDiagonalJacobianProjection << ','
            << mode.avalancheCrossJacobianProjection << ','
            << mode.jacobianProjectionClosure << ','
            << mode.transportRhsProjection << ','
            << mode.recombinationRhsProjection << ','
            << mode.avalancheRhsProjection << ','
            << mode.rhsProjectionClosure << ','
            << mode.noCrossCarrierStepAmplitude << ','
            << mode.noRecombinationStepAmplitude << ','
            << mode.noAvalancheStepAmplitude << ','
            << mode.transportOnlyStepAmplitude << ','
            << mode.rightElectronFraction << ','
            << mode.leftElectronFraction << ',' << mode.topRightCarrier << ','
            << mode.topRightNode << ',' << mode.topRightValue << ','
            << mode.topLeftCarrier << ',' << mode.topLeftNode << ','
            << mode.topLeftValue << '\n';
    }
}

void writeNewtonCarrierBlockSolveVariantsCsv(
    const std::filesystem::path& path,
    const vela::NewtonCarrierBlockDecompositionEvaluation& diagnostics)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write carrier-block solve variants CSV: " + path.string());
    out << std::setprecision(17);
    out << "variant,scaled_step_norm,physical_step_norm_V,"
        << "relative_difference_from_full,cosine_with_full,"
        << "relative_linear_closure\n";
    for (const auto& variant : diagnostics.solveVariants) {
        out << variant.name << ',' << variant.scaledStepNorm << ','
            << variant.physicalStepNorm_V << ','
            << variant.relativeDifferenceFromFull << ','
            << variant.cosineWithFull << ','
            << variant.relativeLinearClosure << '\n';
    }
}

void writeNewtonCarrierBlockSolveNodesCsv(
    const std::filesystem::path& path,
    const vela::DeviceMesh& mesh,
    const vela::NewtonCarrierBlockDecompositionEvaluation& diagnostics)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write carrier-block solve nodes CSV: " + path.string());
    out << std::setprecision(17);
    out << "variant,node_id,x,y,delta_phin_V,delta_phip_V\n";
    for (const auto& variant : diagnostics.solveVariants) {
        for (int nodeId = 0; nodeId < variant.deltaPhin.size(); ++nodeId) {
            const vela::Node& node = mesh.getNode(static_cast<vela::Index>(nodeId));
            out << variant.name << ',' << nodeId << ',' << node.x << ',' << node.y
                << ',' << variant.deltaPhin(nodeId) << ','
                << variant.deltaPhip(nodeId) << '\n';
        }
    }
}

nlohmann::json runNewtonCarrierBlockDecompositionProbe(
    const std::string& configFile,
    const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state =
        readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const auto diagnostics =
        solver.evaluateCarrierBlockDecomposition(state);
    const std::filesystem::path prefix =
        resolvePath(cfgDir, cfg.at("output_prefix").get<std::string>());
    if (!prefix.parent_path().empty())
        std::filesystem::create_directories(prefix.parent_path());
    const auto suffixed = [&](const std::string& suffix) {
        return std::filesystem::path(prefix.string() + suffix);
    };
    writeNewtonCarrierBlockColumnsCsv(
        suffixed("_columns.csv"), problem.mesh, diagnostics);
    writeNewtonCarrierBlockSingularModesCsv(
        suffixed("_singular_modes.csv"), diagnostics);
    writeNewtonCarrierBlockSolveVariantsCsv(
        suffixed("_solve_variants.csv"), diagnostics);
    writeNewtonCarrierBlockSolveNodesCsv(
        suffixed("_solve_nodes.csv"), problem.mesh, diagnostics);

    nlohmann::json summary = {
        {"nodes", problem.mesh.numNodes()},
        {"free_electron_unknowns", diagnostics.freeElectronUnknowns},
        {"free_hole_unknowns", diagnostics.freeHoleUnknowns},
        {"raw_condition", matrixConditionJson(diagnostics.rawCondition)},
        {"row_scaled_condition", matrixConditionJson(diagnostics.rowScaledCondition)},
        {"l2_equilibrated_condition", matrixConditionJson(diagnostics.l2EquilibratedCondition)},
        {"electron_electron_norm", diagnostics.electronElectronNorm},
        {"electron_hole_norm", diagnostics.electronHoleNorm},
        {"hole_electron_norm", diagnostics.holeElectronNorm},
        {"hole_hole_norm", diagnostics.holeHoleNorm},
        {"cross_carrier_norm_fraction", diagnostics.crossCarrierNormFraction},
        {"recombination_cross_norm", diagnostics.recombinationCrossNorm},
        {"avalanche_cross_norm", diagnostics.avalancheCrossNorm},
        {"transport_cross_norm", diagnostics.transportCrossNorm},
        {"free_column_norm_spread", diagnostics.freeColumnNormSpread},
        {"free_row_norm_spread", diagnostics.freeRowNormSpread},
        {"row_weight_spread", diagnostics.rowWeightSpread},
        {"block_residuals", {
            {"psi", diagnostics.residual.blockNorms.psi},
            {"phin", diagnostics.residual.blockNorms.phin},
            {"phip", diagnostics.residual.blockNorms.phip},
            {"combined", diagnostics.residual.blockNorms.combined},
        }},
    };
    std::ofstream summaryOut(suffixed("_summary.json"));
    if (!summaryOut.is_open())
        throw std::runtime_error("Cannot write carrier-block summary JSON: " +
                                 suffixed("_summary.json").string());
    summaryOut << std::setw(2) << summary << '\n';
    return summary;
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
    out << std::setprecision(17);
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
    const std::vector<vela::CoupledDDEdgeFluxDiagnostic>& edges,
    const vela::UnitScalingConfig& scaling)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write SG edge flux probe CSV: " + path.string());
    out << "edge_id,node0,node1,x0,y0,x1,y1,length_m,couple_m,"
        << "net_doping_avg_m3,ni0_m3,ni1_m3,psi0_V,psi1_V,phin0_V,phin1_V,"
        << "phip0_V,phip1_V,electric_field_V_m,electron_mobility_m2_V_s,"
        << "hole_mobility_m2_V_s,electron_flux,hole_flux\n";
    const vela::PhysicalUnitSystem& units = scaling.unitSystem();
    out << std::setprecision(17);
    for (const auto& edge : edges) {
        out << edge.edgeId << ','
            << edge.node0 << ','
            << edge.node1 << ','
            << units.internalLengthToMeters(edge.x0) << ','
            << units.internalLengthToMeters(edge.y0) << ','
            << units.internalLengthToMeters(edge.x1) << ','
            << units.internalLengthToMeters(edge.y1) << ','
            << units.internalLengthToMeters(edge.length_m) << ','
            << units.internalLengthToMeters(edge.couple_m) << ','
            << units.internalConcentrationToM3(edge.netDopingAvg_m3) << ','
            << units.internalConcentrationToM3(edge.ni0_m3) << ','
            << units.internalConcentrationToM3(edge.ni1_m3) << ','
            << edge.psi0_V << ','
            << edge.psi1_V << ','
            << edge.phin0_V << ','
            << edge.phin1_V << ','
            << edge.phip0_V << ','
            << edge.phip1_V << ','
            << units.internalElectricFieldToVPerM(edge.electricField_V_m) << ','
            << units.internalMobilityToM2PerVS(edge.electronMobility_m2_V_s) << ','
            << units.internalMobilityToM2PerVS(edge.holeMobility_m2_V_s) << ','
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
        edges,
        problem.newton.inputScaling);
    return {
        {"nodes", problem.mesh.numNodes()},
        {"edge_count", edges.size()},
    };
}

void writeTransportEdgeJacobianProbeCsv(
    const std::filesystem::path& path,
    const std::vector<vela::CoupledDDTransportEdgeJacobianDiagnostic>& records,
    const vela::UnitScalingConfig& scaling)
{
    const vela::PhysicalUnitSystem& units = scaling.unitSystem();
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error(
            "Cannot write transport edge Jacobian probe CSV: " + path.string());
    out << "edge_id,carrier,node0,node1,row_node,column_node,row_endpoint,column_endpoint,"
        << "length_m,couple_m,qfp_drive_V_m,d_qfp_drive_d_column_per_m,"
        << "mobility_m2_V_s,d_mobility_d_column_m2_V2_s,"
        << "bernoulli_node0,bernoulli_node1,carrier_density_node0_m3,"
        << "carrier_density_node1_m3,flux_physical,flux_scaled,row_sign,"
        << "production_frozen_mobility_derivative_physical,"
        << "frozen_mobility_fd_derivative_physical,"
        << "live_mobility_fd_derivative_physical,"
        << "mobility_response_derivative_physical,"
        << "live_minus_frozen_fd_derivative_physical,"
        << "bernoulli_qfp_derivative_physical,"
        << "carrier_population_derivative_physical,"
        << "production_row_derivative_scaled,live_mobility_row_derivative_scaled,"
        << "row_constrained,column_constrained,"
        << "contact_eliminated_production_edge_derivative,contact_identity_entry,"
        << "continuity_row_weight,solver_production_edge_derivative,"
        << "solver_contact_identity_entry\n";
    out << std::setprecision(17);
    for (const auto& record : records) {
        out << record.edgeId << ',' << record.carrier << ','
            << record.node0 << ',' << record.node1 << ','
            << record.rowNode << ',' << record.columnNode << ','
            << record.rowEndpoint << ',' << record.columnEndpoint << ','
            << units.internalLengthToMeters(record.lengthInternal) << ','
            << units.internalLengthToMeters(record.coupleInternal) << ','
            << units.internalElectricFieldToVPerM(record.qfpDriveInternal) << ','
            << units.internalElectricFieldToVPerM(record.dQfpDriveDColumnInternal) << ','
            << units.internalMobilityToM2PerVS(record.mobilityInternal) << ','
            << units.internalMobilityToM2PerVS(record.dMobilityDColumnInternal) << ','
            << record.bernoulliNode0 << ',' << record.bernoulliNode1 << ','
            << units.internalConcentrationToM3(record.carrierDensityNode0Internal) << ','
            << units.internalConcentrationToM3(record.carrierDensityNode1Internal) << ','
            << record.fluxPhysical << ',' << record.fluxScaled << ','
            << record.rowSign << ','
            << record.productionFrozenMobilityDerivativePhysical << ','
            << record.frozenMobilityFiniteDifferenceDerivativePhysical << ','
            << record.liveMobilityFiniteDifferenceDerivativePhysical << ','
            << record.mobilityResponseDerivativePhysical << ','
            << record.liveMinusFrozenFiniteDifferenceDerivativePhysical << ','
            << record.bernoulliQfpDerivativePhysical << ','
            << record.carrierPopulationDerivativePhysical << ','
            << record.productionRowDerivativeScaled << ','
            << record.liveMobilityRowDerivativeScaled << ','
            << static_cast<int>(record.rowConstrained) << ','
            << static_cast<int>(record.columnConstrained) << ','
            << record.contactEliminatedProductionEdgeDerivative << ','
            << record.contactIdentityEntry << ','
            << record.continuityRowWeight << ','
            << record.solverProductionEdgeDerivative << ','
            << record.solverContactIdentityEntry << '\n';
    }
}

nlohmann::json runTransportEdgeJacobianProbe(
    const std::string& configFile,
    const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution state =
        readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const vela::Real step = cfg.value("physical_finite_difference_step_V", 1.0e-7);
    const auto records =
        solver.evaluateTransportEdgeJacobianDiagnostics(state, step);
    writeTransportEdgeJacobianProbeCsv(
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>()),
        records,
        problem.newton.inputScaling);
    return {
        {"nodes", problem.mesh.numNodes()},
        {"record_count", records.size()},
        {"physical_finite_difference_step_V", step},
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
    const vela::PhysicalUnitSystem& units = newton.inputScaling.unitSystem();

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
        const vela::Real fieldFactor = units.fieldFromCoordinateDeltaFactor();
        const vela::Real electricField = std::abs(state.psi(i1) - state.psi(i0)) / length * fieldFactor;
        const vela::Real electronQfField = std::abs(state.phin(i1) - state.phin(i0)) / length * fieldFactor;
        const vela::Real holeQfField = std::abs(state.phip(i1) - state.phip(i0)) / length * fieldFactor;
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
            << units.internalLengthToMeters(n0.x) << ','
            << units.internalLengthToMeters(n0.y) << ','
            << units.internalLengthToMeters(n1.x) << ','
            << units.internalLengthToMeters(n1.y) << ','
            << units.internalLengthToMeters(length) << ','
            << units.internalLengthToMeters(edge.couple) << ','
            << units.internalConcentrationToM3(netDoping) << ','
            << units.internalElectricFieldToVPerM(electricField) << ','
            << units.internalElectricFieldToVPerM(electronQfField) << ','
            << units.internalElectricFieldToVPerM(holeQfField) << ','
            << units.internalElectricFieldToVPerM(electronMobilityField) << ','
            << units.internalElectricFieldToVPerM(holeMobilityField) << ','
            << units.internalMobilityToM2PerVS(electronLowField) << ','
            << units.internalMobilityToM2PerVS(holeLowField) << ','
            << units.internalMobilityToM2PerVS(electronFinal) << ','
            << units.internalMobilityToM2PerVS(holeFinal) << ','
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
        vela::readDDSolutionStateCsv(
            statePath,
            problem.mesh.numNodes(),
            problem.newton.inputScaling);
    const vela::Real fdStep = cfg.value("finite_difference_step", 1.0e-7);
    const std::string fdMode =
        cfg.value("finite_difference_mode", std::string("double_symmetric"));
    std::vector<std::string> blocks;
    if (cfg.contains("blocks")) {
        if (!cfg.at("blocks").is_array())
            throw std::invalid_argument("newton_jacobian_block_probe blocks must be an array.");
        for (const auto& value : cfg.at("blocks"))
            blocks.push_back(value.get<std::string>());
    }

    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const auto rows =
        solver.evaluateJacobianBlockAudit(state, fdStep, blocks, fdMode);

    const std::filesystem::path outputPath =
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>());
    if (!outputPath.parent_path().empty())
        std::filesystem::create_directories(outputPath.parent_path());
    std::ofstream out(outputPath);
    if (!out.is_open())
        throw std::runtime_error("Cannot write jacobian block probe CSV: " + outputPath.string());
    out << "block,analytic_norm,fd_norm,diff_norm,rel_diff,"
        << "analytic_psi_column_norm,fd_psi_column_norm,diff_psi_column_norm,rel_psi_column_diff,"
        << "analytic_phin_column_norm,fd_phin_column_norm,diff_phin_column_norm,rel_phin_column_diff,"
        << "analytic_phip_column_norm,fd_phip_column_norm,diff_phip_column_norm,rel_phip_column_diff,"
        << "analytic_electron_phin_norm,fd_electron_phin_norm,diff_electron_phin_norm,rel_electron_phin_diff,"
        << "analytic_electron_phip_norm,fd_electron_phip_norm,diff_electron_phip_norm,rel_electron_phip_diff,"
        << "analytic_hole_phin_norm,fd_hole_phin_norm,diff_hole_phin_norm,rel_hole_phin_diff,"
        << "analytic_hole_phip_norm,fd_hole_phip_norm,diff_hole_phip_norm,rel_hole_phip_diff\n";
    out << std::setprecision(17);
    for (const auto& row : rows) {
        out << row.block << ','
            << row.analyticNorm << ','
            << row.fdNorm << ','
            << row.diffNorm << ','
            << row.relDiff << ','
            << row.analyticPsiColumnNorm << ','
            << row.fdPsiColumnNorm << ','
            << row.diffPsiColumnNorm << ','
            << row.relPsiColumnDiff << ','
            << row.analyticPhinColumnNorm << ','
            << row.fdPhinColumnNorm << ','
            << row.diffPhinColumnNorm << ','
            << row.relPhinColumnDiff << ','
            << row.analyticPhipColumnNorm << ','
            << row.fdPhipColumnNorm << ','
            << row.diffPhipColumnNorm << ','
            << row.relPhipColumnDiff << ','
            << row.analyticElectronPhinNorm << ','
            << row.fdElectronPhinNorm << ','
            << row.diffElectronPhinNorm << ','
            << row.relElectronPhinDiff << ','
            << row.analyticElectronPhipNorm << ','
            << row.fdElectronPhipNorm << ','
            << row.diffElectronPhipNorm << ','
            << row.relElectronPhipDiff << ','
            << row.analyticHolePhinNorm << ','
            << row.fdHolePhinNorm << ','
            << row.diffHolePhinNorm << ','
            << row.relHolePhinDiff << ','
            << row.analyticHolePhipNorm << ','
            << row.fdHolePhipNorm << ','
            << row.diffHolePhipNorm << ','
            << row.relHolePhipDiff << '\n';
    }

    nlohmann::json result = {
        {"nodes", problem.mesh.numNodes()},
        {"blocks", rows.size()},
        {"finite_difference_mode", fdMode},
        {"output_csv", outputPath.string()},
    };
    for (const auto& row : rows) {
        if (!row.configurationFingerprint.empty()) {
            result["impact_configuration_fingerprint"] =
                row.configurationFingerprint;
            result["impact_active_branch_fingerprint"] =
                row.activeBranchFingerprint;
            break;
        }
    }
    return result;
}

} // namespace

int main(int argc, char** argv)
{
    std::string configFile;
    bool includeMeshReport = false;
    vela::RuntimeLogCliOverrides logOverrides;
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--config" && i + 1 < argc) {
                configFile = argv[++i];
            } else if (arg == "--mesh-report") {
                includeMeshReport = true;
            } else if (arg == "--log" && i + 1 < argc) {
                const std::string value = argv[++i];
                if (value == "auto") {
                    logOverrides.enabled.reset();
                    logOverrides.file.reset();
                } else if (value == "off") {
                    logOverrides.enabled = false;
                } else {
                    logOverrides.enabled = true;
                    logOverrides.file = value;
                }
            } else if (arg == "--log-profile" && i + 1 < argc) {
                logOverrides.profile = vela::runtimeLogProfileFromString(argv[++i]);
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

        RuntimeLogOverrideGuard logOverrideGuard(logOverrides);
        std::ifstream ifs(configFile);
        if (!ifs.is_open()) {
            std::cerr << "Cannot open config file: " << configFile << '\n';
            return 1;
        }

        nlohmann::json cfg;
        ifs >> cfg;
        const std::string type = cfg.value("simulation_type", cfg.contains("sweep") ? "dc_sweep" : "poisson");
        std::optional<vela::RuntimeLogSession> runnerLogSession;
        if (type != "dc_sweep" && type != "poisson")
            runnerLogSession = vela::RuntimeLogSession::fromConfig(cfg, configFile, type);

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
            if (result.releaseBVConfigAudit.has_value())
                status["avalanche_config"] =
                    releaseBVConfigAuditJson(*result.releaseBVConfigAudit);
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
            status["convergence_reason"] = result.result.convergenceReason;
            status["failure_reason"] = result.result.failureDiagnostics.failureReason;
            status["final_block_residuals"] = blockResidualsJson(result.result.finalBlockNorms);
            status["carrier_row_convergence"] =
                carrierRowConvergenceJson(result.result.finalCarrierRowConvergence);
            status["carrier_row_recovery"] =
                carrierRowRecoveryJson(result.result.carrierRowRecovery);
            if (includeMeshReport)
                status["mesh_report"] = meshReportJson(result.mesh.lastGeometryBuildReport());
        } else if (type == "newton_solve_from_state") {
            status.update(runNewtonSolveFromState(configFile, cfg));
        } else if (type == "newton_residual_probe") {
            status.update(runNewtonResidualProbe(configFile, cfg));
        } else if (type == "newton_step_probe") {
            status.update(runNewtonStepProbe(configFile, cfg));
        } else if (type == "newton_feedback_substitution_probe") {
            status.update(runNewtonFeedbackSubstitutionProbe(configFile, cfg));
        } else if (type == "newton_poisson_qfp_cross_block_probe") {
            status.update(runNewtonPoissonQfpCrossBlockProbe(configFile, cfg));
        } else if (type == "newton_jvp_probe") {
            status.update(runNewtonJvpProbe(configFile, cfg));
        } else if (type == "newton_block_step_probe") {
            status.update(runNewtonBlockStepProbe(configFile, cfg));
        } else if (type == "newton_regularized_carrier_step_probe") {
            status.update(runNewtonRegularizedCarrierStepProbe(configFile, cfg));
        } else if (type == "newton_carrier_row_probe") {
            status.update(runNewtonCarrierRowProbe(configFile, cfg));
        } else if (type == "newton_carrier_block_decomposition_probe") {
            status.update(runNewtonCarrierBlockDecompositionProbe(configFile, cfg));
        } else if (type == "newton_carrier_term_probe") {
            status.update(runNewtonCarrierTermProbe(configFile, cfg));
        } else if (type == "sg_edge_flux_probe") {
            status.update(runSgEdgeFluxProbe(configFile, cfg));
        } else if (type == "transport_edge_jacobian_probe") {
            status.update(runTransportEdgeJacobianProbe(configFile, cfg));
        } else if (type == "edge_mobility_probe") {
            status.update(runEdgeMobilityProbe(configFile, cfg));
        } else if (type == "newton_jacobian_block_probe") {
            status.update(runNewtonJacobianBlockProbe(configFile, cfg));
        } else {
            std::cerr << "Unknown simulation_type: " << type << '\n';
            return 2;
        }

        const bool converged = status.value("converged", false);
        if (runnerLogSession.has_value() && runnerLogSession->active())
            runnerLogSession->finish(converged);
        std::cout << status.dump() << '\n';
        return converged ? 0 : 1;
    } catch (const std::exception& ex) {
        std::cerr << "vela_example_runner failed: " << ex.what() << '\n';
        return 1;
    }
}
