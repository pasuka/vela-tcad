#include "vela/simulation/DCSweep.h"
#include "vela/io/CSVWriter.h"
#include "vela/io/MeshReader.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/post/ContactCurrent.h"
#include "vela/solver/NewtonSolver.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace vela {

namespace {

enum class SolverMethod {
    Gummel,
    Newton,
};

std::string formatReal(Real value)
{
    std::ostringstream oss;
    oss << std::setprecision(17) << value;
    return oss.str();
}

bool isFiniteSolution(const DDSolution& sol)
{
    auto finiteVector = [](const VectorXd& values) {
        for (int i = 0; i < values.size(); ++i) {
            if (!std::isfinite(values(i)))
                return false;
        }
        return true;
    };
    return finiteVector(sol.psi) && finiteVector(sol.phin) &&
           finiteVector(sol.phip) && finiteVector(sol.n) && finiteVector(sol.p);
}

std::string normalizedSolverMethod(const nlohmann::json& cfg)
{
    std::string method = "gummel";
    if (cfg.contains("solver")) {
        const auto& solver = cfg.at("solver");
        if (solver.contains("method"))
            method = solver.at("method").get<std::string>();
        else if (solver.contains("type"))
            method = solver.at("type").get<std::string>();
    }

    std::transform(method.begin(), method.end(), method.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return method;
}

SolverMethod solverMethodFromJson(const nlohmann::json& cfg)
{
    const std::string method = normalizedSolverMethod(cfg);
    if (method == "gummel")
        return SolverMethod::Gummel;
    if (method == "newton")
        return SolverMethod::Newton;
    throw std::invalid_argument(
        "DCSweep: solver.method must be 'gummel' or 'newton'.");
}

DCSweepConfig dcSweepConfigFromJson(const nlohmann::json& cfg,
                                    const std::filesystem::path& cfgDir)
{
    const auto& j = cfg.at("sweep");
    DCSweepConfig sweep;
    sweep.contact = j.at("contact").get<std::string>();
    sweep.start = j.at("start").get<Real>();
    sweep.stop = j.at("stop").get<Real>();
    sweep.step = j.at("step").get<Real>();
    const Real nominalStep = std::abs(sweep.step);
    sweep.shrinkFactor = j.value("shrink_factor", sweep.shrinkFactor);
    sweep.growthFactor = j.value("growth_factor", sweep.growthFactor);
    sweep.maxRetries = j.value("max_retries", sweep.maxRetries);
    sweep.minStep = j.value("min_step", nominalStep * std::pow(sweep.shrinkFactor, sweep.maxRetries));
    sweep.maxStep = j.value("max_step", nominalStep);
    sweep.stopOnFailure = j.value("stop_on_failure", sweep.stopOnFailure);
    sweep.currentContact = j.value("current_contact", sweep.contact);
    sweep.writeVtk = j.value("write_vtk", cfg.value("write_vtk", false));
    sweep.csvFile = j.value("csv_file", cfg.value("output_csv", sweep.csvFile));
    sweep.vtkPrefix = j.value("vtk_prefix", cfg.value("output_vtk_prefix", std::string("dc_sweep")));

    auto resolve = [&](std::string path) {
        std::filesystem::path fp(path);
        if (fp.is_relative())
            fp = cfgDir / fp;
        return fp.string();
    };
    sweep.csvFile = resolve(sweep.csvFile);
    sweep.vtkPrefix = resolve(sweep.vtkPrefix);

    if (sweep.step == 0.0)
        throw std::invalid_argument("DCSweep: sweep.step must be non-zero.");
    if ((sweep.stop - sweep.start) * sweep.step < 0.0)
        throw std::invalid_argument("DCSweep: sweep.step sign must move start toward stop.");
    if (sweep.minStep <= 0.0)
        throw std::invalid_argument("DCSweep: sweep.min_step must be positive.");
    if (sweep.maxStep <= 0.0)
        throw std::invalid_argument("DCSweep: sweep.max_step must be positive.");
    if (sweep.minStep > sweep.maxStep)
        throw std::invalid_argument("DCSweep: sweep.min_step must not exceed sweep.max_step.");
    if (sweep.growthFactor < 1.0)
        throw std::invalid_argument("DCSweep: sweep.growth_factor must be at least 1.");
    if (sweep.shrinkFactor <= 0.0 || sweep.shrinkFactor >= 1.0)
        throw std::invalid_argument("DCSweep: sweep.shrink_factor must be greater than 0 and less than 1.");
    if (sweep.maxRetries < 0)
        throw std::invalid_argument("DCSweep: sweep.max_retries must be non-negative.");
    return sweep;
}

DopingModel dopingFromJson(const DeviceMesh& mesh, const nlohmann::json& cfg)
{
    std::vector<RegionDopingSpec> specs;
    for (const auto& entry : cfg.at("doping")) {
        RegionDopingSpec spec;
        spec.region = entry.at("region").get<std::string>();
        spec.donors = entry.at("donors").get<Real>();
        spec.acceptors = entry.at("acceptors").get<Real>();
        specs.push_back(std::move(spec));
    }
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

std::unordered_map<std::string, Real> contactBiasesFromJson(const nlohmann::json& cfg)
{
    std::unordered_map<std::string, Real> biases;
    for (const auto& contact : cfg.at("contacts")) {
        biases[contact.at("name").get<std::string>()] =
            contact.at("bias").get<Real>();
    }
    return biases;
}

std::string vtkFilename(const std::string& prefix, int index, Real voltage)
{
    std::ostringstream oss;
    oss << prefix << "_" << std::setw(4) << std::setfill('0') << index
        << "_" << std::setprecision(6) << std::defaultfloat << voltage << "V.vtk";
    return oss.str();
}

} // namespace

std::vector<DCSweepPoint> DCSweep::run(const std::string& configFile) const
{
    return runWithResult(configFile).points;
}

DCSweepResult DCSweep::runWithResult(const std::string& configFile) const
{
    std::ifstream ifs(configFile);
    if (!ifs.is_open())
        throw std::runtime_error("DCSweep: cannot open config file: " + configFile);

    nlohmann::json cfg;
    ifs >> cfg;

    const std::filesystem::path cfgDir = std::filesystem::path(configFile).parent_path();
    auto resolve = [&](const std::string& path) {
        std::filesystem::path fp(path);
        if (fp.is_relative())
            fp = cfgDir / fp;
        return fp.string();
    };

    JsonMeshReader reader;
    DeviceMesh mesh = reader.read(resolve(cfg.at("mesh_file").get<std::string>()));
    MaterialDatabase matdb;
    if (cfg.contains("materials_file"))
        matdb.loadJson(resolve(cfg.at("materials_file").get<std::string>()));
    DopingModel doping = dopingFromJson(mesh, cfg);
    std::unordered_map<std::string, Real> baseBiases = contactBiasesFromJson(cfg);
    DCSweepConfig sweep = dcSweepConfigFromJson(cfg, cfgDir);
    const nlohmann::json solverCfg = cfg.value("solver", nlohmann::json::object());
    const SolverMethod solverMethod = solverMethodFromJson(cfg);
    GummelConfig gummel = gummelConfigFromJson(solverCfg);
    NewtonConfig newton = newtonConfigFromJson(solverCfg);
    MobilityModelConfig mobilityConfig = solverMethod == SolverMethod::Newton
        ? mobilityModelConfig(newton.mobility)
        : mobilityModelConfig(gummel.mobility);
    ContactCurrent contactCurrent(mesh, matdb, doping, mobilityConfig);

    CSVWriter csv(sweep.csvFile);
    csv.writeHeader({"voltage", "electron_current", "hole_current",
                     "total_current", "converged", "iterations",
                     "attempted_step", "accepted_step", "retry_count"});

    std::vector<DCSweepPoint> points;
    DDSolution previousSolution;
    Real previousVoltage = sweep.start;
    int vtkIndex = 0;

    auto solvePoint = [&](Real voltage, const DDSolution* initial) -> std::pair<bool, DDSolution> {
        auto biases = baseBiases;
        biases[sweep.contact] = voltage;
        try {
            if (solverMethod == SolverMethod::Newton) {
                NewtonResult result = initial != nullptr
                    ? runNewton(mesh, matdb, doping, biases, *initial, newton)
                    : runNewton(mesh, matdb, doping, biases, newton);
                DDSolution sol = std::move(result.solution);
                return {result.converged && isFiniteSolution(sol), std::move(sol)};
            }

            DDSolution sol = initial != nullptr
                ? runGummel(mesh, matdb, doping, biases, gummel, *initial)
                : runGummel(mesh, matdb, doping, biases, gummel);
            return {sol.converged && isFiniteSolution(sol), std::move(sol)};
        } catch (const std::exception& ex) {
            std::throw_with_nested(std::runtime_error(
                "DCSweep: solver threw at voltage " + formatReal(voltage) +
                " V: " + ex.what()));
        } catch (...) {
            std::throw_with_nested(std::runtime_error(
                "DCSweep: solver threw an unknown exception at voltage " +
                formatReal(voltage) + " V."));
        }
    };

    auto recordPoint = [&](Real voltage, const DDSolution& sol, bool converged,
                           Real attemptedStep, Real acceptedStep, int retryCount) {
        ContactCurrentResult current{};
        if (converged)
            current = contactCurrent.compute(sol, sweep.currentContact);

        DCSweepPoint point;
        point.voltage = voltage;
        point.electronCurrent = current.electronCurrent;
        point.holeCurrent = current.holeCurrent;
        point.totalCurrent = current.totalCurrent;
        point.converged = converged;
        point.iterations = sol.iters;
        point.attemptedStep = attemptedStep;
        point.acceptedStep = acceptedStep;
        point.retryCount = retryCount;
        points.push_back(point);

        csv.writeRow({formatReal(point.voltage),
                      formatReal(point.electronCurrent),
                      formatReal(point.holeCurrent),
                      formatReal(point.totalCurrent),
                      point.converged ? "1" : "0",
                      std::to_string(point.iterations),
                      formatReal(point.attemptedStep),
                      formatReal(point.acceptedStep),
                      std::to_string(point.retryCount)});

        if (converged && sweep.writeVtk)
            writeDDSolutionVTK(vtkFilename(sweep.vtkPrefix, vtkIndex++, voltage), mesh, doping, sol);
    };

    auto [startOk, startSol] = solvePoint(sweep.start, nullptr);
    recordPoint(sweep.start, startSol, startOk, 0.0, 0.0, 0);
    if (!startOk)
        return DCSweepResult{std::move(mesh), std::move(points)};
    previousSolution = std::move(startSol);

    const Real direction = (sweep.step > 0.0) ? 1.0 : -1.0;
    const Real tolerance = 1.0e-12;
    Real adaptiveStep = std::min(std::abs(sweep.step), sweep.maxStep);

    auto limitedTarget = [&](Real target, Real stepMagnitude) {
        const Real remaining = direction * (target - previousVoltage);
        const Real limited = previousVoltage + direction * std::min(stepMagnitude, remaining);
        return limited;
    };

    auto advanceToward = [&](Real target) -> bool {
        int retryCount = 0;
        Real trialStep = std::min(adaptiveStep, sweep.maxStep);
        DDSolution lastSol;
        Real lastAttempted = 0.0;
        Real lastCandidate = previousVoltage;

        while (true) {
            const Real remaining = direction * (target - previousVoltage);
            if (remaining <= tolerance)
                return true;

            const Real stepMagnitude = std::min(trialStep, remaining);
            const Real candidate = limitedTarget(target, stepMagnitude);
            const Real attemptedStep = candidate - previousVoltage;
            auto [ok, sol] = solvePoint(candidate, &previousSolution);
            lastSol = sol;
            lastAttempted = attemptedStep;
            lastCandidate = candidate;

            if (ok) {
                recordPoint(candidate, sol, true, attemptedStep, attemptedStep, retryCount);
                previousSolution = std::move(sol);
                previousVoltage = candidate;
                adaptiveStep = std::min(sweep.maxStep, stepMagnitude * sweep.growthFactor);
                return true;
            }

            if (retryCount >= sweep.maxRetries)
                break;

            const Real shrunken = stepMagnitude * sweep.shrinkFactor;
            if (shrunken < sweep.minStep - std::numeric_limits<Real>::epsilon())
                break;

            trialStep = shrunken;
            ++retryCount;
        }

        recordPoint(lastCandidate, lastSol, false, lastAttempted, 0.0, retryCount);
        adaptiveStep = std::max(sweep.minStep, std::min(sweep.maxStep, std::abs(lastAttempted) * sweep.shrinkFactor));
        return false;
    };

    bool blockedByFailedStep = false;
    Real nominalTarget = sweep.start + sweep.step;
    while (!blockedByFailedStep && direction * (nominalTarget - sweep.stop) <= tolerance) {
        while (direction * (previousVoltage - nominalTarget) < -tolerance) {
            if (!advanceToward(nominalTarget)) {
                if (sweep.stopOnFailure)
                    return DCSweepResult{std::move(mesh), std::move(points)};
                blockedByFailedStep = true;
                break;
            }
        }
        nominalTarget += sweep.step;
    }

    while (!blockedByFailedStep && direction * (previousVoltage - sweep.stop) < -tolerance) {
        if (!advanceToward(sweep.stop)) {
            if (sweep.stopOnFailure)
                return DCSweepResult{std::move(mesh), std::move(points)};
            break;
        }
    }

    return DCSweepResult{std::move(mesh), std::move(points)};
}

} // namespace vela
