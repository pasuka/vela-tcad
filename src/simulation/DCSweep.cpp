#include "vela/simulation/DCSweep.h"
#include "vela/io/CSVWriter.h"
#include "vela/io/MeshReader.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/post/ContactCurrent.h"
#include <nlohmann/json.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace vela {

namespace {

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

DCSweepConfig dcSweepConfigFromJson(const nlohmann::json& cfg,
                                    const std::filesystem::path& cfgDir)
{
    const auto& j = cfg.at("sweep");
    DCSweepConfig sweep;
    sweep.contact = j.at("contact").get<std::string>();
    sweep.start = j.at("start").get<Real>();
    sweep.stop = j.at("stop").get<Real>();
    sweep.step = j.at("step").get<Real>();
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
    DopingModel doping = dopingFromJson(mesh, cfg);
    std::unordered_map<std::string, Real> baseBiases = contactBiasesFromJson(cfg);
    DCSweepConfig sweep = dcSweepConfigFromJson(cfg, cfgDir);
    GummelConfig gummel = cfg.contains("solver") ? gummelConfigFromJson(cfg.at("solver")) : GummelConfig{};
    MobilityModelConfig mobilityConfig = mobilityModelConfig(gummel.mobility);
    ContactCurrent contactCurrent(mesh, matdb, doping, mobilityConfig);

    CSVWriter csv(sweep.csvFile);
    csv.writeHeader({"voltage", "electron_current", "hole_current",
                     "total_current", "converged", "iterations"});

    std::vector<DCSweepPoint> points;
    DDSolution previousSolution;
    Real previousVoltage = sweep.start;
    int vtkIndex = 0;

    auto solvePoint = [&](Real voltage, const DDSolution* initial) -> std::pair<bool, DDSolution> {
        auto biases = baseBiases;
        biases[sweep.contact] = voltage;
        try {
            DDSolution sol = initial != nullptr
                ? runGummel(mesh, matdb, doping, biases, gummel, *initial)
                : runGummel(mesh, matdb, doping, biases, gummel);
            return {sol.converged && isFiniteSolution(sol), std::move(sol)};
        } catch (const std::exception& ex) {
            throw std::runtime_error(
                "DCSweep: solver threw at voltage " + formatReal(voltage) +
                " V: " + ex.what());
        } catch (...) {
            throw std::runtime_error(
                "DCSweep: solver threw an unknown exception at voltage " +
                formatReal(voltage) + " V.");
        }
    };

    auto recordPoint = [&](Real voltage, const DDSolution& sol, bool converged) {
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
        points.push_back(point);

        csv.writeRow({formatReal(point.voltage),
                      formatReal(point.electronCurrent),
                      formatReal(point.holeCurrent),
                      formatReal(point.totalCurrent),
                      point.converged ? "1" : "0",
                      std::to_string(point.iterations)});

        if (converged && sweep.writeVtk)
            writeDDSolutionVTK(vtkFilename(sweep.vtkPrefix, vtkIndex++, voltage), mesh, doping, sol);
    };

    auto [startOk, startSol] = solvePoint(sweep.start, nullptr);
    recordPoint(sweep.start, startSol, startOk);
    if (!startOk)
        return points;
    previousSolution = std::move(startSol);
    auto advanceToward = [&](Real target) -> bool {
        Real candidate = target;
        int depth = 0;

        while (true) {
            auto [ok, sol] = solvePoint(candidate, &previousSolution);
            if (ok) {
                recordPoint(candidate, sol, true);
                previousSolution = std::move(sol);
                previousVoltage = candidate;
                return true;
            }

            if (depth >= 5) {
                recordPoint(candidate, sol, false);
                return false;
            }

            candidate = 0.5 * (previousVoltage + candidate);
            ++depth;
        }
    };

    const Real direction = (sweep.step > 0.0) ? 1.0 : -1.0;
    Real nominalTarget = sweep.start + sweep.step;
    while (direction * (nominalTarget - sweep.stop) <= 1.0e-12) {
        while (direction * (previousVoltage - nominalTarget) < -1.0e-12) {
            if (!advanceToward(nominalTarget))
                return points;
        }
        nominalTarget += sweep.step;
    }

    while (direction * (previousVoltage - sweep.stop) < -1.0e-12) {
        if (!advanceToward(sweep.stop))
            return points;
    }

    return points;
}

} // namespace vela
