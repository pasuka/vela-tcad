#include "vela/simulation/DCSweep.h"
#include "vela/boundary/BoundaryCondition.h"
#include "vela/simulation/DCSweepStepControl.h"
#include "vela/simulation/ConfigParsing.h"
#include "vela/io/CSVWriter.h"
#include "vela/io/MeshReader.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/post/ContactCurrent.h"
#include "vela/post/ElectricFieldDiagnostics.h"
#include "vela/post/TerminalCharge.h"
#include "vela/post/StoredCharge.h"
#include "vela/solver/NewtonSolver.h"
#include "vela/solver/SolutionValidation.h"
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
#include <unordered_set>

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


std::string sanitizedColumnToken(std::string value)
{
    for (char& ch : value) {
        const unsigned char uch = static_cast<unsigned char>(ch);
        if (!std::isalnum(uch))
            ch = '_';
        else
            ch = static_cast<char>(std::tolower(uch));
    }
    if (value.empty())
        value = "terminal";
    return value;
}

std::string terminalChargeName(const TerminalChargeConfig& cfg, std::size_t index)
{
    if (!cfg.name.empty())
        return sanitizedColumnToken(cfg.name);
    if (!cfg.contact.empty())
        return sanitizedColumnToken(cfg.contact);
    return "terminal" + std::to_string(index + 1);
}

std::string capacitanceMnemonic(const std::string& sweepContact, const std::string& terminalName)
{
    return "C" + sanitizedColumnToken(sweepContact) + "_" + sanitizedColumnToken(terminalName);
}

TerminalChargeConfig terminalChargeConfigFromJson(const nlohmann::json& chargeCfg,
                                                  const DCSweepConfig& sweep,
                                                  std::size_t index)
{
    TerminalChargeConfig config;
    config.name = chargeCfg.value("name", std::string{});
    config.contact = chargeCfg.value("contact", sweep.chargeContact.empty() ? sweep.contact : sweep.chargeContact);
    config.regions = chargeCfg.value("regions", sweep.chargeRegions);
    config.contactRadius = chargeCfg.value("contact_radius", sweep.chargeContactRadius);
    config.includeMobileCharge = chargeCfg.value("include_mobile_charge", config.includeMobileCharge);
    config.includeIonizedDopants = chargeCfg.value("include_ionized_dopants", config.includeIonizedDopants);
    config.perMeter = chargeCfg.value("per_meter", sweep.chargePerMeter);
    config.depth_m = chargeCfg.value("depth_m", sweep.chargeDepth_m);
    if (config.name.empty())
        config.name = terminalChargeName(config, index);
    if (!config.perMeter && config.depth_m <= 0.0)
        throw std::invalid_argument("DCSweep: sweep terminal charge depth_m must be positive.");
    return config;
}

DDSolutionValidationOptions validationOptionsFromJson(const nlohmann::json& cfg)
{
    DDSolutionValidationOptions options;
    const auto validation = cfg.value("validation", nlohmann::json::object());
    options.carrierFloor = validation.value("carrier_floor", options.carrierFloor);
    options.enforceMinimumCarrierDensity =
        validation.value("enforce_minimum_carrier_density", options.enforceMinimumCarrierDensity);
    options.minimumCarrierDensity =
        validation.value("minimum_carrier_density", options.minimumCarrierDensity);
    options.checkContactQuasiFermiBias =
        validation.value("check_contact_quasi_fermi_bias", options.checkContactQuasiFermiBias);
    options.contactPotentialAbsTolerance =
        validation.value("contact_potential_abs_tolerance", options.contactPotentialAbsTolerance);
    options.contactPotentialRelTolerance =
        validation.value("contact_potential_rel_tolerance", options.contactPotentialRelTolerance);
    if (options.carrierFloor < 0.0)
        throw std::invalid_argument("DCSweep: validation.carrier_floor must be non-negative.");
    if (options.enforceMinimumCarrierDensity && options.minimumCarrierDensity < 0.0)
        throw std::invalid_argument(
            "DCSweep: validation.minimum_carrier_density must be non-negative.");
    if (options.contactPotentialAbsTolerance < 0.0 ||
        options.contactPotentialRelTolerance < 0.0) {
        throw std::invalid_argument(
            "DCSweep: validation contact tolerances must be non-negative.");
    }
    return options;
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
        "DCSweep: solver.method/type must be 'gummel' or 'newton'.");
}

DCSweepConfig dcSweepConfigFromJson(const nlohmann::json& cfg,
                                    const std::filesystem::path& cfgDir)
{
    const auto& j = cfg.at("sweep");
    DCSweepConfig sweep;
    sweep.mode = curveSweepModeFromString(j.value("mode", std::string("iv")));
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

    const auto chargeCfg = j.value("terminal_charge", nlohmann::json::object());
    sweep.chargeContact = chargeCfg.value("contact", j.value("charge_contact", sweep.contact));
    sweep.chargeRegions = chargeCfg.value("regions", j.value("charge_regions", std::vector<std::string>{}));
    sweep.chargeContactRadius = chargeCfg.value("contact_radius", j.value("charge_contact_radius", 0.0));
    sweep.chargePerMeter = chargeCfg.value("per_meter", j.value("charge_per_meter", true));
    sweep.chargeDepth_m = chargeCfg.value("depth_m", j.value("charge_depth_m", 1.0));

    sweep.storedChargeEnabled = j.contains("stored_charge");
    const auto storedCfg = j.value("stored_charge", nlohmann::json::object());
    sweep.storedCharge.regions = storedCfg.value("regions", std::vector<std::string>{});
    sweep.storedCharge.perMeter = storedCfg.value("per_meter", true);
    sweep.storedCharge.depth_m = storedCfg.value("depth_m", 1.0);

    if (j.contains("terminal_charges")) {
        const auto& charges = j.at("terminal_charges");
        if (!charges.is_array())
            throw std::invalid_argument("DCSweep: sweep.terminal_charges must be an array.");
        if (charges.empty())
            throw std::invalid_argument("DCSweep: sweep.terminal_charges must not be empty.");
        std::unordered_set<std::string> names;
        for (std::size_t i = 0; i < charges.size(); ++i) {
            TerminalChargeConfig config = terminalChargeConfigFromJson(charges.at(i), sweep, i);
            const std::string name = terminalChargeName(config, i);
            if (!names.insert(name).second)
                throw std::invalid_argument("DCSweep: duplicate terminal_charges name '" + name + "'.");
            config.name = name;
            sweep.terminalCharges.push_back(std::move(config));
        }
    } else {
        sweep.terminalCharges.push_back(terminalChargeConfigFromJson(chargeCfg, sweep, 0));
    }

    const auto bvCfg = j.value("breakdown", nlohmann::json::object());
    sweep.breakdown.maxElectricField_V_per_m = bvCfg.value("max_electric_field_V_per_m", j.value("breakdown_max_electric_field_V_per_m", 0.0));
    sweep.breakdown.currentJumpRatio = bvCfg.value("current_jump_ratio", j.value("breakdown_current_jump_ratio", 0.0));
    sweep.breakdown.nonConvergenceBreakdown = bvCfg.value("non_convergence", j.value("breakdown_on_non_convergence", true));

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
    if (sweep.mode == CurveSweepMode::BVReverse && sweep.start * sweep.stop < 0.0)
        throw std::invalid_argument("DCSweep: bv_reverse sweeps must stay on one reverse-bias polarity side.");
    if (!sweep.chargePerMeter && sweep.chargeDepth_m <= 0.0)
        throw std::invalid_argument("DCSweep: sweep terminal charge depth_m must be positive.");
    if (sweep.storedChargeEnabled && !sweep.storedCharge.perMeter && sweep.storedCharge.depth_m <= 0.0)
        throw std::invalid_argument("DCSweep: sweep stored_charge depth_m must be positive.");
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

struct ContactConfig {
    std::unordered_map<std::string, Real> biases;
    ContactSpecsMap specs;
};

ContactConfig contactConfigFromJson(const nlohmann::json& cfg)
{
    // Route legacy ``contacts[]`` parsing through the unified boundary parser
    // so the optional ``type`` field is recognised and normalised.  For the
    // DD/Gummel paths the bias map preserves the historical semantics:
    // ohmic contacts apply the raw bias, while Dirichlet/MetalGate use the
    // same effective potential expression as the Poisson driver.  Schottky
    // contacts route through a barrier-aware ``ContactSpecsMap`` and use the
    // raw bias so the swept value reaches the metal Fermi level.  Floating
    // contacts have no DD physics yet and are rejected so a misconfigured
    // deck fails loudly instead of silently downgrading.
    ContactConfig out;
    for (const auto& spec : parseContactBoundarySpecs(cfg)) {
        switch (spec.type) {
            case ContactType::Ohmic:
                out.biases[spec.name] = spec.bias;
                break;
            case ContactType::Dirichlet:
            case ContactType::MetalGate:
                out.biases[spec.name] = effectivePoissonDirichletPotential(spec);
                break;
            case ContactType::Schottky:
                out.biases[spec.name] = spec.bias;
                out.specs[spec.name] = spec;
                break;
            case ContactType::Floating:
                throw std::runtime_error(
                    "DCSweep: contact '" + spec.name + "' has type '" +
                    toString(spec.type) + "' which is not yet implemented "
                    "for drift-diffusion sweeps. Use ohmic or schottky for now.");
        }
    }
    return out;
}


std::string vtkFilename(const std::string& prefix, int index, Real voltage)
{
    std::ostringstream oss;
    oss << prefix << "_" << std::setw(4) << std::setfill('0') << index
        << "_" << std::setprecision(6) << std::defaultfloat << voltage << "V.vtk";
    return oss.str();
}


std::string stepDiagnostics(const DCSweepPoint& point)
{
    return "attempted_step=" + formatReal(point.attemptedStep) +
           ";accepted_step=" + formatReal(point.acceptedStep) +
           ";retry_count=" + std::to_string(point.retryCount);
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
    std::vector<RegionFixedChargeSpec> fixedChargeSpecs = parseRegionFixedChargeSpecs(cfg);
    std::vector<InterfaceSheetChargeSpec> sheetChargeSpecs = parseInterfaceSheetChargeSpecs(cfg);
    ContactConfig contactConfig = contactConfigFromJson(cfg);
    std::unordered_map<std::string, Real>& baseBiases = contactConfig.biases;
    ContactSpecsMap& contactSpecs = contactConfig.specs;
    DCSweepConfig sweep = dcSweepConfigFromJson(cfg, cfgDir);

    const nlohmann::json solverCfg = cfg.value("solver", nlohmann::json::object());
    const SolverMethod solverMethod = solverMethodFromJson(cfg);
    const DDSolutionValidationOptions validationOptions = validationOptionsFromJson(cfg);
    GummelConfig gummel;
    NewtonConfig newton;
    MobilityModelConfig mobilityConfig;
    if (solverMethod == SolverMethod::Newton) {
        newton = newtonConfigFromJson(solverCfg);
        mobilityConfig = newton.mobility;
    } else {
        gummel = gummelConfigFromJson(solverCfg);
        mobilityConfig = gummel.mobility;
    }
    const Real temperature_K = (solverMethod == SolverMethod::Newton)
        ? newton.temperature_K
        : gummel.temperature_K;
    ContactCurrent contactCurrent(mesh, matdb, doping, mobilityConfig, temperature_K);
    TerminalCharge terminalCharge(mesh, doping);
    StoredCharge storedCharge(mesh);
    const bool hasMultiTerminalCharges = cfg.at("sweep").contains("terminal_charges");
    const TerminalChargeConfig& legacyChargeConfig = sweep.terminalCharges.front();

    CSVWriter csv(sweep.csvFile);
    std::vector<std::string> header = {"mode", "bias_contact", "bias_V",
        "current_contact", "current_electron", "current_hole", "current_total",
        "converged", "iterations", "step_diagnostics", "validation_diagnostics"};
    std::vector<std::pair<std::string, std::string>> chargeColumns;
    std::vector<std::pair<std::string, std::string>> capacitanceColumns;
    if (sweep.mode == CurveSweepMode::IV && sweep.storedChargeEnabled)
        header.push_back(sweep.storedCharge.perMeter ? "stored_charge_C_per_m" : "stored_charge_C");
    if (sweep.mode == CurveSweepMode::CVQuasistatic) {
        header.push_back(legacyChargeConfig.perMeter ? "charge_C_per_m" : "charge_C");
        header.push_back(legacyChargeConfig.perMeter ? "capacitance_F_per_m" : "capacitance_F");
        if (hasMultiTerminalCharges) {
            for (std::size_t i = 0; i < sweep.terminalCharges.size(); ++i) {
                const TerminalChargeConfig& config = sweep.terminalCharges[i];
                const std::string name = terminalChargeName(config, i);
                const std::string chargeColumn = "charge_" + name + (config.perMeter ? "_C_per_m" : "_C");
                const std::string capColumn = "capacitance_" + capacitanceMnemonic(sweep.contact, name) +
                    (config.perMeter ? "_F_per_m" : "_F");
                chargeColumns.emplace_back(name, chargeColumn);
                capacitanceColumns.emplace_back(name, capColumn);
                header.push_back(chargeColumn);
                header.push_back(capColumn);
            }
        }
    }
    if (sweep.mode == CurveSweepMode::BVReverse) {
        header.push_back("max_electric_field_V_per_m");
        header.push_back("current_jump_ratio");
        header.push_back("breakdown_detected");
        header.push_back("breakdown_voltage");
        header.push_back("criterion");
        header.push_back("last_stable_bias");
        header.push_back("failed_bias");
        header.push_back("failure_reason");
    }
    csv.writeHeader(header);

    std::vector<DCSweepPoint> points;
    DDSolution previousSolution;
    int vtkIndex = 0;

    struct SolvePointAttempt {
        bool ok = false;
        DDSolution solution;
        std::string failureReason;
        std::string validationDiagnostics;
    };

    if (solverMethod == SolverMethod::Newton && !contactSpecs.empty()) {
        // The Newton coupled-DD path does not yet construct Schottky-aware
        // boundary conditions.  Fail loudly so a user mis-configuring the
        // deck does not silently get the legacy Ohmic interpretation.
        std::string firstSchottky;
        for (const auto& [name, _] : contactSpecs) {
            firstSchottky = name;
            break;
        }
        throw std::runtime_error(
            "DCSweep: solver.method='newton' does not yet support Schottky contacts "
            "(contact '" + firstSchottky + "'). Use solver.method='gummel' for the "
            "Schottky prototype, or switch the contact to type='ohmic'.");
    }

    auto solvePoint = [&](Real voltage, const DDSolution* initial) -> SolvePointAttempt {
        auto biases = baseBiases;
        biases[sweep.contact] = voltage;
        try {
            bool solverConverged = false;
            DDSolution sol;
            if (solverMethod == SolverMethod::Newton) {
                NewtonResult result = initial != nullptr
                    ? runNewton(mesh, matdb, doping, biases, *initial, newton, fixedChargeSpecs, sheetChargeSpecs)
                    : runNewton(mesh, matdb, doping, biases, newton, fixedChargeSpecs, sheetChargeSpecs);
                solverConverged = result.converged;
                sol = std::move(result.solution);
            } else {
                sol = initial != nullptr
                    ? runGummel(mesh, matdb, doping, biases, contactSpecs, gummel, *initial, fixedChargeSpecs, sheetChargeSpecs)
                    : runGummel(mesh, matdb, doping, biases, contactSpecs, gummel, fixedChargeSpecs, sheetChargeSpecs);
                solverConverged = sol.converged;
            }


            const DDSolutionValidationResult validation =
                validateDDSolution(sol, mesh, biases, validationOptions);
            SolvePointAttempt attempt;
            attempt.ok = solverConverged && validation.valid;
            attempt.solution = std::move(sol);
            attempt.validationDiagnostics = validation.diagnosticsString();
            if (!solverConverged)
                attempt.failureReason = "non_convergence";
            else if (!validation.valid)
                attempt.failureReason = "validation_failed";
            return attempt;
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

    auto lastConvergedPoint = [&]() -> const DCSweepPoint* {
        for (auto it = points.rbegin(); it != points.rend(); ++it) {
            if (it->converged)
                return &(*it);
        }
        return nullptr;
    };

    auto recordPoint = [&](Real voltage, const DDSolution& sol, bool converged,
                           Real attemptedStep, Real acceptedStep, int retryCount,
                           const std::string& failureReason = std::string(),
                           const std::string& validationDiagnostics = std::string()) {
        ContactCurrentResult current{};
        if (converged)
            current = contactCurrent.compute(sol, sweep.currentContact);

        DCSweepPoint point;
        point.voltage = voltage;
        point.bias = voltage;
        point.outputCsv = sweep.csvFile;
        point.electronCurrent = current.electronCurrent;
        point.holeCurrent = current.holeCurrent;
        point.totalCurrent = current.totalCurrent;
        point.converged = converged;
        point.iterations = sol.iters;
        point.attemptedStep = attemptedStep;
        point.acceptedStep = acceptedStep;
        point.retryCount = retryCount;
        point.validationDiagnostics = validationDiagnostics;
        if (converged && sweep.mode == CurveSweepMode::IV && sweep.storedChargeEnabled) {
            point.extraFields.emplace_back(
                sweep.storedCharge.perMeter ? "stored_charge_C_per_m" : "stored_charge_C",
                storedCharge.compute(sol, sweep.storedCharge).charge);
        }
        if (!converged && sweep.mode != CurveSweepMode::BVReverse)
            point.failureReason = failureReason;
        if (converged && sweep.mode == CurveSweepMode::CVQuasistatic) {
            for (std::size_t i = 0; i < sweep.terminalCharges.size(); ++i) {
                const TerminalChargeConfig& config = sweep.terminalCharges[i];
                const std::string name = terminalChargeName(config, i);
                const Real charge = terminalCharge.compute(sol, config).charge;
                point.terminalChargeValues.emplace_back(name, charge);
                Real capacitance = 0.0;
                if (!points.empty()) {
                    const DCSweepPoint& prev = points.back();
                    const Real dV = point.bias - prev.bias;
                    auto prevIt = std::find_if(prev.terminalChargeValues.begin(), prev.terminalChargeValues.end(),
                        [&](const auto& entry) { return entry.first == name; });
                    if (dV != 0.0 && prevIt != prev.terminalChargeValues.end())
                        capacitance = (charge - prevIt->second) / dV;
                }
                point.terminalCapacitanceValues.emplace_back(name, capacitance);
            }
            if (!point.terminalChargeValues.empty())
                point.terminalCharge = point.terminalChargeValues.front().second;
            if (!point.terminalCapacitanceValues.empty())
                point.capacitance = point.terminalCapacitanceValues.front().second;
            if (hasMultiTerminalCharges) {
                for (std::size_t i = 0; i < point.terminalChargeValues.size(); ++i) {
                    point.extraFields.emplace_back(chargeColumns.at(i).second, point.terminalChargeValues.at(i).second);
                    point.extraFields.emplace_back(capacitanceColumns.at(i).second, point.terminalCapacitanceValues.at(i).second);
                }
            }
        }
        if (converged && sweep.mode == CurveSweepMode::BVReverse) {
            point.maxElectricField = maxEdgeElectricFieldMagnitude(mesh, sol.psi);
            if (!points.empty()) {
                const Real previous = std::abs(points.back().totalCurrent);
                const Real currentAbs = std::abs(point.totalCurrent);
                if (previous > 0.0)
                    point.currentJumpRatio = currentAbs / previous;
                else if (currentAbs > 0.0)
                    point.currentJumpRatio = std::numeric_limits<Real>::infinity();
            }
            if (sweep.breakdown.maxElectricField_V_per_m > 0.0 &&
                point.maxElectricField >= sweep.breakdown.maxElectricField_V_per_m) {
                point.breakdownDetected = true;
                point.breakdownVoltage = point.bias;
                point.breakdownCriterion = "max_electric_field";
            } else if (sweep.breakdown.currentJumpRatio > 0.0 &&
                       point.currentJumpRatio >= sweep.breakdown.currentJumpRatio) {
                point.breakdownDetected = true;
                point.breakdownVoltage = point.bias;
                point.breakdownCriterion = "current_jump";
            }
        } else if (!converged && sweep.mode == CurveSweepMode::BVReverse &&
                   sweep.breakdown.nonConvergenceBreakdown) {
            point.failed = true;
            point.failedBias = voltage;
            point.failureReason = failureReason.empty() ? "non_convergence" : failureReason;
            if (const DCSweepPoint* stable = lastConvergedPoint()) {
                point.lastStableBias = stable->bias;
                point.breakdownDetected = true;
                point.breakdownVoltage = point.lastStableBias;
                point.breakdownCriterion = "last_stable_before_nonconvergence";
            }
        }
        std::vector<std::string> row = {
            toString(sweep.mode),
            sweep.contact,
            formatReal(point.bias),
            sweep.currentContact,
            formatReal(point.electronCurrent),
            formatReal(point.holeCurrent),
            formatReal(point.totalCurrent),
            point.converged ? "1" : "0",
            std::to_string(point.iterations),
            stepDiagnostics(point),
            point.validationDiagnostics};
        if (sweep.mode == CurveSweepMode::IV && sweep.storedChargeEnabled) {
            const char* storedColumn = sweep.storedCharge.perMeter ? "stored_charge_C_per_m" : "stored_charge_C";
            Real storedValue = 0.0;
            for (const auto& [name, value] : point.extraFields) {
                if (name == storedColumn) {
                    storedValue = value;
                    break;
                }
            }
            row.push_back(formatReal(storedValue));
        }
        if (sweep.mode == CurveSweepMode::CVQuasistatic) {
            row.push_back(formatReal(point.terminalCharge));
            row.push_back(formatReal(point.capacitance));
            if (hasMultiTerminalCharges) {
                std::unordered_map<std::string, Real> extraFieldValues;
                for (const auto& [name, value] : point.extraFields)
                    extraFieldValues[name] = value;
                for (std::size_t i = 0; i < chargeColumns.size(); ++i) {
                    const auto chargeIt = extraFieldValues.find(chargeColumns.at(i).second);
                    row.push_back(formatReal(chargeIt != extraFieldValues.end() ? chargeIt->second : 0.0));
                    const auto capacitanceIt = extraFieldValues.find(capacitanceColumns.at(i).second);
                    row.push_back(formatReal(capacitanceIt != extraFieldValues.end() ? capacitanceIt->second : 0.0));
                }
            }
        }
        if (sweep.mode == CurveSweepMode::BVReverse) {
            row.push_back(formatReal(point.maxElectricField));
            row.push_back(formatReal(point.currentJumpRatio));
            row.push_back(point.breakdownDetected ? "1" : "0");
            row.push_back(formatReal(point.breakdownVoltage));
            row.push_back(point.breakdownCriterion);
            row.push_back(formatReal(point.lastStableBias));
            row.push_back(formatReal(point.failedBias));
            row.push_back(point.failureReason);
        }
        csv.writeRow(row);

        if (converged && sweep.writeVtk) {
            point.outputVtk = vtkFilename(sweep.vtkPrefix, vtkIndex++, voltage);
            writeDDSolutionVTK(point.outputVtk, mesh, doping, sol);
        }

        points.push_back(std::move(point));
    };

    bool startOk = false;
    DDSolution startSol;
    std::string startFailureReason;
    std::string startValidationDiagnostics;
    try {
        auto startAttempt = solvePoint(sweep.start, nullptr);
        startOk = startAttempt.ok;
        startSol = std::move(startAttempt.solution);
        startFailureReason = startAttempt.failureReason;
        startValidationDiagnostics = startAttempt.validationDiagnostics;
        if (!startOk && startFailureReason.empty())
            startFailureReason = "non_convergence";
    } catch (const std::exception&) {
        if (sweep.mode == CurveSweepMode::BVReverse && sweep.breakdown.nonConvergenceBreakdown)
            startFailureReason = "solver_exception";
        else
            throw;
    }
    recordPoint(sweep.start, startSol, startOk, 0.0, 0.0, 0, startFailureReason,
                startValidationDiagnostics);
    if (!startOk)
        return DCSweepResult{std::move(mesh), std::move(points)};
    previousSolution = std::move(startSol);

    DDSolution lastStepSolution;
    std::string lastStepFailureReason;
    std::string lastStepValidationDiagnostics;
    detail::DCSweepStepControlConfig stepControl;
    stepControl.start = sweep.start;
    stepControl.stop = sweep.stop;
    stepControl.step = sweep.step;
    stepControl.minStep = sweep.minStep;
    stepControl.maxStep = sweep.maxStep;
    stepControl.growthFactor = sweep.growthFactor;
    stepControl.shrinkFactor = sweep.shrinkFactor;
    stepControl.maxRetries = sweep.maxRetries;
    stepControl.stopOnFailure = sweep.stopOnFailure;

    detail::runDCSweepStepControl(
        stepControl,
        [&](Real voltage, Real, int) {
            try {
                SolvePointAttempt attempt = solvePoint(voltage, &previousSolution);
                lastStepSolution = std::move(attempt.solution);
                lastStepFailureReason = attempt.ok ? std::string() : attempt.failureReason;
                lastStepValidationDiagnostics = attempt.validationDiagnostics;
                return attempt.ok;
            } catch (const std::exception&) {
                if (sweep.mode == CurveSweepMode::BVReverse &&
                    sweep.breakdown.nonConvergenceBreakdown) {
                    lastStepSolution = DDSolution{};
                    lastStepFailureReason = "solver_exception";
                    lastStepValidationDiagnostics.clear();
                    return false;
                }
                throw;
            }
        },
        [&](const detail::DCSweepStepControlEvent& event) {
            std::string failureReason;
            if (!event.converged)
                failureReason = !lastStepFailureReason.empty()
                    ? lastStepFailureReason
                    : event.failureReason;
            recordPoint(event.voltage, lastStepSolution, event.converged,
                        event.attemptedStep, event.acceptedStep, event.retryCount, failureReason,
                        lastStepValidationDiagnostics);
            if (event.converged) {
                previousSolution = std::move(lastStepSolution);
            }
        });

    return DCSweepResult{std::move(mesh), std::move(points)};
}

namespace detail {

namespace {

void validateDCSweepStepControlConfig(const DCSweepStepControlConfig& cfg)
{
    const auto requireFinite = [](Real value, const char* name) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(std::string("DCSweep step control: ") + name +
                                        " must be finite.");
        }
    };

    requireFinite(cfg.start, "start");
    requireFinite(cfg.stop, "stop");
    requireFinite(cfg.step, "step");
    requireFinite(cfg.minStep, "minStep");
    requireFinite(cfg.maxStep, "maxStep");
    requireFinite(cfg.growthFactor, "growthFactor");
    requireFinite(cfg.shrinkFactor, "shrinkFactor");

    if (cfg.step == 0.0)
        throw std::invalid_argument("DCSweep step control: step must be non-zero.");
    if ((cfg.stop - cfg.start) * cfg.step < 0.0) {
        throw std::invalid_argument(
            "DCSweep step control: step sign must move start toward stop.");
    }
    if (cfg.minStep <= 0.0)
        throw std::invalid_argument("DCSweep step control: minStep must be positive.");
    if (cfg.maxStep <= 0.0)
        throw std::invalid_argument("DCSweep step control: maxStep must be positive.");
    if (cfg.minStep > cfg.maxStep) {
        throw std::invalid_argument(
            "DCSweep step control: minStep must not exceed maxStep.");
    }
    if (cfg.growthFactor < 1.0) {
        throw std::invalid_argument(
            "DCSweep step control: growthFactor must be at least 1.");
    }
    if (cfg.shrinkFactor <= 0.0 || cfg.shrinkFactor >= 1.0) {
        throw std::invalid_argument(
            "DCSweep step control: shrinkFactor must be greater than 0 and less than 1.");
    }
    if (cfg.maxRetries < 0) {
        throw std::invalid_argument(
            "DCSweep step control: maxRetries must be non-negative.");
    }
}

} // namespace

void runDCSweepStepControl(const DCSweepStepControlConfig& cfg,
                           const DCSweepStepAttempt& attempt,
                           const DCSweepStepRecorder& record)
{
    validateDCSweepStepControlConfig(cfg);

    Real previousVoltage = cfg.start;
    const Real direction = (cfg.step > 0.0) ? 1.0 : -1.0;
    const Real tolerance = 1.0e-12;
    Real adaptiveStep = std::min(std::abs(cfg.step), cfg.maxStep);

    auto limitedTarget = [&](Real target, Real stepMagnitude) {
        const Real remaining = direction * (target - previousVoltage);
        const Real limited = previousVoltage + direction * std::min(stepMagnitude, remaining);
        return limited;
    };

    auto advanceToward = [&](Real target) -> bool {
        int retryCount = 0;
        Real trialStep = std::min(adaptiveStep, cfg.maxStep);
        Real lastAttempted = 0.0;
        Real lastCandidate = previousVoltage;

        while (true) {
            const Real remaining = direction * (target - previousVoltage);
            if (remaining <= tolerance)
                return true;

            const Real stepMagnitude = std::min(trialStep, remaining);
            const Real candidate = limitedTarget(target, stepMagnitude);
            const Real attemptedStep = candidate - previousVoltage;
            const bool ok = attempt(candidate, attemptedStep, retryCount);
            lastAttempted = attemptedStep;
            lastCandidate = candidate;

            if (ok) {
                record({candidate, true, attemptedStep, attemptedStep, retryCount});
                previousVoltage = candidate;
                adaptiveStep = std::min(cfg.maxStep, stepMagnitude * cfg.growthFactor);
                return true;
            }

            std::string failureReason;
            if (retryCount >= cfg.maxRetries) {
                failureReason = "non_convergence";
                record({lastCandidate, false, lastAttempted, 0.0, retryCount, failureReason});
                adaptiveStep = std::max(cfg.minStep,
                                        std::min(cfg.maxStep, std::abs(lastAttempted) * cfg.shrinkFactor));
                return false;
            }

            const Real shrunken = stepMagnitude * cfg.shrinkFactor;
            if (shrunken < cfg.minStep - std::numeric_limits<Real>::epsilon()) {
                failureReason = "min_step_exhausted";
                record({lastCandidate, false, lastAttempted, 0.0, retryCount, failureReason});
                adaptiveStep = std::max(cfg.minStep,
                                        std::min(cfg.maxStep, std::abs(lastAttempted) * cfg.shrinkFactor));
                return false;
            }

            trialStep = shrunken;
            ++retryCount;
        }
    };

    bool blockedByFailedStep = false;
    Real nominalTarget = cfg.start + cfg.step;
    while (!blockedByFailedStep && direction * (nominalTarget - cfg.stop) <= tolerance) {
        while (direction * (previousVoltage - nominalTarget) < -tolerance) {
            if (!advanceToward(nominalTarget)) {
                if (cfg.stopOnFailure)
                    return;
                blockedByFailedStep = true;
                break;
            }
        }
        nominalTarget += cfg.step;
    }

    while (!blockedByFailedStep && direction * (previousVoltage - cfg.stop) < -tolerance) {
        if (!advanceToward(cfg.stop)) {
            if (cfg.stopOnFailure)
                return;
            break;
        }
    }
}

} // namespace detail

} // namespace vela
