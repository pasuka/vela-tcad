#include "vela/io/MeshReader.h"
#include "vela/io/VTKWriter.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/simulation/CurveSweep.h"
#include "vela/simulation/DCSweep.h"
#include "vela/simulation/PoissonSimulation.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace {

using MaterialLoadJsonLegacy = void (vela::MaterialDatabase::*)(const std::string&);

struct LastResult {
    std::shared_ptr<vela::DeviceMesh> mesh;
    std::vector<std::pair<std::string, std::vector<vela::Real>>> nodeScalars;
};

struct SweepPythonMetadata {
    std::string curveType = "iv";
    std::string biasContact;
    std::string currentContact;
    std::string chargeContact;
    bool unitScaling = false;
    bool writeVtk = false;
    std::string outputCsv;
    std::string vtkPrefix;
};

LastResult& lastResult()
{
    static LastResult result;
    return result;
}

std::vector<vela::Real> eigenToVector(const vela::VectorXd& values)
{
    std::vector<vela::Real> out(static_cast<std::size_t>(values.size()));
    for (int i = 0; i < values.size(); ++i)
        out[static_cast<std::size_t>(i)] = values(i);
    return out;
}

std::shared_ptr<vela::DeviceMesh> loadMeshImpl(const std::string& jsonFile)
{
    vela::JsonMeshReader reader;
    return std::make_shared<vela::DeviceMesh>(reader.read(jsonFile));
}

void rememberMesh(std::shared_ptr<vela::DeviceMesh> mesh)
{
    LastResult& result = lastResult();
    result.mesh = std::move(mesh);
    result.nodeScalars.clear();
}

std::filesystem::path configDirectory(const std::string& configFile)
{
    const std::filesystem::path configPath(configFile);
    const std::filesystem::path parent = configPath.parent_path();
    if (!parent.empty())
        return parent;
    return std::filesystem::current_path();
}

std::string resolvePath(const std::filesystem::path& baseDir, const std::string& path)
{
    std::filesystem::path resolved(path);
    if (resolved.is_relative())
        resolved = baseDir / resolved;
    return resolved.string();
}

nlohmann::json readJsonFile(const std::string& configFile)
{
    std::ifstream input(configFile);
    if (!input)
        throw std::runtime_error("Python API: unable to open sweep config file: " + configFile);

    nlohmann::json cfg;
    input >> cfg;
    return cfg;
}

SweepPythonMetadata sweepMetadataFromConfig(const std::string& configFile)
{
    const nlohmann::json cfg = readJsonFile(configFile);
    const auto& sweep = cfg.at("sweep");
    const std::filesystem::path cfgDir = configDirectory(configFile);

    SweepPythonMetadata metadata;
    metadata.curveType = vela::toString(vela::curveSweepModeFromString(sweep.value("mode", std::string("iv"))));
    metadata.biasContact = sweep.at("contact").get<std::string>();
    metadata.currentContact = sweep.value("current_contact", metadata.biasContact);
    if (cfg.contains("scaling")) {
        const auto& scaling = cfg.at("scaling");
        if (!scaling.is_object())
            throw std::invalid_argument("scaling must be an object.");
        const std::string mode = scaling.value("mode", std::string{});
        if (mode == "unit_scaling")
            metadata.unitScaling = true;
        else if (!mode.empty())
            throw std::invalid_argument("Unsupported scaling.mode '" + mode + "'.");
    }
    const nlohmann::json chargeCfg = sweep.value("terminal_charge", nlohmann::json::object());
    metadata.chargeContact = chargeCfg.value("contact", sweep.value("charge_contact", metadata.biasContact));
    metadata.writeVtk = sweep.value("write_vtk", cfg.value("write_vtk", false));
    metadata.outputCsv = resolvePath(
        cfgDir, sweep.value("csv_file", cfg.value("output_csv", std::string("dc_sweep.csv"))));
    metadata.vtkPrefix = resolvePath(
        cfgDir, sweep.value("vtk_prefix", cfg.value("output_vtk_prefix", std::string("dc_sweep"))));
    return metadata;
}

py::dict convergenceDiagnostics(const vela::DCSweepPoint& point)
{
    py::dict diagnostics;
    diagnostics["converged"] = point.converged;
    diagnostics["iterations"] = point.iterations;
    diagnostics["attempted_step"] = point.attemptedStep;
    diagnostics["accepted_step"] = point.acceptedStep;
    diagnostics["retry_count"] = point.retryCount;
    diagnostics["failed"] = point.failed;
    diagnostics["last_stable_bias"] = point.lastStableBias;
    diagnostics["failed_bias"] = point.failedBias;
    diagnostics["failure_reason"] = point.failureReason;
    diagnostics["validation_diagnostics"] = point.validationDiagnostics;
    return diagnostics;
}

std::vector<py::dict> sweepPointsToPython(const std::vector<vela::DCSweepPoint>& points,
                                          const SweepPythonMetadata& metadata)
{
    std::vector<py::dict> out;
    out.reserve(points.size());
    for (const vela::DCSweepPoint& point : points) {
        py::dict row;
        row["curve_type"] = metadata.curveType;
        row["scaling_mode"] = metadata.unitScaling ? "unit_scaling" : "legacy";
        row["bias_contact"] = metadata.biasContact;
        row["current_contact"] = metadata.currentContact;
        row["charge_contact"] = metadata.chargeContact;
        row["voltage"] = point.voltage;
        row["bias"] = point.bias;
        row["electron_current"] = point.electronCurrent;
        row["hole_current"] = point.holeCurrent;
        row["total_current"] = point.totalCurrent;
        row["converged"] = point.converged;
        row["iterations"] = point.iterations;
        row["attempted_step"] = point.attemptedStep;
        row["accepted_step"] = point.acceptedStep;
        row["retry_count"] = point.retryCount;
        row["convergence_diagnostics"] = convergenceDiagnostics(point);
        row["terminal_charge"] = point.terminalCharge;
        row["capacitance"] = point.capacitance;
        py::dict terminalCharges;
        for (const auto& [name, value] : point.terminalChargeValues)
            terminalCharges[py::str(name)] = value;
        row["terminal_charges"] = std::move(terminalCharges);
        py::dict terminalCapacitances;
        for (const auto& [name, value] : point.terminalCapacitanceValues)
            terminalCapacitances[py::str(name)] = value;
        row["terminal_capacitances"] = std::move(terminalCapacitances);
        for (const auto& [name, value] : point.extraFields)
            row[py::str(name)] = value;
        row["max_electric_field"] = point.maxElectricField;
        row["current_jump_ratio"] = point.currentJumpRatio;
        row["breakdown_detected"] = point.breakdownDetected;
        row["breakdown_voltage"] = point.breakdownVoltage;
        row["breakdown_criterion"] = point.breakdownCriterion;
        row["failed"] = point.failed;
        row["last_stable_bias"] = point.lastStableBias;
        row["failed_bias"] = point.failedBias;
        row["failure_reason"] = point.failureReason;
        row["validation_diagnostics"] = point.validationDiagnostics;
        row["output_csv"] = point.outputCsv.empty() ? metadata.outputCsv : point.outputCsv;
        row["output_vtk"] = point.outputVtk;
        py::list outputFiles;
        outputFiles.append(row["output_csv"]);
        if (!point.outputVtk.empty())
            outputFiles.append(point.outputVtk);
        row["output_files"] = std::move(outputFiles);
        out.push_back(std::move(row));
    }
    return out;
}

std::vector<py::dict> sweepPointsToPython(const std::string& configFile,
                                          const std::vector<vela::DCSweepPoint>& points)
{
    return sweepPointsToPython(points, sweepMetadataFromConfig(configFile));
}

std::shared_ptr<vela::DeviceMesh> loadMesh(const std::string& jsonFile)
{
    auto mesh = loadMeshImpl(jsonFile);
    rememberMesh(mesh);
    return mesh;
}

std::vector<vela::Real> runPoisson(const std::string& configFile)
{
    vela::PoissonSimulation simulation;
    vela::PoissonResult poissonResult = simulation.runWithResult(configFile);

    LastResult& result = lastResult();
    result.mesh = std::make_shared<vela::DeviceMesh>(std::move(poissonResult.mesh));
    result.nodeScalars.clear();
    result.nodeScalars.emplace_back("potential_V", eigenToVector(poissonResult.potential));
    result.nodeScalars.emplace_back("net_doping_m3", std::move(poissonResult.netDoping));

    return result.nodeScalars.front().second;
}

std::vector<py::dict> runDCSweep(const std::string& configFile)
{
    vela::DCSweep sweep;
    const std::vector<vela::DCSweepPoint> points = sweep.run(configFile);
    return sweepPointsToPython(configFile, points);
}

void writeVtk(const std::string& outputFile)
{
    const LastResult& result = lastResult();
    if (!result.mesh)
        throw std::runtime_error("write_vtk requires a prior load_mesh() or run_poisson() call.");

    const std::filesystem::path outputPath(outputFile);
    if (const std::filesystem::path parent = outputPath.parent_path(); !parent.empty())
        std::filesystem::create_directories(parent);

    vela::VTKWriter writer(outputFile, *result.mesh);
    writer.write();
    for (const auto& [name, values] : result.nodeScalars)
        writer.addNodeScalar(name, values);
}

} // namespace

PYBIND11_MODULE(_core, m)
{
    m.doc() = "Python bindings for the Vela TCAD C++ core";

    py::class_<vela::DeviceMesh, std::shared_ptr<vela::DeviceMesh>>(m, "DeviceMesh")
        .def("num_nodes", &vela::DeviceMesh::numNodes)
        .def("num_edges", &vela::DeviceMesh::numEdges)
        .def("num_cells", &vela::DeviceMesh::numCells)
        .def("num_regions", &vela::DeviceMesh::numRegions)
        .def("num_contacts", &vela::DeviceMesh::numContacts);

    py::class_<vela::MaterialDatabase>(m, "MaterialDatabase")
        .def(py::init<>())
        .def(py::init<const std::string&>(), py::arg("materials_file"))
        .def("load_json",
             static_cast<MaterialLoadJsonLegacy>(&vela::MaterialDatabase::loadJson),
             py::arg("materials_file"))
        .def("has_material", &vela::MaterialDatabase::hasMaterial);

    py::class_<vela::PoissonSimulation>(m, "PoissonSimulation")
        .def(py::init<>())
        .def("run", [](vela::PoissonSimulation& self, const std::string& configFile) {
            return eigenToVector(self.run(configFile));
        });

    py::class_<vela::DCSweep>(m, "DCSweep")
        .def(py::init<>())
        .def("run", [](const vela::DCSweep& self, const std::string& configFile) {
            return sweepPointsToPython(configFile, self.run(configFile));
        });

    m.def("load_mesh", &loadMesh, py::arg("json_file"),
          "Read a JSON mesh file with the C++ mesh reader.");
    m.def("run_poisson", &runPoisson, py::arg("config_file"),
          "Run a Poisson simulation from a Vela JSON config file.");
    m.def("run_dc_sweep", &runDCSweep, py::arg("config_file"),
          "Run a DC sweep from a Vela JSON config file.");
    m.def("write_vtk", &writeVtk, py::arg("output_file"),
          "Export the most recently loaded mesh or Poisson result to VTK.");
}
