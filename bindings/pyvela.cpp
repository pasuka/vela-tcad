#include "vela/io/MeshReader.h"
#include "vela/io/VTKWriter.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/simulation/DCSweep.h"
#include "vela/simulation/PoissonSimulation.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace {

struct LastResult {
    std::shared_ptr<vela::DeviceMesh> mesh;
    std::vector<std::pair<std::string, std::vector<vela::Real>>> nodeScalars;
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

std::vector<py::dict> sweepPointsToPython(const std::vector<vela::DCSweepPoint>& points)
{
    std::vector<py::dict> out;
    out.reserve(points.size());
    for (const vela::DCSweepPoint& point : points) {
        py::dict row;
        row["voltage"] = point.voltage;
        row["electron_current"] = point.electronCurrent;
        row["hole_current"] = point.holeCurrent;
        row["total_current"] = point.totalCurrent;
        row["converged"] = point.converged;
        row["iterations"] = point.iterations;
        row["attempted_step"] = point.attemptedStep;
        row["accepted_step"] = point.acceptedStep;
        row["retry_count"] = point.retryCount;
        out.push_back(std::move(row));
    }
    return out;
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
    return sweepPointsToPython(points);
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
        .def("load_json", &vela::MaterialDatabase::loadJson, py::arg("materials_file"))
        .def("has_material", &vela::MaterialDatabase::hasMaterial);

    py::class_<vela::PoissonSimulation>(m, "PoissonSimulation")
        .def(py::init<>())
        .def("run", [](vela::PoissonSimulation& self, const std::string& configFile) {
            return eigenToVector(self.run(configFile));
        });

    py::class_<vela::DCSweep>(m, "DCSweep")
        .def(py::init<>())
        .def("run", [](const vela::DCSweep& self, const std::string& configFile) {
            return sweepPointsToPython(self.run(configFile));
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
