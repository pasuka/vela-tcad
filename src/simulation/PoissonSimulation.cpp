#include "vela/simulation/PoissonSimulation.h"
#include "vela/equation/PoissonAssembler.h"
#include "vela/io/MeshReader.h"
#include "vela/io/VTKWriter.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/LinearSolver.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <unordered_map>

namespace vela {

namespace {

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

} // namespace

VectorXd PoissonSimulation::run(const std::string& configFile)
{
    return runWithResult(configFile).potential;
}

PoissonResult PoissonSimulation::runWithResult(const std::string& configFile)
{
    // ------------------------------------------------------------------
    // Load config JSON
    // ------------------------------------------------------------------
    std::ifstream ifs(configFile);
    if (!ifs.is_open())
        throw std::runtime_error(
            "PoissonSimulation: cannot open config file: " + configFile);

    nlohmann::json cfg;
    ifs >> cfg;

    // Resolve paths relative to the config file's directory
    const std::filesystem::path cfgDir = configDirectory(configFile);
    const std::string meshFile = resolvePath(cfgDir, cfg.at("mesh_file").get<std::string>());
    const std::string outputVtk = resolvePath(cfgDir, cfg.at("output_vtk").get<std::string>());

    // ------------------------------------------------------------------
    // Build mesh
    // ------------------------------------------------------------------
    JsonMeshReader reader;
    DeviceMesh mesh = reader.read(meshFile);

    // ------------------------------------------------------------------
    // Material database (built-in Si, SiO2)
    // ------------------------------------------------------------------
    MaterialDatabase matdb;

    // ------------------------------------------------------------------
    // Doping model
    // ------------------------------------------------------------------
    std::vector<RegionDopingSpec> dopingSpecs;
    for (const auto& entry : cfg.at("doping")) {
        RegionDopingSpec spec;
        spec.region    = entry.at("region").get<std::string>();
        spec.donors    = entry.at("donors").get<Real>();
        spec.acceptors = entry.at("acceptors").get<Real>();
        dopingSpecs.push_back(std::move(spec));
    }
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, dopingSpecs);

    // ------------------------------------------------------------------
    // Assemble Poisson equation
    // ------------------------------------------------------------------
    PoissonAssembler assembler(mesh, matdb, doping);
    assembler.assemble();

    // ------------------------------------------------------------------
    // Dirichlet boundary conditions from contacts
    // ------------------------------------------------------------------
    // Build name → bias map from config
    std::unordered_map<std::string, Real> contactBias;
    for (const auto& ct : cfg.at("contacts")) {
        contactBias[ct.at("name").get<std::string>()] =
            ct.at("bias").get<Real>();
    }

    // Map contact nodes → prescribed potential
    std::unordered_map<Index, Real> dirichletBCs;
    for (Index c = 0; c < mesh.numContacts(); ++c) {
        const Contact& contact = mesh.getContact(c);
        auto it = contactBias.find(contact.name);
        if (it == contactBias.end()) continue;
        for (Index nid : contact.node_ids)
            dirichletBCs[nid] = it->second;
    }

    assembler.applyDirichlet(dirichletBCs);

    // ------------------------------------------------------------------
    // Solve
    // ------------------------------------------------------------------
    LinearSolver solver;
    VectorXd psi = solver.solve(assembler.matrix(), assembler.rhs());

    std::vector<Real> psiVec(mesh.numNodes());
    for (Index i = 0; i < mesh.numNodes(); ++i)
        psiVec[i] = psi(static_cast<int>(i));

    // Also write net doping for visualisation
    std::vector<Real> dopingVec(mesh.numNodes());
    for (Index i = 0; i < mesh.numNodes(); ++i)
        dopingVec[i] = doping.netDoping(i);

    // ------------------------------------------------------------------
    // Write VTK output
    // ------------------------------------------------------------------
    VTKWriter writer(outputVtk, mesh);
    writer.write();
    writer.addNodeScalar("potential_V", psiVec);
    writer.addNodeScalar("net_doping_m3", dopingVec);

    return PoissonResult{std::move(mesh), std::move(psi), std::move(dopingVec)};
}

} // namespace vela
