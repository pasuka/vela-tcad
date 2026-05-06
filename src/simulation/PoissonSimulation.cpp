#include "vela/simulation/PoissonSimulation.h"
#include "vela/io/MeshReader.h"
#include "vela/io/VTKWriter.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/equation/PoissonAssembler.h"
#include "vela/solver/LinearSolver.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include <filesystem>
#include <unordered_map>

namespace vela {

VectorXd PoissonSimulation::run(const std::string& configFile)
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
    std::filesystem::path cfgDir =
        std::filesystem::path(configFile).parent_path();

    auto resolvePath = [&](const std::string& p) -> std::string {
        std::filesystem::path fp(p);
        if (fp.is_relative())
            return (cfgDir / fp).string();
        return p;
    };

    const std::string meshFile   = resolvePath(cfg.at("mesh_file").get<std::string>());
    const std::string outputVtk  = resolvePath(cfg.at("output_vtk").get<std::string>());

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
    DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, cfg.at("doping"));

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

    // ------------------------------------------------------------------
    // Write VTK output
    // ------------------------------------------------------------------
    VTKWriter writer(outputVtk, mesh);
    writer.write();

    std::vector<Real> psiVec(mesh.numNodes());
    for (Index i = 0; i < mesh.numNodes(); ++i)
        psiVec[i] = psi(static_cast<int>(i));
    writer.addNodeScalar("potential_V", psiVec);

    // Also write net doping for visualisation
    std::vector<Real> dopingVec(mesh.numNodes());
    for (Index i = 0; i < mesh.numNodes(); ++i)
        dopingVec[i] = doping.netDoping(i);
    writer.addNodeScalar("net_doping_m3", dopingVec);

    return psi;
}

} // namespace vela
