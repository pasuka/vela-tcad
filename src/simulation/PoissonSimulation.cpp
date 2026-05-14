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
#include <utility>
#include <vector>

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

void appendFixedChargeSpec(std::vector<RegionFixedChargeSpec>& specs,
                           std::unordered_map<std::string, std::string>& sourcesByRegion,
                           std::string region,
                           Real fixedCharge,
                           std::string source)
{
    const auto [_, inserted] = sourcesByRegion.emplace(region, source);
    if (!inserted)
        throw std::runtime_error(
            "PoissonSimulation: duplicate fixed_charge_m3 for region '" +
            region + "' from " + sourcesByRegion.at(region) + " and " +
            source + ". Specify fixed charge for each region only once.");

    specs.push_back(RegionFixedChargeSpec{std::move(region), fixedCharge});
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
    const std::string materialsFile = cfg.contains("materials_file")
        ? resolvePath(cfgDir, cfg.at("materials_file").get<std::string>())
        : std::string{};

    // ------------------------------------------------------------------
    // Build mesh
    // ------------------------------------------------------------------
    JsonMeshReader reader;
    DeviceMesh mesh = reader.read(meshFile);

    // ------------------------------------------------------------------
    // Material database (built-in Si, SiO2 plus optional config override)
    // ------------------------------------------------------------------
    MaterialDatabase matdb;
    if (!materialsFile.empty())
        matdb.loadJson(materialsFile);

    // ------------------------------------------------------------------
    // Doping model
    // ------------------------------------------------------------------
    std::vector<RegionDopingSpec> dopingSpecs;
    std::vector<RegionFixedChargeSpec> fixedChargeSpecs;
    std::unordered_map<std::string, std::string> fixedChargeSourcesByRegion;
    for (const auto& entry : cfg.at("doping")) {
        const auto region = entry.at("region").get<std::string>();

        RegionDopingSpec spec;
        spec.region    = region;
        spec.donors    = entry.at("donors").get<Real>();
        spec.acceptors = entry.at("acceptors").get<Real>();
        dopingSpecs.push_back(std::move(spec));

        if (entry.contains("fixed_charge_m3")) {
            appendFixedChargeSpec(
                fixedChargeSpecs,
                fixedChargeSourcesByRegion,
                region,
                entry.at("fixed_charge_m3").get<Real>(),
                "doping entry");
        }
    }
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, dopingSpecs);

    if (cfg.contains("regions")) {
        for (const auto& entry : cfg.at("regions")) {
            if (!entry.contains("fixed_charge_m3")) continue;
            appendFixedChargeSpec(
                fixedChargeSpecs,
                fixedChargeSourcesByRegion,
                entry.at("name").get<std::string>(),
                entry.at("fixed_charge_m3").get<Real>(),
                "regions entry");
        }
    }

    std::vector<InterfaceSheetChargeSpec> sheetChargeSpecs;
    if (cfg.contains("interfaces")) {
        for (const auto& entry : cfg.at("interfaces")) {
            if (!entry.contains("sheet_charge_m2") &&
                !entry.contains("fixed_charge_m2") &&
                !entry.contains("trap_density_m2")) continue;

            const Real sheetCharge = entry.value("sheet_charge_m2", 0.0);
            const Real fixedCharge = entry.value("fixed_charge_m2", 0.0);
            const Real trapDensity = entry.value("trap_density_m2", 0.0);
            const Real trapOccupancy = entry.value("trap_occupancy", 0.0);

            if (entry.contains("regions")) {
                const auto regions = entry.at("regions").get<std::vector<std::string>>();
                if (regions.size() != 2)
                    throw std::runtime_error(
                        "PoissonSimulation: interface regions must contain exactly two names.");
                sheetChargeSpecs.push_back(InterfaceSheetChargeSpec{
                    regions[0], regions[1], sheetCharge, fixedCharge, trapDensity, trapOccupancy});
            } else {
                sheetChargeSpecs.push_back(InterfaceSheetChargeSpec{
                    entry.at("region0").get<std::string>(),
                    entry.at("region1").get<std::string>(),
                    sheetCharge, fixedCharge, trapDensity, trapOccupancy});
            }
        }
    }

    // ------------------------------------------------------------------
    // Assemble Poisson equation
    // ------------------------------------------------------------------
    PoissonAssembler assembler(mesh, matdb, doping,
                               std::move(fixedChargeSpecs),
                               std::move(sheetChargeSpecs));
    assembler.assemble();

    // ------------------------------------------------------------------
    // Dirichlet boundary conditions from contacts
    // ------------------------------------------------------------------
    // Build name -> effective contact potential map from config. A configured
    // flatband voltage or work function shifts the electrostatic Dirichlet
    // potential as psi_contact = bias - offset. work_function_eV is interpreted
    // as its equivalent voltage because 1 eV/q is 1 V.
    std::unordered_map<std::string, Real> contactBias;
    for (const auto& ct : cfg.at("contacts")) {
        const bool hasFlatband = ct.contains("flatband_voltage");
        const bool hasWorkFunction = ct.contains("work_function_eV");
        if (hasFlatband && hasWorkFunction)
            throw std::runtime_error(
                "PoissonSimulation: contact cannot set both flatband_voltage and work_function_eV.");

        Real value = ct.at("bias").get<Real>();
        if (hasFlatband)
            value -= ct.at("flatband_voltage").get<Real>();
        if (hasWorkFunction)
            value -= ct.at("work_function_eV").get<Real>();

        contactBias[ct.at("name").get<std::string>()] = value;
    }

    // Map contact nodes -> prescribed potential
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
