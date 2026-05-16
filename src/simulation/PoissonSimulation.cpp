#include "vela/simulation/PoissonSimulation.h"
#include "vela/boundary/BoundaryCondition.h"
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
            if (entry.contains("trap_density_m2") &&
                (trapOccupancy < 0.0 || trapOccupancy > 1.0)) {
                throw std::runtime_error(
                    "PoissonSimulation: trap_occupancy must be in [0, 1] when trap_density_m2 is provided.");
            }

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
    // Parse explicit boundary conditions
    // ------------------------------------------------------------------
    const std::vector<BoundarySegmentSpec> boundarySpecs =
        parseBoundarySegmentSpecs(cfg);

    std::vector<PoissonNeumannBoundarySpec> neumannBoundarySpecs;
    for (const auto& spec : boundarySpecs) {
        switch (spec.type) {
            case BoundaryType::Neumann:
            case BoundaryType::Insulating:
            case BoundaryType::Symmetry: {
                // Insulating and symmetry are zero Neumann
                const Real displacement = (spec.type == BoundaryType::Neumann)
                    ? spec.value : 0.0;
                neumannBoundarySpecs.push_back(
                    PoissonNeumannBoundarySpec{spec.node_ids, displacement});
                break;
            }
            case BoundaryType::Dirichlet:
                // Already rejected by parser
                throw std::runtime_error(
                    "PoissonSimulation: boundary '" + spec.name +
                    "' has type 'dirichlet' which should have been rejected by parser.");
        }
    }

    // ------------------------------------------------------------------
    // Assemble Poisson equation
    // ------------------------------------------------------------------
    PoissonAssembler assembler(mesh, matdb, doping,
                               std::move(fixedChargeSpecs),
                               std::move(sheetChargeSpecs),
                               std::move(neumannBoundarySpecs));
    assembler.assemble();

    // ------------------------------------------------------------------
    // Dirichlet boundary conditions from contacts
    // ------------------------------------------------------------------
    // The unified boundary parser interprets each ``contacts[]`` entry and
    // normalises optional ``type`` / ``flatband_voltage`` / ``work_function_eV``
    // fields.  Legacy decks without ``type`` are treated as ohmic; ohmic,
    // dirichlet, and metal_gate contacts all map to an effective Poisson
    // Dirichlet potential.  Schottky and floating contacts do not yet have a
    // dedicated Poisson model and are rejected here so misconfiguration
    // surfaces early.
    const std::vector<ContactBoundarySpec> contactSpecs =
        parseContactBoundarySpecs(cfg);

    std::unordered_map<std::string, Real> contactBias;
    for (const auto& spec : contactSpecs) {
        switch (spec.type) {
            case ContactType::Ohmic:
            case ContactType::Dirichlet:
            case ContactType::MetalGate:
                contactBias[spec.name] = effectivePoissonDirichletPotential(spec);
                break;
            case ContactType::Schottky:
            case ContactType::Floating:
                throw std::runtime_error(
                    "PoissonSimulation: contact '" + spec.name +
                    "' has type '" + toString(spec.type) +
                    "' which is not yet implemented for the Poisson driver. "
                    "Use ohmic, dirichlet, or metal_gate for now.");
        }
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
