#include "vela/simulation/PoissonSimulation.h"
#include "vela/boundary/BoundaryCondition.h"
#include "vela/equation/PoissonAssembler.h"
#include "vela/io/MeshReader.h"
#include "vela/io/VTKWriter.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/LinearSolver.h"
#include "vela/simulation/ConfigParsing.h"
#include <nlohmann/json.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
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
    const UnitScalingConfig scaling = parseUnitScalingConfig(cfg);

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
    DeviceMesh mesh = reader.read(meshFile, scaling);

    // ------------------------------------------------------------------
    // Material database (built-in Si, SiO2 plus optional config override)
    // ------------------------------------------------------------------
    MaterialDatabase matdb;
    if (!materialsFile.empty())
        matdb.loadJson(materialsFile, scaling);

    // ------------------------------------------------------------------
    // Doping model
    // ------------------------------------------------------------------
    std::vector<RegionDopingSpec> dopingSpecs = parseDopingSpecs(cfg, scaling);
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, dopingSpecs);

    std::vector<RegionFixedChargeSpec> fixedChargeSpecs =
        parseRegionFixedChargeSpecs(cfg, scaling);
    std::vector<InterfaceSheetChargeSpec> sheetChargeSpecs =
        parseInterfaceSheetChargeSpecs(cfg, scaling);

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

    // Helper used for Schottky contacts to look up the contact-region material.
    auto contactMaterialByName = [&](const std::string& name) -> std::optional<Material> {
        for (Index c = 0; c < mesh.numContacts(); ++c) {
            const Contact& ct = mesh.getContact(c);
            if (ct.name != name) continue;
            if (ct.region_id < mesh.numRegions()) {
                const Region& region = mesh.getRegion(ct.region_id);
                if (matdb.hasMaterial(region.material))
                    return matdb.getMaterial(region.material);
            }
            for (Index ci = 0; ci < mesh.numCells(); ++ci) {
                const Cell& cell = mesh.getCell(ci);
                for (Index nid : cell.node_ids) {
                    for (Index cnid : ct.node_ids) {
                        if (nid == cnid) {
                            const Region& region = mesh.getRegion(cell.region_id);
                            if (matdb.hasMaterial(region.material))
                                return matdb.getMaterial(region.material);
                        }
                    }
                }
            }
        }
        return std::nullopt;
    };

    std::unordered_map<std::string, Real> contactBias;
    for (const auto& spec : contactSpecs) {
        switch (spec.type) {
            case ContactType::Ohmic:
            case ContactType::Dirichlet:
            case ContactType::MetalGate:
                contactBias[spec.name] = effectivePoissonDirichletPotential(spec);
                break;
            case ContactType::Schottky: {
                // Schottky contacts: use the same Dirichlet-barrier prototype
                // as the DD path so the Poisson-only field plot is consistent
                // with the DD result.  psi_contact = bias - (phi_Bn - Eg/2)
                // when bandgap+affinity are known, else psi_contact = bias -
                // phi_Bn (1 eV/q == 1 V).
                const auto materialOpt = contactMaterialByName(spec.name);
                const Real bandgap_eV = materialOpt && materialOpt->bandgap_eV
                    ? *materialOpt->bandgap_eV
                    : std::numeric_limits<Real>::quiet_NaN();
                const Real affinity_eV = materialOpt && materialOpt->electron_affinity_eV
                    ? *materialOpt->electron_affinity_eV
                    : std::numeric_limits<Real>::quiet_NaN();
                const Real phiBn = schottkyElectronBarrier_eV(spec, affinity_eV);
                Real psiOffset = phiBn;
                if (std::isfinite(bandgap_eV) && bandgap_eV > 0.0 &&
                    std::isfinite(affinity_eV)) {
                    psiOffset = phiBn - 0.5 * bandgap_eV;
                }
                if (spec.workFunction_eV && std::isfinite(affinity_eV) &&
                    !spec.barrier_eV && !spec.electronBarrier_eV &&
                    std::isfinite(bandgap_eV) && bandgap_eV > 0.0) {
                    psiOffset = *spec.workFunction_eV - affinity_eV - 0.5 * bandgap_eV;
                }
                contactBias[spec.name] = spec.bias - psiOffset;
                break;
            }
            case ContactType::Floating:
                throw std::runtime_error(
                    "PoissonSimulation: contact '" + spec.name +
                    "' has type '" + toString(spec.type) +
                    "' which is not yet implemented for the Poisson driver. "
                    "Use ohmic, dirichlet, metal_gate, or schottky for now.");
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
