#pragma once

/**
 * @file AssemblerUtils.h
 * @brief Shared FVM assembly helpers used by PoissonAssembler and DDAssembler.
 *
 * All functions are free (non-member) to avoid duplication between assemblers.
 * Include this header from assembler .cpp files; do not expose it as part of
 * the public library API.
 */

#include "vela/core/Types.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/equation/ChargeSpec.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/Material.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/physics/BandgapNarrowing.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <string>
#include <utility>

namespace vela::detail {


// ---------------------------------------------------------------------------
// Fixed and interface charge helpers
// ---------------------------------------------------------------------------

struct RegionPairKey {
    std::string first;
    std::string second;

    bool operator==(const RegionPairKey& other) const
    {
        return first == other.first && second == other.second;
    }
};

struct RegionPairKeyHash {
    std::size_t operator()(const RegionPairKey& key) const
    {
        const std::hash<std::string> hash;
        std::size_t seed = hash(key.first);
        seed ^= hash(key.second) + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        return seed;
    }
};

inline RegionPairKey makeRegionPairKey(std::string a, std::string b)
{
    if (b < a)
        std::swap(a, b);
    return RegionPairKey{std::move(a), std::move(b)};
}

inline Real triangleArea(const DeviceMesh& mesh, const Cell& cell)
{
    if (cell.node_ids.size() < 3) return 0.0;

    const Node& a = mesh.getNode(cell.node_ids[0]);
    const Node& b = mesh.getNode(cell.node_ids[1]);
    const Node& c = mesh.getNode(cell.node_ids[2]);

    return 0.5 * std::abs((b.x - a.x) * (c.y - a.y) -
                          (c.x - a.x) * (b.y - a.y));
}

inline std::unordered_map<std::string, Real> fixedChargeByRegion(
    const std::vector<RegionFixedChargeSpec>& fixedCharges,
    const std::string& context)
{
    std::unordered_map<std::string, Real> fixedByRegion;
    for (const auto& spec : fixedCharges) {
        const auto [_, inserted] = fixedByRegion.emplace(spec.region, spec.fixedCharge);
        if (!inserted)
            throw std::invalid_argument(
                context + ": duplicate fixed_charge_m3 for region '" + spec.region + "'.");
    }
    return fixedByRegion;
}

inline std::unordered_map<RegionPairKey, Real, RegionPairKeyHash> sheetChargeByRegionPair(
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges)
{
    std::unordered_map<RegionPairKey, Real, RegionPairKeyHash> sheetByRegionPair;
    for (const auto& spec : sheetCharges)
        sheetByRegionPair[makeRegionPairKey(spec.region0, spec.region1)] += spec.totalSheetCharge();
    return sheetByRegionPair;
}

inline VectorXd computeFixedAndInterfaceChargeRhs(
    const DeviceMesh& mesh,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<RegionFixedChargeSpec>& fixedCharges,
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges,
    const std::string& context)
{
    VectorXd contribution = VectorXd::Zero(static_cast<int>(mesh.numNodes()));

    const auto fixedByRegion = fixedChargeByRegion(fixedCharges, context);
    if (!fixedByRegion.empty()) {
        for (Index c = 0; c < mesh.numCells(); ++c) {
            const Cell& cell = mesh.getCell(c);
            const Region& region = mesh.getRegion(cell.region_id);
            auto it = fixedByRegion.find(region.name);
            if (it == fixedByRegion.end()) continue;

            const Real nodeCharge = constants::q * it->second * triangleArea(mesh, cell) / 3.0;
            for (Index nid : cell.node_ids)
                contribution(static_cast<int>(nid)) += nodeCharge;
        }
    }

    const auto sheetByRegionPair = sheetChargeByRegionPair(sheetCharges);
    if (!sheetByRegionPair.empty()) {
        for (Index e = 0; e < mesh.numEdges(); ++e) {
            const auto& cells = edgeCells[e];
            if (cells.size() != 2) continue;

            const Region& r0 = mesh.getRegion(mesh.getCell(cells[0]).region_id);
            const Region& r1 = mesh.getRegion(mesh.getCell(cells[1]).region_id);
            const auto it = sheetByRegionPair.find(makeRegionPairKey(r0.name, r1.name));
            if (it == sheetByRegionPair.end()) continue;

            const Edge& edge = mesh.getEdge(e);
            const Real endpointCharge = constants::q * it->second * edge.length * 0.5;
            contribution(static_cast<int>(edge.n0)) += endpointCharge;
            contribution(static_cast<int>(edge.n1)) += endpointCharge;
        }
    }

    return contribution;
}

inline void addFixedAndInterfaceChargeToRhs(
    const DeviceMesh& mesh,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<RegionFixedChargeSpec>& fixedCharges,
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges,
    VectorXd& rhs,
    const std::string& context)
{
    rhs += computeFixedAndInterfaceChargeRhs(mesh, edgeCells, fixedCharges, sheetCharges, context);
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

/// Return precomputed per-node control-volume areas.
inline std::vector<Real> computeNodeVolumes(const DeviceMesh& mesh)
{
    const Index N = mesh.numNodes();
    std::vector<Real> vol(N, 0.0);
    for (Index i = 0; i < N; ++i)
        vol[i] = mesh.getNode(i).volume;
    return vol;
}

/// Return precomputed per-edge box coupling lengths.
inline std::vector<Real> computeEdgeCouplings(const DeviceMesh& mesh)
{
    const Index E = mesh.numEdges();
    std::vector<Real> couple(E, 0.0);
    for (Index e = 0; e < E; ++e)
        couple[e] = mesh.getEdge(e).couple;
    return couple;
}

/// Return a per-node max adjacent-edge scalar-gradient magnitude [scalar unit/m].
inline std::vector<Real> computeNodeScalarGradientMagnitudes(const VectorXd& value,
                                                            const DeviceMesh& mesh)
{
    std::vector<Real> maxField(mesh.numNodes(), 0.0);
    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        if (edge.length <= 1.0e-30)
            continue;
        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real edgeField = std::abs((value(j) - value(i)) / edge.length);
        maxField[edge.n0] = std::max(maxField[edge.n0], edgeField);
        maxField[edge.n1] = std::max(maxField[edge.n1], edgeField);
    }
    return maxField;
}

/// Return a per-node max adjacent-edge electric-field magnitude [V/m].
inline std::vector<Real> computeNodeElectricFields(const VectorXd& psi, const DeviceMesh& mesh)
{
    return computeNodeScalarGradientMagnitudes(psi, mesh);
}

/// Build edge -> adjacent cell ids map.
inline std::vector<std::vector<Index>> buildEdgeCellMap(const DeviceMesh& mesh)
{
    std::vector<std::vector<Index>> edgeCells(mesh.numEdges());
    std::unordered_map<Index, Index> pairToEdge;
    const Index N = mesh.numNodes();
    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        pairToEdge[edge.n0 * N + edge.n1] = e;
    }
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& cell = mesh.getCell(c);
        if (cell.node_ids.size() < 3) continue;
        const Index nids[3] = {
            cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};
        for (int k = 0; k < 3; ++k) {
            Index a = nids[k];
            Index b = nids[(k + 1) % 3];
            if (a > b) std::swap(a, b);
            auto it = pairToEdge.find(a * N + b);
            if (it != pairToEdge.end())
                edgeCells[it->second].push_back(c);
        }
    }
    return edgeCells;
}

/// Build node -> adjacent cell ids map.
inline std::vector<std::vector<Index>> buildNodeCellMap(const DeviceMesh& mesh)
{
    std::vector<std::vector<Index>> nodeCells(mesh.numNodes());
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const Cell& cell = mesh.getCell(c);
        for (Index nodeId : cell.node_ids)
            nodeCells[nodeId].push_back(c);
    }
    return nodeCells;
}

// ---------------------------------------------------------------------------
// Material helpers
// ---------------------------------------------------------------------------

/// Return the average value of a material property over edge-adjacent cells.
/// Falls back to @p fallback only when the edge has no adjacent cells.
/// Throws if any adjacent cell references an unknown material.
inline Real edgeAvgMaterialProp(
    const std::vector<Index>& cells,
    const DeviceMesh&          mesh,
    const MaterialDatabase&    matdb,
    Real Material::*           prop,
    Real                       fallback)
{
    if (cells.empty()) return fallback;
    Real sum = 0.0;
    for (Index c : cells) {
        const auto& region = mesh.getRegion(mesh.getCell(c).region_id);
        sum += matdb.getMaterial(region.material).*prop;
    }
    return sum / static_cast<Real>(cells.size());
}


/// Build one temperature-adjusted material per mesh cell for hot-path reuse.
inline std::vector<Material> buildCellMaterials(const DeviceMesh&       mesh,
                                                const MaterialDatabase& matdb,
                                                Real                    temperature_K)
{
    std::vector<Material> materials;
    materials.reserve(mesh.numCells());
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& region = mesh.getRegion(mesh.getCell(c).region_id);
        materials.push_back(matdb.getMaterial(region.material, temperature_K));
    }
    return materials;
}


inline Real cellCentroidPotential(const DeviceMesh& mesh, const VectorXd& psi, Index cellId)
{
    const Cell& cell = mesh.getCell(cellId);
    if (cell.node_ids.empty())
        return std::numeric_limits<Real>::quiet_NaN();

    Real sum = 0.0;
    for (Index nodeId : cell.node_ids)
        sum += psi(static_cast<int>(nodeId));
    return sum / static_cast<Real>(cell.node_ids.size());
}

inline std::pair<Real, Real> cellCentroid(const DeviceMesh& mesh, Index cellId)
{
    const Cell& cell = mesh.getCell(cellId);
    if (cell.node_ids.empty())
        return {std::numeric_limits<Real>::quiet_NaN(),
                std::numeric_limits<Real>::quiet_NaN()};

    Real x = 0.0;
    Real y = 0.0;
    for (Index nodeId : cell.node_ids) {
        const Node& node = mesh.getNode(nodeId);
        x += node.x;
        y += node.y;
    }
    const Real invCount = 1.0 / static_cast<Real>(cell.node_ids.size());
    return {x * invCount, y * invCount};
}

inline Real estimateSurfaceNormalField(const std::vector<Index>& cells,
                                       const DeviceMesh& mesh,
                                       const VectorXd& psi,
                                       Index edgeId,
                                       Index cellId)
{
    const Edge& edge = mesh.getEdge(edgeId);
    if (edge.length <= 1.0e-30)
        return std::numeric_limits<Real>::quiet_NaN();

    const Node& n0 = mesh.getNode(edge.n0);
    const Node& n1 = mesh.getNode(edge.n1);
    const Real normalX = -(n1.y - n0.y) / edge.length;
    const Real normalY =  (n1.x - n0.x) / edge.length;

    const auto [cx, cy] = cellCentroid(mesh, cellId);
    const Real cellPhi = cellCentroidPotential(mesh, psi, cellId);
    if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(cellPhi))
        return std::numeric_limits<Real>::quiet_NaN();

    Real maxField = std::numeric_limits<Real>::quiet_NaN();
    for (Index otherCellId : cells) {
        if (otherCellId == cellId)
            continue;
        const auto [ox, oy] = cellCentroid(mesh, otherCellId);
        const Real otherPhi = cellCentroidPotential(mesh, psi, otherCellId);
        if (!std::isfinite(ox) || !std::isfinite(oy) || !std::isfinite(otherPhi))
            continue;
        const Real normalDistance = std::abs((ox - cx) * normalX + (oy - cy) * normalY);
        if (normalDistance <= 1.0e-30)
            continue;
        const Real field = std::abs((otherPhi - cellPhi) / normalDistance);
        if (!std::isfinite(maxField) || field > maxField)
            maxField = field;
    }
    if (std::isfinite(maxField))
        return maxField;

    const Real edgePhi = 0.5 * (psi(static_cast<int>(edge.n0)) + psi(static_cast<int>(edge.n1)));
    const Real mx = 0.5 * (n0.x + n1.x);
    const Real my = 0.5 * (n0.y + n1.y);
    const Real normalDistance = std::abs((cx - mx) * normalX + (cy - my) * normalY);
    if (normalDistance <= 1.0e-30)
        return std::numeric_limits<Real>::quiet_NaN();
    return std::abs((cellPhi - edgePhi) / normalDistance);
}

/// Return average model mobility [m^2/V/s] for edge @p edgeId.
inline Real edgeMobility(const std::vector<std::vector<Index>>& edgeCells,
                         const DeviceMesh&                       mesh,
                         const DopingModel&                      doping,
                         const MobilityModel&                    mobility,
                         const std::vector<Material>&            cellMaterials,
                         Index                                   edgeId,
                         CarrierType                             carrier,
                         Real                                    electricField,
                         const MobilityModelConfig*              mobilityConfig = nullptr,
                         const VectorXd*                         psi = nullptr)
{
    const auto& cells = edgeCells[edgeId];
    if (cells.empty()) return 0.0;

    const Edge& edge = mesh.getEdge(edgeId);
    const Real netDoping = 0.5 * (doping.netDoping(edge.n0) +
                                  doping.netDoping(edge.n1));

    const bool surfaceEnabled =
        mobilityConfig != nullptr && isSurfaceMobilityModel(*mobilityConfig);
    std::vector<std::string> adjacentRegionNames;
    if (surfaceEnabled) {
        adjacentRegionNames.reserve(cells.size());
        for (Index c : cells)
            adjacentRegionNames.push_back(mesh.getRegion(mesh.getCell(c).region_id).name);
    }

    Real sum = 0.0;
    Index contributingCells = 0;
    for (Index c : cells) {
        const Material& material = cellMaterials.at(static_cast<std::size_t>(c));
        const Real baseMobility = (carrier == CarrierType::Electron) ? material.mun : material.mup;
        if (baseMobility <= 0.0)
            continue;

        // Average only transport-capable cells.  This keeps oxide-only
        // edges pinned while preserving lateral semiconductor transport on
        // edges that lie along a semiconductor/oxide interface. Surface
        // mobility is enabled only on configured regions/interfaces; when no
        // normal-field estimate is available the NaN field disables the surface
        // factor while preserving any high-field velocity saturation.
        const Region& region = mesh.getRegion(mesh.getCell(c).region_id);
        const bool surfaceApplies = surfaceEnabled &&
            surfaceMobilityAppliesToRegionPair(*mobilityConfig, region.name, adjacentRegionNames);
        const Real surfaceNormalField = (surfaceApplies && psi != nullptr)
            ? estimateSurfaceNormalField(cells, mesh, *psi, edgeId, c)
            : std::numeric_limits<Real>::quiet_NaN();
        const Real modelMobility = (carrier == CarrierType::Electron)
            ? mobility.electronMobility(
                material, netDoping, 0.0, 0.0, electricField, surfaceNormalField)
            : mobility.holeMobility(
                material, netDoping, 0.0, 0.0, electricField, surfaceNormalField);
        if (modelMobility <= 0.0)
            continue;

        sum += modelMobility;
        ++contributingCells;
    }
    if (contributingCells == 0)
        return 0.0;
    return sum / static_cast<Real>(contributingCells);
}

/// Return average model mobility [m^2/V/s] over semiconductor cells adjacent to a node.
inline Real nodeMobility(const std::vector<std::vector<Index>>& nodeCells,
                         const DeviceMesh&                     mesh,
                         const DopingModel&                    doping,
                         const MobilityModel&                  mobility,
                         const std::vector<Material>&          cellMaterials,
                         Index                                 nodeId,
                         CarrierType                           carrier,
                         Real                                  drivingField)
{
    const auto& cells = nodeCells[nodeId];
    if (cells.empty())
        return 0.0;

    Real sum = 0.0;
    Index contributingCells = 0;
    for (Index c : cells) {
        const Material& material = cellMaterials.at(static_cast<std::size_t>(c));
        const Real baseMobility = (carrier == CarrierType::Electron) ? material.mun : material.mup;
        if (baseMobility <= 0.0)
            continue;
        const Real modelMobility = (carrier == CarrierType::Electron)
            ? mobility.electronMobility(
                material, doping.netDoping(nodeId), 0.0, 0.0, drivingField)
            : mobility.holeMobility(
                material, doping.netDoping(nodeId), 0.0, 0.0, drivingField);
        if (modelMobility <= 0.0)
            continue;
        sum += modelMobility;
        ++contributingCells;
    }
    if (contributingCells == 0)
        return 0.0;
    return sum / static_cast<Real>(contributingCells);
}

inline Real interpolatedAvalancheDrivingField(const ImpactIonizationModelConfig& config,
                                              Real                               drivingField,
                                              Real                               electricField,
                                              Real                               carrierDensity,
                                              Real                               referenceDensity)
{
    if (config.drivingForceInterpolation != "quasi_fermi_to_electric_field" ||
        referenceDensity <= 0.0) {
        return drivingField;
    }
    const Real carrier = std::max(carrierDensity, 0.0);
    const Real weight = carrier / (carrier + referenceDensity);
    return weight * drivingField + (1.0 - weight) * electricField;
}

inline Real electronAvalancheDrivingField(const ImpactIonizationModelConfig& config,
                                          Real                               drivingField,
                                          Real                               electricField,
                                          Real                               electronDensity)
{
    return interpolatedAvalancheDrivingField(
        config,
        drivingField,
        electricField,
        electronDensity,
        config.electronDrivingForceRefDensity);
}

inline Real holeAvalancheDrivingField(const ImpactIonizationModelConfig& config,
                                      Real                               drivingField,
                                      Real                               electricField,
                                      Real                               holeDensity)
{
    return interpolatedAvalancheDrivingField(
        config,
        drivingField,
        electricField,
        holeDensity,
        config.holeDrivingForceRefDensity);
}

inline bool usesCurrentAlignedAvalancheDrivingForce(
    const ImpactIonizationModelConfig& config)
{
    return config.drivingForce == "grad_potential_parallel_j" ||
           config.drivingForce == "effective_field_parallel_j";
}

inline Real parallelCurrentAvalancheDrivingField(Real signedDrivingField,
                                                Real signedCurrentProxy)
{
    if (!std::isfinite(signedDrivingField) || !std::isfinite(signedCurrentProxy) ||
        std::abs(signedCurrentProxy) <= 0.0) {
        return 0.0;
    }
    const Real currentSign = signedCurrentProxy > 0.0 ? 1.0 : -1.0;
    return std::max(signedDrivingField * currentSign, 0.0);
}

/// Resolves the SG edge-current avalanche source-volume factor used in
/// `factor * h * edge.couple`. A finite `source_volume_factor` overrides the
/// named `source_volume_policy` preset; `0` falls back to the preset.
inline Real avalancheSourceVolumeFactor(const ImpactIonizationModelConfig& config)
{
    if (config.sourceVolumeFactor > 0.0)
        return config.sourceVolumeFactor;
    return config.sourceVolumePolicy == "edge_box" ? 1.0 : 0.5;
}

/// Validates the impact-ionization configuration shared by the Gummel and
/// Newton solver config loaders. `context` is prefixed to any thrown message.
inline void validateImpactIonizationDrivingForce(const ImpactIonizationModelConfig& config,
                                                 const char* context)
{
    const bool currentAlignedDrivingForce = usesCurrentAlignedAvalancheDrivingForce(config);
    if (config.drivingForce != "electric_field" &&
        config.drivingForce != "quasi_fermi_gradient" &&
        !currentAlignedDrivingForce) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.driving_force must be 'electric_field', "
            "'quasi_fermi_gradient', 'grad_potential_parallel_j', or "
            "'effective_field_parallel_j'.");
    }
    if (config.generation != "carrier_density" &&
        config.generation != "current_density") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.generation must be 'carrier_density' or "
            "'current_density'.");
    }
    if (config.currentApproximation != "mobility_density_gradient" &&
        config.currentApproximation != "density_gradient" &&
        config.currentApproximation != "grad_qf") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.current_approximation must be "
            "'mobility_density_gradient', 'density_gradient', or 'grad_qf'.");
    }
    if (config.drivingForceInterpolation != "none" &&
        config.drivingForceInterpolation != "quasi_fermi_to_electric_field") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.driving_force_interpolation.mode must be "
            "'none' or 'quasi_fermi_to_electric_field'.");
    }
    if (currentAlignedDrivingForce &&
        (config.generation != "current_density" ||
         (config.currentApproximation != "density_gradient" &&
          config.currentApproximation != "grad_qf"))) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization current-aligned driving forces require "
            "generation='current_density' with current_approximation='density_gradient' "
            "or 'grad_qf'.");
    }
    if (config.drivingForceInterpolation != "none" &&
        config.drivingForce != "quasi_fermi_gradient") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.driving_force_interpolation requires "
            "driving_force='quasi_fermi_gradient'.");
    }
    if (!std::isfinite(config.electronDrivingForceRefDensity) ||
        !std::isfinite(config.holeDrivingForceRefDensity) ||
        config.electronDrivingForceRefDensity < 0.0 ||
        config.holeDrivingForceRefDensity < 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization driving-force reference densities must be "
            "finite and non-negative.");
    }
    if (!std::isfinite(config.sourceGeometryScale) ||
        config.sourceGeometryScale <= 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.source_geometry_scale must be positive and finite.");
    }
    if (config.sourceVolumePolicy != "edge_half_box" &&
        config.sourceVolumePolicy != "edge_box") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.source_volume_policy must be 'edge_half_box' or 'edge_box'.");
    }
    if (config.sourceVolumeFactor != 0.0 &&
        (!std::isfinite(config.sourceVolumeFactor) ||
         config.sourceVolumeFactor < 0.5 ||
         config.sourceVolumeFactor > 1.0)) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.source_volume_factor must be 0 or within [0.5, 1.0].");
    }
    if (!std::isfinite(config.quasiFermiCarrierTruncation) ||
        config.quasiFermiCarrierTruncation < 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.quasi_fermi_carrier_truncation must be non-negative and finite.");
    }
    if (!std::isfinite(config.minimumField) || config.minimumField < 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.minimum_field_V_m must be non-negative and finite.");
    }
}

inline bool usesDensityGradientAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "density_gradient";
}

inline bool usesEdgeCurrentAvalancheSource(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           (config.currentApproximation == "density_gradient" ||
            config.currentApproximation == "grad_qf");
}

inline bool usesQuasiFermiCarrierTruncation(const ImpactIonizationModelConfig& config)
{
    return config.quasiFermiCarrierTruncation > 0.0;
}

inline std::vector<bool> contactNodeMask(const DeviceMesh& mesh)
{
    std::vector<bool> mask(static_cast<std::size_t>(mesh.numNodes()), false);
    for (const Contact& contact : mesh.contacts()) {
        for (Index nodeId : contact.node_ids) {
            if (nodeId < mesh.numNodes())
                mask[static_cast<std::size_t>(nodeId)] = true;
        }
    }
    return mask;
}

inline bool edgeTouchesContactElement(const DeviceMesh& mesh,
                                      const std::vector<std::vector<Index>>& edgeCells,
                                      Index edgeId,
                                      const std::vector<bool>& contactNodes)
{
    const Edge& edge = mesh.getEdge(edgeId);
    if (contactNodes[static_cast<std::size_t>(edge.n0)] ||
        contactNodes[static_cast<std::size_t>(edge.n1)]) {
        return true;
    }
    if (edgeId >= edgeCells.size())
        return false;
    for (Index cellId : edgeCells[edgeId]) {
        const Cell& cell = mesh.getCell(cellId);
        for (Index nodeId : cell.node_ids) {
            if (nodeId < mesh.numNodes() &&
                contactNodes[static_cast<std::size_t>(nodeId)]) {
                return true;
            }
        }
    }
    return false;
}

inline Real edgeHighFieldDrivingField(bool qfDrivingForce,
                                      Real qfField,
                                      Real electricField,
                                      const std::vector<std::vector<Index>>& edgeCells,
                                      const DeviceMesh& mesh,
                                      Index edgeId,
                                      const std::vector<bool>& contactNodes)
{
    if (!qfDrivingForce)
        return electricField;
    if (edgeTouchesContactElement(mesh, edgeCells, edgeId, contactNodes))
        return electricField;
    return qfField;
}

inline Real electronQfForAvalancheGradient(Real psi,
                                           Real phin,
                                           Real electronDensity,
                                           Real intrinsicDensity,
                                           Real Vt,
                                           const ImpactIonizationModelConfig& config)
{
    if (!usesQuasiFermiCarrierTruncation(config) || intrinsicDensity <= 0.0)
        return phin;
    const Real carrier = std::max(
        std::max(electronDensity, 0.0),
        config.quasiFermiCarrierTruncation * intrinsicDensity);
    return psi - Vt * std::log(carrier / intrinsicDensity);
}

inline Real holeQfForAvalancheGradient(Real psi,
                                       Real phip,
                                       Real holeDensity,
                                       Real intrinsicDensity,
                                       Real Vt,
                                       const ImpactIonizationModelConfig& config)
{
    if (!usesQuasiFermiCarrierTruncation(config) || intrinsicDensity <= 0.0)
        return phip;
    const Real carrier = std::max(
        std::max(holeDensity, 0.0),
        config.quasiFermiCarrierTruncation * intrinsicDensity);
    return psi + Vt * std::log(carrier / intrinsicDensity);
}

struct EdgeAvalancheDirectionalWeights {
    Real electronNode0 = 0.5;
    Real electronNode1 = 0.5;
    Real holeNode0 = 0.5;
    Real holeNode1 = 0.5;
};

struct CellScalarGradientCache {
    std::vector<Point2> gradients;
    std::vector<Real> areas;
    std::vector<bool> valid;
};

template <typename ValueAt>
inline Point2 cellScalarGradient(
    const DeviceMesh& mesh,
    const Cell&       cell,
    ValueAt&&         valueAt,
    bool&             valid,
    Real&             area)
{
    valid = false;
    area = 0.0;
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
        return Point2::Zero();

    const Index n0 = cell.node_ids[0];
    const Index n1 = cell.node_ids[1];
    const Index n2 = cell.node_ids[2];
    const Node& p0 = mesh.getNode(n0);
    const Node& p1 = mesh.getNode(n1);
    const Node& p2 = mesh.getNode(n2);
    const Real dx10 = p1.x - p0.x;
    const Real dy10 = p1.y - p0.y;
    const Real dx20 = p2.x - p0.x;
    const Real dy20 = p2.y - p0.y;
    const Real det = dx10 * dy20 - dy10 * dx20;
    if (std::abs(det) <= 1.0e-300)
        return Point2::Zero();

    const Real dv10 = valueAt(n1) - valueAt(n0);
    const Real dv20 = valueAt(n2) - valueAt(n0);
    valid = true;
    area = 0.5 * std::abs(det);
    return Point2{
        (dv10 * dy20 - dv20 * dy10) / det,
        (dx10 * dv20 - dx20 * dv10) / det,
    };
}

template <typename ValueAt>
inline CellScalarGradientCache computeCellScalarGradientCache(
    const DeviceMesh& mesh,
    ValueAt&&         valueAt)
{
    CellScalarGradientCache cache;
    cache.gradients.assign(static_cast<std::size_t>(mesh.numCells()), Point2::Zero());
    cache.areas.assign(static_cast<std::size_t>(mesh.numCells()), 0.0);
    cache.valid.assign(static_cast<std::size_t>(mesh.numCells()), false);

    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        bool valid = false;
        Real area = 0.0;
        cache.gradients[cellId] = cellScalarGradient(
            mesh, mesh.getCell(cellId), valueAt, valid, area);
        cache.areas[cellId] = area;
        cache.valid[cellId] = valid;
    }
    return cache;
}

template <typename ValueAt>
inline Point2 edgeAveragedCellScalarGradient(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    ValueAt&&                              valueAt,
    bool&                                  valid)
{
    valid = false;
    if (edgeId >= edgeCells.size())
        return Point2::Zero();

    Point2 weightedGradient = Point2::Zero();
    Real totalArea = 0.0;
    for (const Index cellId : edgeCells[edgeId]) {
        bool cellValid = false;
        Real area = 0.0;
        const Point2 gradient = cellScalarGradient(
            mesh, mesh.getCell(cellId), valueAt, cellValid, area);
        if (!cellValid || area <= 0.0)
            continue;
        weightedGradient += area * gradient;
        totalArea += area;
    }

    if (totalArea <= 0.0)
        return Point2::Zero();
    valid = true;
    return weightedGradient / totalArea;
}

inline Point2 edgeAveragedCellScalarGradient(
    const std::vector<std::vector<Index>>& edgeCells,
    Index                                  edgeId,
    const CellScalarGradientCache&         cache,
    bool&                                  valid)
{
    valid = false;
    if (edgeId >= edgeCells.size())
        return Point2::Zero();

    Point2 weightedGradient = Point2::Zero();
    Real totalArea = 0.0;
    for (const Index cellId : edgeCells[edgeId]) {
        if (cellId >= cache.valid.size() || !cache.valid[cellId])
            continue;
        const Real area = cache.areas[cellId];
        if (area <= 0.0)
            continue;
        weightedGradient += area * cache.gradients[cellId];
        totalArea += area;
    }

    if (totalArea <= 0.0)
        return Point2::Zero();
    valid = true;
    return weightedGradient / totalArea;
}

template <typename ValueAt>
inline Real edgeMinusGradientUnitDot(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    ValueAt&&                              valueAt)
{
    const Edge& edge = mesh.getEdge(edgeId);
    const Node& n0 = mesh.getNode(edge.n0);
    const Node& n1 = mesh.getNode(edge.n1);
    if (edge.length <= 1.0e-30)
        return 0.0;

    bool validGradient = false;
    const Point2 gradient = edgeAveragedCellScalarGradient(
        edgeCells, mesh, edgeId, valueAt, validGradient);
    if (!validGradient)
        return 0.0;

    const Point2 minusGradient = -gradient;
    const Real gradientNorm = minusGradient.norm();
    if (gradientNorm <= 1.0e-300)
        return 0.0;

    const Point2 edgeUnit{(n1.x - n0.x) / edge.length, (n1.y - n0.y) / edge.length};
    return std::clamp(edgeUnit.dot(minusGradient / gradientNorm), -1.0, 1.0);
}

inline Real edgeMinusGradientUnitDot(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    const CellScalarGradientCache&         cache)
{
    const Edge& edge = mesh.getEdge(edgeId);
    const Node& n0 = mesh.getNode(edge.n0);
    const Node& n1 = mesh.getNode(edge.n1);
    if (edge.length <= 1.0e-30)
        return 0.0;

    bool validGradient = false;
    const Point2 gradient = edgeAveragedCellScalarGradient(
        edgeCells, edgeId, cache, validGradient);
    if (!validGradient)
        return 0.0;

    const Point2 minusGradient = -gradient;
    const Real gradientNorm = minusGradient.norm();
    if (gradientNorm <= 1.0e-300)
        return 0.0;

    const Point2 edgeUnit{(n1.x - n0.x) / edge.length, (n1.y - n0.y) / edge.length};
    return std::clamp(edgeUnit.dot(minusGradient / gradientNorm), -1.0, 1.0);
}

inline EdgeAvalancheDirectionalWeights edgeAvalancheDirectionalWeights(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    const CellScalarGradientCache&         electronGradientCache,
    const CellScalarGradientCache&         holeGradientCache)
{
    EdgeAvalancheDirectionalWeights weights;
    const Real electronDot = edgeMinusGradientUnitDot(
        edgeCells, mesh, edgeId, electronGradientCache);
    const Real holeDot = edgeMinusGradientUnitDot(
        edgeCells, mesh, edgeId, holeGradientCache);

    weights.electronNode0 = std::clamp(0.5 + 0.5 * electronDot, 0.0, 1.0);
    weights.electronNode1 = 1.0 - weights.electronNode0;
    weights.holeNode1 = std::clamp(0.5 + 0.5 * holeDot, 0.0, 1.0);
    weights.holeNode0 = 1.0 - weights.holeNode1;
    return weights;
}

template <typename ElectronQfAt, typename HoleQfAt>
inline EdgeAvalancheDirectionalWeights edgeAvalancheDirectionalWeights(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    ElectronQfAt&&                         electronQfAt,
    HoleQfAt&&                             holeQfAt)
{
    EdgeAvalancheDirectionalWeights weights;
    const Real electronDot = edgeMinusGradientUnitDot(
        edgeCells, mesh, edgeId, electronQfAt);
    const Real holeDot = edgeMinusGradientUnitDot(
        edgeCells, mesh, edgeId, holeQfAt);

    weights.electronNode0 = std::clamp(0.5 + 0.5 * electronDot, 0.0, 1.0);
    weights.electronNode1 = 1.0 - weights.electronNode0;
    weights.holeNode1 = std::clamp(0.5 + 0.5 * holeDot, 0.0, 1.0);
    weights.holeNode0 = 1.0 - weights.holeNode1;
    return weights;
}

struct SgEdgeCurrentAvalancheSourceRecord {
    Index edgeId = 0;
    Index node0 = 0;
    Index node1 = 0;
    Real edgeLength = 0.0;
    Real edgeCouple = 0.0;
    Real edgeAreaProxy = 0.0;
    Real electricField = 0.0;
    Real electronImpactField = 0.0;
    Real holeImpactField = 0.0;
    Real electronAlpha = 0.0;
    Real holeAlpha = 0.0;
    Real electronMobility = 0.0;
    Real holeMobility = 0.0;
    Real electronFluxProxy = 0.0;
    Real holeFluxProxy = 0.0;
    Real electronSourceIntegral = 0.0;
    Real holeSourceIntegral = 0.0;
    Real edgeSourceIntegral = 0.0;
    Real electronNode0SourceIntegral = 0.0;
    Real electronNode1SourceIntegral = 0.0;
    Real holeNode0SourceIntegral = 0.0;
    Real holeNode1SourceIntegral = 0.0;
    Real node0SourceIntegral = 0.0;
    Real node1SourceIntegral = 0.0;
};

struct SgAvalancheSourceComponentIntegrals {
    std::vector<Real> electron;
    std::vector<Real> hole;
    std::vector<Real> combined;
};

inline std::vector<SgEdgeCurrentAvalancheSourceRecord> sgEdgeCurrentAvalancheSourceRecords(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&           ni,
    Real                               Vt)
{
    std::vector<SgEdgeCurrentAvalancheSourceRecord> records;
    records.reserve(mesh.numEdges());
    const bool qfImpact = config.drivingForce == "quasi_fermi_gradient";
    const bool currentAlignedImpact = usesCurrentAlignedAvalancheDrivingForce(config);
    const bool qfMobility = mobilityConfig.highFieldDrivingForce == "quasi_fermi_gradient";
    const std::vector<bool> contactNodes = contactNodeMask(mesh);
    const CellScalarGradientCache electronQfGradientCache = computeCellScalarGradientCache(
        mesh, [&](Index node) {
            const int idx = static_cast<int>(node);
            return electronQfForAvalancheGradient(
                psi(idx), phin(idx), n(idx), ni[node], Vt, config);
        });
    const CellScalarGradientCache holeQfGradientCache = computeCellScalarGradientCache(
        mesh, [&](Index node) {
            const int idx = static_cast<int>(node);
            return holeQfForAvalancheGradient(
                psi(idx), phip(idx), p(idx), ni[node], Vt, config);
        });

    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        const Real h = edge.length;
        if (h <= 1.0e-30 || edge.couple <= 0.0)
            continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real psi_i = psi(i);
        const Real psi_j = psi(j);
        const Real phin_i = phin(i);
        const Real phin_j = phin(j);
        const Real phip_i = phip(i);
        const Real phip_j = phip(j);

        const Real electronQf_i = electronQfForAvalancheGradient(
            psi_i, phin_i, n(i), ni[edge.n0], Vt, config);
        const Real electronQf_j = electronQfForAvalancheGradient(
            psi_j, phin_j, n(j), ni[edge.n1], Vt, config);
        const Real holeQf_i = holeQfForAvalancheGradient(
            psi_i, phip_i, p(i), ni[edge.n0], Vt, config);
        const Real holeQf_j = holeQfForAvalancheGradient(
            psi_j, phip_j, p(j), ni[edge.n1], Vt, config);
        const Real electricField = std::abs((psi_j - psi_i) / h);
        const Real electronQfField = std::abs((electronQf_j - electronQf_i) / h);
        const Real holeQfField = std::abs((holeQf_j - holeQf_i) / h);
        const Real electronCoefficientField = edgeHighFieldDrivingField(
            qfImpact, electronQfField, electricField, edgeCells, mesh, e, contactNodes);
        const Real holeCoefficientField = edgeHighFieldDrivingField(
            qfImpact, holeQfField, electricField, edgeCells, mesh, e, contactNodes);
        const Real electronMobilityField = qfMobility ? electronQfField : electricField;
        const Real holeMobilityField = qfMobility ? holeQfField : electricField;

        const Real nAvg = 0.5 * (n(i) + n(j));
        const Real pAvg = 0.5 * (p(i) + p(j));
        const Real signedElectricField01 = -(psi_j - psi_i) / h;

        const Real sourceVolumeFactor = avalancheSourceVolumeFactor(config);
        const Real edgeArea = sourceVolumeFactor * h * edge.couple * config.sourceGeometryScale;
        SgEdgeCurrentAvalancheSourceRecord record;
        record.edgeId = e;
        record.node0 = edge.n0;
        record.node1 = edge.n1;
        record.edgeLength = h;
        record.edgeCouple = edge.couple;
        record.edgeAreaProxy = edgeArea;
        record.electricField = electricField;

        const Real mun = edgeMobility(
            edgeCells, mesh, doping, mobility, cellMaterials, e, CarrierType::Electron,
            electronMobilityField, &mobilityConfig, &psi);
        record.electronMobility = mun;
        if (mun > 0.0) {
            const Real electronContinuityFlux01 =
                sgElectronContinuityFluxFromQuasiFermiVariableNi(
                    ni[edge.n0],
                    ni[edge.n1],
                    psi_i,
                    psi_j,
                    phin_i,
                    phin_j,
                    Vt,
                    mun * Vt / h);
            record.electronFluxProxy = std::abs(electronContinuityFlux01);
            record.electronImpactField = currentAlignedImpact
                ? parallelCurrentAvalancheDrivingField(
                    signedElectricField01, electronContinuityFlux01)
                : electronAvalancheDrivingField(
                    config, electronCoefficientField, electricField, nAvg);
            record.electronAlpha = impact.electronCoefficient(record.electronImpactField);
            record.electronSourceIntegral =
                record.electronAlpha * record.electronFluxProxy * edgeArea;
            record.edgeSourceIntegral += record.electronSourceIntegral;
        }

        const Real mup = edgeMobility(
            edgeCells, mesh, doping, mobility, cellMaterials, e, CarrierType::Hole,
            holeMobilityField, &mobilityConfig, &psi);
        record.holeMobility = mup;
        if (mup > 0.0) {
            const Real holeContinuityFlux01 =
                sgHoleContinuityFluxFromQuasiFermiVariableNi(
                    ni[edge.n0],
                    ni[edge.n1],
                    psi_i,
                    psi_j,
                    phip_i,
                    phip_j,
                    Vt,
                    mup * Vt / h);
            record.holeFluxProxy = std::abs(holeContinuityFlux01);
            record.holeImpactField = currentAlignedImpact
                ? parallelCurrentAvalancheDrivingField(
                    signedElectricField01, holeContinuityFlux01)
                : holeAvalancheDrivingField(
                    config, holeCoefficientField, electricField, pAvg);
            record.holeAlpha = impact.holeCoefficient(record.holeImpactField);
            record.holeSourceIntegral =
                record.holeAlpha * record.holeFluxProxy * edgeArea;
            record.edgeSourceIntegral += record.holeSourceIntegral;
        }

        const EdgeAvalancheDirectionalWeights weights = edgeAvalancheDirectionalWeights(
            edgeCells,
            mesh,
            e,
            electronQfGradientCache,
            holeQfGradientCache);
        record.electronNode0SourceIntegral =
            weights.electronNode0 * record.electronSourceIntegral;
        record.electronNode1SourceIntegral =
            weights.electronNode1 * record.electronSourceIntegral;
        record.holeNode0SourceIntegral =
            weights.holeNode0 * record.holeSourceIntegral;
        record.holeNode1SourceIntegral =
            weights.holeNode1 * record.holeSourceIntegral;
        record.node0SourceIntegral =
            record.electronNode0SourceIntegral + record.holeNode0SourceIntegral;
        record.node1SourceIntegral =
            record.electronNode1SourceIntegral + record.holeNode1SourceIntegral;
        records.push_back(record);
    }
    return records;
}

inline SgAvalancheSourceComponentIntegrals sgEdgeCurrentAvalancheSourceComponentIntegrals(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&           ni,
    Real                               Vt)
{
    SgAvalancheSourceComponentIntegrals source;
    source.electron.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    source.hole.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    source.combined.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    const auto records = sgEdgeCurrentAvalancheSourceRecords(
        config,
        impact,
        mobilityConfig,
        mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        psi,
        phin,
        phip,
        n,
        p,
        ni,
        Vt);
    for (const auto& record : records) {
        source.electron[record.node0] += record.electronNode0SourceIntegral;
        source.electron[record.node1] += record.electronNode1SourceIntegral;
        source.hole[record.node0] += record.holeNode0SourceIntegral;
        source.hole[record.node1] += record.holeNode1SourceIntegral;
        source.combined[record.node0] += record.node0SourceIntegral;
        source.combined[record.node1] += record.node1SourceIntegral;
    }
    return source;
}

inline std::vector<Real> sgEdgeCurrentAvalancheSourceIntegrals(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&           ni,
    Real                               Vt)
{
    std::vector<Real> source(mesh.numNodes(), 0.0);
    const auto records = sgEdgeCurrentAvalancheSourceRecords(
        config,
        impact,
        mobilityConfig,
        mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        psi,
        phin,
        phip,
        n,
        p,
        ni,
        Vt);
    for (const auto& record : records) {
        source[record.node0] += record.node0SourceIntegral;
        source[record.node1] += record.node1SourceIntegral;
    }
    return source;
}

inline Real impactIonizationGenerationRate(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& nodeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    Index                              nodeId,
    Real                               electricField,
    Real                               electronDrivingField,
    Real                               holeDrivingField,
    Real                               n,
    Real                               p)
{
    if (config.generation != "current_density")
        return impact.generationRate(electricField, n, p);

    const Real electronImpactField = electronAvalancheDrivingField(
        config, electronDrivingField, electricField, n);
    const Real holeImpactField = holeAvalancheDrivingField(
        config, holeDrivingField, electricField, p);
    const Real alphaN = impact.electronCoefficient(electronImpactField);
    const Real alphaP = impact.holeCoefficient(holeImpactField);
    const Real mun = nodeMobility(
        nodeCells, mesh, doping, mobility, cellMaterials, nodeId, CarrierType::Electron,
        electronImpactField);
    const Real mup = nodeMobility(
        nodeCells, mesh, doping, mobility, cellMaterials, nodeId, CarrierType::Hole,
        holeImpactField);
    return alphaN * mun * std::max(n, 0.0) * std::abs(electronImpactField) +
           alphaP * mup * std::max(p, 0.0) * std::abs(holeImpactField);
}

/// Return average dielectric constant [F/m] for edge @p edgeId.
inline Real edgeEpsilon(const std::vector<std::vector<Index>>& edgeCells,
                        const DeviceMesh&                       mesh,
                        const MaterialDatabase&                 matdb,
                        Index                                   edgeId)
{
    return edgeAvgMaterialProp(edgeCells[edgeId], mesh, matdb,
                               &Material::eps_r, 1.0) * constants::eps0;
}

// ---------------------------------------------------------------------------
// Per-node ni vector
// ---------------------------------------------------------------------------

/// Build per-node intrinsic concentration vector from the material database.
/// For interface nodes shared by multiple regions, uses the first-found value.
inline std::vector<Real> buildNodeNi(const DeviceMesh&       mesh,
                                     const MaterialDatabase& matdb,
                                     Real                    temperature_K = constants::T0)
{
    const Index N = mesh.numNodes();
    std::vector<Real> ni_v(N, 0.0);
    std::vector<bool> found(N, false);
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& cell   = mesh.getCell(c);
        const auto& region = mesh.getRegion(cell.region_id);
        const Real ni_mat = matdb.getMaterial(region.material, temperature_K).ni;
        for (Index nid : cell.node_ids) {
            if (!found[nid]) {
                ni_v[nid]  = ni_mat;
                found[nid] = true;
            }
        }
    }
    return ni_v;
}

/// Validate that the doping model has one entry per mesh node.
inline void validateDopingMeshSize(const DeviceMesh& mesh,
                                   const DopingModel& doping,
                                   const std::string& context)
{
    if (doping.numNodes() != mesh.numNodes())
        throw std::invalid_argument(
            context + ": doping model size does not match mesh node count.");
}

/// Build per-node effective intrinsic concentration including bandgap narrowing.
inline std::vector<Real> buildEffectiveNodeNi(const DeviceMesh&       mesh,
                                              const MaterialDatabase& matdb,
                                              const DopingModel&      doping,
                                              const BandgapNarrowing& bgn,
                                              Real                    thermalVoltage)
{
    const Real temperature_K = thermalVoltage * constants::q / constants::kb;
    std::vector<Real> ni_v = buildNodeNi(mesh, matdb, temperature_K);
    for (Index i = 0; i < mesh.numNodes(); ++i) {
        const Real delta = bgn.deltaEg(doping.totalImpurity(i), 0.0, 0.0);
        ni_v[i] = effectiveIntrinsicDensity(ni_v[i], thermalVoltage, delta);
    }
    return ni_v;
}

/// Validate inputs before building effective intrinsic concentrations.
inline std::vector<Real> buildValidatedEffectiveNodeNi(
    const std::string&             context,
    const DeviceMesh&              mesh,
    const MaterialDatabase&        matdb,
    const DopingModel&             doping,
    const BandgapNarrowingConfig&  bandgapNarrowingConfig,
    Real                           thermalVoltage)
{
    validateDopingMeshSize(mesh, doping, context);
    return buildEffectiveNodeNi(
        mesh,
        matdb,
        doping,
        *makeBandgapNarrowingModel(bandgapNarrowingConfig),
        thermalVoltage);
}

// ---------------------------------------------------------------------------
// Dirichlet boundary conditions
// ---------------------------------------------------------------------------

/**
 * @brief Apply strong Dirichlet BCs (row-replacement) to a sparse system.
 *
 * For each constrained node i with prescribed value v:
 *   1. Reduce rhs(k) by A(k,i)*v for all free rows k.
 *   2. Zero column i and row i.
 *   3. Set A(i,i) = 1, rhs(i) = v.
 *
 * The diagonal A(i,i) = 1 is always set explicitly via coeffRef so that
 * the constraint is enforced even if node i had no prior stiffness entries
 * (e.g. insulator nodes in continuity assemblies where all adjacent edges
 * were skipped due to zero mobility).
 */
inline void applyDirichletBC(SparseMatrixd&                         A,
                              VectorXd&                              b,
                              const std::unordered_map<Index, Real>& bcs)
{
    A.makeCompressed();

    // Step 1: propagate prescribed values into free-node RHS
    for (const auto& [nodeId, value] : bcs) {
        const int i = static_cast<int>(nodeId);
        for (SparseMatrixd::InnerIterator it(A, i); it; ++it) {
            const int k = static_cast<int>(it.row());
            if (k == i) continue;
            if (bcs.count(static_cast<Index>(k)) == 0)
                b(k) -= it.value() * value;
        }
    }

    // Step 2 & 3: zero Dirichlet rows/cols
    std::unordered_set<int> dirichletSet;
    for (const auto& [nodeId, _] : bcs)
        dirichletSet.insert(static_cast<int>(nodeId));

    for (int col = 0; col < A.outerSize(); ++col) {
        const bool colIsDirichlet = dirichletSet.count(col) > 0;
        for (SparseMatrixd::InnerIterator it(A, col); it; ++it) {
            const int  row          = static_cast<int>(it.row());
            const bool rowIsDirichlet = dirichletSet.count(row) > 0;
            if (rowIsDirichlet || colIsDirichlet)
                it.valueRef() = (row == col) ? 1.0 : 0.0;
        }
    }
    A.prune(0.0);

    // Explicitly ensure A(i,i) = 1 even if the node had no prior entries.
    for (const auto& [nodeId, value] : bcs) {
        const int i = static_cast<int>(nodeId);
        A.coeffRef(i, i) = 1.0;
        b(i) = value;
    }
}

} // namespace vela::detail
