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
#include "vela/core/UnitScaling.h"
#include "vela/equation/ChargeSpec.h"
#include "vela/equation/Tri3LocalForwardAD.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/Material.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/BandToBandTunnelingModel.h"
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
#include <type_traits>

namespace vela::detail {

inline Real physicalPotentialCentralDifferenceStep(
    Real physicalValue, Real potentialScale, Real relativeStep)
{
    return relativeStep * std::max(potentialScale, std::abs(physicalValue));
}


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
    const std::string& context,
    Real chargeAreaFactor = 1.0,
    Real chargeLineFactor = 1.0)
{
    VectorXd contribution = VectorXd::Zero(static_cast<int>(mesh.numNodes()));

    const auto fixedByRegion = fixedChargeByRegion(fixedCharges, context);
    if (!fixedByRegion.empty()) {
        for (Index c = 0; c < mesh.numCells(); ++c) {
            const Cell& cell = mesh.getCell(c);
            const Region& region = mesh.getRegion(cell.region_id);
            auto it = fixedByRegion.find(region.name);
            if (it == fixedByRegion.end()) continue;

            const Real nodeCharge = constants::q * it->second * triangleArea(mesh, cell) *
                chargeAreaFactor / 3.0;
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
            const Real endpointCharge = constants::q * it->second * edge.length *
                chargeLineFactor * 0.5;
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
    const std::string& context,
    Real chargeAreaFactor = 1.0,
    Real chargeLineFactor = 1.0)
{
    rhs += computeFixedAndInterfaceChargeRhs(
        mesh, edgeCells, fixedCharges, sheetCharges, context,
        chargeAreaFactor, chargeLineFactor);
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

/// Return carrier-transport couplings, honoring explicit per-edge diagnostic
/// overrides while leaving the mesh electrostatic geometry untouched.
inline Real transportEdgeCouple(const Edge& edge)
{
    return edge.transport_couple >= 0.0
        ? edge.transport_couple
        : edge.couple;
}

inline std::vector<Real> computeTransportEdgeCouplings(const DeviceMesh& mesh)
{
    const Index E = mesh.numEdges();
    std::vector<Real> couple(E, 0.0);
    for (Index e = 0; e < E; ++e)
        couple[e] = transportEdgeCouple(mesh.getEdge(e));
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

/// Return a per-node Sentaurus-like cell-gradient electric-field magnitude [V/m].
inline std::vector<Real> computeNodeElectricFields(const VectorXd& psi, const DeviceMesh& mesh, Real fieldFactor = 1.0);

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

inline std::vector<std::vector<Index>> buildCellEdgeMap(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh)
{
    std::vector<std::vector<Index>> cellEdges(static_cast<std::size_t>(mesh.numCells()));
    for (Index edgeId = 0; edgeId < edgeCells.size(); ++edgeId) {
        for (Index cellId : edgeCells[edgeId]) {
            if (cellId < mesh.numCells())
                cellEdges[static_cast<std::size_t>(cellId)].push_back(edgeId);
        }
    }
    return cellEdges;
}

inline Real cellSmoothedEdgeFluxMagnitude(
    Index                                  edgeId,
    const std::vector<Real>&               rawEdgeFlux,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<std::vector<Index>>& cellEdges)
{
    if (edgeId >= edgeCells.size())
        return 0.0;
    Real edgeSum = 0.0;
    int adjacentCellCount = 0;
    for (Index cellId : edgeCells[edgeId]) {
        if (cellId >= cellEdges.size())
            continue;
        Real cellSum = 0.0;
        int cellEdgeCount = 0;
        for (Index otherEdgeId : cellEdges[cellId]) {
            if (otherEdgeId >= rawEdgeFlux.size())
                continue;
            cellSum += rawEdgeFlux[otherEdgeId];
            ++cellEdgeCount;
        }
        if (cellEdgeCount <= 0)
            continue;
        edgeSum += cellSum / static_cast<Real>(cellEdgeCount);
        ++adjacentCellCount;
    }
    if (adjacentCellCount <= 0)
        return edgeId < rawEdgeFlux.size() ? rawEdgeFlux[edgeId] : 0.0;
    return edgeSum / static_cast<Real>(adjacentCellCount);
}

inline Point2 cellVectorCurrent(
    Index                                  cellId,
    const std::vector<Real>&               signedEdgeFlux,
    const std::vector<std::vector<Index>>& cellEdges,
    const DeviceMesh&                      mesh)
{
    if (cellId >= cellEdges.size())
        return Point2::Zero();

    Real a00 = 0.0;
    Real a01 = 0.0;
    Real a11 = 0.0;
    Real b0 = 0.0;
    Real b1 = 0.0;
    int used = 0;
    for (Index edgeId : cellEdges[cellId]) {
        if (edgeId >= signedEdgeFlux.size())
            continue;
        const Edge& edge = mesh.getEdge(edgeId);
        if (edge.length <= 1.0e-30)
            continue;
        const Node& n0 = mesh.getNode(edge.n0);
        const Node& n1 = mesh.getNode(edge.n1);
        const Real tx = (n1.x - n0.x) / edge.length;
        const Real ty = (n1.y - n0.y) / edge.length;
        const Real flux = signedEdgeFlux[edgeId];
        a00 += tx * tx;
        a01 += tx * ty;
        a11 += ty * ty;
        b0 += tx * flux;
        b1 += ty * flux;
        ++used;
    }

    const Real det = a00 * a11 - a01 * a01;
    const Real scale = std::max({std::abs(a00 * a11), std::abs(a01 * a01), Real{1.0}});
    if (used < 2 || std::abs(det) <= 1.0e-24 * scale)
        return Point2::Zero();

    const Real jx = (b0 * a11 - b1 * a01) / det;
    const Real jy = (a00 * b1 - a01 * b0) / det;
    return Point2{jx, jy};
}

inline Real cellVectorCurrentMagnitude(
    Index                                  cellId,
    const std::vector<Real>&               signedEdgeFlux,
    const std::vector<std::vector<Index>>& cellEdges,
    const DeviceMesh&                      mesh)
{
    return cellVectorCurrent(cellId, signedEdgeFlux, cellEdges, mesh).norm();
}

/// Average the constant current vectors reconstructed in cells adjacent to an
/// edge. SG fluxes use the particle-flux convention; callers negate the
/// electron vector when a conventional electron-current direction is needed.
inline Point2 edgeAveragedCellVectorCurrent(
    Index                                  edgeId,
    const std::vector<Real>&               signedEdgeFlux,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<std::vector<Index>>& cellEdges,
    const DeviceMesh&                      mesh)
{
    if (edgeId >= edgeCells.size())
        return Point2::Zero();
    Point2 sum = Point2::Zero();
    int count = 0;
    for (Index cellId : edgeCells[edgeId]) {
        const Point2 current = cellVectorCurrent(
            cellId, signedEdgeFlux, cellEdges, mesh);
        if (current.squaredNorm() <= 0.0)
            continue;
        sum += current;
        ++count;
    }
    if (count > 0)
        return sum / static_cast<Real>(count);
    if (edgeId >= signedEdgeFlux.size() || edgeId >= mesh.numEdges())
        return Point2::Zero();
    const Edge& edge = mesh.getEdge(edgeId);
    if (edge.length <= 1.0e-30)
        return Point2::Zero();
    const Node& n0 = mesh.getNode(edge.n0);
    const Node& n1 = mesh.getNode(edge.n1);
    return signedEdgeFlux[edgeId] *
        Point2{(n1.x - n0.x) / edge.length, (n1.y - n0.y) / edge.length};
}

inline Real cellVectorReconstructedEdgeFluxMagnitude(
    Index                                  edgeId,
    const std::vector<Real>&               signedEdgeFlux,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<std::vector<Index>>& cellEdges,
    const DeviceMesh&                      mesh)
{
    if (edgeId >= edgeCells.size())
        return 0.0;
    Real edgeSum = 0.0;
    int adjacentCellCount = 0;
    for (Index cellId : edgeCells[edgeId]) {
        const Real cellMagnitude = cellVectorCurrentMagnitude(
            cellId, signedEdgeFlux, cellEdges, mesh);
        if (cellMagnitude <= 0.0)
            continue;
        edgeSum += cellMagnitude;
        ++adjacentCellCount;
    }
    if (adjacentCellCount <= 0)
        return edgeId < signedEdgeFlux.size() ? std::abs(signedEdgeFlux[edgeId]) : 0.0;
    return edgeSum / static_cast<Real>(adjacentCellCount);
}

inline std::vector<std::vector<Index>> buildNodeEdgeMap(
    const DeviceMesh& mesh)
{
    std::vector<std::vector<Index>> nodeEdges(
        static_cast<std::size_t>(mesh.numNodes()));
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        nodeEdges[edge.n0].push_back(edgeId);
        nodeEdges[edge.n1].push_back(edgeId);
    }
    return nodeEdges;
}

template <typename FluxAccessor, typename ActiveAccessor, typename WeightAccessor>
inline Point2 nodalWeightedLeastSquaresCurrentVector(
    Index                                  node,
    const std::vector<std::vector<Index>>& nodeEdges,
    const DeviceMesh&                      mesh,
    FluxAccessor&&                         signedFlux,
    ActiveAccessor&&                       activeEdge,
    WeightAccessor&&                       edgeWeight)
{
    if (node >= nodeEdges.size())
        return Point2::Zero();
    Real a00 = 0.0;
    Real a01 = 0.0;
    Real a11 = 0.0;
    Real b0 = 0.0;
    Real b1 = 0.0;
    Real fallbackWeight = 0.0;
    Point2 fallback = Point2::Zero();
    int used = 0;
    for (Index edgeId : nodeEdges[node]) {
        if (!activeEdge(edgeId))
            continue;
        const Edge& edge = mesh.getEdge(edgeId);
        if (edge.length <= 1.0e-30)
            continue;
        const Node& node0 = mesh.getNode(edge.n0);
        const Node& node1 = mesh.getNode(edge.n1);
        const Point2 tangent{
            (node1.x - node0.x) / edge.length,
            (node1.y - node0.y) / edge.length};
        const Real weight = edgeWeight(edgeId);
        if (!(weight > 0.0) || !std::isfinite(weight))
            continue;
        const Real flux = signedFlux(edgeId);
        a00 += weight * tangent.x() * tangent.x();
        a01 += weight * tangent.x() * tangent.y();
        a11 += weight * tangent.y() * tangent.y();
        b0 += weight * tangent.x() * flux;
        b1 += weight * tangent.y() * flux;
        fallback += weight * flux * tangent;
        fallbackWeight += weight;
        ++used;
    }
    const Real determinant = a00 * a11 - a01 * a01;
    const Real scale = std::max({
        std::abs(a00 * a11), std::abs(a01 * a01), Real{1.0e-300}});
    if (used >= 2 && std::abs(determinant) > 1.0e-24 * scale) {
        return Point2{
            (b0 * a11 - b1 * a01) / determinant,
            (a00 * b1 - a01 * b0) / determinant};
    }
    if (fallbackWeight > 0.0)
        return fallback / fallbackWeight;
    return Point2::Zero();
}

/// Recover a nodal current vector from all active incident SG edge-flux
/// projections.  The edge-length weight is the discrete analogue of the
/// element support used when Sentaurus interpolates current density to mesh
/// vertices.  A one-dimensional boundary stencil falls back to its resolved
/// tangential projection instead of returning zero.
template <typename FluxAccessor, typename ActiveAccessor>
inline Point2 nodalLeastSquaresCurrentVector(
    Index                                  node,
    const std::vector<std::vector<Index>>& nodeEdges,
    const DeviceMesh&                      mesh,
    FluxAccessor&&                         signedFlux,
    ActiveAccessor&&                       activeEdge)
{
    return nodalWeightedLeastSquaresCurrentVector(
        node, nodeEdges, mesh,
        std::forward<FluxAccessor>(signedFlux),
        std::forward<ActiveAccessor>(activeEdge),
        [&](Index edgeId) { return mesh.getEdge(edgeId).length; });
}

/// Recover a diagnostic nodal current from the conservative SG fluxes crossing
/// the box-method dual faces.  Weighting the edge-normal projection by the
/// corresponding dual-face length is the boundary-L2 fit of the finite-volume
/// current trace.  The authoritative conservative quantities remain the SG
/// line fluxes; this vector is their vertex representation for field output.
template <typename FluxAccessor, typename ActiveAccessor>
inline Point2 nodalDualFaceSgCurrentVector(
    Index                                  node,
    const std::vector<std::vector<Index>>& nodeEdges,
    const DeviceMesh&                      mesh,
    FluxAccessor&&                         signedFlux,
    ActiveAccessor&&                       activeEdge)
{
    return nodalWeightedLeastSquaresCurrentVector(
        node, nodeEdges, mesh,
        std::forward<FluxAccessor>(signedFlux),
        std::forward<ActiveAccessor>(activeEdge),
        [&](Index edgeId) { return mesh.getEdge(edgeId).couple; });
}

/// Recover a diagnostic nodal current by first fitting one constant vector in
/// every triangle from its three SG edge projections, then area-averaging the
/// incident cell vectors at each vertex.  This deliberately uses the edge
/// opposite a target vertex and therefore has a wider support than the direct
/// nodal fit above.  The SG edge fluxes remain the conservative oracle; this
/// helper changes only their point-field representation.
template <typename FluxAccessor, typename ActiveAccessor>
inline std::vector<Point2> cellFirstAreaWeightedSgCurrentVectors(
    const DeviceMesh& mesh,
    FluxAccessor&&    signedFlux,
    ActiveAccessor&&  activeEdge)
{
    const auto edgeCells = buildEdgeCellMap(mesh);
    const auto cellEdges = buildCellEdgeMap(edgeCells, mesh);
    std::vector<Point2> sums(
        static_cast<std::size_t>(mesh.numNodes()), Point2::Zero());
    std::vector<Real> weights(
        static_cast<std::size_t>(mesh.numNodes()), 0.0);

    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        const Cell& cell = mesh.getCell(cellId);
        if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
            continue;

        Real a00 = 0.0;
        Real a01 = 0.0;
        Real a11 = 0.0;
        Real b0 = 0.0;
        Real b1 = 0.0;
        int used = 0;
        for (Index edgeId : cellEdges[static_cast<std::size_t>(cellId)]) {
            if (!activeEdge(edgeId))
                continue;
            const Edge& edge = mesh.getEdge(edgeId);
            if (!(edge.length > 1.0e-30) ||
                !(edge.couple > 0.0) || !std::isfinite(edge.couple)) {
                continue;
            }
            const Node& node0 = mesh.getNode(edge.n0);
            const Node& node1 = mesh.getNode(edge.n1);
            const Real tx = (node1.x - node0.x) / edge.length;
            const Real ty = (node1.y - node0.y) / edge.length;
            const Real projection = signedFlux(edgeId);
            if (!std::isfinite(projection))
                continue;
            const Real weight = edge.couple;
            a00 += weight * tx * tx;
            a01 += weight * tx * ty;
            a11 += weight * ty * ty;
            b0 += weight * tx * projection;
            b1 += weight * ty * projection;
            ++used;
        }

        const Real determinant = a00 * a11 - a01 * a01;
        const Real matrixScale = std::max({
            std::abs(a00 * a11), std::abs(a01 * a01), Real{1.0e-300}});
        if (used < 2 ||
            std::abs(determinant) <= 1.0e-24 * matrixScale) {
            continue;
        }
        const Point2 cellCurrent{
            (b0 * a11 - b1 * a01) / determinant,
            (a00 * b1 - a01 * b0) / determinant};

        const Node& first = mesh.getNode(cell.node_ids[0]);
        const Node& second = mesh.getNode(cell.node_ids[1]);
        const Node& third = mesh.getNode(cell.node_ids[2]);
        const Real area = 0.5 * std::abs(
            (second.x - first.x) * (third.y - first.y) -
            (second.y - first.y) * (third.x - first.x));
        if (!(area > 0.0) || !std::isfinite(area))
            continue;
        for (Index nodeId : cell.node_ids) {
            sums[static_cast<std::size_t>(nodeId)] += area * cellCurrent;
            weights[static_cast<std::size_t>(nodeId)] += area;
        }
    }

    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const Real weight = weights[static_cast<std::size_t>(node)];
        if (weight > 0.0)
            sums[static_cast<std::size_t>(node)] /= weight;
    }
    return sums;
}

struct EdgeAveragedNodalCurrent {
    Point2 vector = Point2::Zero();
    Real magnitude = 0.0;
};

template <typename FluxAccessor, typename ActiveAccessor>
inline EdgeAveragedNodalCurrent edgeAveragedNodalCurrent(
    Index                                  edgeId,
    const std::vector<std::vector<Index>>& nodeEdges,
    const DeviceMesh&                      mesh,
    FluxAccessor&&                         signedFlux,
    ActiveAccessor&&                       activeEdge)
{
    if (edgeId >= mesh.numEdges())
        return {};
    const Edge& edge = mesh.getEdge(edgeId);
    const Point2 current0 = nodalLeastSquaresCurrentVector(
        edge.n0, nodeEdges, mesh, signedFlux, activeEdge);
    const Point2 current1 = nodalLeastSquaresCurrentVector(
        edge.n1, nodeEdges, mesh, signedFlux, activeEdge);
    return {
        0.5 * (current0 + current1),
        0.5 * (current0.norm() + current1.norm())};
}

template <typename FluxAccessor, typename ActiveAccessor>
inline Point2 edgeAveragedNodalCurrentVector(
    Index                                  edgeId,
    const std::vector<std::vector<Index>>& nodeEdges,
    const DeviceMesh&                      mesh,
    FluxAccessor&&                         signedFlux,
    ActiveAccessor&&                       activeEdge)
{
    return edgeAveragedNodalCurrent(
        edgeId, nodeEdges, mesh,
        std::forward<FluxAccessor>(signedFlux),
        std::forward<ActiveAccessor>(activeEdge)).vector;
}

template <typename FluxAccessor, typename ActiveAccessor>
inline Real edgeAveragedNodalCurrentMagnitude(
    Index                                  edgeId,
    const std::vector<std::vector<Index>>& nodeEdges,
    const DeviceMesh&                      mesh,
    FluxAccessor&&                         signedFlux,
    ActiveAccessor&&                       activeEdge)
{
    return edgeAveragedNodalCurrent(
        edgeId, nodeEdges, mesh,
        std::forward<FluxAccessor>(signedFlux),
        std::forward<ActiveAccessor>(activeEdge)).magnitude;
}


inline Index edgeIdForNodePair(
    const DeviceMesh&          mesh,
    const std::vector<Index>&  candidateEdges,
    Index                      a,
    Index                      b)
{
    for (Index edgeId : candidateEdges) {
        if (edgeId >= mesh.numEdges())
            continue;
        const Edge& edge = mesh.getEdge(edgeId);
        if ((edge.n0 == a && edge.n1 == b) || (edge.n0 == b && edge.n1 == a))
            return edgeId;
    }
    return mesh.numEdges();
}

inline Point2 cellCentroid(const DeviceMesh& mesh, const Cell& cell)
{
    Point2 centroid = Point2::Zero();
    if (cell.node_ids.empty())
        return centroid;
    for (Index nodeId : cell.node_ids) {
        const Node& node = mesh.getNode(nodeId);
        centroid += Point2{node.x, node.y};
    }
    return centroid / static_cast<Real>(cell.node_ids.size());
}

inline Point2 medianDualFaceNormal(
    const DeviceMesh& mesh,
    const Cell&       cell,
    Index             ownerNode,
    Index             neighborNode)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() < 3)
        return Point2::Zero();
    const Node& owner = mesh.getNode(ownerNode);
    const Node& neighbor = mesh.getNode(neighborNode);
    const Point2 centroid = cellCentroid(mesh, cell);
    const Point2 midpoint{0.5 * (owner.x + neighbor.x), 0.5 * (owner.y + neighbor.y)};
    const Point2 segment = midpoint - centroid;
    const Real length = segment.norm();
    if (length <= 1.0e-30)
        return Point2::Zero();
    Point2 normal{segment.y() / length, -segment.x() / length};
    const Point2 towardNeighbor{neighbor.x - owner.x, neighbor.y - owner.y};
    if (normal.dot(towardNeighbor) < 0.0)
        normal = -normal;
    return normal;
}

inline Real medianDualFaceLength(const DeviceMesh& mesh, const Cell& cell, Index a, Index b)
{
    const Node& na = mesh.getNode(a);
    const Node& nb = mesh.getNode(b);
    const Point2 centroid = cellCentroid(mesh, cell);
    const Point2 midpoint{0.5 * (na.x + nb.x), 0.5 * (na.y + nb.y)};
    return (midpoint - centroid).norm();
}

inline Real medianDualCellVectorCurrentMagnitude(
    Index                                  cellId,
    const std::vector<Real>&               signedEdgeFlux,
    const std::vector<std::vector<Index>>& cellEdges,
    const DeviceMesh&                      mesh)
{
    if (cellId >= mesh.numCells() || cellId >= cellEdges.size())
        return 0.0;
    const Cell& cell = mesh.getCell(cellId);
    if (cell.type != CellType::Tri3 || cell.node_ids.size() < 3)
        return 0.0;

    Real a00 = 0.0;
    Real a01 = 0.0;
    Real a11 = 0.0;
    Real b0 = 0.0;
    Real b1 = 0.0;
    Real absSum = 0.0;
    int used = 0;
    for (int k = 0; k < 3; ++k) {
        const Index owner = cell.node_ids[static_cast<std::size_t>(k)];
        const Index neighbor = cell.node_ids[static_cast<std::size_t>((k + 1) % 3)];
        const Index edgeId = edgeIdForNodePair(mesh, cellEdges[cellId], owner, neighbor);
        if (edgeId >= signedEdgeFlux.size())
            continue;
        const Point2 normal = medianDualFaceNormal(mesh, cell, owner, neighbor);
        const Real normalNorm = normal.norm();
        if (normalNorm <= 1.0e-30)
            continue;
        const Edge& edge = mesh.getEdge(edgeId);
        const Real orientation = (edge.n0 == owner && edge.n1 == neighbor) ? 1.0 : -1.0;
        const Real flux = orientation * signedEdgeFlux[edgeId];
        const Real weight = std::max(medianDualFaceLength(mesh, cell, owner, neighbor), Real{1.0e-300});
        const Real nx = normal.x() / normalNorm;
        const Real ny = normal.y() / normalNorm;
        a00 += weight * nx * nx;
        a01 += weight * nx * ny;
        a11 += weight * ny * ny;
        b0 += weight * nx * flux;
        b1 += weight * ny * flux;
        absSum += std::abs(flux);
        ++used;
    }

    const Real det = a00 * a11 - a01 * a01;
    const Real scale = std::max({std::abs(a00 * a11), std::abs(a01 * a01), Real{1.0}});
    if (used < 2 || std::abs(det) <= 1.0e-24 * scale)
        return used > 0 ? absSum / static_cast<Real>(used) : 0.0;
    const Real jx = (b0 * a11 - b1 * a01) / det;
    const Real jy = (a00 * b1 - a01 * b0) / det;
    return std::sqrt(jx * jx + jy * jy);
}

inline Real medianDualFaceVectorReconstructedEdgeFluxMagnitude(
    Index                                  edgeId,
    const std::vector<Real>&               signedEdgeFlux,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<std::vector<Index>>& cellEdges,
    const DeviceMesh&                      mesh)
{
    if (edgeId >= edgeCells.size())
        return 0.0;
    Real edgeSum = 0.0;
    int adjacentCellCount = 0;
    for (Index cellId : edgeCells[edgeId]) {
        const Real cellMagnitude = medianDualCellVectorCurrentMagnitude(
            cellId, signedEdgeFlux, cellEdges, mesh);
        if (cellMagnitude <= 0.0)
            continue;
        edgeSum += cellMagnitude;
        ++adjacentCellCount;
    }
    if (adjacentCellCount <= 0)
        return edgeId < signedEdgeFlux.size() ? std::abs(signedEdgeFlux[edgeId]) : 0.0;
    return edgeSum / static_cast<Real>(adjacentCellCount);
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

inline Real estimateSurfaceDistance(const DeviceMesh& mesh,
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
    const Real mx = 0.5 * (n0.x + n1.x);
    const Real my = 0.5 * (n0.y + n1.y);
    return std::abs((cx - mx) * normalX + (cy - my) * normalY);
}

inline std::pair<Real, Real> nearestSurfaceFieldAndDistanceForCell(
    const DeviceMesh& mesh,
    const std::vector<std::vector<Index>>& edgeCells,
    const VectorXd& psi,
    Index cellId,
    const MobilityModelConfig& mobilityConfig,
    Real fieldFactor)
{
    const Cell& cell = mesh.getCell(cellId);
    const Region& region = mesh.getRegion(cell.region_id);
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
        return {std::numeric_limits<Real>::quiet_NaN(),
                std::numeric_limits<Real>::quiet_NaN()};
    const Node& p0 = mesh.getNode(cell.node_ids[0]);
    const Node& p1 = mesh.getNode(cell.node_ids[1]);
    const Node& p2 = mesh.getNode(cell.node_ids[2]);
    const Real dx10 = p1.x - p0.x;
    const Real dy10 = p1.y - p0.y;
    const Real dx20 = p2.x - p0.x;
    const Real dy20 = p2.y - p0.y;
    const Real det = dx10 * dy20 - dy10 * dx20;
    if (std::abs(det) <= 1.0e-300)
        return {std::numeric_limits<Real>::quiet_NaN(),
                std::numeric_limits<Real>::quiet_NaN()};
    const Real dv10 = psi(static_cast<int>(cell.node_ids[1])) -
        psi(static_cast<int>(cell.node_ids[0]));
    const Real dv20 = psi(static_cast<int>(cell.node_ids[2])) -
        psi(static_cast<int>(cell.node_ids[0]));
    const Point2 gradient{
        (dv10 * dy20 - dv20 * dy10) / det,
        (dx10 * dv20 - dx20 * dv10) / det};

    Real nearestDistance = std::numeric_limits<Real>::infinity();
    Real nearestField = std::numeric_limits<Real>::quiet_NaN();
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        std::vector<std::string> adjacentRegions;
        for (Index candidateCell : edgeCells.at(edgeId))
            adjacentRegions.push_back(
                mesh.getRegion(mesh.getCell(candidateCell).region_id).name);
        if (!surfaceMobilityAppliesToRegionPair(
                mobilityConfig, region.name, adjacentRegions))
            continue;
        const Node& a = mesh.getNode(edge.n0);
        const Node& b = mesh.getNode(edge.n1);
        if (edge.length <= 1.0e-30)
            continue;
        const Point2 normal{-(b.y - a.y) / edge.length,
                             (b.x - a.x) / edge.length};
        const auto [cx, cy] = cellCentroid(mesh, cellId);
        const Real abx = b.x - a.x;
        const Real aby = b.y - a.y;
        const Real projection = std::clamp(
            ((cx - a.x) * abx + (cy - a.y) * aby) /
                (edge.length * edge.length),
            0.0, 1.0);
        const Real nearestX = a.x + projection * abx;
        const Real nearestY = a.y + projection * aby;
        const Real distance = std::hypot(cx - nearestX, cy - nearestY);
        if (distance < nearestDistance) {
            nearestDistance = distance;
            nearestField = std::abs(gradient.dot(normal)) * fieldFactor;
        }
    }
    return {nearestField, nearestDistance};
}

inline void updateSurfaceMobilityCellGeometry(
    MobilityModelConfig& config,
    const DeviceMesh& mesh,
    const std::vector<std::vector<Index>>& edgeCells,
    const VectorXd& psi,
    Real fieldFactor)
{
    if (!isSurfaceMobilityModel(config))
        return;
    if (config.surface.cellNormalX.size() != mesh.numCells() ||
        config.surface.cellNormalY.size() != mesh.numCells() ||
        config.surface.cellDistances.size() != mesh.numCells()) {
        config.surface.cellNormalX.assign(
            mesh.numCells(), std::numeric_limits<Real>::quiet_NaN());
        config.surface.cellNormalY.assign(
            mesh.numCells(), std::numeric_limits<Real>::quiet_NaN());
        config.surface.cellDistances.assign(
            mesh.numCells(), std::numeric_limits<Real>::quiet_NaN());
        struct InterfaceSegment {
            Real ax;
            Real ay;
            Real bx;
            Real by;
            Real nx;
            Real ny;
            std::vector<std::string> regions;
        };
        std::vector<InterfaceSegment> interfaces;
        for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
            const Edge& edge = mesh.getEdge(edgeId);
            std::vector<std::string> regions;
            for (Index adjacentCell : edgeCells.at(edgeId))
                regions.push_back(mesh.getRegion(
                    mesh.getCell(adjacentCell).region_id).name);
            if (regions.size() < 2 || edge.length <= 1.0e-30)
                continue;
            // A configured interface is an exact physical selector, not only
            // a per-cell applicability test.  Filter the geometry candidates
            // here so internal same-region edges do not participate in the
            // nearest-interface search.  This also changes the preprocessing
            // cost from O(cells * all_edges) to O(cells * selected_interface).
            if (!config.surface.surfaceInterface.empty()) {
                if (config.surface.surfaceInterface.size() != 2)
                    throw std::invalid_argument(
                        "surface mobility surface_interface must contain exactly two region names.");
                const std::string& interfaceA =
                    config.surface.surfaceInterface.at(0);
                const std::string& interfaceB =
                    config.surface.surfaceInterface.at(1);
                const bool hasA = std::find(
                    regions.begin(), regions.end(), interfaceA) != regions.end();
                const bool hasB = std::find(
                    regions.begin(), regions.end(), interfaceB) != regions.end();
                if (!hasA || !hasB)
                    continue;
            }
            const Node& a = mesh.getNode(edge.n0);
            const Node& b = mesh.getNode(edge.n1);
            interfaces.push_back({
                a.x, a.y, b.x, b.y,
                -(b.y - a.y) / edge.length,
                 (b.x - a.x) / edge.length,
                std::move(regions)});
        }
        for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
            const Region& region = mesh.getRegion(mesh.getCell(cellId).region_id);
            const auto [cx, cy] = cellCentroid(mesh, cellId);
            Real nearest = std::numeric_limits<Real>::infinity();
            for (const InterfaceSegment& segment : interfaces) {
                if (!surfaceMobilityAppliesToRegionPair(
                        config, region.name, segment.regions))
                    continue;
                const Real abx = segment.bx - segment.ax;
                const Real aby = segment.by - segment.ay;
                const Real length2 = abx * abx + aby * aby;
                const Real projection = std::clamp(
                    ((cx - segment.ax) * abx + (cy - segment.ay) * aby) /
                        length2,
                    0.0, 1.0);
                const Real nearestX = segment.ax + projection * abx;
                const Real nearestY = segment.ay + projection * aby;
                const Real distance = std::hypot(cx - nearestX, cy - nearestY);
                if (distance < nearest) {
                    nearest = distance;
                    config.surface.cellNormalX[cellId] = segment.nx;
                    config.surface.cellNormalY[cellId] = segment.ny;
                    config.surface.cellDistances[cellId] = distance;
                }
            }
        }
    }
    config.surface.cellNormalFields.assign(
        mesh.numCells(), std::numeric_limits<Real>::quiet_NaN());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        if (!std::isfinite(config.surface.cellNormalX[cellId]) ||
            !std::isfinite(config.surface.cellNormalY[cellId]))
            continue;
        const Cell& cell = mesh.getCell(cellId);
        if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
            continue;
        const Node& p0 = mesh.getNode(cell.node_ids[0]);
        const Node& p1 = mesh.getNode(cell.node_ids[1]);
        const Node& p2 = mesh.getNode(cell.node_ids[2]);
        const Real dx10 = p1.x - p0.x;
        const Real dy10 = p1.y - p0.y;
        const Real dx20 = p2.x - p0.x;
        const Real dy20 = p2.y - p0.y;
        const Real det = dx10 * dy20 - dy10 * dx20;
        if (std::abs(det) <= 1.0e-300)
            continue;
        const Real dv10 = psi(static_cast<int>(cell.node_ids[1])) -
            psi(static_cast<int>(cell.node_ids[0]));
        const Real dv20 = psi(static_cast<int>(cell.node_ids[2])) -
            psi(static_cast<int>(cell.node_ids[0]));
        const Real gradientX = (dv10 * dy20 - dv20 * dy10) / det;
        const Real gradientY = (dx10 * dv20 - dx20 * dv10) / det;
        config.surface.cellNormalFields[cellId] = std::abs(
            gradientX * config.surface.cellNormalX[cellId] +
            gradientY * config.surface.cellNormalY[cellId]) * fieldFactor;
    }
}

inline Real cellAverageTotalImpurity(const DeviceMesh& mesh,
                                     const DopingModel& doping,
                                     Index cellId)
{
    const Cell& cell = mesh.getCell(cellId);
    if (cell.node_ids.empty())
        return 0.0;
    Real sum = 0.0;
    for (Index nodeId : cell.node_ids)
        sum += doping.totalImpurity(nodeId);
    return sum / static_cast<Real>(cell.node_ids.size());
}

inline Real edgeMobilityDopingConcentration(
    const DeviceMesh& mesh,
    const DopingModel& doping,
    const Edge& edge,
    Index cellId,
    const MobilityModelConfig* config)
{
    const std::string basis = config != nullptr
        ? config->dopingConcentrationBasis
        : "net_doping";
    if (basis == "total_impurity") {
        return 0.5 * (doping.totalImpurity(edge.n0) +
                      doping.totalImpurity(edge.n1));
    }
    if (basis == "cell_reconstructed_total_impurity")
        return cellAverageTotalImpurity(mesh, doping, cellId);
    return 0.5 * (doping.netDoping(edge.n0) + doping.netDoping(edge.n1));
}

inline Real nodeMobilityDopingConcentration(
    const DeviceMesh& mesh,
    const DopingModel& doping,
    Index nodeId,
    Index cellId,
    const MobilityModelConfig* config)
{
    const std::string basis = config != nullptr
        ? config->dopingConcentrationBasis
        : "net_doping";
    if (basis == "total_impurity")
        return doping.totalImpurity(nodeId);
    if (basis == "cell_reconstructed_total_impurity")
        return cellAverageTotalImpurity(mesh, doping, cellId);
    return doping.netDoping(nodeId);
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
        const bool lombardi = mobilityConfig != nullptr &&
            (mobilityConfig->model == "masetti_lombardi" ||
             mobilityConfig->model == "masetti_field_lombardi");
        const bool surfaceApplies = surfaceEnabled &&
            (lombardi
                ? (mobilityConfig->surface.surfaceRegion.empty() ||
                   mobilityConfig->surface.surfaceRegion == region.name)
                : surfaceMobilityAppliesToRegionPair(
                    *mobilityConfig, region.name, adjacentRegionNames));
        const Real surfaceNormalField = (surfaceApplies && psi != nullptr)
            ? (c < mobilityConfig->surface.cellNormalFields.size()
                ? mobilityConfig->surface.cellNormalFields[c]
                : estimateSurfaceNormalField(cells, mesh, *psi, edgeId, c) *
                    mobilityConfig->surface.coordinateFieldFactor)
            : std::numeric_limits<Real>::quiet_NaN();
        const Real surfaceDistance = surfaceApplies
            ? (c < mobilityConfig->surface.cellDistances.size()
                ? mobilityConfig->surface.cellDistances[c]
                : estimateSurfaceDistance(mesh, edgeId, c))
            : std::numeric_limits<Real>::quiet_NaN();
        const Real mobilityDoping = edgeMobilityDopingConcentration(
            mesh, doping, edge, c, mobilityConfig);
        const Real modelMobility = (carrier == CarrierType::Electron)
            ? mobility.electronMobility(
                material, mobilityDoping, 0.0, 0.0, electricField,
                surfaceNormalField, surfaceDistance)
            : mobility.holeMobility(
                material, mobilityDoping, 0.0, 0.0, electricField,
                surfaceNormalField, surfaceDistance);
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
                         Real                                  drivingField,
                         const MobilityModelConfig*            mobilityConfig = nullptr)
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
        const bool hasSurfaceGeometry =
            mobilityConfig != nullptr &&
            c < mobilityConfig->surface.cellNormalFields.size() &&
            c < mobilityConfig->surface.cellDistances.size() &&
            std::isfinite(mobilityConfig->surface.cellNormalFields[c]) &&
            std::isfinite(mobilityConfig->surface.cellDistances[c]);
        const Real surfaceNormalField = hasSurfaceGeometry
            ? mobilityConfig->surface.cellNormalFields[c]
            : std::numeric_limits<Real>::quiet_NaN();
        const Real surfaceDistance = hasSurfaceGeometry
            ? mobilityConfig->surface.cellDistances[c]
            : std::numeric_limits<Real>::quiet_NaN();
        const Real mobilityDoping = nodeMobilityDopingConcentration(
            mesh, doping, nodeId, c, mobilityConfig);
        const Real modelMobility = (carrier == CarrierType::Electron)
            ? mobility.electronMobility(
                material, mobilityDoping, 0.0, 0.0, drivingField,
                surfaceNormalField, surfaceDistance)
            : mobility.holeMobility(
                material, mobilityDoping, 0.0, 0.0, drivingField,
                surfaceNormalField, surfaceDistance);
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
    if (config.debugRawVanOverstraeten ||
        config.drivingForceInterpolation != "quasi_fermi_to_electric_field" ||
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
           config.drivingForce == "effective_field_parallel_j" ||
           config.drivingForce == "eparallel";
}

inline bool usesSentaurusEparallelAvalancheDrivingForce(
    const ImpactIonizationModelConfig& config)
{
    return config.drivingForce == "eparallel";
}

inline bool usesQuasiFermiAvalancheDrivingForce(
    const ImpactIonizationModelConfig& config)
{
    return config.debugRawVanOverstraeten ||
           config.drivingForce == "quasi_fermi_gradient";
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

/// Sentaurus Eparallel drive: the non-negative electric-field component in
/// the direction of the conventional carrier current, E dot J / |J|.
inline Real sentaurusEparallelAvalancheDrivingField(
    const Point2& electricField,
    const Point2& conventionalCurrent)
{
    if (!electricField.allFinite() || !conventionalCurrent.allFinite())
        return 0.0;
    const Real currentMagnitude = conventionalCurrent.norm();
    if (currentMagnitude <= 1.0e-300)
        return 0.0;
    return std::max(
        electricField.dot(conventionalCurrent) / currentMagnitude, 0.0);
}

/// Resolves the legacy SG edge-current avalanche source-volume factor used in
/// `factor * h * edge.couple`. A finite `source_volume_factor` overrides the
/// named `source_volume_policy` preset; `0` falls back to the preset.
inline Real avalancheSourceVolumeFactor(const ImpactIonizationModelConfig& config)
{
    if (config.sourceVolumeFactor > 0.0)
        return config.sourceVolumeFactor;
    return config.sourceVolumePolicy == "edge_box" ? 1.0 : 0.5;
}

inline Real triangleSignedDoubleArea(const Point2& a, const Point2& b, const Point2& c)
{
    return (b.x() - a.x()) * (c.y() - a.y()) -
           (c.x() - a.x()) * (b.y() - a.y());
}

inline Point2 meshPoint(const DeviceMesh& mesh, Index node)
{
    const Node& n = mesh.getNode(node);
    return Point2{n.x, n.y};
}

inline int tri3LocalEdgeIndex(const Cell& cell, Index edgeNode0, Index edgeNode1)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
        return -1;
    for (int local = 0; local < 3; ++local) {
        const Index a = cell.node_ids[static_cast<std::size_t>(local)];
        const Index b = cell.node_ids[static_cast<std::size_t>((local + 1) % 3)];
        if ((a == edgeNode0 && b == edgeNode1) ||
            (a == edgeNode1 && b == edgeNode0)) {
            return local;
        }
    }
    return -1;
}

inline Real angleBetween(const Point2& a, const Point2& b)
{
    const Real denom = a.norm() * b.norm();
    if (denom <= 1.0e-300)
        return 0.0;
    const Real cosTheta = std::clamp(a.dot(b) / denom, -1.0, 1.0);
    return std::acos(cosTheta);
}

inline Real geniusTri3TruncatedPartialVolumeWithEdge(
    const DeviceMesh& mesh,
    const Cell&       cell,
    Index             edgeNode0,
    Index             edgeNode1)
{
    const int localEdge = tri3LocalEdgeIndex(cell, edgeNode0, edgeNode1);
    if (localEdge < 0)
        return 0.0;

    const std::array<Index, 3> ids = {
        cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};
    const std::array<Point2, 3> p = {
        meshPoint(mesh, ids[0]), meshPoint(mesh, ids[1]), meshPoint(mesh, ids[2])};
    const Real lengthScale = std::max({
        (p[0] - p[1]).norm(), (p[1] - p[2]).norm(), (p[2] - p[0]).norm()});
    const Real distanceTolerance =
        64.0 * std::numeric_limits<Real>::epsilon() * lengthScale;
    const Real det = triangleSignedDoubleArea(p[0], p[1], p[2]);
    if (std::abs(det) <= 1.0e-300)
        return 0.0;

    const Real a2 = p[0].squaredNorm();
    const Real b2 = p[1].squaredNorm();
    const Real c2 = p[2].squaredNorm();
    const Real invDenom = 1.0 / (2.0 * det);
    const Point2 circumcenter{
        (a2 * (p[1].y() - p[2].y()) +
         b2 * (p[2].y() - p[0].y()) +
         c2 * (p[0].y() - p[1].y())) * invDenom,
        (a2 * (p[2].x() - p[1].x()) +
         b2 * (p[0].x() - p[2].x()) +
         c2 * (p[1].x() - p[0].x())) * invDenom};

    constexpr int sideNodes[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    std::array<Real, 3> lengths = {0.0, 0.0, 0.0};
    std::array<Real, 3> dt = {0.0, 0.0, 0.0};
    int obtuseEdge = -1;
    for (int local = 0; local < 3; ++local) {
        const Point2& p1 = p[sideNodes[local][0]];
        const Point2& p2 = p[sideNodes[local][1]];
        const Point2& p3 = p[(2 + local) % 3];
        const Point2 sideCenter = 0.5 * (p1 + p2);
        lengths[static_cast<std::size_t>(local)] = (p1 - p2).norm();
        const Real rawDistance = (sideCenter - circumcenter).norm();
        const Real distance =
            rawDistance <= distanceTolerance ? 0.0 : rawDistance;
        if ((p1 - p3).dot(p2 - p3) < 0.0) {
            dt[static_cast<std::size_t>(local)] = -distance;
            obtuseEdge = local;
        } else {
            dt[static_cast<std::size_t>(local)] = distance;
        }
    }

    if (obtuseEdge >= 0) {
        const int obtuseNode = (2 + obtuseEdge) % 3;
        const Point2& p1 = p[sideNodes[obtuseEdge][0]];
        const Point2& p2 = p[sideNodes[obtuseEdge][1]];
        const Point2& p3 = p[obtuseNode];
        const Real theta1 = angleBetween(p2 - p1, p3 - p1);
        const Real theta2 = angleBetween(p1 - p2, p3 - p2);
        const Real cos1 = std::cos(theta1);
        const Real cos2 = std::cos(theta2);
        const Point2 preEdgeCenter = 0.5 * (p1 + p3);
        const Point2 posEdgeCenter = 0.5 * (p2 + p3);
        dt[static_cast<std::size_t>(obtuseEdge)] = 0.0;
        if (std::abs(cos1) > 1.0e-300 && (p2 - p1).norm() > 1.0e-300) {
            const Point2 m1 = p1 + (p2 - p1).normalized() *
                ((preEdgeCenter - p1).norm() / cos1);
            dt[static_cast<std::size_t>((obtuseEdge + 2) % 3)] =
                (preEdgeCenter - m1).norm();
        }
        if (std::abs(cos2) > 1.0e-300 && (p1 - p2).norm() > 1.0e-300) {
            const Point2 m2 = p2 + (p1 - p2).normalized() *
                ((posEdgeCenter - p2).norm() / cos2);
            dt[static_cast<std::size_t>((obtuseEdge + 1) % 3)] =
                (posEdgeCenter - m2).norm();
        }
    }

    return 0.5 * lengths[static_cast<std::size_t>(localEdge)] *
           std::max(0.0, dt[static_cast<std::size_t>(localEdge)]);
}

inline std::array<Real, 3> tri3ElementEdgeBoxPartialVolumes(
    const DeviceMesh& mesh,
    const Cell&       cell)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3) {
        throw std::invalid_argument(
            "element-edge GSS/Laux reconstruction requires a Tri3 cell");
    }
    std::array<Real, 3> partialVolumes{};
    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const Index node0 =
            cell.node_ids[static_cast<std::size_t>(localEdge)];
        const Index node1 =
            cell.node_ids[static_cast<std::size_t>((localEdge + 1) % 3)];
        partialVolumes[static_cast<std::size_t>(localEdge)] =
            geniusTri3TruncatedPartialVolumeWithEdge(
                mesh, cell, node0, node1);
    }
    const Real exactArea = 0.5 * std::abs(triangleSignedDoubleArea(
        meshPoint(mesh, cell.node_ids[0]),
        meshPoint(mesh, cell.node_ids[1]),
        meshPoint(mesh, cell.node_ids[2])));
    const Real truncatedArea =
        partialVolumes[0] + partialVolumes[1] + partialVolumes[2];
    const Real closureTolerance =
        128.0 * std::numeric_limits<Real>::epsilon() * exactArea;
    if (exactArea > 0.0 && truncatedArea > 0.0 &&
        std::abs(truncatedArea - exactArea) > closureTolerance) {
        const Real conservativeScale = exactArea / truncatedArea;
        for (Real& partialVolume : partialVolumes)
            partialVolume *= conservativeScale;
    }
    return partialVolumes;
}

inline std::array<Real, 3> tri3ElementVertexBoxMeasures(
    const DeviceMesh& mesh,
    const Cell&       cell)
{
    const auto edgePartialVolumes =
        tri3ElementEdgeBoxPartialVolumes(mesh, cell);
    std::array<Real, 3> vertexMeasures{};
    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const Real halfMeasure =
            0.5 * edgePartialVolumes[static_cast<std::size_t>(localEdge)];
        vertexMeasures[static_cast<std::size_t>(localEdge)] += halfMeasure;
        vertexMeasures[static_cast<std::size_t>((localEdge + 1) % 3)] +=
            halfMeasure;
    }
    return vertexMeasures;
}

inline Point2 solveTri3EdgeProjectionPair(
    const Point2& tangentA,
    Real          valueA,
    const Point2& tangentB,
    Real          valueB)
{
    const Real determinant =
        tangentA.x() * tangentB.y() - tangentA.y() * tangentB.x();
    if (std::abs(determinant) <= 1.0e-14) {
        throw std::invalid_argument(
            "parallel triangle edges cannot reconstruct a current vector");
    }
    return Point2{
        (valueA * tangentB.y() - tangentA.y() * valueB) / determinant,
        (tangentA.x() * valueB - valueA * tangentB.x()) / determinant};
}

inline Point2 gssLauxTri3CurrentVector(
    const DeviceMesh&          mesh,
    const Cell&                cell,
    const std::array<Real, 3>& signedEdgeCurrent)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3) {
        throw std::invalid_argument(
            "element-edge GSS/Laux reconstruction requires a Tri3 cell");
    }

    std::array<Point2, 3> tangents{};
    std::array<Real, 3> lengths{};
    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const Point2 delta =
            meshPoint(
                mesh,
                cell.node_ids[static_cast<std::size_t>(
                    (localEdge + 1) % 3)]) -
            meshPoint(
                mesh,
                cell.node_ids[static_cast<std::size_t>(localEdge)]);
        lengths[static_cast<std::size_t>(localEdge)] = delta.norm();
        if (lengths[static_cast<std::size_t>(localEdge)] <= 1.0e-30) {
            throw std::invalid_argument(
                "degenerate triangle edge cannot reconstruct a current vector");
        }
        tangents[static_cast<std::size_t>(localEdge)] =
            delta / lengths[static_cast<std::size_t>(localEdge)];
    }

    const std::array<Point2, 3> pairVectors = {
        solveTri3EdgeProjectionPair(
            tangents[0], signedEdgeCurrent[0],
            tangents[1], signedEdgeCurrent[1]),
        solveTri3EdgeProjectionPair(
            tangents[0], signedEdgeCurrent[0],
            tangents[2], signedEdgeCurrent[2]),
        solveTri3EdgeProjectionPair(
            tangents[1], signedEdgeCurrent[1],
            tangents[2], signedEdgeCurrent[2])};
    const auto pairIndex = [](int first, int second) {
        const int low = std::min(first, second);
        const int high = std::max(first, second);
        return low == 0 ? high - 1 : 2;
    };

    const auto partialVolumes =
        tri3ElementEdgeBoxPartialVolumes(mesh, cell);
    std::array<Point2, 3> edgeVectors{};
    for (int target = 0; target < 3; ++target) {
        std::array<int, 2> others{};
        int next = 0;
        for (int candidate = 0; candidate < 3; ++candidate) {
            if (candidate != target)
                others[static_cast<std::size_t>(next++)] = candidate;
        }
        const int first = others[0];
        const int second = others[1];
        const Real firstWeight =
            2.0 * partialVolumes[static_cast<std::size_t>(first)] /
            lengths[static_cast<std::size_t>(first)];
        const Real secondWeight =
            2.0 * partialVolumes[static_cast<std::size_t>(second)] /
            lengths[static_cast<std::size_t>(second)];
        const Point2& firstPair =
            pairVectors[static_cast<std::size_t>(pairIndex(target, first))];
        const Point2& secondPair =
            pairVectors[static_cast<std::size_t>(pairIndex(target, second))];
        const Real weightSum = firstWeight + secondWeight;
        if (weightSum > 0.0) {
            edgeVectors[static_cast<std::size_t>(target)] =
                (firstWeight * firstPair + secondWeight * secondPair) /
                weightSum;
        } else {
            edgeVectors[static_cast<std::size_t>(target)] =
                0.5 * (firstPair + secondPair);
        }
    }

    const Real totalVolume =
        partialVolumes[0] + partialVolumes[1] + partialVolumes[2];
    if (totalVolume <= 1.0e-300) {
        throw std::invalid_argument(
            "triangle has no positive box partial volume");
    }
    return (
        partialVolumes[0] * edgeVectors[0] +
        partialVolumes[1] * edgeVectors[1] +
        partialVolumes[2] * edgeVectors[2]) / totalVolume;
}

inline Real geniusTruncatedEdgeSourceVolume(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId)
{
    if (edgeId >= mesh.numEdges() || edgeId >= edgeCells.size())
        return 0.0;
    const Edge& edge = mesh.getEdge(edgeId);
    Real volume = 0.0;
    for (Index cellId : edgeCells[edgeId]) {
        if (cellId >= mesh.numCells())
            continue;
        volume += geniusTri3TruncatedPartialVolumeWithEdge(
            mesh, mesh.getCell(cellId), edge.n0, edge.n1);
    }
    return volume;
}

/// Conservative counterpart of the legacy Genius truncated edge support.
/// Each transport triangle first normalizes its three truncated edge pieces
/// to the exact triangle area, then contributes the piece belonging to the
/// queried edge. Non-transport cells are excluded so a Si/oxide interface
/// cannot add oxide area to IntegrSemiconductor AvalancheGeneration.
inline Real geniusConservativeEdgeSourceVolume(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    const std::vector<Material>&           cellMaterials,
    Index                                  edgeId)
{
    if (edgeId >= mesh.numEdges() || edgeId >= edgeCells.size())
        return 0.0;
    const Edge& edge = mesh.getEdge(edgeId);
    Real volume = 0.0;
    for (Index cellId : edgeCells[edgeId]) {
        if (cellId >= mesh.numCells() || cellId >= cellMaterials.size())
            continue;
        const Material& material = cellMaterials[cellId];
        if (!(material.ni > 0.0 || material.mun > 0.0 || material.mup > 0.0))
            continue;
        const Cell& cell = mesh.getCell(cellId);
        const int localEdge = tri3LocalEdgeIndex(cell, edge.n0, edge.n1);
        if (localEdge < 0)
            continue;
        const auto partialVolumes = tri3ElementEdgeBoxPartialVolumes(mesh, cell);
        volume += partialVolumes[static_cast<std::size_t>(localEdge)];
    }
    return volume;
}

inline Real avalancheSourceEdgeArea(
    const ImpactIonizationModelConfig&     config,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    const std::vector<Material>*           cellMaterials = nullptr)
{
    if (edgeId >= mesh.numEdges())
        return 0.0;
    const Edge& edge = mesh.getEdge(edgeId);
    Real area = 0.0;
    if (config.sourceVolumeFactor > 0.0) {
        area = config.sourceVolumeFactor * edge.length * edge.couple;
    } else if (config.sourceVolumePolicy == "genius_truncated") {
        area = geniusTruncatedEdgeSourceVolume(edgeCells, mesh, edgeId);
    } else if (config.sourceVolumePolicy == "genius_conservative") {
        if (cellMaterials == nullptr) {
            throw std::invalid_argument(
                "genius_conservative avalanche source volume requires cell materials");
        }
        area = geniusConservativeEdgeSourceVolume(
            edgeCells, mesh, *cellMaterials, edgeId);
    } else {
        area = avalancheSourceVolumeFactor(config) * edge.length * edge.couple;
    }
    return area * config.sourceGeometryScale;
}

/// Validates the impact-ionization configuration shared by the Gummel and
/// Newton solver config loaders. `context` is prefixed to any thrown message.
inline void validateImpactIonizationDrivingForce(const ImpactIonizationModelConfig& config,
                                                 const char* context)
{
    if (config.couplingMode != "self_consistent" &&
        config.couplingMode != "postprocess_only") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.coupling_mode must be 'self_consistent' or "
            "'postprocess_only'.");
    }
    const bool configuredCurrentAlignedDrivingForce =
        usesCurrentAlignedAvalancheDrivingForce(config);
    const bool currentAlignedDrivingForce =
        !config.debugRawVanOverstraeten && configuredCurrentAlignedDrivingForce;
    if (config.contactElectricFieldFallbackScope != "contact_node_cell" &&
        config.contactElectricFieldFallbackScope != "contact_boundary_face") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.contact_electric_field_fallback_scope must be "
            "'contact_node_cell' or 'contact_boundary_face'.");
    }

    const bool validContactFallbackMode =
        config.contactElectricFieldFallbackMode == "cell_gradient_magnitude" ||
        config.contactElectricFieldFallbackMode == "face_normal" ||
        config.contactElectricFieldFallbackMode == "one_sided" ||
        config.contactElectricFieldFallbackMode == "distance_weighted_blend";
    if (!validContactFallbackMode) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.contact_electric_field_fallback_mode must be "
            "'cell_gradient_magnitude', 'face_normal', 'one_sided', or "
            "'distance_weighted_blend'.");
    }
    if (config.contactElectricFieldFallbackMode != "cell_gradient_magnitude" &&
        config.contactElectricFieldFallbackScope != "contact_boundary_face") {
        throw std::invalid_argument(
            std::string(context) +
            ": contact-face fallback modes require scope 'contact_boundary_face'.");
    }

    if (config.contactElectricFieldFallback &&
        (config.drivingForce != "quasi_fermi_gradient" ||
         config.quasiFermiGradientDiscretization != "cell_gradient" ||
         config.generation != "current_density" ||
         (config.sourceMappingMode != "triangle_gss_gradqf_truncated" &&
          config.sourceMappingMode != "element_vertex_box_measure"))) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.contact_electric_field_fallback requires "
            "current-density GradQF with cell_gradient and a triangle or "
            "element-vertex-box source mapping.");
    }

    if (config.drivingForce != "electric_field" &&
        config.drivingForce != "quasi_fermi_gradient" &&
        !configuredCurrentAlignedDrivingForce) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.driving_force must be 'electric_field', "
            "'quasi_fermi_gradient', 'grad_potential_parallel_j', "
            "'effective_field_parallel_j', or 'eparallel'.");
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
        config.currentApproximation != "grad_qf" &&
        config.currentApproximation != "cell_reconstructed" &&
        config.currentApproximation != "psi_gradient_proxy" &&
        config.currentApproximation != "cell_current_reconstructed" &&
        config.currentApproximation != "cell_vector_current_reconstructed" &&
        config.currentApproximation != "nodal_vector_current_reconstructed" &&
        config.currentApproximation != "element_edge_sg_gss_laux" &&
        config.currentApproximation != "conserved_total_current") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.current_approximation must be "
            "'mobility_density_gradient', 'density_gradient', 'grad_qf', "
            "'cell_reconstructed', 'psi_gradient_proxy', 'cell_current_reconstructed', "
            "'cell_vector_current_reconstructed', 'nodal_vector_current_reconstructed', "
            "'element_edge_sg_gss_laux', "
            "or 'conserved_total_current'.");
    }
    if (config.currentMagnitudeMode != "edge_scalar_abs" &&
        config.currentMagnitudeMode != "dual_face_vector_mag") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.current_magnitude_mode must be "
            "'edge_scalar_abs' or 'dual_face_vector_mag'.");
    }
    if (config.eparallelFieldRecovery != "edge_adjacent_cells" &&
        config.eparallelFieldRecovery != "nodal_vertex_star") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.eparallel_field_recovery must be "
            "'edge_adjacent_cells' or 'nodal_vertex_star'.");
    }
    if (config.eparallelFieldRecovery == "nodal_vertex_star" &&
        !usesSentaurusEparallelAvalancheDrivingForce(config)) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.eparallel_field_recovery='nodal_vertex_star' "
            "requires driving_force='eparallel'.");
    }
    if (config.cellReconstructedMidpointDensity != "bernoulli" &&
        config.cellReconstructedMidpointDensity != "arithmetic" &&
        config.cellReconstructedMidpointDensity != "gss_logistic") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.cell_reconstructed_midpoint_density must be "
            "'bernoulli', 'arithmetic', or 'gss_logistic'.");
    }
    if (config.drivingForceInterpolation != "none" &&
        config.drivingForceInterpolation != "quasi_fermi_to_electric_field") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.driving_force_interpolation.mode must be "
            "'none' or 'quasi_fermi_to_electric_field'.");
    }
    const bool eparallelVectorCurrent =
        usesSentaurusEparallelAvalancheDrivingForce(config) &&
        (config.currentApproximation == "cell_vector_current_reconstructed" ||
         config.currentApproximation == "nodal_vector_current_reconstructed");
    if (currentAlignedDrivingForce &&
        (config.generation != "current_density" ||
         (config.currentApproximation != "density_gradient" &&
          config.currentApproximation != "grad_qf" &&
          !eparallelVectorCurrent))) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization current-aligned driving forces require "
            "generation='current_density' with current_approximation='density_gradient', "
            "'grad_qf', or Sentaurus Eparallel with "
            "'cell_vector_current_reconstructed' or "
            "'nodal_vector_current_reconstructed'.");
    }
    if (config.drivingForceInterpolation != "none" &&
        config.drivingForce != "quasi_fermi_gradient" &&
        !config.debugRawVanOverstraeten) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.driving_force_interpolation requires "
            "driving_force='quasi_fermi_gradient'.");
    }
    if (config.quasiFermiGradientDiscretization != "edge_difference" &&
        config.quasiFermiGradientDiscretization != "cell_gradient") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.quasi_fermi_gradient_discretization must be "
            "'edge_difference' or 'cell_gradient'.");
    }
    if (config.quasiFermiGradientDiscretization == "cell_gradient" &&
        config.drivingForce != "quasi_fermi_gradient" &&
        !config.debugRawVanOverstraeten) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.quasi_fermi_gradient_discretization='cell_gradient' "
            "requires driving_force='quasi_fermi_gradient'.");
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
    if (config.sourceVolumePolicy != "genius_truncated" &&
        config.sourceVolumePolicy != "genius_conservative" &&
        config.sourceVolumePolicy != "edge_half_box" &&
        config.sourceVolumePolicy != "edge_box") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.source_volume_policy must be 'genius_truncated', "
            "'genius_conservative', 'edge_half_box', or 'edge_box'.");
    }
    if (config.sourceVolumeFactor != 0.0 &&
        (!std::isfinite(config.sourceVolumeFactor) ||
         config.sourceVolumeFactor < 0.5 ||
         config.sourceVolumeFactor > 1.0)) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.source_volume_factor must be 0 or within [0.5, 1.0].");
    }
    if (config.sourceMappingMode != "node_F_node_alpha_node_G" &&
        config.sourceMappingMode != "edge_F_edge_alpha_edge_G_to_node" &&
        config.sourceMappingMode != "cell_F_cell_alpha_cell_G_to_node" &&
        config.sourceMappingMode != "nodal_eparallel_p1" &&
        config.sourceMappingMode != "triangle_gss_gradqf_truncated" &&
        config.sourceMappingMode != "element_vertex_box_measure") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.source_mapping_mode must be 'node_F_node_alpha_node_G', "
            "'edge_F_edge_alpha_edge_G_to_node', 'cell_F_cell_alpha_cell_G_to_node', "
            "'nodal_eparallel_p1', 'triangle_gss_gradqf_truncated', or "
            "'element_vertex_box_measure'.");
    }
    if (config.sourceJacobianMode != "local_ad" &&
        config.sourceJacobianMode != "finite_difference" &&
        config.sourceJacobianMode != "frozen") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.source_jacobian must be "
            "'local_ad', 'finite_difference', or 'frozen'.");
    }
    if (config.sourceMappingMode == "nodal_eparallel_p1" &&
        (config.couplingMode != "postprocess_only" ||
         config.drivingForce != "eparallel" ||
         config.generation != "current_density" ||
         config.currentApproximation != "nodal_vector_current_reconstructed" ||
         config.eparallelFieldRecovery != "nodal_vertex_star")) {
        throw std::invalid_argument(
            std::string(context) +
            ": nodal_eparallel_p1 requires postprocess_only current-density "
            "eparallel with nodal_vector_current_reconstructed and "
            "eparallel_field_recovery='nodal_vertex_star'.");
    }
    const bool elementEdgeGssLauxRequested =
        config.currentApproximation == "element_edge_sg_gss_laux" ||
        config.sourceMappingMode == "element_vertex_box_measure";
    if (elementEdgeGssLauxRequested &&
        (config.generation != "current_density" ||
         (config.drivingForce != "quasi_fermi_gradient" &&
          config.drivingForce != "electric_field") ||
         config.currentApproximation != "element_edge_sg_gss_laux" ||
         (config.drivingForce == "quasi_fermi_gradient" &&
          config.quasiFermiGradientDiscretization != "cell_gradient") ||
         config.sourceMappingMode != "element_vertex_box_measure")) {
        throw std::invalid_argument(
            std::string(context) +
            ": element_edge_sg_gss_laux requires the canonical element-box configuration.");
    }
    if (config.sourceMappingMode == "triangle_gss_gradqf_truncated" &&
        (config.generation != "current_density" ||
         config.drivingForce != "quasi_fermi_gradient" ||
         config.currentApproximation != "cell_reconstructed" ||
         config.currentMagnitudeMode != "edge_scalar_abs" ||
         config.cellReconstructedMidpointDensity != "gss_logistic" ||
         config.quasiFermiGradientDiscretization != "cell_gradient" ||
         config.sourceVolumePolicy != "genius_truncated" ||
         config.sourceVolumeFactor != 0.0 ||
         config.sourceGeometryScale != 1.0 ||
         config.edgeSourcePartition != "symmetric" ||
         config.drivingForceInterpolation != "none" ||
         config.minimumField != 0.0 ||
         config.electronDrivingForceRefDensity != 0.0 ||
         config.holeDrivingForceRefDensity != 0.0)) {
        throw std::invalid_argument(
            std::string(context) +
            ": triangle_gss_gradqf_truncated requires the canonical GSS GradQf configuration.");
    }
    if (config.cellReconstructedMidpointDensity == "gss_logistic" &&
        config.sourceMappingMode != "triangle_gss_gradqf_truncated") {
        throw std::invalid_argument(
            std::string(context) +
            ": cell_reconstructed_midpoint_density='gss_logistic' requires "
            "source_mapping_mode='triangle_gss_gradqf_truncated'.");
    }
    if (config.edgeSourcePartition != "symmetric" &&
        config.edgeSourcePartition != "qf_gradient") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.edge_source_partition must be 'symmetric' or 'qf_gradient'.");
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
    if (!std::isfinite(config.aScale) || config.aScale <= 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.A_scale must be positive and finite.");
    }
    if (config.aScale != 1.0 && config.model != "van_overstraeten") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.A_scale requires model='van_overstraeten'.");
    }
    if (!std::isfinite(config.bScale) || config.bScale <= 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.B_scale must be positive and finite.");
    }
    if (config.bScale != 1.0 && config.model != "van_overstraeten") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.B_scale requires model='van_overstraeten'.");
    }
    if (config.debugRawVanOverstraeten && config.model != "van_overstraeten") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.debug_raw_vanoverstraeten requires model='van_overstraeten'.");
    }
}

inline bool usesTriangleGssAvalancheSource(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.sourceMappingMode == "triangle_gss_gradqf_truncated";
}

inline bool usesElementEdgeGssLauxAvalancheSource(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "element_edge_sg_gss_laux" &&
           config.sourceMappingMode == "element_vertex_box_measure";
}

inline bool usesDensityGradientAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "density_gradient";
}

inline bool usesCellReconstructedAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "cell_reconstructed";
}

inline bool usesPsiGradientProxyAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "psi_gradient_proxy";
}

inline bool usesCellCurrentReconstructedAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "cell_current_reconstructed";
}

inline bool usesCellVectorCurrentReconstructedAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "cell_vector_current_reconstructed";
}

inline bool usesNodalVectorCurrentReconstructedAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "nodal_vector_current_reconstructed";
}

/// Avalanche source driven by the conserved total-current magnitude |F_p-F_n|
/// on each edge instead of the per-carrier local-density SG flux. The total
/// charge current is divergence-free in the converged state, so it does not
/// collapse on the depleted side of a reverse-biased junction where the
/// per-carrier flux (and hence the generation seed) otherwise vanishes.
inline bool usesConservedTotalCurrentAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "conserved_total_current";
}

inline Real conservedTotalCurrentFluxMagnitude(Real electronContinuityFlux,
                                               Real holeContinuityFlux)
{
    return std::abs(holeContinuityFlux - electronContinuityFlux);
}

inline Real reconstructedAvalancheCurrentDensityMagnitude(Real mobility,
                                                         Real carrierDensity,
                                                         Real drivingField)
{
    if (mobility <= 0.0)
        return 0.0;
    return mobility * std::max(carrierDensity, 0.0) * std::abs(drivingField);
}

inline Real selectAvalancheCurrentFluxProxy(
    const ImpactIonizationModelConfig& config,
    Real rawFluxMagnitude,
    Real reconstructedFluxMagnitude,
    Real mobility,
    Real midpointDensity,
    Real impactField,
    Real electricField,
    Real conservedTotalFluxMagnitude)
{
    const bool reconstructedCurrent =
        config.currentMagnitudeMode == "dual_face_vector_mag" ||
        usesCellCurrentReconstructedAvalancheCurrent(config) ||
        usesCellVectorCurrentReconstructedAvalancheCurrent(config) ||
        usesNodalVectorCurrentReconstructedAvalancheCurrent(config);
    if (reconstructedCurrent)
        return reconstructedFluxMagnitude;
    if (usesPsiGradientProxyAvalancheCurrent(config)) {
        return reconstructedAvalancheCurrentDensityMagnitude(
            mobility, midpointDensity, electricField);
    }
    if (usesCellReconstructedAvalancheCurrent(config)) {
        return reconstructedAvalancheCurrentDensityMagnitude(
            mobility, midpointDensity, impactField);
    }
    if (usesConservedTotalCurrentAvalancheCurrent(config))
        return conservedTotalFluxMagnitude;
    return rawFluxMagnitude;
}

/// Fermi/logistic weight aux2(x) = 1 / (1 + exp(x)); numerically stable.
inline Real avalancheMidpointAux2(Real x)
{
    if (x >= 0.0)
        return 1.0 / (1.0 + std::exp(x));
    const Real ex = std::exp(x);
    return 1.0 / (1.0 + ex);
}

/// Bernoulli/exponentially weighted edge-midpoint carrier density:
///   n_mid = n_i * aux2((V_i - V_j) / (2 Vt)) + n_j * aux2((V_j - V_i) / (2 Vt))
/// with aux2(x) = 1/(1+exp(x)); the two weights sum to 1. For electrons pass the
/// electrostatic potentials as (V_i, V_j); for holes swap them so the potential
/// enters with the opposite sign.
inline Real bernoulliWeightedMidpointDensity(Real density_i,
                                             Real density_j,
                                             Real potential_i,
                                             Real potential_j,
                                             Real Vt)
{
    if (Vt <= 0.0)
        return 0.5 * (density_i + density_j);
    const Real arg = (potential_i - potential_j) / (2.0 * Vt);
    const Real weight_i = avalancheMidpointAux2(arg);
    const Real weight_j = avalancheMidpointAux2(-arg);
    return density_i * weight_i + density_j * weight_j;
}

inline Real cellReconstructedAvalancheMidpointDensity(
    const ImpactIonizationModelConfig& config,
    Real density_i,
    Real density_j,
    Real potential_i,
    Real potential_j,
    Real Vt)
{
    if (config.cellReconstructedMidpointDensity == "arithmetic")
        return 0.5 * (density_i + density_j);
    return bernoulliWeightedMidpointDensity(density_i, density_j, potential_i, potential_j, Vt);
}

inline Real gssElectronAvalancheMidpointDensity(
    Real density_i,
    Real density_j,
    Real potential_i,
    Real potential_j,
    Real Vt)
{
    if (Vt <= 0.0)
        return 0.5 * (density_i + density_j);
    const Real arg = (potential_j - potential_i) / (2.0 * Vt);
    return density_i * avalancheMidpointAux2(arg) +
           density_j * avalancheMidpointAux2(-arg);
}

inline Real gssHoleAvalancheMidpointDensity(
    Real density_i,
    Real density_j,
    Real potential_i,
    Real potential_j,
    Real Vt)
{
    if (Vt <= 0.0)
        return 0.5 * (density_i + density_j);
    const Real arg = (potential_i - potential_j) / (2.0 * Vt);
    return density_i * avalancheMidpointAux2(arg) +
           density_j * avalancheMidpointAux2(-arg);
}

inline bool usesEdgeCurrentAvalancheSource(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           (config.currentApproximation == "density_gradient" ||
            config.currentApproximation == "grad_qf" ||
            config.currentApproximation == "cell_reconstructed" ||
            config.currentApproximation == "psi_gradient_proxy" ||
            config.currentApproximation == "cell_current_reconstructed" ||
            config.currentApproximation == "cell_vector_current_reconstructed" ||
            config.currentApproximation == "nodal_vector_current_reconstructed" ||
            config.currentApproximation == "element_edge_sg_gss_laux" ||
            config.currentApproximation == "conserved_total_current");
}

inline bool usesDirectionalEdgeAvalancheSourcePartition(
    const ImpactIonizationModelConfig& config)
{
    return usesEdgeCurrentAvalancheSource(config) &&
           (config.currentApproximation == "grad_qf" ||
            config.edgeSourcePartition == "qf_gradient");
}

inline bool usesQuasiFermiCarrierTruncation(const ImpactIonizationModelConfig& config)
{
    return !config.debugRawVanOverstraeten &&
           config.quasiFermiCarrierTruncation > 0.0;
}

inline bool usesCellGradientQuasiFermiAvalancheDrive(
    const ImpactIonizationModelConfig& config)
{
    return usesQuasiFermiAvalancheDrivingForce(config) &&
           config.quasiFermiGradientDiscretization == "cell_gradient";
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
inline bool cellTouchesContact(const DeviceMesh& mesh, const Cell& cell)
{
    for (const Contact& contact : mesh.contacts()) {
        for (Index contactNode : contact.node_ids) {
            for (Index cellNode : cell.node_ids) {
                if (contactNode == cellNode)
                    return true;
            }
        }
    }
    return false;
}

inline bool cellTouchesContact(
    const Cell& cell,
    const std::vector<bool>& contactNodes)
{
    return std::any_of(
        cell.node_ids.begin(), cell.node_ids.end(),
        [&](Index node) {
            return node < contactNodes.size() && contactNodes[node];
        });
}

inline bool cellTouchesContactBoundaryFace(
    const DeviceMesh& mesh,
    const Cell& cell)
{
    for (const Contact& contact : mesh.contacts()) {
        std::size_t contactVertexCount = 0;
        for (Index cellNode : cell.node_ids) {
            for (Index contactNode : contact.node_ids) {
                if (cellNode == contactNode) {
                    ++contactVertexCount;
                    break;
                }
            }
        }
        if (contactVertexCount >= 2)
            return true;
    }
    return false;
}

struct ContactBoundaryFaceGeometry {
    std::array<std::size_t, 2> faceLocalNodes{};
    std::size_t oppositeLocalNode = 0;
    Point2 unitNormal = Point2::Zero();
    Point2 faceMidpoint = Point2::Zero();
    Real oppositeToFaceMidpointDistance = 0.0;
    Real normalHeight = 0.0;
    Real centroidElectricWeight = 0.0;
};

inline bool contactBoundaryFaceGeometry(
    const DeviceMesh& mesh,
    const Cell& cell,
    ContactBoundaryFaceGeometry& geometry)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
        return false;
    for (const Contact& contact : mesh.contacts()) {
        std::array<bool, 3> onContact{};
        for (std::size_t local = 0; local < 3; ++local) {
            onContact[local] = std::find(
                contact.node_ids.begin(), contact.node_ids.end(),
                cell.node_ids[local]) != contact.node_ids.end();
        }
        for (std::size_t local = 0; local < 3; ++local) {
            const std::size_t next = (local + 1) % 3;
            if (!onContact[local] || !onContact[next])
                continue;
            if (!contact.edge_node_ids.empty()) {
                const Index first = cell.node_ids[local];
                const Index second = cell.node_ids[next];
                const bool exactEdge = std::any_of(
                    contact.edge_node_ids.begin(), contact.edge_node_ids.end(),
                    [&](const std::array<Index, 2>& edge) {
                        return (edge[0] == first && edge[1] == second) ||
                            (edge[0] == second && edge[1] == first);
                    });
                if (!exactEdge)
                    continue;
            }
            const std::size_t opposite = (local + 2) % 3;
            const Point2 face0 = meshPoint(mesh, cell.node_ids[local]);
            const Point2 face1 = meshPoint(mesh, cell.node_ids[next]);
            const Point2 oppositePoint =
                meshPoint(mesh, cell.node_ids[opposite]);
            const Point2 tangent = face1 - face0;
            const Real faceLength = tangent.norm();
            if (faceLength <= 1.0e-30)
                return false;
            geometry.faceLocalNodes = {local, next};
            geometry.oppositeLocalNode = opposite;
            geometry.unitNormal =
                Point2{-tangent.y() / faceLength, tangent.x() / faceLength};
            geometry.faceMidpoint = 0.5 * (face0 + face1);
            const Point2 oppositeOffset =
                oppositePoint - geometry.faceMidpoint;
            geometry.oppositeToFaceMidpointDistance = oppositeOffset.norm();
            geometry.normalHeight =
                std::abs(oppositeOffset.dot(geometry.unitNormal));
            if (geometry.oppositeToFaceMidpointDistance <= 1.0e-30 ||
                geometry.normalHeight <= 1.0e-30) {
                return false;
            }
            const Point2 centroid = (face0 + face1 + oppositePoint) / 3.0;
            const Real centroidDistance =
                std::abs((centroid - geometry.faceMidpoint)
                             .dot(geometry.unitNormal));
            geometry.centroidElectricWeight = std::clamp(
                1.0 - centroidDistance / geometry.normalHeight, 0.0, 1.0);
            return true;
        }
    }
    return false;
}

inline bool cellUsesContactElectricFieldFallback(
    const ImpactIonizationModelConfig& config,
    const DeviceMesh& mesh,
    const Cell& cell)
{
    return config.contactElectricFieldFallback &&
           config.drivingForce == "quasi_fermi_gradient" &&
           (config.contactElectricFieldFallbackScope == "contact_boundary_face"
                ? cellTouchesContactBoundaryFace(mesh, cell)
                : cellTouchesContact(mesh, cell));
}

inline Real contactElectricFallbackImpactField(
    const ImpactIonizationModelConfig& config,
    const DeviceMesh& mesh,
    const Cell& cell,
    const Point2& electricGradient,
    Real quasiFermiField,
    const std::function<Real(Index)>& potential,
    Real fieldFactor)
{
    const Real electricMagnitude = electricGradient.norm() * fieldFactor;
    if (config.contactElectricFieldFallbackMode == "cell_gradient_magnitude")
        return electricMagnitude;

    ContactBoundaryFaceGeometry geometry;
    if (!contactBoundaryFaceGeometry(mesh, cell, geometry))
        return quasiFermiField;
    const Real faceNormalField =
        std::abs(electricGradient.dot(geometry.unitNormal)) * fieldFactor;
    if (config.contactElectricFieldFallbackMode == "face_normal")
        return faceNormalField;

    const Index face0 = cell.node_ids[geometry.faceLocalNodes[0]];
    const Index face1 = cell.node_ids[geometry.faceLocalNodes[1]];
    const Index opposite = cell.node_ids[geometry.oppositeLocalNode];
    const Real oneSidedField = std::abs(
        potential(opposite) - 0.5 * (potential(face0) + potential(face1))) /
        geometry.oppositeToFaceMidpointDistance * fieldFactor;
    if (config.contactElectricFieldFallbackMode == "one_sided")
        return oneSidedField;

    return geometry.centroidElectricWeight * faceNormalField +
           (1.0 - geometry.centroidElectricWeight) * quasiFermiField;
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

struct NodalScalarGradientCache {
    std::vector<Point2> gradients;
    std::vector<bool> valid;
};

inline bool isTransportMaterial(const Material& material);

/// Recover a P1 nodal gradient over the complete transport-material vertex
/// star.  Sentaurus exposes vertex fields after this wider recovery rather
/// than retaining the two-cell stencil local to an individual edge.
inline NodalScalarGradientCache computeTransportNodalScalarGradientCache(
    const DeviceMesh&                      mesh,
    const std::vector<std::vector<Index>>& nodeCells,
    const std::vector<Material>&           cellMaterials,
    const CellScalarGradientCache&         cellCache)
{
    NodalScalarGradientCache output;
    output.gradients.assign(mesh.numNodes(), Point2::Zero());
    output.valid.assign(mesh.numNodes(), false);
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        if (node >= nodeCells.size())
            continue;
        Point2 weightedGradient = Point2::Zero();
        Real totalArea = 0.0;
        for (const Index cellId : nodeCells[node]) {
            if (cellId >= cellMaterials.size() ||
                cellId >= cellCache.valid.size() ||
                !cellCache.valid[cellId] ||
                !isTransportMaterial(cellMaterials[cellId])) {
                continue;
            }
            const Real area = cellCache.areas[cellId];
            if (area <= 0.0)
                continue;
            weightedGradient += area * cellCache.gradients[cellId];
            totalArea += area;
        }
        if (totalArea > 0.0) {
            output.gradients[node] = weightedGradient / totalArea;
            output.valid[node] = true;
        }
    }
    return output;
}

inline Point2 edgeAveragedNodalScalarGradient(
    const DeviceMesh&                mesh,
    Index                            edgeId,
    const NodalScalarGradientCache&  cache,
    bool&                            valid)
{
    valid = false;
    if (edgeId >= mesh.numEdges())
        return Point2::Zero();
    const Edge& edge = mesh.getEdge(edgeId);
    Point2 gradient = Point2::Zero();
    Real count = 0.0;
    for (const Index node : {edge.n0, edge.n1}) {
        if (node >= cache.valid.size() || !cache.valid[node])
            continue;
        gradient += cache.gradients[node];
        count += 1.0;
    }
    if (count <= 0.0)
        return Point2::Zero();
    valid = true;
    return gradient / count;
}

inline Point2 eparallelElectricGradientForEdge(
    const ImpactIonizationModelConfig&       config,
    const std::vector<std::vector<Index>>&   edgeCells,
    const DeviceMesh&                        mesh,
    Index                                    edgeId,
    const CellScalarGradientCache&           cellCache,
    const NodalScalarGradientCache&          nodalCache,
    bool&                                    valid)
{
    if (config.eparallelFieldRecovery == "nodal_vertex_star") {
        return edgeAveragedNodalScalarGradient(
            mesh, edgeId, nodalCache, valid);
    }
    return edgeAveragedCellScalarGradient(
        edgeCells, edgeId, cellCache, valid);
}

inline std::vector<Real> computeNodeCellGradientMagnitudes(
    const std::vector<std::vector<Index>>& nodeCells,
    const CellScalarGradientCache&         cache)
{
    std::vector<Point2> gradients(nodeCells.size(), Point2::Zero());
    for (std::size_t node = 0; node < nodeCells.size(); ++node) {
        Point2 weightedGradient = Point2::Zero();
        Real totalArea = 0.0;
        for (const Index cellId : nodeCells[node]) {
            if (cellId >= cache.valid.size() || !cache.valid[cellId])
                continue;
            const Real area = cache.areas[cellId];
            if (area <= 0.0)
                continue;
            weightedGradient += area * cache.gradients[cellId];
            totalArea += area;
        }
        if (totalArea > 0.0)
            gradients[node] = weightedGradient / totalArea;
    }
    std::vector<Real> fields(nodeCells.size(), 0.0);
    for (std::size_t node = 0; node < nodeCells.size(); ++node)
        fields[node] = gradients[node].norm();
    return fields;
}

template <typename ValueAt>
inline std::vector<Real> computeNodeCellGradientMagnitudes(
    const DeviceMesh&                      mesh,
    const std::vector<std::vector<Index>>& nodeCells,
    ValueAt&&                              valueAt)
{
    return computeNodeCellGradientMagnitudes(
        nodeCells, computeCellScalarGradientCache(mesh, valueAt));
}

template <typename ValueAt>
inline std::vector<Point2> computeNodeWeightedLeastSquaresGradients(
    const DeviceMesh&                      mesh,
    const std::vector<std::vector<Index>>& nodeCells,
    ValueAt&&                              valueAt)
{
    std::vector<std::unordered_set<Index>> nodeNeighbors(mesh.numNodes());
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        nodeNeighbors[edge.n0].insert(edge.n1);
        nodeNeighbors[edge.n1].insert(edge.n0);
    }

    std::vector<Point2> fields(mesh.numNodes(), Point2::Zero());
    const CellScalarGradientCache fallbackCache = computeCellScalarGradientCache(
        mesh, [&](Index node) { return valueAt(node); });

    for (Index nodeId = 0; nodeId < mesh.numNodes(); ++nodeId) {
        const Node& center = mesh.getNode(nodeId);
        const Real centerValue = valueAt(nodeId);
        Real sxx = 0.0;
        Real sxy = 0.0;
        Real syy = 0.0;
        Real sxv = 0.0;
        Real syv = 0.0;

        for (const Index neighborId : nodeNeighbors[nodeId]) {
            const Node& neighbor = mesh.getNode(neighborId);
            const Real dx = neighbor.x - center.x;
            const Real dy = neighbor.y - center.y;
            const Real distance = std::hypot(dx, dy);
            if (distance <= 1.0e-30)
                continue;
            const Real weight = 1.0 / distance;
            const Real dv = valueAt(neighborId) - centerValue;
            sxx += weight * dx * dx;
            sxy += weight * dx * dy;
            syy += weight * dy * dy;
            sxv += weight * dx * dv;
            syv += weight * dy * dv;
        }

        const Real det = sxx * syy - sxy * sxy;
        if (std::abs(det) <= 1.0e-60) {
            Point2 weightedGradient = Point2::Zero();
            Real totalArea = 0.0;
            for (const Index cellId : nodeCells[nodeId]) {
                if (cellId >= fallbackCache.valid.size() || !fallbackCache.valid[cellId])
                    continue;
                const Real area = fallbackCache.areas[cellId];
                if (area <= 0.0)
                    continue;
                weightedGradient += area * fallbackCache.gradients[cellId];
                totalArea += area;
            }
            if (totalArea > 0.0)
                fields[nodeId] = weightedGradient / totalArea;
            continue;
        }

        const Real gradX = (sxv * syy - syv * sxy) / det;
        const Real gradY = (sxx * syv - sxy * sxv) / det;
        fields[nodeId] = Point2{gradX, gradY};
    }
    return fields;
}

template <typename ValueAt>
inline std::vector<Real> computeNodeWeightedLeastSquaresGradientMagnitudes(
    const DeviceMesh&                      mesh,
    const std::vector<std::vector<Index>>& nodeCells,
    ValueAt&&                              valueAt)
{
    const std::vector<Point2> gradients = computeNodeWeightedLeastSquaresGradients(
        mesh, nodeCells, std::forward<ValueAt>(valueAt));
    std::vector<Real> fields(gradients.size(), 0.0);
    for (std::size_t node = 0; node < gradients.size(); ++node)
        fields[node] = gradients[node].norm();
    return fields;
}

inline std::vector<Real> computeNodeElectricFields(const VectorXd& psi, const DeviceMesh& mesh, Real fieldFactor)
{
    const std::vector<std::vector<Index>> nodeCells = buildNodeCellMap(mesh);
    std::vector<Real> fields = computeNodeWeightedLeastSquaresGradientMagnitudes(
        mesh, nodeCells, [&](Index node) { return psi(static_cast<int>(node)); });
    for (Real& field : fields)
        field *= fieldFactor;
    return fields;
}

inline Real bandToBandGenerationRateInternal(
    const BandToBandTunnelingModel& model,
    const PhysicalUnitSystem& unitSystem,
    Real electricFieldInternal)
{
    const Real field_V_per_m =
        unitSystem.internalElectricFieldToVPerM(electricFieldInternal);
    const Real generation_m3_per_s = model.generationRate(field_V_per_m);
    return unitSystem.m3ToInternalConcentration(generation_m3_per_s);
}

inline bool isTransportMaterial(const Material& material)
{
    return material.ni > 0.0 || material.mun > 0.0 || material.mup > 0.0;
}

inline std::vector<Real> transportCellVectorEdgeGradientMagnitudes(
    const DeviceMesh& mesh,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<Material>& cellMaterials,
    const VectorXd& values,
    Real fieldFactor)
{
    const CellScalarGradientCache gradients = computeCellScalarGradientCache(
        mesh, [&](Index node) { return values(static_cast<int>(node)); });
    std::vector<Real> fields(mesh.numEdges(), 0.0);
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        if (edgeId >= edgeCells.size())
            continue;
        Point2 weightedGradient = Point2::Zero();
        Real totalArea = 0.0;
        for (const Index cellId : edgeCells[edgeId]) {
            if (cellId >= cellMaterials.size() ||
                !isTransportMaterial(cellMaterials[cellId]) ||
                cellId >= gradients.valid.size() || !gradients.valid[cellId] ||
                gradients.areas[cellId] <= 0.0) {
                continue;
            }
            const Real area = gradients.areas[cellId];
            weightedGradient += area * gradients.gradients[cellId];
            totalArea += area;
        }
        if (totalArea > 0.0)
            fields[edgeId] = (weightedGradient / totalArea).norm() * fieldFactor;
    }
    return fields;
}

template <typename ValueAt>
inline Real contactCellElectricFieldMagnitudeForMobility(
    const MobilityModelConfig&               config,
    const DeviceMesh&                        mesh,
    const std::vector<std::vector<Index>>&   edgeCells,
    const std::vector<Material>&             cellMaterials,
    Index                                    edgeId,
    ValueAt&&                                potentialAt,
    Real                                     fieldFactor,
    Real                                     fallbackField,
    bool*                                    applied = nullptr,
    const std::vector<bool>*                 contactNodes = nullptr)
{
    if (applied != nullptr)
        *applied = false;
    if (!config.contactElectricFieldFallback ||
        config.highFieldDrivingForce != "quasi_fermi_gradient" ||
        edgeId >= edgeCells.size()) {
        return fallbackField;
    }

    Point2 weightedGradient = Point2::Zero();
    Real totalArea = 0.0;
    for (const Index cellId : edgeCells[edgeId]) {
        if (cellId >= cellMaterials.size() ||
            !isTransportMaterial(cellMaterials[cellId])) {
            continue;
        }
        const Cell& cell = mesh.getCell(cellId);
        const bool touchesContact = contactNodes != nullptr
            ? cellTouchesContact(cell, *contactNodes)
            : cellTouchesContact(mesh, cell);
        if (!touchesContact)
            continue;
        bool valid = false;
        Real area = 0.0;
        const Point2 gradient = cellScalarGradient(
            mesh, cell, potentialAt, valid, area);
        if (!valid || area <= 0.0)
            continue;
        weightedGradient += area * gradient;
        totalArea += area;
    }
    if (totalArea <= 0.0)
        return fallbackField;
    if (applied != nullptr)
        *applied = true;
    return (weightedGradient / totalArea).norm() * fieldFactor;
}

inline std::vector<Real> contactCellElectricFieldEdgeMagnitudesForMobility(
    const MobilityModelConfig&               config,
    const DeviceMesh&                        mesh,
    const std::vector<std::vector<Index>>&   edgeCells,
    const std::vector<Material>&             cellMaterials,
    const VectorXd&                          potential,
    Real                                     fieldFactor)
{
    std::vector<Real> fields(mesh.numEdges(), 0.0);
    if (!config.contactElectricFieldFallback ||
        config.highFieldDrivingForce != "quasi_fermi_gradient") {
        return fields;
    }
    const std::vector<bool> contactNodes = contactNodeMask(mesh);
    const CellScalarGradientCache gradients = computeCellScalarGradientCache(
        mesh, [&](Index node) { return potential(static_cast<int>(node)); });
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        if (edgeId >= edgeCells.size())
            continue;
        Point2 weightedGradient = Point2::Zero();
        Real totalArea = 0.0;
        for (const Index cellId : edgeCells[edgeId]) {
            if (cellId >= cellMaterials.size() ||
                !isTransportMaterial(cellMaterials[cellId]) ||
                cellId >= gradients.valid.size() || !gradients.valid[cellId] ||
                gradients.areas[cellId] <= 0.0 ||
                !cellTouchesContact(mesh.getCell(cellId), contactNodes)) {
                continue;
            }
            const Real area = gradients.areas[cellId];
            weightedGradient += area * gradients.gradients[cellId];
            totalArea += area;
        }
        if (totalArea > 0.0)
            fields[edgeId] = (weightedGradient / totalArea).norm() * fieldFactor;
    }
    return fields;
}

inline Real mobilityHighFieldDrivingField(
    const MobilityModelConfig& config,
    Index edgeId,
    Real quasiFermiField,
    Real electricField,
    const std::vector<Real>& contactElectricFields)
{
    if (config.highFieldDrivingForce != "quasi_fermi_gradient")
        return electricField;
    if (config.contactElectricFieldFallback &&
        edgeId < contactElectricFields.size() &&
        contactElectricFields[edgeId] > 0.0) {
        return contactElectricFields[edgeId];
    }
    return quasiFermiField;
}

/**
 * Integrate local BTBT generation over semiconductor triangles and lump one
 * third of each cell source to its vertices.  Computing the field per
 * semiconductor cell is important at shared Si/oxide nodes: an unrestricted
 * nodal least-squares gradient mixes the oxide field into a silicon-only
 * material model and can exponentially over-predict E2 generation.
 */
inline std::vector<Real> bandToBandGenerationNodeSourceIntegrals(
    const BandToBandTunnelingModel& model,
    const PhysicalUnitSystem&       unitSystem,
    const DeviceMesh&               mesh,
    const std::vector<Material>&    cellMaterials,
    const VectorXd&                 psi,
    Real                            fieldFactor)
{
    std::vector<Real> sources(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    if (!model.enabled())
        return sources;

    const CellScalarGradientCache gradients = computeCellScalarGradientCache(
        mesh, [&](Index node) { return psi(static_cast<int>(node)); });
    if (model.config().sourceIntegration == "transport_node_lumped") {
        std::vector<std::unordered_set<Index>> transportNeighbors(mesh.numNodes());
        std::vector<std::vector<Index>> transportNodeCells(mesh.numNodes());
        std::vector<Real> transportNodeAreas(mesh.numNodes(), 0.0);

        for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
            if (cellId >= cellMaterials.size() ||
                !isTransportMaterial(cellMaterials[cellId]) ||
                cellId >= gradients.valid.size() || !gradients.valid[cellId] ||
                gradients.areas[cellId] <= 0.0) {
                continue;
            }
            const Cell& cell = mesh.getCell(cellId);
            if (cell.node_ids.empty())
                continue;
            const Real lumpedArea = gradients.areas[cellId] /
                static_cast<Real>(cell.node_ids.size());
            for (const Index node : cell.node_ids) {
                transportNodeCells[node].push_back(cellId);
                transportNodeAreas[node] += lumpedArea;
                for (const Index neighbor : cell.node_ids) {
                    if (neighbor != node)
                        transportNeighbors[node].insert(neighbor);
                }
            }
        }

        for (Index nodeId = 0; nodeId < mesh.numNodes(); ++nodeId) {
            if (transportNodeAreas[nodeId] <= 0.0)
                continue;
            const Node& center = mesh.getNode(nodeId);
            const Real centerValue = psi(static_cast<int>(nodeId));
            Real sxx = 0.0;
            Real sxy = 0.0;
            Real syy = 0.0;
            Real sxv = 0.0;
            Real syv = 0.0;
            for (const Index neighborId : transportNeighbors[nodeId]) {
                const Node& neighbor = mesh.getNode(neighborId);
                const Real dx = neighbor.x - center.x;
                const Real dy = neighbor.y - center.y;
                const Real distance = std::hypot(dx, dy);
                if (distance <= 1.0e-30)
                    continue;
                const Real weight = 1.0 / distance;
                const Real dv = psi(static_cast<int>(neighborId)) - centerValue;
                sxx += weight * dx * dx;
                sxy += weight * dx * dy;
                syy += weight * dy * dy;
                sxv += weight * dx * dv;
                syv += weight * dy * dv;
            }

            Point2 gradient = Point2::Zero();
            const Real det = sxx * syy - sxy * sxy;
            if (std::abs(det) > 1.0e-60) {
                gradient.x() = (sxv * syy - syv * sxy) / det;
                gradient.y() = (sxx * syv - sxy * sxv) / det;
            } else {
                Real totalArea = 0.0;
                for (const Index cellId : transportNodeCells[nodeId]) {
                    const Real area = gradients.areas[cellId];
                    gradient += area * gradients.gradients[cellId];
                    totalArea += area;
                }
                if (totalArea > 0.0)
                    gradient /= totalArea;
            }
            const Real generation = bandToBandGenerationRateInternal(
                model, unitSystem, gradient.norm() * fieldFactor);
            sources[nodeId] = generation * transportNodeAreas[nodeId];
        }
        return sources;
    }

    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        if (cellId >= cellMaterials.size() ||
            !isTransportMaterial(cellMaterials[cellId]) ||
            cellId >= gradients.valid.size() || !gradients.valid[cellId]) {
            continue;
        }
        const Cell& cell = mesh.getCell(cellId);
        if (cell.node_ids.empty() || gradients.areas[cellId] <= 0.0)
            continue;
        const Real generation = bandToBandGenerationRateInternal(
            model, unitSystem, gradients.gradients[cellId].norm() * fieldFactor);
        const Real lumped = generation * gradients.areas[cellId] /
            static_cast<Real>(cell.node_ids.size());
        for (Index node : cell.node_ids)
            sources[node] += lumped;
    }
    return sources;
}

inline std::vector<Real> transportNodeLumpedAreas(
    const DeviceMesh& mesh,
    const std::vector<Material>& cellMaterials)
{
    std::vector<Real> areas(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    const CellScalarGradientCache geometry = computeCellScalarGradientCache(
        mesh, [](Index) { return 0.0; });
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        if (cellId >= cellMaterials.size() ||
            !isTransportMaterial(cellMaterials[cellId]) ||
            cellId >= geometry.valid.size() || !geometry.valid[cellId]) {
            continue;
        }
        const Cell& cell = mesh.getCell(cellId);
        if (cell.node_ids.empty())
            continue;
        const Real lumped = geometry.areas[cellId] /
            static_cast<Real>(cell.node_ids.size());
        for (Index node : cell.node_ids)
            areas[node] += lumped;
    }
    return areas;
}

inline Real edgeQuasiFermiCoefficientField(
    const ImpactIonizationModelConfig&     config,
    Real                                   edgeQfField,
    Real                                   electricField,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    const std::vector<bool>&               contactNodes,
    const CellScalarGradientCache&         qfGradientCache,
    Real                                   fieldFactor)
{
    if (usesCellGradientQuasiFermiAvalancheDrive(config)) {
        bool validGradient = false;
        const Point2 gradient = edgeAveragedCellScalarGradient(
            edgeCells, edgeId, qfGradientCache, validGradient);
        return validGradient ? gradient.norm() * fieldFactor : edgeQfField;
    }
    if (config.debugRawVanOverstraeten)
        return edgeQfField;
    return edgeHighFieldDrivingField(
        true, edgeQfField, electricField, edgeCells, mesh, edgeId, contactNodes);
}

inline std::vector<Real> computeElectronAvalancheNodeQuasiFermiDrivingFields(
    const ImpactIonizationModelConfig&     config,
    const DeviceMesh&                      mesh,
    const std::vector<std::vector<Index>>& nodeCells,
    const VectorXd&                        psi,
    const VectorXd&                        phin,
    const VectorXd&                        n,
    const std::vector<Real>&               ni,
    Real                                   Vt,
    Real                                   fieldFactor = 1.0)
{
    std::vector<Real> fields = !usesCellGradientQuasiFermiAvalancheDrive(config)
        ? computeNodeScalarGradientMagnitudes(phin, mesh)
        : computeNodeCellGradientMagnitudes(
            mesh, nodeCells, [&](Index node) {
                const int idx = static_cast<int>(node);
                return electronQfForAvalancheGradient(
                    psi(idx), phin(idx), n(idx), ni[node], Vt, config);
            });
    for (Real& field : fields)
        field *= fieldFactor;
    return fields;
}

inline std::vector<Real> computeHoleAvalancheNodeQuasiFermiDrivingFields(
    const ImpactIonizationModelConfig&     config,
    const DeviceMesh&                      mesh,
    const std::vector<std::vector<Index>>& nodeCells,
    const VectorXd&                        psi,
    const VectorXd&                        phip,
    const VectorXd&                        p,
    const std::vector<Real>&               ni,
    Real                                   Vt,
    Real                                   fieldFactor = 1.0)
{
    std::vector<Real> fields = !usesCellGradientQuasiFermiAvalancheDrive(config)
        ? computeNodeScalarGradientMagnitudes(phip, mesh)
        : computeNodeCellGradientMagnitudes(
            mesh, nodeCells, [&](Index node) {
                const int idx = static_cast<int>(node);
                return holeQfForAvalancheGradient(
                    psi(idx), phip(idx), p(idx), ni[node], Vt, config);
            });
    for (Real& field : fields)
        field *= fieldFactor;
    return fields;
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

struct TriangleGssAvalancheSourceRecord {
    Index cellId = 0;
    int localEdge = -1;
    Index edgeId = 0;
    Index node0 = 0;
    Index node1 = 0;
    Real edgeLength = 0.0;
    Real truncatedPartialVolume = 0.0;
    Real electronCellQfField = 0.0;
    Real holeCellQfField = 0.0;
    Real electronImpactField = 0.0;
    Real holeImpactField = 0.0;
    Real electronEdgeQfField = 0.0;
    Real holeEdgeQfField = 0.0;
    Real electronMidpointDensity = 0.0;
    Real holeMidpointDensity = 0.0;
    Real electronMobilityDrivingField = 0.0;
    Real holeMobilityDrivingField = 0.0;
    Real electronLowFieldMobility = 0.0;
    Real holeLowFieldMobility = 0.0;
    Real electronMobility = 0.0;
    Real holeMobility = 0.0;
    Real electronAlpha = 0.0;
    Real holeAlpha = 0.0;
    Real electronFluxProxy = 0.0;
    Real holeFluxProxy = 0.0;
    Real electronSourceIntegral = 0.0;
    Real holeSourceIntegral = 0.0;
    Real edgeSourceIntegral = 0.0;
    Real node0SourceIntegral = 0.0;
    Real node1SourceIntegral = 0.0;
};

inline Real triangleGssEndpointAveragedMobility(
    const MobilityModelConfig&    mobilityConfig,
    const MobilityModel&          mobility,
    const DeviceMesh&             mesh,
    const DopingModel&            doping,
    const std::vector<Material>&  cellMaterials,
    const VectorXd&               n,
    const VectorXd&               p,
    Index                         cellId,
    Index                         node0,
    Index                         node1,
    CarrierType                   carrier,
    Real                          drivingField,
    const VectorXd*               psi = nullptr,
    const std::vector<Index>*     cellEdgeIds = nullptr,
    Real                          fieldFactor = 1.0)
{
    const Material& material = cellMaterials.at(static_cast<std::size_t>(cellId));
    const auto [surfaceField, surfaceDistance] =
        isSurfaceMobilityModel(mobilityConfig) &&
            cellId < mobilityConfig.surface.cellNormalFields.size() &&
            cellId < mobilityConfig.surface.cellDistances.size()
        ? std::pair<Real, Real>{
              mobilityConfig.surface.cellNormalFields[cellId],
              mobilityConfig.surface.cellDistances[cellId]}
        : isSurfaceMobilityModel(mobilityConfig) && psi != nullptr && cellEdgeIds != nullptr
        ? nearestSurfaceFieldAndDistanceForCell(
              mesh, buildEdgeCellMap(mesh), *psi, cellId, mobilityConfig,
              fieldFactor)
        : std::pair<Real, Real>{0.0, 0.0};
    const auto atNode = [&](Index node) {
        const Real mobilityDoping = nodeMobilityDopingConcentration(
            mesh, doping, node, cellId, &mobilityConfig);
        return carrier == CarrierType::Electron
            ? mobility.electronMobility(
                material, mobilityDoping, n(static_cast<int>(node)),
                p(static_cast<int>(node)), drivingField,
                surfaceField, surfaceDistance)
            : mobility.holeMobility(
                material, mobilityDoping, n(static_cast<int>(node)),
                p(static_cast<int>(node)), drivingField,
                surfaceField, surfaceDistance);
    };
    return 0.5 * (atNode(node0) + atNode(node1));
}

inline std::vector<TriangleGssAvalancheSourceRecord>
triangleGssAvalancheSourceRecordsForCell(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<Index>&          cellEdgeIds,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    Index                              cellId,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    const Cell& cell = mesh.getCell(cellId);
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3) {
        throw std::invalid_argument(
            "triangle GSS avalanche source requires Tri3 cells; unsupported cell " +
            std::to_string(cellId));
    }

    // Keep the primal triangle-GSS source on exactly the same material
    // support as its local-AD/finite-difference Jacobian.  Doping-based
    // mobility models can otherwise return a finite mobility for an
    // insulator cell whose vertices are shared with silicon, which creates a
    // nonzero oxide-side avalanche residual with no corresponding Jacobian
    // contribution.
    const Material& material =
        cellMaterials.at(static_cast<std::size_t>(cellId));
    if (material.ni <= 0.0 ||
        (material.mun <= 0.0 && material.mup <= 0.0)) {
        return {};
    }

    bool electronGradientValid = false;
    bool holeGradientValid = false;
    Real electronDoubleArea = 0.0;
    Real holeDoubleArea = 0.0;
    const Point2 electronGradient =
        cellScalarGradient(mesh, cell, [&](Index node) {
            return phin(static_cast<int>(node));
        }, electronGradientValid, electronDoubleArea);
    const Point2 holeGradient =
        cellScalarGradient(mesh, cell, [&](Index node) {
            return phip(static_cast<int>(node));
        }, holeGradientValid, holeDoubleArea);
    if (!electronGradientValid || !holeGradientValid)
        return {};

    const Real electronCellField = electronGradient.norm() * fieldFactor;
    const Real holeCellField = holeGradient.norm() * fieldFactor;
    Real electronImpactField = electronCellField;
    Real holeImpactField = holeCellField;
    if (cellUsesContactElectricFieldFallback(config, mesh, cell)) {
        bool electricGradientValid = false;
        Real electricDoubleArea = 0.0;
        const Point2 electricGradient =
            cellScalarGradient(mesh, cell, [&](Index node) {
                return psi(static_cast<int>(node));
            }, electricGradientValid, electricDoubleArea);
        if (!electricGradientValid)
            return {};
        const auto potential = [&](Index node) {
            return psi(static_cast<int>(node));
        };
        electronImpactField = contactElectricFallbackImpactField(
            config, mesh, cell, electricGradient, electronCellField,
            potential, fieldFactor);
        holeImpactField = contactElectricFallbackImpactField(
            config, mesh, cell, electricGradient, holeCellField,
            potential, fieldFactor);
    }
    const Real electronAlpha = impact.electronCoefficient(electronImpactField);
    const Real holeAlpha = impact.holeCoefficient(holeImpactField);
    std::vector<TriangleGssAvalancheSourceRecord> records;
    records.reserve(3);

    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const Index node0 = cell.node_ids[static_cast<std::size_t>(localEdge)];
        const Index node1 =
            cell.node_ids[static_cast<std::size_t>((localEdge + 1) % 3)];
        Index edgeId = mesh.numEdges();
        for (Index candidate : cellEdgeIds) {
            const Edge& edge = mesh.getEdge(candidate);
            if ((edge.n0 == node0 && edge.n1 == node1) ||
                (edge.n0 == node1 && edge.n1 == node0)) {
                edgeId = candidate;
                break;
            }
        }
        if (edgeId >= mesh.numEdges())
            throw std::runtime_error("triangle GSS avalanche source could not map a cell edge");
        const Edge& edge = mesh.getEdge(edgeId);
        if (edge.length <= 1.0e-30)
            continue;

        TriangleGssAvalancheSourceRecord record;
        record.cellId = cellId;
        record.localEdge = localEdge;
        record.edgeId = edgeId;
        record.node0 = node0;
        record.node1 = node1;
        record.edgeLength = edge.length;
        record.truncatedPartialVolume =
            geniusTri3TruncatedPartialVolumeWithEdge(
                mesh, cell, node0, node1) * config.sourceGeometryScale;
        record.electronCellQfField = electronCellField;
        record.holeCellQfField = holeCellField;
        record.electronImpactField = electronImpactField;
        record.holeImpactField = holeImpactField;
        record.electronEdgeQfField =
            std::abs(phin(static_cast<int>(node1)) - phin(static_cast<int>(node0))) /
            edge.length * fieldFactor;
        record.holeEdgeQfField =
            std::abs(phip(static_cast<int>(node1)) - phip(static_cast<int>(node0))) /
            edge.length * fieldFactor;
        record.electronMidpointDensity = gssElectronAvalancheMidpointDensity(
            n(static_cast<int>(node0)), n(static_cast<int>(node1)),
            psi(static_cast<int>(node0)), psi(static_cast<int>(node1)), Vt);
        record.holeMidpointDensity = gssHoleAvalancheMidpointDensity(
            p(static_cast<int>(node0)), p(static_cast<int>(node1)),
            psi(static_cast<int>(node0)), psi(static_cast<int>(node1)), Vt);
        record.electronMobilityDrivingField = record.electronEdgeQfField;
        record.holeMobilityDrivingField = record.holeEdgeQfField;
        record.electronLowFieldMobility = triangleGssEndpointAveragedMobility(
            mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId, node0, node1,
            CarrierType::Electron, 0.0, &psi, &cellEdgeIds, fieldFactor);
        record.holeLowFieldMobility = triangleGssEndpointAveragedMobility(
            mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId, node0, node1,
            CarrierType::Hole, 0.0, &psi, &cellEdgeIds, fieldFactor);
        record.electronMobility = triangleGssEndpointAveragedMobility(
            mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId, node0, node1,
            CarrierType::Electron, record.electronEdgeQfField,
            &psi, &cellEdgeIds, fieldFactor);
        record.holeMobility = triangleGssEndpointAveragedMobility(
            mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId, node0, node1,
            CarrierType::Hole, record.holeEdgeQfField,
            &psi, &cellEdgeIds, fieldFactor);
        record.electronAlpha = electronAlpha;
        record.holeAlpha = holeAlpha;
        record.electronFluxProxy = record.electronMobility *
            record.electronMidpointDensity * record.electronEdgeQfField;
        record.holeFluxProxy = record.holeMobility *
            record.holeMidpointDensity * record.holeEdgeQfField;
        record.electronSourceIntegral = electronAlpha *
            record.electronFluxProxy * record.truncatedPartialVolume;
        record.holeSourceIntegral = holeAlpha *
            record.holeFluxProxy * record.truncatedPartialVolume;
        record.edgeSourceIntegral =
            record.electronSourceIntegral + record.holeSourceIntegral;
        record.node0SourceIntegral = 0.5 * record.edgeSourceIntegral;
        record.node1SourceIntegral = 0.5 * record.edgeSourceIntegral;
        records.push_back(record);
    }
    return records;
}

inline std::vector<TriangleGssAvalancheSourceRecord>
triangleGssAvalancheSourceRecords(
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
    const std::vector<Real>&,
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    std::vector<TriangleGssAvalancheSourceRecord> records;
    const auto cellEdges = buildCellEdgeMap(edgeCells, mesh);
    records.reserve(static_cast<std::size_t>(mesh.numCells()) * 3);

    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        auto cellRecords = triangleGssAvalancheSourceRecordsForCell(
            config, impact, mobilityConfig, mobility,
            cellEdges.at(static_cast<std::size_t>(cellId)),
            mesh, doping, cellMaterials, cellId, psi, phin, phip, n, p, Vt,
            fieldFactor);
        records.insert(records.end(), cellRecords.begin(), cellRecords.end());
    }
    return records;
}

#include "vela/equation/ElementEdgeGssLauxAD.inl"

struct ElementEdgeGssLauxAvalancheSourceRecord {
    Index cellId = 0;
    std::array<Index, 3> edgeIds{};
    std::array<Real, 3> edgeLengths{};
    std::array<Real, 3> edgePartialVolumes{};
    std::array<Real, 3> vertexMeasures{};
    std::array<Real, 3> electronMobilities{};
    std::array<Real, 3> holeMobilities{};
    std::array<Real, 3> electronLowFieldMobilities{};
    std::array<Real, 3> holeLowFieldMobilities{};
    std::array<Real, 3> electronMobilityDrivingFields{};
    std::array<Real, 3> holeMobilityDrivingFields{};
    std::array<Real, 3> electronSignedEdgeFlux{};
    std::array<Real, 3> holeSignedEdgeFlux{};
    Point2 electronCurrentVector = Point2::Zero();
    Point2 holeCurrentVector = Point2::Zero();
    Real electronImpactField = 0.0;
    Real holeImpactField = 0.0;
    Real electronAlpha = 0.0;
    Real holeAlpha = 0.0;
    std::array<Real, 3> electronSourceIntegrals{};
    std::array<Real, 3> holeSourceIntegrals{};
    std::array<Real, 3> combinedSourceIntegrals{};
};

inline ElementEdgeGssLauxAvalancheSourceRecord
elementEdgeGssLauxAvalancheSourceRecordForCell(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<Index>&          cellEdgeIds,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    Index                              cellId,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&           ni,
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    if (!usesElementEdgeGssLauxAvalancheSource(config)) {
        throw std::invalid_argument(
            "element-edge GSS/Laux record requires its canonical configuration");
    }
    const Cell& cell = mesh.getCell(cellId);
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3) {
        throw std::invalid_argument(
            "element-edge GSS/Laux avalanche source requires Tri3 cells");
    }

    bool electricGradientValid = false;
    Real electricDoubleArea = 0.0;
    const Point2 electricGradient =
        cellScalarGradient(mesh, cell, [&](Index node) {
            return psi(static_cast<int>(node));
        }, electricGradientValid, electricDoubleArea);
    bool electronGradientValid = false;
    bool holeGradientValid = false;
    Real electronDoubleArea = 0.0;
    Real holeDoubleArea = 0.0;
    const Point2 electronGradient =
        cellScalarGradient(mesh, cell, [&](Index node) {
            return phin(static_cast<int>(node));
        }, electronGradientValid, electronDoubleArea);
    const Point2 holeGradient =
        cellScalarGradient(mesh, cell, [&](Index node) {
            return phip(static_cast<int>(node));
        }, holeGradientValid, holeDoubleArea);
    if (!electricGradientValid || !electronGradientValid ||
        !holeGradientValid) {
        throw std::invalid_argument(
            "degenerate triangle cannot evaluate avalanche driving fields");
    }

    ElementEdgeGssLauxAvalancheSourceRecord record;
    record.cellId = cellId;
    record.edgePartialVolumes =
        tri3ElementEdgeBoxPartialVolumes(mesh, cell);
    record.vertexMeasures =
        tri3ElementVertexBoxMeasures(mesh, cell);
    const bool useContactFallback =
        cellUsesContactElectricFieldFallback(config, mesh, cell);
    const Real electronQfField = electronGradient.norm() * fieldFactor;
    const Real holeQfField = holeGradient.norm() * fieldFactor;
    if (config.drivingForce == "electric_field") {
        record.electronImpactField = electricGradient.norm() * fieldFactor;
        record.holeImpactField = record.electronImpactField;
    } else if (useContactFallback) {
        const auto potential = [&](Index node) {
            return psi(static_cast<int>(node));
        };
        record.electronImpactField = contactElectricFallbackImpactField(
            config, mesh, cell, electricGradient, electronQfField,
            potential, fieldFactor);
        record.holeImpactField = contactElectricFallbackImpactField(
            config, mesh, cell, electricGradient, holeQfField,
            potential, fieldFactor);
    } else {
        record.electronImpactField = electronQfField;
        record.holeImpactField = holeQfField;
    }
    record.electronAlpha =
        impact.electronCoefficient(record.electronImpactField);
    record.holeAlpha =
        impact.holeCoefficient(record.holeImpactField);

    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const std::size_t local = static_cast<std::size_t>(localEdge);
        const Index node0 = cell.node_ids[local];
        const Index node1 =
            cell.node_ids[static_cast<std::size_t>((localEdge + 1) % 3)];
        const Index edgeId =
            edgeIdForNodePair(mesh, cellEdgeIds, node0, node1);
        if (edgeId >= mesh.numEdges()) {
            throw std::runtime_error(
                "element-edge GSS/Laux source could not map a cell edge");
        }
        const Real edgeLength =
            (meshPoint(mesh, node1) - meshPoint(mesh, node0)).norm();
        if (edgeLength <= 1.0e-30) {
            throw std::invalid_argument(
                "degenerate triangle edge cannot evaluate an SG current");
        }
        const Real electricEdgeField =
            std::abs(
                psi(static_cast<int>(node1)) -
                psi(static_cast<int>(node0))) /
            edgeLength * fieldFactor;
        const bool qfMobility =
            mobilityConfig.highFieldDrivingForce == "quasi_fermi_gradient";
        record.edgeIds[local] = edgeId;
        record.edgeLengths[local] = edgeLength;

        const Real electronEdgeField =
            std::abs(
                phin(static_cast<int>(node1)) -
                phin(static_cast<int>(node0))) /
            edgeLength * fieldFactor;
        const Real holeEdgeField =
            std::abs(
                phip(static_cast<int>(node1)) -
                phip(static_cast<int>(node0))) /
            edgeLength * fieldFactor;
        const Real electronMobilityDrive =
            qfMobility
            ? (mobilityConfig.highFieldGradientDiscretization == "transport_cell_vector"
                ? electronGradient.norm() * fieldFactor
                : electronEdgeField)
            : electricEdgeField;
        const Real holeMobilityDrive =
            qfMobility
            ? (mobilityConfig.highFieldGradientDiscretization == "transport_cell_vector"
                ? holeGradient.norm() * fieldFactor
                : holeEdgeField)
            : electricEdgeField;
        const Real electronLowFieldMobility =
            triangleGssEndpointAveragedMobility(
                mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId,
                node0, node1, CarrierType::Electron, 0.0,
                &psi, &cellEdgeIds, fieldFactor);
        const Real holeLowFieldMobility =
            triangleGssEndpointAveragedMobility(
                mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId,
                node0, node1, CarrierType::Hole, 0.0,
                &psi, &cellEdgeIds, fieldFactor);
        const Real electronMobility =
            triangleGssEndpointAveragedMobility(
                mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId,
                node0, node1, CarrierType::Electron,
                electronMobilityDrive, &psi, &cellEdgeIds, fieldFactor);
        const Real holeMobility =
            triangleGssEndpointAveragedMobility(
                mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId,
                node0, node1, CarrierType::Hole,
                holeMobilityDrive, &psi, &cellEdgeIds, fieldFactor);
        record.electronLowFieldMobilities[local] = electronLowFieldMobility;
        record.holeLowFieldMobilities[local] = holeLowFieldMobility;
        record.electronMobilityDrivingFields[local] = electronMobilityDrive;
        record.holeMobilityDrivingFields[local] = holeMobilityDrive;
        record.electronMobilities[local] = electronMobility;
        record.holeMobilities[local] = holeMobility;
        record.electronSignedEdgeFlux[local] =
            electronMobility > 0.0
            ? sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni.at(node0), ni.at(node1),
                psi(static_cast<int>(node0)), psi(static_cast<int>(node1)),
                phin(static_cast<int>(node0)), phin(static_cast<int>(node1)),
                Vt, electronMobility * Vt * fieldFactor / edgeLength)
            : 0.0;
        record.holeSignedEdgeFlux[local] =
            holeMobility > 0.0
            ? sgHoleContinuityFluxFromQuasiFermiVariableNi(
                ni.at(node0), ni.at(node1),
                psi(static_cast<int>(node0)), psi(static_cast<int>(node1)),
                phip(static_cast<int>(node0)), phip(static_cast<int>(node1)),
                Vt, holeMobility * Vt * fieldFactor / edgeLength)
            : 0.0;
    }

    record.electronCurrentVector = gssLauxTri3CurrentVector(
        mesh, cell, record.electronSignedEdgeFlux);
    record.holeCurrentVector = gssLauxTri3CurrentVector(
        mesh, cell, record.holeSignedEdgeFlux);
    std::array<Real, 3> localPsi{};
    std::array<Real, 3> localPhin{};
    std::array<Real, 3> localPhip{};
    std::array<Real, 3> localElectronDensity{};
    std::array<Real, 3> localHoleDensity{};
    std::array<Real, 3> localIntrinsicDensity{};
    for (std::size_t localNode = 0; localNode < 3; ++localNode) {
        const Index node = cell.node_ids[localNode];
        localPsi[localNode] = psi(static_cast<int>(node));
        localPhin[localNode] = phin(static_cast<int>(node));
        localPhip[localNode] = phip(static_cast<int>(node));
        localElectronDensity[localNode] = n(static_cast<int>(node));
        localHoleDensity[localNode] = p(static_cast<int>(node));
        localIntrinsicDensity[localNode] = ni.at(node);
    }
    if (isSurfaceMobilityModel(mobilityConfig)) {
        const Real electronGeneration = record.electronAlpha *
            record.electronCurrentVector.norm();
        const Real holeGeneration = record.holeAlpha *
            record.holeCurrentVector.norm();
        const auto vertexMeasures = tri3ElementVertexBoxMeasures(mesh, cell);
        for (std::size_t localNode = 0; localNode < 3; ++localNode) {
            const Real measure = vertexMeasures[localNode] *
                config.sourceGeometryScale;
            record.electronSourceIntegrals[localNode] =
                electronGeneration * measure;
            record.holeSourceIntegrals[localNode] =
                holeGeneration * measure;
            record.combinedSourceIntegrals[localNode] =
                record.electronSourceIntegrals[localNode] +
                record.holeSourceIntegrals[localNode];
        }
    } else {
        const auto sharedSourceIntegrals =
            elementEdgeGssLauxAvalancheSourceIntegralsLocal<Real>(
                config, mobilityConfig, mobility, cellEdgeIds, mesh, doping,
                cellMaterials, cellId, localPsi, localPhin, localPhip,
                localElectronDensity, localHoleDensity, localIntrinsicDensity,
                Vt, fieldFactor);
        for (std::size_t localNode = 0; localNode < 3; ++localNode) {
            record.electronSourceIntegrals[localNode] =
                sharedSourceIntegrals.electron[localNode];
            record.holeSourceIntegrals[localNode] =
                sharedSourceIntegrals.hole[localNode];
            record.combinedSourceIntegrals[localNode] =
                sharedSourceIntegrals.combined[localNode];
        }
    }
    return record;
}

inline std::vector<ElementEdgeGssLauxAvalancheSourceRecord>
elementEdgeGssLauxAvalancheSourceRecords(
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
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    const auto cellEdges = buildCellEdgeMap(edgeCells, mesh);
    std::vector<ElementEdgeGssLauxAvalancheSourceRecord> records;
    records.reserve(static_cast<std::size_t>(mesh.numCells()));
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        records.push_back(
            elementEdgeGssLauxAvalancheSourceRecordForCell(
                config, impact, mobilityConfig, mobility,
                cellEdges.at(static_cast<std::size_t>(cellId)),
                mesh, doping, cellMaterials, cellId, psi, phin, phip, n, p,
                ni, Vt, fieldFactor));
    }
    return records;
}

struct SgEdgeCurrentAvalancheSourceRecord {
    Index edgeId = 0;
    Index node0 = 0;
    Index node1 = 0;
    Real edgeLength = 0.0;
    Real edgeCouple = 0.0;
    Real edgeAreaProxy = 0.0;
    Real electricField = 0.0;
    Point2 electricFieldVector = Point2::Zero();
    Point2 electronCurrentVector = Point2::Zero();
    Point2 holeCurrentVector = Point2::Zero();
    Real electronImpactField = 0.0;
    Real holeImpactField = 0.0;
    Real electronAlpha = 0.0;
    Real holeAlpha = 0.0;
    Real electronMobility = 0.0;
    Real holeMobility = 0.0;
    Real electronLowFieldMobility = 0.0;
    Real holeLowFieldMobility = 0.0;
    Real electronMobilityDrivingField = 0.0;
    Real holeMobilityDrivingField = 0.0;
    Real electronRawFluxProxy = 0.0;
    Real holeRawFluxProxy = 0.0;
    Real electronRawSignedFluxProxy = 0.0;
    Real holeRawSignedFluxProxy = 0.0;
    Real electronReconstructedFluxProxy = 0.0;
    Real holeReconstructedFluxProxy = 0.0;
    Real electronFluxProxy = 0.0;
    Real holeFluxProxy = 0.0;
    Real electronFinalOverRawFluxProxy = 0.0;
    Real holeFinalOverRawFluxProxy = 0.0;
    Real electronSourceIntegral = 0.0;
    Real holeSourceIntegral = 0.0;
    Real edgeSourceIntegral = 0.0;
    Real electronNode0SourceIntegral = 0.0;
    Real electronNode1SourceIntegral = 0.0;
    Real holeNode0SourceIntegral = 0.0;
    Real holeNode1SourceIntegral = 0.0;
    Point2 node0ElectricFieldVector = Point2::Zero();
    Point2 node1ElectricFieldVector = Point2::Zero();
    Point2 electronNode0CurrentVector = Point2::Zero();
    Point2 electronNode1CurrentVector = Point2::Zero();
    Point2 holeNode0CurrentVector = Point2::Zero();
    Point2 holeNode1CurrentVector = Point2::Zero();
    Real electronNode0ImpactField = 0.0;
    Real electronNode1ImpactField = 0.0;
    Real holeNode0ImpactField = 0.0;
    Real holeNode1ImpactField = 0.0;
    Real electronNode0Alpha = 0.0;
    Real electronNode1Alpha = 0.0;
    Real holeNode0Alpha = 0.0;
    Real holeNode1Alpha = 0.0;
    Real node0SourceIntegral = 0.0;
    Real node1SourceIntegral = 0.0;
    bool electronSgUsesFermiDirac = false;
    Real electronSgGeneralizedEinsteinFactor = 1.0;
    Real electronSgGeneralizedBernoulliArgument = 0.0;
    bool holeSgUsesFermiDirac = false;
    Real holeSgGeneralizedEinsteinFactor = 1.0;
    Real holeSgGeneralizedBernoulliArgument = 0.0;
    SgElectronVariableNiFluxDecomposition electronSgFluxDecomposition;
    Real electronSgProductionSignedFluxNative = 0.0;
    Real electronSgReconstructionRelativeError = 0.0;
    Real electronSgProductionVsReferenceRelativeError = 0.0;
    Real electronSgStableVsReferenceRelativeError = 0.0;
    bool electronSgDiagnosticsCollected = false;
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
    Real                               Vt,
    Real                               fieldFactor = 1.0,
    bool                               collectElectronSgDiagnostics = false,
    const CarrierStatisticsConfig&     carrierStatistics = {},
    const std::vector<Real>&           Nc = {},
    const std::vector<Real>&           Nv = {})
{
    std::vector<SgEdgeCurrentAvalancheSourceRecord> records;
    records.reserve(mesh.numEdges());
    const bool qfImpact = usesQuasiFermiAvalancheDrivingForce(config);
    const bool currentAlignedImpact =
        !config.debugRawVanOverstraeten && usesCurrentAlignedAvalancheDrivingForce(config);
    const bool sentaurusEparallelImpact =
        !config.debugRawVanOverstraeten &&
        usesSentaurusEparallelAvalancheDrivingForce(config);
    const bool cellCurrentReconstructedCurrent = usesCellCurrentReconstructedAvalancheCurrent(config);
    const bool cellVectorCurrentReconstructedCurrent = usesCellVectorCurrentReconstructedAvalancheCurrent(config);
    const bool nodalVectorCurrentReconstructedCurrent =
        usesNodalVectorCurrentReconstructedAvalancheCurrent(config);
    const bool dualFaceVectorCurrentMagnitude = config.currentMagnitudeMode == "dual_face_vector_mag";
    const bool usesReconstructedSgCurrent =
        cellCurrentReconstructedCurrent ||
        cellVectorCurrentReconstructedCurrent ||
        nodalVectorCurrentReconstructedCurrent ||
        dualFaceVectorCurrentMagnitude;
    const bool needsFullEdgeFlux =
        usesReconstructedSgCurrent || sentaurusEparallelImpact;
    const bool directionalEdgePartition = usesDirectionalEdgeAvalancheSourcePartition(config);
    const bool qfMobility = mobilityConfig.highFieldDrivingForce == "quasi_fermi_gradient";
    const bool vectorQfMobility = qfMobility &&
        mobilityConfig.highFieldGradientDiscretization == "transport_cell_vector";
    const std::vector<Real> electronVectorMobilityFields = vectorQfMobility
        ? transportCellVectorEdgeGradientMagnitudes(
              mesh, edgeCells, cellMaterials, phin, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> holeVectorMobilityFields = vectorQfMobility
        ? transportCellVectorEdgeGradientMagnitudes(
              mesh, edgeCells, cellMaterials, phip, fieldFactor)
        : std::vector<Real>{};
    const bool fermiDirac = usesFermiDirac(carrierStatistics);
    if (fermiDirac && (Nc.size() != mesh.numNodes() || Nv.size() != mesh.numNodes())) {
        throw std::invalid_argument(
            "sgEdgeCurrentAvalancheSourceRecords: Fermi-Dirac statistics require "
            "per-node Nc and Nv vectors.");
    }
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
    const CellScalarGradientCache electricGradientCache =
        computeCellScalarGradientCache(mesh, [&](Index node) {
            return psi(static_cast<int>(node));
        });
    const bool nodalVertexStarEparallel = sentaurusEparallelImpact &&
        config.eparallelFieldRecovery == "nodal_vertex_star";
    const NodalScalarGradientCache nodalElectricGradientCache =
        nodalVertexStarEparallel
        ? computeTransportNodalScalarGradientCache(
              mesh, buildNodeCellMap(mesh), cellMaterials, electricGradientCache)
        : NodalScalarGradientCache{};
    const std::vector<std::vector<Index>> cellEdges =
        buildCellEdgeMap(edgeCells, mesh);
    const std::vector<std::vector<Index>> nodeEdges = buildNodeEdgeMap(mesh);

    struct EdgeSgFluxEvaluation {
        Real continuityFlux = 0.0;
        Real generalizedEinsteinFactor = 1.0;
        Real bernoulliArgument = 0.0;
    };
    const auto electronEdgeSgFlux = [&](const Edge& edge, Real mobilityValue) {
        EdgeSgFluxEvaluation evaluation;
        if (!(mobilityValue > 0.0) || !(edge.length > 1.0e-30))
            return evaluation;
        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real coefficient = mobilityValue * Vt * fieldFactor / edge.length;
        if (!fermiDirac) {
            evaluation.bernoulliArgument = (psi(j) - psi(i)) / Vt
                + std::log(ni[edge.n0] / ni[edge.n1]);
            evaluation.continuityFlux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni[edge.n0], ni[edge.n1], psi(i), psi(j), phin(i), phin(j),
                Vt, coefficient);
            return evaluation;
        }
        if (!(ni[edge.n0] > 0.0) || !(ni[edge.n1] > 0.0) ||
            !(Nc[edge.n0] > 0.0) || !(Nc[edge.n1] > 0.0))
            return evaluation;
        const Real eta0 = (psi(i) - phin(i)) / Vt
            + std::log(ni[edge.n0] / Nc[edge.n0]);
        const Real eta1 = (psi(j) - phin(j)) / Vt
            + std::log(ni[edge.n1] / Nc[edge.n1]);
        const Real driftPotential = psi(j) - psi(i) + Vt * std::log(
            (ni[edge.n1] / Nc[edge.n1]) / (ni[edge.n0] / Nc[edge.n0]));
        evaluation.generalizedEinsteinFactor =
            sgGeneralizedEinsteinFactor(n(i), n(j), eta0, eta1);
        evaluation.bernoulliArgument = driftPotential /
            (Vt * evaluation.generalizedEinsteinFactor);
        evaluation.continuityFlux = sgElectronFermiDiracContinuityFlux(
            n(i), n(j), eta0, eta1, driftPotential, phin(i), phin(j),
            Vt, coefficient);
        return evaluation;
    };
    const auto holeEdgeSgFlux = [&](const Edge& edge, Real mobilityValue) {
        EdgeSgFluxEvaluation evaluation;
        if (!(mobilityValue > 0.0) || !(edge.length > 1.0e-30))
            return evaluation;
        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real coefficient = mobilityValue * Vt * fieldFactor / edge.length;
        if (!fermiDirac) {
            evaluation.bernoulliArgument = (psi(j) - psi(i)) / Vt
                + std::log(ni[edge.n0] / ni[edge.n1]);
            evaluation.continuityFlux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                ni[edge.n0], ni[edge.n1], psi(i), psi(j), phip(i), phip(j),
                Vt, coefficient);
            return evaluation;
        }
        if (!(ni[edge.n0] > 0.0) || !(ni[edge.n1] > 0.0) ||
            !(Nv[edge.n0] > 0.0) || !(Nv[edge.n1] > 0.0))
            return evaluation;
        const Real eta0 = (phip(i) - psi(i)) / Vt
            + std::log(ni[edge.n0] / Nv[edge.n0]);
        const Real eta1 = (phip(j) - psi(j)) / Vt
            + std::log(ni[edge.n1] / Nv[edge.n1]);
        const Real driftPotential = psi(j) - psi(i) + Vt * std::log(
            (ni[edge.n0] / Nv[edge.n0]) / (ni[edge.n1] / Nv[edge.n1]));
        evaluation.generalizedEinsteinFactor =
            sgGeneralizedEinsteinFactor(p(i), p(j), eta0, eta1);
        evaluation.bernoulliArgument = driftPotential /
            (Vt * evaluation.generalizedEinsteinFactor);
        evaluation.continuityFlux = sgHoleFermiDiracContinuityFlux(
            p(i), p(j), eta0, eta1, driftPotential, phip(i), phip(j),
            Vt, coefficient);
        return evaluation;
    };

    std::vector<Real> rawElectronFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> rawHoleFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> rawSignedElectronFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> rawSignedHoleFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> reconstructedElectronFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> reconstructedHoleFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<bool> activeElectronEdge(static_cast<std::size_t>(mesh.numEdges()), false);
    std::vector<bool> activeHoleEdge(static_cast<std::size_t>(mesh.numEdges()), false);
    std::vector<Point2> nodalElectronCurrent;
    std::vector<Point2> nodalHoleCurrent;
    if (needsFullEdgeFlux) {
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
            const Real electricField = std::abs((psi_j - psi_i) / h) * fieldFactor;
            const Real electronQfField = std::abs((electronQf_j - electronQf_i) / h) * fieldFactor;
            const Real holeQfField = std::abs((holeQf_j - holeQf_i) / h) * fieldFactor;
            const Real electronMobilityField = vectorQfMobility
                ? electronVectorMobilityFields[e]
                : (qfMobility ? electronQfField : electricField);
            const Real holeMobilityField = vectorQfMobility
                ? holeVectorMobilityFields[e]
                : (qfMobility ? holeQfField : electricField);
            const Real mun = edgeMobility(
                edgeCells, mesh, doping, mobility, cellMaterials, e, CarrierType::Electron,
                electronMobilityField, &mobilityConfig, &psi);
            if (mun > 0.0) {
                const Real signedFlux = electronEdgeSgFlux(edge, mun).continuityFlux;
                rawSignedElectronFlux[static_cast<std::size_t>(e)] = signedFlux;
                rawElectronFlux[static_cast<std::size_t>(e)] = std::abs(signedFlux);
                activeElectronEdge[static_cast<std::size_t>(e)] = true;
            }
            const Real mup = edgeMobility(
                edgeCells, mesh, doping, mobility, cellMaterials, e, CarrierType::Hole,
                holeMobilityField, &mobilityConfig, &psi);
            if (mup > 0.0) {
                const Real signedFlux = holeEdgeSgFlux(edge, mup).continuityFlux;
                rawSignedHoleFlux[static_cast<std::size_t>(e)] = signedFlux;
                rawHoleFlux[static_cast<std::size_t>(e)] = std::abs(signedFlux);
                activeHoleEdge[static_cast<std::size_t>(e)] = true;
            }
        }

        if (nodalVectorCurrentReconstructedCurrent) {
            nodalElectronCurrent.resize(
                static_cast<std::size_t>(mesh.numNodes()), Point2::Zero());
            nodalHoleCurrent.resize(
                static_cast<std::size_t>(mesh.numNodes()), Point2::Zero());
            const auto electronFlux = [&](Index edgeId) {
                return rawSignedElectronFlux[static_cast<std::size_t>(edgeId)];
            };
            const auto holeFlux = [&](Index edgeId) {
                return rawSignedHoleFlux[static_cast<std::size_t>(edgeId)];
            };
            const auto electronActive = [&](Index edgeId) {
                return activeElectronEdge[static_cast<std::size_t>(edgeId)];
            };
            const auto holeActive = [&](Index edgeId) {
                return activeHoleEdge[static_cast<std::size_t>(edgeId)];
            };
            for (Index node = 0; node < mesh.numNodes(); ++node) {
                nodalElectronCurrent[static_cast<std::size_t>(node)] =
                    nodalLeastSquaresCurrentVector(
                        node, nodeEdges, mesh, electronFlux, electronActive);
                nodalHoleCurrent[static_cast<std::size_t>(node)] =
                    nodalLeastSquaresCurrentVector(
                        node, nodeEdges, mesh, holeFlux, holeActive);
            }
        }

        if (usesReconstructedSgCurrent) {
            for (Index e = 0; e < mesh.numEdges(); ++e) {
                reconstructedElectronFlux[static_cast<std::size_t>(e)] = dualFaceVectorCurrentMagnitude
                    ? medianDualFaceVectorReconstructedEdgeFluxMagnitude(
                        e, rawSignedElectronFlux, edgeCells, cellEdges, mesh)
                    : (nodalVectorCurrentReconstructedCurrent
                        ? 0.5 * (
                            nodalElectronCurrent[mesh.getEdge(e).n0].norm() +
                            nodalElectronCurrent[mesh.getEdge(e).n1].norm())
                        : (cellVectorCurrentReconstructedCurrent
                            ? cellVectorReconstructedEdgeFluxMagnitude(
                                e, rawSignedElectronFlux, edgeCells, cellEdges, mesh)
                            : cellSmoothedEdgeFluxMagnitude(
                                e, rawElectronFlux, edgeCells, cellEdges)));
                reconstructedHoleFlux[static_cast<std::size_t>(e)] = dualFaceVectorCurrentMagnitude
                    ? medianDualFaceVectorReconstructedEdgeFluxMagnitude(
                        e, rawSignedHoleFlux, edgeCells, cellEdges, mesh)
                    : (nodalVectorCurrentReconstructedCurrent
                        ? 0.5 * (
                            nodalHoleCurrent[mesh.getEdge(e).n0].norm() +
                            nodalHoleCurrent[mesh.getEdge(e).n1].norm())
                        : (cellVectorCurrentReconstructedCurrent
                            ? cellVectorReconstructedEdgeFluxMagnitude(
                                e, rawSignedHoleFlux, edgeCells, cellEdges, mesh)
                            : cellSmoothedEdgeFluxMagnitude(
                                e, rawHoleFlux, edgeCells, cellEdges)));
            }
        }
    }

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
        const Real electricField = std::abs((psi_j - psi_i) / h) * fieldFactor;
        const Real electronQfField = std::abs((electronQf_j - electronQf_i) / h) * fieldFactor;
        const Real holeQfField = std::abs((holeQf_j - holeQf_i) / h) * fieldFactor;
        const Real electronCoefficientField = qfImpact
            ? edgeQuasiFermiCoefficientField(
                config, electronQfField, electricField, edgeCells, mesh, e,
                contactNodes, electronQfGradientCache, fieldFactor)
            : electricField;
        const Real holeCoefficientField = qfImpact
            ? edgeQuasiFermiCoefficientField(
                config, holeQfField, electricField, edgeCells, mesh, e,
                contactNodes, holeQfGradientCache, fieldFactor)
            : electricField;
        const Real electronMobilityField = vectorQfMobility
            ? electronVectorMobilityFields[e]
            : (qfMobility ? electronQfField : electricField);
        const Real holeMobilityField = vectorQfMobility
            ? holeVectorMobilityFields[e]
            : (qfMobility ? holeQfField : electricField);

        const Real nAvg = 0.5 * (n(i) + n(j));
        const Real pAvg = 0.5 * (p(i) + p(j));
        const Real nMid = cellReconstructedAvalancheMidpointDensity(
            config, n(i), n(j), psi_i, psi_j, Vt);
        const Real pMid = cellReconstructedAvalancheMidpointDensity(
            config, p(i), p(j), psi_j, psi_i, Vt);
        const Real signedElectricField01 = -(psi_j - psi_i) / h * fieldFactor;

        const Real edgeArea = avalancheSourceEdgeArea(
            config, edgeCells, mesh, e, &cellMaterials);
        SgEdgeCurrentAvalancheSourceRecord record;
        record.edgeId = e;
        record.node0 = edge.n0;
        record.node1 = edge.n1;
        record.edgeLength = h;
        record.edgeCouple = edge.couple;
        record.edgeAreaProxy = edgeArea;
        record.electricField = electricField;
        if (sentaurusEparallelImpact) {
            bool validElectricGradient = false;
            const Point2 electricGradient = eparallelElectricGradientForEdge(
                config,
                edgeCells,
                mesh,
                e,
                electricGradientCache,
                nodalElectricGradientCache,
                validElectricGradient);
            if (validElectricGradient) {
                record.electricFieldVector = -fieldFactor * electricGradient;
            } else {
                const Node& node0 = mesh.getNode(edge.n0);
                const Node& node1 = mesh.getNode(edge.n1);
                record.electricFieldVector = signedElectricField01 *
                    Point2{(node1.x - node0.x) / h, (node1.y - node0.y) / h};
            }
            // SG electron flux is a particle flux; conventional Jn points in
            // the opposite direction. Hole particle flux and Jp are aligned.
            if (nodalVectorCurrentReconstructedCurrent) {
                record.electronCurrentVector = -0.5 * (
                    nodalElectronCurrent[edge.n0] +
                    nodalElectronCurrent[edge.n1]);
                record.holeCurrentVector = 0.5 * (
                    nodalHoleCurrent[edge.n0] +
                    nodalHoleCurrent[edge.n1]);
            } else {
                record.electronCurrentVector = -edgeAveragedCellVectorCurrent(
                    e, rawSignedElectronFlux, edgeCells, cellEdges, mesh);
                record.holeCurrentVector = edgeAveragedCellVectorCurrent(
                    e, rawSignedHoleFlux, edgeCells, cellEdges, mesh);
            }
        }
        record.electronMobilityDrivingField = electronMobilityField;
        record.holeMobilityDrivingField = holeMobilityField;

        const Real mun = edgeMobility(
            edgeCells, mesh, doping, mobility, cellMaterials, e, CarrierType::Electron,
            electronMobilityField, &mobilityConfig, &psi);
        record.electronMobility = mun;
        record.electronLowFieldMobility = edgeMobility(
            edgeCells, mesh, doping, mobility, cellMaterials, e,
            CarrierType::Electron, 0.0, &mobilityConfig, &psi);
        const Real mup = edgeMobility(
            edgeCells, mesh, doping, mobility, cellMaterials, e, CarrierType::Hole,
            holeMobilityField, &mobilityConfig, &psi);
        record.holeMobility = mup;
        record.holeLowFieldMobility = edgeMobility(
            edgeCells, mesh, doping, mobility, cellMaterials, e,
            CarrierType::Hole, 0.0, &mobilityConfig, &psi);
        constexpr bool IncludeElectronNiGradientDrift = true;
        const Real electronContinuityCoefficient =
            mun > 0.0 ? mun * Vt * fieldFactor / h : 0.0;
        const EdgeSgFluxEvaluation electronSg = electronEdgeSgFlux(edge, mun);
        const EdgeSgFluxEvaluation holeSg = holeEdgeSgFlux(edge, mup);
        const Real electronContinuityFlux01 = electronSg.continuityFlux;
        const Real holeContinuityFlux01 = holeSg.continuityFlux;
        record.electronSgUsesFermiDirac = fermiDirac;
        record.electronSgGeneralizedEinsteinFactor =
            electronSg.generalizedEinsteinFactor;
        record.electronSgGeneralizedBernoulliArgument =
            electronSg.bernoulliArgument;
        record.holeSgUsesFermiDirac = fermiDirac;
        record.holeSgGeneralizedEinsteinFactor =
            holeSg.generalizedEinsteinFactor;
        record.holeSgGeneralizedBernoulliArgument = holeSg.bernoulliArgument;
        // SG continuity fluxes use the particle-flux convention, so the
        // physical charge-current magnitude is |F_p - F_n|.
        const Real conservedTotalFluxMagnitude =
            conservedTotalCurrentFluxMagnitude(
                electronContinuityFlux01, holeContinuityFlux01);
        if (collectElectronSgDiagnostics && !fermiDirac) {
            record.electronSgDiagnosticsCollected = true;
            record.electronSgFluxDecomposition =
                sgElectronContinuityFluxFromQuasiFermiVariableNiDecomposition(
                    ni[edge.n0],
                    ni[edge.n1],
                    psi_i,
                    psi_j,
                    phin_i,
                    phin_j,
                    Vt,
                    electronContinuityCoefficient,
                    IncludeElectronNiGradientDrift);
            record.electronSgProductionSignedFluxNative = electronContinuityFlux01;
            const Real electronSgReconstructionScale = std::max({
                std::abs(record.electronSgProductionSignedFluxNative),
                std::abs(record.electronSgFluxDecomposition.reconstructedFlux),
                std::abs(record.electronSgFluxDecomposition.coef)
                    * (std::abs(record.electronSgFluxDecomposition.leftTerm)
                       + std::abs(record.electronSgFluxDecomposition.rightTerm)),
                Real{1.0e-300}});
            record.electronSgReconstructionRelativeError =
                std::abs(record.electronSgProductionSignedFluxNative
                         - record.electronSgFluxDecomposition.reconstructedFlux)
                / electronSgReconstructionScale;

            const Real referenceScale = std::max({
                std::abs(record.electronSgProductionSignedFluxNative),
                record.electronSgFluxDecomposition.highPrecisionReferenceTermScale,
                Real{1.0e-300}});
            const long double referenceScaleLong =
                static_cast<long double>(referenceScale);
            const long double referenceFluxLong = static_cast<long double>(
                record.electronSgFluxDecomposition.highPrecisionReferenceFlux);
            record.electronSgProductionVsReferenceRelativeError =
                static_cast<Real>(std::abs(
                    static_cast<long double>(
                        record.electronSgProductionSignedFluxNative)
                    - referenceFluxLong) / referenceScaleLong);
            record.electronSgStableVsReferenceRelativeError =
                static_cast<Real>(std::abs(
                    static_cast<long double>(
                        record.electronSgFluxDecomposition.stableFactorizedFlux)
                    - referenceFluxLong) / referenceScaleLong);
        } else if (collectElectronSgDiagnostics) {
            // The legacy left/right decomposition is specific to the
            // Boltzmann variable-ni formula.  Fermi-Dirac diagnostics expose
            // the generalized Einstein factor and Bernoulli argument.  Keep
            // the endpoint and algebraic fields populated so edge-audit CSVs
            // remain useful across both statistics modes.
            record.electronSgDiagnosticsCollected = true;
            record.electronSgProductionSignedFluxNative =
                electronContinuityFlux01;
            auto& decomposition = record.electronSgFluxDecomposition;
            decomposition.ni0 = ni[edge.n0];
            decomposition.ni1 = ni[edge.n1];
            decomposition.n0 = n(i);
            decomposition.n1 = n(j);
            decomposition.psi0 = psi_i;
            decomposition.psi1 = psi_j;
            decomposition.phin0 = phin_i;
            decomposition.phin1 = phin_j;
            decomposition.eta = electronSg.bernoulliArgument;
            const SGEdgeWeights generalizedWeights =
                sgEdgeWeights(electronSg.bernoulliArgument, 1.0);
            decomposition.bernoulliMinusEta = generalizedWeights.b_minus;
            decomposition.bernoulliEta = generalizedWeights.b_plus;
            decomposition.coef = electronContinuityCoefficient
                * electronSg.generalizedEinsteinFactor;
            decomposition.leftTerm = generalizedWeights.b_minus * n(i);
            decomposition.rightTerm = generalizedWeights.b_plus * n(j);
            decomposition.signedDifference =
                decomposition.leftTerm - decomposition.rightTerm;
            decomposition.reconstructedFlux =
                decomposition.coef * decomposition.signedDifference;
            decomposition.stableFactorizedFlux = electronContinuityFlux01;
            decomposition.highPrecisionReferenceFlux = electronContinuityFlux01;
            decomposition.includeNiGradientDrift = true;
            decomposition.flatQuasiFermiShortCircuit = phin_i == phin_j;
        }

        if (mun > 0.0) {
            record.electronImpactField = sentaurusEparallelImpact
                ? sentaurusEparallelAvalancheDrivingField(
                    record.electricFieldVector, record.electronCurrentVector)
                : (currentAlignedImpact
                    ? parallelCurrentAvalancheDrivingField(
                        signedElectricField01, electronContinuityFlux01)
                    : electronAvalancheDrivingField(
                        config, electronCoefficientField, electricField, nAvg));
            record.electronRawSignedFluxProxy = electronContinuityFlux01;
            record.electronRawFluxProxy = std::abs(electronContinuityFlux01);
            record.electronReconstructedFluxProxy = usesReconstructedSgCurrent
                ? reconstructedElectronFlux[static_cast<std::size_t>(e)]
                : record.electronRawFluxProxy;
            record.electronFluxProxy = selectAvalancheCurrentFluxProxy(
                config,
                record.electronRawFluxProxy,
                record.electronReconstructedFluxProxy,
                mun,
                nMid,
                record.electronImpactField,
                electricField,
                conservedTotalFluxMagnitude);
            record.electronFinalOverRawFluxProxy = record.electronRawFluxProxy > 0.0
                ? record.electronFluxProxy / record.electronRawFluxProxy
                : 0.0;
            record.electronAlpha = impact.electronCoefficient(record.electronImpactField);
            record.electronSourceIntegral =
                record.electronAlpha * record.electronFluxProxy * edgeArea;
            record.edgeSourceIntegral += record.electronSourceIntegral;
        }

        if (mup > 0.0) {
            record.holeImpactField = sentaurusEparallelImpact
                ? sentaurusEparallelAvalancheDrivingField(
                    record.electricFieldVector, record.holeCurrentVector)
                : (currentAlignedImpact
                    ? parallelCurrentAvalancheDrivingField(
                        signedElectricField01, holeContinuityFlux01)
                    : holeAvalancheDrivingField(
                        config, holeCoefficientField, electricField, pAvg));
            record.holeRawSignedFluxProxy = holeContinuityFlux01;
            record.holeRawFluxProxy = std::abs(holeContinuityFlux01);
            record.holeReconstructedFluxProxy = usesReconstructedSgCurrent
                ? reconstructedHoleFlux[static_cast<std::size_t>(e)]
                : record.holeRawFluxProxy;
            record.holeFluxProxy = selectAvalancheCurrentFluxProxy(
                config,
                record.holeRawFluxProxy,
                record.holeReconstructedFluxProxy,
                mup,
                pMid,
                record.holeImpactField,
                electricField,
                conservedTotalFluxMagnitude);
            record.holeFinalOverRawFluxProxy = record.holeRawFluxProxy > 0.0
                ? record.holeFluxProxy / record.holeRawFluxProxy
                : 0.0;
            record.holeAlpha = impact.holeCoefficient(record.holeImpactField);
            record.holeSourceIntegral =
                record.holeAlpha * record.holeFluxProxy * edgeArea;
            record.edgeSourceIntegral += record.holeSourceIntegral;
        }

        EdgeAvalancheDirectionalWeights weights;
        if (directionalEdgePartition) {
            weights = edgeAvalancheDirectionalWeights(
                edgeCells,
                mesh,
                e,
                electronQfGradientCache,
                holeQfGradientCache);
        }
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
    if (config.sourceMappingMode == "nodal_eparallel_p1") {
        std::vector<Point2> electricGradient(mesh.numNodes(), Point2::Zero());
        for (Index node = 0; node < mesh.numNodes(); ++node) {
            if (node < nodalElectricGradientCache.valid.size() &&
                nodalElectricGradientCache.valid[node]) {
                electricGradient[node] =
                    fieldFactor * nodalElectricGradientCache.gradients[node];
            }
        }
        std::vector<Real> nodeMeasure(mesh.numNodes(), 0.0);
        for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
            if (cellId >= cellMaterials.size() ||
                !isTransportMaterial(cellMaterials[cellId])) {
                continue;
            }
            const Cell& cell = mesh.getCell(cellId);
            if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
                continue;
            const Real share = triangleArea(mesh, cell) / 3.0;
            for (const Index node : cell.node_ids)
                nodeMeasure[node] += share;
        }

        std::vector<Real> electronNodeSource(mesh.numNodes(), 0.0);
        std::vector<Real> holeNodeSource(mesh.numNodes(), 0.0);
        std::vector<Real> electronNodeDrive(mesh.numNodes(), 0.0);
        std::vector<Real> holeNodeDrive(mesh.numNodes(), 0.0);
        std::vector<Real> electronNodeAlpha(mesh.numNodes(), 0.0);
        std::vector<Real> holeNodeAlpha(mesh.numNodes(), 0.0);
        for (Index node = 0; node < mesh.numNodes(); ++node) {
            const Point2 electricField = -electricGradient[node];
            const Point2 electronCurrent = -nodalElectronCurrent[node];
            const Point2 holeCurrent = nodalHoleCurrent[node];
            electronNodeDrive[node] = sentaurusEparallelAvalancheDrivingField(
                electricField, electronCurrent);
            holeNodeDrive[node] = sentaurusEparallelAvalancheDrivingField(
                electricField, holeCurrent);
            electronNodeAlpha[node] =
                impact.electronCoefficient(electronNodeDrive[node]);
            holeNodeAlpha[node] =
                impact.holeCoefficient(holeNodeDrive[node]);
            electronNodeSource[node] = electronNodeAlpha[node] *
                electronCurrent.norm() * nodeMeasure[node];
            holeNodeSource[node] = holeNodeAlpha[node] *
                holeCurrent.norm() * nodeMeasure[node];
        }

        std::vector<Real> nodeEdgeMeasure(mesh.numNodes(), 0.0);
        for (const auto& record : records) {
            const Real halfMeasure = 0.5 * record.edgeAreaProxy;
            nodeEdgeMeasure[record.node0] += halfMeasure;
            nodeEdgeMeasure[record.node1] += halfMeasure;
        }
        for (auto& record : records) {
            const Real halfMeasure = 0.5 * record.edgeAreaProxy;
            const Real weight0 = nodeEdgeMeasure[record.node0] > 0.0
                ? halfMeasure / nodeEdgeMeasure[record.node0] : 0.0;
            const Real weight1 = nodeEdgeMeasure[record.node1] > 0.0
                ? halfMeasure / nodeEdgeMeasure[record.node1] : 0.0;
            record.electronNode0SourceIntegral =
                weight0 * electronNodeSource[record.node0];
            record.electronNode1SourceIntegral =
                weight1 * electronNodeSource[record.node1];
            record.holeNode0SourceIntegral = weight0 * holeNodeSource[record.node0];
            record.holeNode1SourceIntegral = weight1 * holeNodeSource[record.node1];
            record.node0ElectricFieldVector = -electricGradient[record.node0];
            record.node1ElectricFieldVector = -electricGradient[record.node1];
            record.electronNode0CurrentVector = -nodalElectronCurrent[record.node0];
            record.electronNode1CurrentVector = -nodalElectronCurrent[record.node1];
            record.holeNode0CurrentVector = nodalHoleCurrent[record.node0];
            record.holeNode1CurrentVector = nodalHoleCurrent[record.node1];
            record.electronNode0ImpactField = electronNodeDrive[record.node0];
            record.electronNode1ImpactField = electronNodeDrive[record.node1];
            record.holeNode0ImpactField = holeNodeDrive[record.node0];
            record.holeNode1ImpactField = holeNodeDrive[record.node1];
            record.electronNode0Alpha = electronNodeAlpha[record.node0];
            record.electronNode1Alpha = electronNodeAlpha[record.node1];
            record.holeNode0Alpha = holeNodeAlpha[record.node0];
            record.holeNode1Alpha = holeNodeAlpha[record.node1];
            record.electronSourceIntegral =
                record.electronNode0SourceIntegral + record.electronNode1SourceIntegral;
            record.holeSourceIntegral =
                record.holeNode0SourceIntegral + record.holeNode1SourceIntegral;
            record.node0SourceIntegral =
                record.electronNode0SourceIntegral + record.holeNode0SourceIntegral;
            record.node1SourceIntegral =
                record.electronNode1SourceIntegral + record.holeNode1SourceIntegral;
            record.edgeSourceIntegral =
                record.electronSourceIntegral + record.holeSourceIntegral;
        }
    }
    return records;
}

inline void addCellMappedEdgeSourceToNodes(
    std::vector<Real>&                       target,
    const std::vector<std::vector<Index>>&   edgeCells,
    const DeviceMesh&                        mesh,
    const SgEdgeCurrentAvalancheSourceRecord& record,
    Real                                     sourceIntegral)
{
    if (record.edgeId >= edgeCells.size() || sourceIntegral == 0.0)
        return;
    const auto& cells = edgeCells[record.edgeId];
    Real areaSum = 0.0;
    for (Index cellId : cells) {
        if (cellId < mesh.numCells())
            areaSum += triangleArea(mesh, mesh.getCell(cellId));
    }
    if (areaSum <= 0.0) {
        if (record.node0 < target.size()) target[record.node0] += 0.5 * sourceIntegral;
        if (record.node1 < target.size()) target[record.node1] += 0.5 * sourceIntegral;
        return;
    }
    for (Index cellId : cells) {
        if (cellId >= mesh.numCells())
            continue;
        const Cell& cell = mesh.getCell(cellId);
        if (cell.node_ids.empty())
            continue;
        const Real cellShare = sourceIntegral * triangleArea(mesh, cell) / areaSum;
        const Real nodeShare = cellShare / static_cast<Real>(cell.node_ids.size());
        for (Index nodeId : cell.node_ids) {
            if (nodeId < target.size())
                target[nodeId] += nodeShare;
        }
    }
}

inline void addMappedEdgeSourceToNodes(
    const ImpactIonizationModelConfig&       config,
    std::vector<Real>&                       target,
    const std::vector<std::vector<Index>>&   edgeCells,
    const DeviceMesh&                        mesh,
    const SgEdgeCurrentAvalancheSourceRecord& record,
    Real                                     node0DirectionalSource,
    Real                                     node1DirectionalSource,
    Real                                     sourceIntegral)
{
    if (config.sourceMappingMode == "cell_F_cell_alpha_cell_G_to_node") {
        addCellMappedEdgeSourceToNodes(target, edgeCells, mesh, record, sourceIntegral);
        return;
    }
    if (config.sourceMappingMode == "edge_F_edge_alpha_edge_G_to_node") {
        if (record.node0 < target.size()) target[record.node0] += 0.5 * sourceIntegral;
        if (record.node1 < target.size()) target[record.node1] += 0.5 * sourceIntegral;
        return;
    }
    if (record.node0 < target.size()) target[record.node0] += node0DirectionalSource;
    if (record.node1 < target.size()) target[record.node1] += node1DirectionalSource;
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
    Real                               Vt,
    Real                               fieldFactor = 1.0,
    const CarrierStatisticsConfig&     carrierStatistics = {},
    const std::vector<Real>&           Nc = {},
    const std::vector<Real>&           Nv = {})
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
        Vt,
        fieldFactor,
        false,
        carrierStatistics,
        Nc,
        Nv);
    for (const auto& record : records) {
        addMappedEdgeSourceToNodes(
            config, source.electron, edgeCells, mesh, record,
            record.electronNode0SourceIntegral,
            record.electronNode1SourceIntegral,
            record.electronSourceIntegral);
        addMappedEdgeSourceToNodes(
            config, source.hole, edgeCells, mesh, record,
            record.holeNode0SourceIntegral,
            record.holeNode1SourceIntegral,
            record.holeSourceIntegral);
        addMappedEdgeSourceToNodes(
            config, source.combined, edgeCells, mesh, record,
            record.node0SourceIntegral,
            record.node1SourceIntegral,
            record.edgeSourceIntegral);
    }
    return source;
}

inline SgAvalancheSourceComponentIntegrals
currentDensityAvalancheSourceComponentIntegrals(
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
    Real                               Vt,
    Real                               fieldFactor = 1.0,
    const CarrierStatisticsConfig&     carrierStatistics = {},
    const std::vector<Real>&           Nc = {},
    const std::vector<Real>&           Nv = {})
{
    if (usesElementEdgeGssLauxAvalancheSource(config)) {
        SgAvalancheSourceComponentIntegrals source;
        source.electron.assign(
            static_cast<std::size_t>(mesh.numNodes()), 0.0);
        source.hole.assign(
            static_cast<std::size_t>(mesh.numNodes()), 0.0);
        source.combined.assign(
            static_cast<std::size_t>(mesh.numNodes()), 0.0);
        const auto records = elementEdgeGssLauxAvalancheSourceRecords(
            config, impact, mobilityConfig, mobility, edgeCells, mesh, doping,
            cellMaterials, psi, phin, phip, n, p, ni, Vt, fieldFactor);
        for (const auto& record : records) {
            const Cell& cell = mesh.getCell(record.cellId);
            for (std::size_t localNode = 0; localNode < 3; ++localNode) {
                const Index node = cell.node_ids[localNode];
                source.electron[node] +=
                    record.electronSourceIntegrals[localNode];
                source.hole[node] +=
                    record.holeSourceIntegrals[localNode];
                source.combined[node] +=
                    record.combinedSourceIntegrals[localNode];
            }
        }
        return source;
    }

    if (!usesTriangleGssAvalancheSource(config)) {
        return sgEdgeCurrentAvalancheSourceComponentIntegrals(
            config, impact, mobilityConfig, mobility, edgeCells, mesh, doping,
            cellMaterials, psi, phin, phip, n, p, ni, Vt, fieldFactor,
            carrierStatistics, Nc, Nv);
    }

    SgAvalancheSourceComponentIntegrals source;
    source.electron.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    source.hole.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    source.combined.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    const auto records = triangleGssAvalancheSourceRecords(
        config, impact, mobilityConfig, mobility, edgeCells, mesh, doping,
        cellMaterials, psi, phin, phip, n, p, ni, Vt, fieldFactor);
    for (const auto& record : records) {
        const Real electronNodeSource = 0.5 * record.electronSourceIntegral;
        const Real holeNodeSource = 0.5 * record.holeSourceIntegral;
        source.electron[record.node0] += electronNodeSource;
        source.electron[record.node1] += electronNodeSource;
        source.hole[record.node0] += holeNodeSource;
        source.hole[record.node1] += holeNodeSource;
        source.combined[record.node0] += record.node0SourceIntegral;
        source.combined[record.node1] += record.node1SourceIntegral;
    }
    return source;
}

inline std::vector<Real> currentDensityAvalancheSourceIntegrals(
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
    Real                               Vt,
    Real                               fieldFactor = 1.0,
    const CarrierStatisticsConfig&     carrierStatistics = {},
    const std::vector<Real>&           Nc = {},
    const std::vector<Real>&           Nv = {})
{
    return currentDensityAvalancheSourceComponentIntegrals(
        config, impact, mobilityConfig, mobility, edgeCells, mesh, doping,
        cellMaterials, psi, phin, phip, n, p, ni, Vt, fieldFactor,
        carrierStatistics, Nc, Nv).combined;
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
    Real                               Vt,
    Real                               fieldFactor = 1.0,
    const CarrierStatisticsConfig&     carrierStatistics = {},
    const std::vector<Real>&           Nc = {},
    const std::vector<Real>&           Nv = {})
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
        Vt,
        fieldFactor,
        false,
        carrierStatistics,
        Nc,
        Nv);
    for (const auto& record : records) {
        addMappedEdgeSourceToNodes(
            config, source, edgeCells, mesh, record,
            record.node0SourceIntegral,
            record.node1SourceIntegral,
            record.edgeSourceIntegral);
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
        electronImpactField, &mobilityConfig);
    const Real mup = nodeMobility(
        nodeCells, mesh, doping, mobility, cellMaterials, nodeId, CarrierType::Hole,
        holeImpactField, &mobilityConfig);
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
/// At semiconductor/insulator interfaces, prefer the transport-capable material
/// so an ordering where oxide cells precede silicon cells cannot assign ni=0 to
/// an Ohmic contact or semiconductor interface node.  Preserve first-found
/// semantics between materials of the same transport class.
inline std::vector<Real> buildNodeNi(const DeviceMesh&       mesh,
                                     const MaterialDatabase& matdb,
                                     Real                    temperature_K = constants::T0)
{
    const Index N = mesh.numNodes();
    std::vector<Real> ni_v(N, 0.0);
    std::vector<bool> found(N, false);
    std::vector<bool> transportFound(N, false);
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& cell   = mesh.getCell(c);
        const auto& region = mesh.getRegion(cell.region_id);
        const Material material = matdb.getMaterial(region.material, temperature_K);
        const Real ni_mat = material.ni;
        const bool isTransportMaterial =
            material.ni > 0.0 || material.mun > 0.0 || material.mup > 0.0;
        for (Index nid : cell.node_ids) {
            if (!found[nid] || (!transportFound[nid] && isTransportMaterial)) {
                ni_v[nid]  = ni_mat;
                found[nid] = true;
                transportFound[nid] = isTransportMaterial;
            }
        }
    }
    return ni_v;
}

/// Build per-node effective conduction/valence density of states.  Shared
/// semiconductor/insulator nodes use the transport-capable material, matching
/// buildNodeNi ownership semantics.
inline std::vector<Real> buildNodeDensityOfStates(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    Real temperature_K,
    bool electrons)
{
    const Index N = mesh.numNodes();
    std::vector<Real> values(N, 0.0);
    std::vector<bool> found(N, false);
    std::vector<bool> transportFound(N, false);
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& cell = mesh.getCell(c);
        const auto& region = mesh.getRegion(cell.region_id);
        const Material material = matdb.getMaterial(region.material, temperature_K);
        const bool transport = material.ni > 0.0 || material.mun > 0.0 || material.mup > 0.0;
        const auto& property = electrons ? material.Nc_m3 : material.Nv_m3;
        const Real value = property.value_or(0.0);
        for (Index node : cell.node_ids) {
            if (!found[node] || (!transportFound[node] && transport)) {
                values[node] = value;
                found[node] = true;
                transportFound[node] = transport;
            }
        }
    }
    return values;
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
                                              const BandgapNarrowingConfig& config,
                                              Real                    thermalVoltage)
{
    const Real temperature_K = thermalVoltage * constants::q / constants::kb;
    std::vector<Real> ni_v = buildNodeNi(mesh, matdb, temperature_K);
    const auto bgn = makeBandgapNarrowingModel(config);
    std::vector<Real> Nc300;
    std::vector<Real> Nv300;
    if (config.fermiStatisticsCorrection) {
        Nc300 = buildNodeDensityOfStates(mesh, matdb, constants::T0, true);
        Nv300 = buildNodeDensityOfStates(mesh, matdb, constants::T0, false);
    }
    const Real thermalVoltage300 = constants::kb * constants::T0 / constants::q;
    for (Index i = 0; i < mesh.numNodes(); ++i) {
        Real delta = bgn->deltaEg(doping.totalImpurity(i), 0.0, 0.0);
        if (config.fermiStatisticsCorrection && ni_v[i] > 0.0) {
            delta += fermiStatisticsBandgapCorrection(
                doping.donors(i), doping.acceptors(i), Nc300[i], Nv300[i],
                thermalVoltage300);
        }
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
        bandgapNarrowingConfig,
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
