#pragma once

#include "vela/core/Types.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/mesh/DeviceMesh.h"
#include <array>
#include <utility>

namespace vela {

/** P1 cell-current support reconstructed from quasi-Fermi gradients.
 *
 * The returned continuity contributions are integrated particle line fluxes
 * in the same native units and sign convention as CoupledDDAssembler's SG
 * edge fluxes.  The current vectors omit q; ContactCurrent applies it after
 * the boundary line integration.
 */
struct ElementQfGradientCellCurrent {
    bool valid = false;
    Real area = 0.0;
    Point2 electronParticleCurrent = Point2::Zero();
    Point2 holeParticleCurrent = Point2::Zero();
    std::array<Real, 3> electronResidual{};
    std::array<Real, 3> holeResidual{};
};

template <typename ElectronQfAt, typename HoleQfAt,
          typename ElectronDensityAt, typename HoleDensityAt>
ElementQfGradientCellCurrent elementQfGradientCellCurrent(
    const DeviceMesh& mesh,
    const Cell& cell,
    ElectronQfAt&& electronQfAt,
    HoleQfAt&& holeQfAt,
    ElectronDensityAt&& electronDensityAt,
    HoleDensityAt&& holeDensityAt,
    Real electronMobility,
    Real holeMobility,
    Real fieldFromCoordinateDeltaFactor)
{
    ElementQfGradientCellCurrent result;
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
        return result;

    bool electronValid = false;
    bool holeValid = false;
    Real electronArea = 0.0;
    Real holeArea = 0.0;
    const Point2 electronGradient = detail::cellScalarGradient(
        mesh, cell, std::forward<ElectronQfAt>(electronQfAt),
        electronValid, electronArea);
    const Point2 holeGradient = detail::cellScalarGradient(
        mesh, cell, std::forward<HoleQfAt>(holeQfAt),
        holeValid, holeArea);
    if (!electronValid || !holeValid || electronArea <= 0.0)
        return result;

    Real electronDensity = 0.0;
    Real holeDensity = 0.0;
    for (const Index node : cell.node_ids) {
        electronDensity += electronDensityAt(node);
        holeDensity += holeDensityAt(node);
    }
    electronDensity /= 3.0;
    holeDensity /= 3.0;

    result.valid = true;
    result.area = electronArea;
    result.electronParticleCurrent = electronMobility * electronDensity *
        fieldFromCoordinateDeltaFactor * electronGradient;
    result.holeParticleCurrent = holeMobility * holeDensity *
        fieldFromCoordinateDeltaFactor * holeGradient;

    for (int localNode = 0; localNode < 3; ++localNode) {
        bool basisValid = false;
        Real basisArea = 0.0;
        const Point2 basisGradient = detail::cellScalarGradient(
            mesh, cell,
            [&](Index node) {
                return node == cell.node_ids[static_cast<std::size_t>(localNode)]
                    ? 1.0 : 0.0;
            },
            basisValid, basisArea);
        if (!basisValid)
            return {};
        result.electronResidual[static_cast<std::size_t>(localNode)] =
            -result.area * result.electronParticleCurrent.dot(basisGradient);
        result.holeResidual[static_cast<std::size_t>(localNode)] =
            result.area * result.holeParticleCurrent.dot(basisGradient);
    }
    return result;
}

} // namespace vela
