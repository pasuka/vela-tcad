#pragma once

#include "vela/core/Types.h"
#include "vela/core/ArithmeticPrecision.h"
#include "vela/mesh/DeviceMesh.h"
#include <vector>
#include <string>

namespace vela {

/**
 * @brief Region-based doping specification for a single material region.
 *
 * Concentrations are stored in the active internal unit system: legacy SI [m^-3]
 * or `unit_scaling` TCAD units [cm^-3].
 */
struct RegionDopingSpec {
    std::string region;    ///< Region name (must match a Region::name in the mesh)
    Real        donors;    ///< Donor concentration Nd in active internal units
    Real        acceptors; ///< Acceptor concentration Na in active internal units
};

/**
 * @brief Per-node doping concentrations for the device.
 *
 * Stores donor (Nd) and acceptor (Na) concentrations at each mesh node
 * in active internal units. Net doping is defined as Nd - Na.
 *
 * Populated either directly via setNodeDoping() or via the static
 * factory fromMeshAndRegions() which maps region-based doping specs
 * onto the mesh nodes.
 */
class DopingModel {
public:
    /// Construct an empty model sized for @p numNodes nodes (all zeroes).
    explicit DopingModel(Index numNodes);

    /// Set doping at a single node.
    void setNodeDoping(Index nodeId, Real donors, Real acceptors);
    /// Preserve both source operands until their subtraction. Units are the
    /// same as for the nodal CSV path.
    void setNodeDopingFromDecimal(Index nodeId, std::string_view donors, std::string_view acceptors,
                                  ArithmeticReal concentrationFactor = 1);

    Real donors   (Index nodeId) const;
    Real acceptors(Index nodeId) const;

    /// Net doping = donors - acceptors in active internal units.
    Real netDoping(Index nodeId) const;
    ArithmeticReal netDopingForAssembly(Index nodeId) const;
    std::array<ArithmeticReal,2> dopingOperandsForAssembly(Index nodeId) const;
    bool hasDecimalNodeDoping() const noexcept { return !decimalDoping_.empty(); }

    /// Total ionized impurity concentration = donors + acceptors in active internal units.
    Real totalImpurity(Index nodeId) const;

    Index numNodes() const { return static_cast<Index>(donors_.size()); }

    /**
     * @brief Build a DopingModel from a vector of region-based doping specs.
     *
     * @param mesh   The device mesh supplying region and cell information.
     * @param specs  One entry per material region with donor/acceptor values.
     *
     * Nodes that belong to multiple regions (interface nodes) receive the
     * arithmetic average of the neighbouring regions' doping values.
     */
    static DopingModel fromMeshAndRegions(
        const DeviceMesh&                 mesh,
        const std::vector<RegionDopingSpec>& specs);

private:
    std::vector<Real> donors_;
    std::vector<Real> acceptors_;
    struct DecimalDoping {
        ArithmeticReal donors, acceptors, net;
    };
    std::vector<DecimalDoping> decimalDoping_;
};

} // namespace vela
