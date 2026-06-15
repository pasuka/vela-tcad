#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include <vector>
#include <string>

namespace vela {

/**
 * @brief Region-based doping specification for a single material region.
 *
 * All concentrations are in SI units [m^-3].
 */
struct RegionDopingSpec {
    std::string region;    ///< Region name (must match a Region::name in the mesh)
    Real        donors;    ///< Donor concentration Nd [m^-3]
    Real        acceptors; ///< Acceptor concentration Na [m^-3]
};

/**
 * @brief Per-node doping concentrations for the device.
 *
 * Stores donor (Nd) and acceptor (Na) concentrations at each mesh node
 * in units of m^-3.  Net doping is defined as Nd - Na.
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

    Real donors   (Index nodeId) const;
    Real acceptors(Index nodeId) const;

    /// Net doping = donors - acceptors [m^-3].
    Real netDoping(Index nodeId) const;

    /// Total ionized impurity concentration = donors + acceptors [m^-3].
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
};

} // namespace vela
