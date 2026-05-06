#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include <nlohmann/json.hpp>
#include <vector>
#include <string>
#include <unordered_map>

namespace vela {

/**
 * @brief Per-node doping concentrations for the device.
 *
 * Stores donor (Nd) and acceptor (Na) concentrations at each mesh node
 * in units of m^-3.  Net doping is defined as Nd - Na.
 *
 * Populated either directly via setNodeDoping() or via the static
 * factory fromMeshAndRegions() which maps region-based JSON doping specs
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

    Index numNodes() const { return static_cast<Index>(donors_.size()); }

    /**
     * @brief Build a DopingModel from a region-based JSON doping array.
     *
     * Expected JSON array format:
     * @code
     * [
     *   { "region": "n_region", "donors": 1e23, "acceptors": 0.0 },
     *   { "region": "p_region", "donors": 0.0,  "acceptors": 1e23 }
     * ]
     * @endcode
     *
     * Nodes that belong to multiple regions (interface nodes) receive the
     * arithmetic average of the neighbouring regions' doping values.
     */
    static DopingModel fromMeshAndRegions(const DeviceMesh& mesh,
                                          const nlohmann::json& dopingArray);

private:
    std::vector<Real> donors_;
    std::vector<Real> acceptors_;
};

} // namespace vela
