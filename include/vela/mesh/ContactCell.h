#pragma once

#include "vela/mesh/DeviceMesh.h"

#include <algorithm>
#include <vector>

namespace vela::detail {

inline std::vector<bool> contactNodeMask(const DeviceMesh& mesh)
{
    std::vector<bool> mask(mesh.numNodes(), false);
    for (const auto& contact : mesh.contacts())
        for (Index node : contact.node_ids)
            if (node < mask.size())
                mask[node] = true;
    return mask;
}

inline bool cellTouchesContact(const Cell& cell, const std::vector<bool>& contactNodes)
{
    return std::any_of(cell.node_ids.begin(), cell.node_ids.end(),
        [&](Index node) { return node < contactNodes.size() && contactNodes[node]; });
}

} // namespace vela::detail
