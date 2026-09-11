#pragma once

#include "vela/core/Types.h"
#include <array>
#include <vector>

namespace vela {
struct IalInterfaceSegment {
    std::size_t node0, node1;
    std::array<Real,2> semiconductorPoint_m;
};
struct IalInterfaceNode {
    Real distance_m;
    int orientationFamily;
    std::size_t nearestInterfaceVertex;
    bool onInterface;
};
/// Explicit 2D geometry policy qualified on the LDMOS T-2022.03 mesh:
/// distance to the actual segments; orientation from the nearest interface
/// vertex's equally weighted incident unit normals. Distance and orientation
/// supports are intentionally different. No reference labels are inputs.
/// Ties use the lowest vertex index. X/Y map mesh axes into crystal axes.
std::vector<IalInterfaceNode> buildIalInterfaceGeometry(
    const std::vector<std::array<Real,2>>& coordinates_m,
    const std::vector<IalInterfaceSegment>& segments,
    const std::array<Real,3>& crystalX = {1.,0.,0.},
    const std::array<Real,3>& crystalY = {0.,1.,0.});
} // namespace vela
