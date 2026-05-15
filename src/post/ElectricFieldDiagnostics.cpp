#include "vela/post/ElectricFieldDiagnostics.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace vela {

Real maxEdgeElectricFieldMagnitude(const DeviceMesh& mesh, const VectorXd& potential_V)
{
    if (potential_V.size() < static_cast<int>(mesh.numNodes())) {
        throw std::invalid_argument(
            "maxEdgeElectricFieldMagnitude: potential vector has fewer entries than mesh nodes");
    }

    Real maxField = 0.0;
    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        if (!(edge.length > 0.0) || !std::isfinite(edge.length))
            continue;
        const int n0 = static_cast<int>(edge.n0);
        const int n1 = static_cast<int>(edge.n1);
        const Real dpsi = potential_V(n1) - potential_V(n0);
        const Real field = std::abs(dpsi) / edge.length;
        if (std::isfinite(field))
            maxField = std::max(maxField, field);
    }
    return maxField;
}

} // namespace vela
