#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"

namespace vela {

/**
 * @brief Return the maximum edge-projected electric field magnitude from nodal potential.
 *
 * This is an engineering diagnostic computed as max(|delta psi| / edge length)
 * over valid mesh edges. It is intentionally edge-based and should not be
 * interpreted as a higher-order reconstructed electric-field solution.
 * Degenerate edges with non-positive or non-finite length are ignored.
 */
Real maxEdgeElectricFieldMagnitude(const DeviceMesh& mesh, const VectorXd& potential_V);

} // namespace vela
