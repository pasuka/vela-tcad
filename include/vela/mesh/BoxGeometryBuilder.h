#pragma once

#include "vela/core/Types.h"

namespace vela {

class DeviceMesh;
struct Cell;
struct Node;

/**
 * @brief Builds 2-D box-method geometric coefficients for Tri3 meshes.
 *
 * Node control-volume areas are computed with a barycentric control volume
 * (one third of each adjacent triangle area). Edge couplings use the
 * cotangent formula with an optional positive barycentric fallback for
 * obtuse/invalid local contributions.
 */
class BoxGeometryBuilder {
public:
    struct Options {
        bool fallbackNegativeCotangent = true;
        bool warnOnNegativeCotangent = false;
    };

    static Real triangleArea(const Node& a, const Node& b, const Node& c);
    static void build(DeviceMesh& mesh);
    static void build(DeviceMesh& mesh, const Options& options);

private:
    static Real triangleArea(const DeviceMesh& mesh, const Cell& cell);
};

} // namespace vela
