#include "vela/mesh/GmshMeshGenerator.h"

#include <stdexcept>

namespace vela {

MeshBundle2D GmshMeshGenerator::generate(const DeviceIr2D&, const GmshMeshingOptions& options) const
{
    if (!options.delaunayTriOnly) {
        throw std::invalid_argument("GmshMeshGenerator: phase-1 only supports delaunayTriOnly=true.");
    }
    throw std::runtime_error(
        "GmshMeshGenerator: backend not wired yet. Planned implementation uses Gmsh C++ API "
        "for 2D Delaunay tri meshing with region/boundary tags.");
}

} // namespace vela

