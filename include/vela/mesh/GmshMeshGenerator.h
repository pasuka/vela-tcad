#pragma once

#include "vela/mesh/DeviceMesh.h"
#include "vela/preprocess/DeviceIr2D.h"
#include <string>
#include <vector>

namespace vela {

struct BoundaryEdgeTag {
    Index node0 = 0;
    Index node1 = 0;
    std::string contactName;
    std::string ownerRegion;
};

struct MeshBundle2D {
    DeviceMesh mesh;
    std::vector<BoundaryEdgeTag> boundaryEdges;
};

struct GmshMeshingOptions {
    bool delaunayTriOnly = true;
};

class GmshMeshGenerator {
public:
    MeshBundle2D generate(const DeviceIr2D& ir, const GmshMeshingOptions& options = {}) const;
};

} // namespace vela

