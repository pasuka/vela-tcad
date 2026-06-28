#pragma once

#include "vela/mesh/DeviceMesh.h"
#include "vela/core/Types.h"
#include <string>
#include <vector>

namespace vela {

/**
 * @brief Writer for VTK legacy ASCII format (.vtk).
 *
 * Outputs the mesh and optional node/cell scalar fields that can be
 * visualized directly in ParaView or VisIt.
 */
class VTKWriter {
public:
    /**
     * @param filename  Output file path (e.g. "result.vtk").
     * @param mesh      The device mesh to write.
     */
    VTKWriter(const std::string& filename, const DeviceMesh& mesh);

    /**
     * @brief Write the mesh geometry (points + cells + region ids).
     *
     * Overwrites any existing file at the path given in the constructor.
     */
    void write();

    /**
     * @brief Append a node-centred scalar field (POINT_DATA) to the file.
     *
     * @param fieldName  Name of the field (used as the VTK dataset name).
     * @param values     One value per node, indexed by node id.
     */
    void addNodeScalar(const std::string& fieldName,
                       const std::vector<Real>& values);

    /**
     * @brief Append a node-centred vector field (POINT_DATA) to the file.
     *
     * @param fieldName  Name of the field (used as the VTK dataset name).
     * @param values     One 3-D vector per node, indexed by node id.
     */
    void addNodeVector(const std::string& fieldName,
                       const std::vector<Point3>& values);

private:
    std::string      filename_;
    const DeviceMesh& mesh_;
    mutable bool     pointDataHeaderWritten_ = false;
};

} // namespace vela
