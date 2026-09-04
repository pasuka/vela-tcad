#pragma once

#include "vela/core/UnitScaling.h"
#include "vela/mesh/DeviceMesh.h"
#include <string>

namespace vela {

/**
 * @brief Abstract interface for mesh readers.
 *
 * Concrete implementations may read JSON meshes (current stage),
 * Gmsh .msh files, or other formats.
 */
class MeshReader {
public:
    virtual ~MeshReader() = default;

    /**
     * @brief Read a mesh from the given file and return a DeviceMesh.
     * @param filename  Path to the mesh file.
     * @return Populated DeviceMesh.
     */
    virtual DeviceMesh read(const std::string& filename) = 0;
};

/**
 * @brief JSON-based mesh reader for the current prototype stage.
 *
 * Reads a simple JSON mesh file produced by the example scripts.
 * The expected JSON schema is documented in docs/config_schema.md.
 */
class JsonMeshReader : public MeshReader {
public:
    DeviceMesh read(const std::string& filename) override;
    DeviceMesh read(const std::string& filename, UnitScalingConfig scaling);
};

} // namespace vela
