#pragma once

#include "vela/mesh/GmshMeshGenerator.h"
#include "vela/physics/DopingModel.h"
#include <filesystem>

namespace vela {

class NeutralMeshWriter {
public:
    static void write(const MeshBundle2D& meshBundle,
                      const std::filesystem::path& outputDirectory,
                      const DopingModel* doping = nullptr);
};

} // namespace vela

