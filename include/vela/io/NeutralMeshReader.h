#pragma once

#include "vela/mesh/DeviceMesh.h"
#include "vela/core/UnitScaling.h"
#include "vela/physics/DopingModel.h"
#include <filesystem>
#include <string>

namespace vela {

class NeutralMeshReader {
public:
    DeviceMesh readDirectory(const std::filesystem::path& directory,
                             UnitScalingConfig scaling) const;
    DopingModel readDopingCsv(const std::filesystem::path& path,
                              Index nodeCount,
                              UnitScalingConfig scaling) const;
    DopingModel readNamedDopingCsv(const std::filesystem::path& path,
                                   Index nodeCount,
                                   UnitScalingConfig scaling,
                                   const std::string& fileLabel) const;
};

} // namespace vela
