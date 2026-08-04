#pragma once

#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/preprocess/DeviceIr2D.h"

namespace vela {

class DopingProfileEvaluator {
public:
    static DopingModel evaluate(const DeviceMesh& mesh, const DeviceIr2D& ir);
};

} // namespace vela

