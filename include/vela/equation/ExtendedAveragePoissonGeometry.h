#pragma once

#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"

namespace vela {

// Physical edge stiffness [F/m], and material-local box area in internal
// length units squared. Both are derived from the same AverageBox.
struct ExtendedAveragePoissonGeometry {
    std::vector<long double> edgeStiffness, chargeVolumes;
    long double permittivityReference_F_per_m;
};

using ExtendedPermittivities = std::unordered_map<std::string, long double>;
ExtendedPermittivities readExtendedPermittivities(
    const std::string& materialsFile, const MaterialDatabase& database);

ExtendedAveragePoissonGeometry buildExtendedAveragePoissonGeometry(
    const DeviceMesh& mesh, const MaterialDatabase& database, long double temperature_K);
ExtendedAveragePoissonGeometry buildExtendedAveragePoissonGeometry(
    const DeviceMesh& mesh, const MaterialDatabase& database, long double temperature_K,
    const ExtendedPermittivities& permittivities);

} // namespace vela
