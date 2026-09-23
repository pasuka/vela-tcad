#pragma once
#include "vela/io/StateArchive.h"
#include "vela/core/UnitScaling.h"
#include "vela/solver/GummelSolver.h"

namespace vela {
class DDStateArchiveScope {
public:
    explicit DDStateArchiveScope(nlohmann::json metadata);
    ~DDStateArchiveScope();
    DDStateArchiveScope(const DDStateArchiveScope&) = delete;
    DDStateArchiveScope& operator=(const DDStateArchiveScope&) = delete;
    nlohmann::json metadata;
private:
    const DDStateArchiveScope* previous_;
};
const nlohmann::json* activeDDStateArchiveMetadata();
// File values are SI; internal densities follow the explicit unit system.
StateArchive archiveDDSolution(const DDSolution& solution, nlohmann::json metadata,
                               UnitScalingConfig scaling = UnitScalingConfig{});
DDSolution restoreDDSolution(const StateArchive& archive,
                              UnitScalingConfig scaling = UnitScalingConfig{});
} // namespace vela
