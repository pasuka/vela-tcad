#pragma once
#include <array>
#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace vela {
// Opt-in diagnostic wall intervals. Nested stages are not additive.
struct ElectrothermalCostProfile {
    enum Mode { Off, Heat, Screening, Linear } mode=Off;
    enum Stage { HeatTotal, HeatMatrixPrepare, HeatMatrixFill, HeatMatrixFinish,
                 ScreeningRoot, SolverCreate, Symbolic, Numeric, Count };
    bool residualHeat=false;
    std::array<double,Count> seconds{};
    std::array<std::uint64_t,Count> calls{};
};
inline thread_local ElectrothermalCostProfile electrothermalCostProfile;
class ElectrothermalCostScope {
    ElectrothermalCostProfile previous_;
public:
    explicit ElectrothermalCostScope(const std::string& mode):previous_(electrothermalCostProfile) {
        ElectrothermalCostProfile next;
        if(mode=="heat")next.mode=ElectrothermalCostProfile::Heat;
        else if(mode=="screening")next.mode=ElectrothermalCostProfile::Screening;
        else if(mode=="linear")next.mode=ElectrothermalCostProfile::Linear;
        else if(mode!="off")throw std::invalid_argument("Unknown diagnostic_electrothermal_cost mode");
        electrothermalCostProfile=next;
    }
    ~ElectrothermalCostScope(){electrothermalCostProfile=previous_;}
    ElectrothermalCostScope(const ElectrothermalCostScope&)=delete;
    ElectrothermalCostScope& operator=(const ElectrothermalCostScope&)=delete;
};
class ElectrothermalCostTimer {
    using Clock=std::chrono::steady_clock;
    ElectrothermalCostProfile::Stage stage_;
    bool enabled_;
    Clock::time_point start_;
public:
    ElectrothermalCostTimer(ElectrothermalCostProfile::Stage stage,bool enabled):stage_(stage),enabled_(enabled) {
        if(enabled_){++electrothermalCostProfile.calls[stage_];start_=Clock::now();}
    }
    ~ElectrothermalCostTimer(){if(enabled_)electrothermalCostProfile.seconds[stage_]+=
        std::chrono::duration<double>(Clock::now()-start_).count();}
    ElectrothermalCostTimer(const ElectrothermalCostTimer&)=delete;
    ElectrothermalCostTimer& operator=(const ElectrothermalCostTimer&)=delete;
};
class ResidualHeatCostScope {
    bool previous_;
public:
    explicit ResidualHeatCostScope(bool residual):previous_(electrothermalCostProfile.residualHeat){electrothermalCostProfile.residualHeat=residual;}
    ~ResidualHeatCostScope(){electrothermalCostProfile.residualHeat=previous_;}
};
} // namespace vela
