#pragma once
#include <array>
#include <cstdint>

namespace vela {
/// Point-local diagnostic counters; nested seconds overlap assembly time.
/// Timing is opt-in because clocks in every element perturb a benchmark.
struct IalKernelProfile {
    bool timingEnabled=false;
    std::array<std::uint64_t,3> passes{}; // values, spatial, temperature
    std::array<double,3> seconds{};
    std::uint64_t highFieldEvaluations=0,highFieldReuses=0;
    std::uint64_t screeningRequests=0,screeningHits=0;
    std::uint64_t screeningCandidateCalls=0,screeningFunctionEvaluations=0,screeningFallbacks=0;
    std::uint64_t localPreparationHits=0,localPreparationBuilds=0;
};
inline thread_local IalKernelProfile ialKernelProfile;
class IalKernelProfilingScope {
    IalKernelProfile previous_;
public:
    explicit IalKernelProfilingScope(bool timing):previous_(ialKernelProfile) {
        ialKernelProfile={};ialKernelProfile.timingEnabled=timing;
    }
    ~IalKernelProfilingScope(){ialKernelProfile=previous_;}
    IalKernelProfilingScope(const IalKernelProfilingScope&)=delete;
    IalKernelProfilingScope& operator=(const IalKernelProfilingScope&)=delete;
};
}
