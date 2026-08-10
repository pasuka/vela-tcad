#pragma once

#include <cstdint>

namespace vela {

/**
 * @brief Thread-local call counters for the carrier-statistics primitives.
 *
 * `gprof` reports a single self-time entry for `fermiDiracHalf`, which merges
 * the edge-flux, node-source, residual, and continuity-diagnostic call sites.
 * These counters exist so that the profiling stages can attribute the calls to
 * the stage that issued them instead of inferring the split.
 *
 * The counters are pure instrumentation: they never influence a numerical
 * result, an assembly order, or a convergence decision.
 */
struct PhysicsCallCounters {
    std::uint64_t fermiDiracHalf = 0;
    std::uint64_t fermiDiracHalfDerivative = 0;
    std::uint64_t inverseFermiDiracHalf = 0;
    std::uint64_t equilibriumStateSolves = 0;
    std::uint64_t equilibriumStateIterations = 0;
};

/**
 * @brief Mutable thread-local counter block.
 *
 * Declared in the header so the increments inside the carrier-statistics
 * primitives compile down to a thread-local load/add pair.
 */
extern thread_local PhysicsCallCounters physicsCallCounters;

} // namespace vela
