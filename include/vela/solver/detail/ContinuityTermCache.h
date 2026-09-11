#pragma once

#include "vela/core/PerformanceProfiler.h"
#include "vela/equation/CoupledDDAssembler.h"
#include <cstring>

namespace vela::detail {

// One Newton solve owns this cache, after configuring its assembler and
// boundary conditions. Changes to assembler references, quantum potential,
// mesh or models require clear() before the next evaluation.
// A new solve/quantum outer iteration creates a new cache. Observable and
// feedback-substituted diagnostics deliberately do not use this cache.
class ContinuityTermCache {
public:
    ContinuityTermCache(const CoupledDDAssembler& assembler,
                        const CoupledDDBoundaryConditions& boundaries)
        : assembler_(assembler), boundaries_(boundaries) {}

    // The returned reference is valid until the next evaluation or clear().
    const std::vector<CoupledDDCarrierTermDiagnostic>& evaluate(const VectorXd& state)
    {
        if (valid_ && state.size() == state_.size() &&
            std::memcmp(state.data(), state_.data(),
                        static_cast<std::size_t>(state.size()) * sizeof(Real)) == 0) {
            incrementPerformanceCounter("newton.continuity_terms_cache_hits");
            return terms_;
        }
        incrementPerformanceCounter("newton.continuity_terms_cache_misses");
        // Invalidate before evaluating: an exception cannot leave a stale hit.
        valid_ = false;
        terms_ = assembler_.carrierContinuityEquationTermDiagnostics(state, boundaries_);
        state_ = state;
        valid_ = state.allFinite();
        return terms_;
    }

    void clear() { valid_ = false; }

private:
    const CoupledDDAssembler& assembler_;
    const CoupledDDBoundaryConditions boundaries_;
    VectorXd state_;
    std::vector<CoupledDDCarrierTermDiagnostic> terms_;
    bool valid_ = false;
};

} // namespace vela::detail
