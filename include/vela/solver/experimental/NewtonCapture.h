#pragma once

#include "vela/core/Types.h"
#include "vela/numerics/LineSearch.h"
#include "vela/equation/experimental/PopulationCapture.h"
#include <cstddef>
#include <vector>

namespace vela::experimental {

// Private, bounded capture of values already evaluated by classical Newton.
// No observer callbacks, population evaluation, serialization or acceptance
// predicates run on the production path. Consume copies after solve returns.
struct NewtonCapturedIteration {
    int iteration = 0;
    VectorXd x, residual, rowWeights;
    SparseMatrixd jacobian;
    AssemblyPopulationCapture population;
    VectorXd linearDirection; // before exact-Dirichlet and step protections
    VectorXd protectedDirection;
    LineSearchResult searched;
    bool hasLinearDirection = false;
    bool hasSearch = false;
};

struct NewtonCapture {
    std::size_t limit = 32;
    bool truncated = false;
    bool failed = false;
    bool unsupported = false;
    bool initialized = false;
    Real potentialScale = 1.;
    VectorXd electronReference, holeReference;
    std::vector<NewtonCapturedIteration> iterations;
};

// Scope is thread-local and never owns or exposes mutable solver state.
// Nested scopes restore the previous sink. Allocation errors are contained.
class NewtonCaptureScope {
public:
    explicit NewtonCaptureScope(NewtonCapture* capture) noexcept
        : previous_(active_) { active_ = capture; }
    ~NewtonCaptureScope() { active_ = previous_; }
    NewtonCaptureScope(const NewtonCaptureScope&) = delete;
    NewtonCaptureScope& operator=(const NewtonCaptureScope&) = delete;

    static bool hasCapacity() noexcept {
        return active_ && !active_->failed && !active_->unsupported &&
            active_->iterations.size() < active_->limit;
    }

    template<class Copy>
    static void record(Copy&& copy) noexcept {
        if (!active_ || active_->failed || active_->unsupported) return;
        try { copy(*active_); }
        catch (...) { active_->failed = true; }
    }
private:
    NewtonCapture* previous_;
    static inline thread_local NewtonCapture* active_ = nullptr;
};

} // namespace vela::experimental
