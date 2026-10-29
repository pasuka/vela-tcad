#pragma once
#include "vela/core/Types.h"

namespace vela::experimental {
struct AssemblyPopulationCapture {
    VectorXd x, electrons, holes;
    bool present = false;
    bool failed = false;
};

// Scoped by the diagnostic caller around one analytic Jacobian evaluation.
// Copies already computed inputs; never invokes population or cache evaluation.
class AssemblyPopulationCaptureScope {
public:
    explicit AssemblyPopulationCaptureScope(AssemblyPopulationCapture* sink) noexcept
        : previous_(active_) { active_ = sink; }
    ~AssemblyPopulationCaptureScope() { active_ = previous_; }
    AssemblyPopulationCaptureScope(const AssemblyPopulationCaptureScope&) = delete;
    AssemblyPopulationCaptureScope& operator=(const AssemblyPopulationCaptureScope&) = delete;
    static void record(const VectorXd& x, const VectorXd& n, const VectorXd& p) noexcept {
        if (!active_ || active_->failed) return;
        try {
            active_->x = x;
            active_->electrons = n;
            active_->holes = p;
            active_->present = true;
        } catch (...) { active_->failed = true; active_->present = false; }
    }
private:
    AssemblyPopulationCapture* previous_;
    static inline thread_local AssemblyPopulationCapture* active_ = nullptr;
};
} // namespace vela::experimental
