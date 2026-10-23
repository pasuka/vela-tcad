#pragma once
#include "vela/core/Types.h"
#include "vela/numerics/LineSearch.h"
#include <cstddef>
#include <vector>

namespace vela::experimental {
struct NewtonAuditIteration {
    int iteration=0;
    VectorXd x,residual,rowWeights,linearDirection,protectedDirection;
    SparseMatrixd jacobian;
    LineSearchResult searched;
    bool hasLinearDirection=false,hasSearch=false;
};
struct NewtonAudit {
    std::size_t limit=32;
    bool started=false,unsupported=false,failed=false,truncated=false;
    std::vector<NewtonAuditIteration> iterations;
};

// Independent witness for report-off runs. Owns separate copies, never reads
// the report collector, and never computes physics, hashes or floating norms.
class NewtonAuditScope {
public:
    explicit NewtonAuditScope(NewtonAudit* audit) noexcept : previous_(active_) { active_=audit; }
    ~NewtonAuditScope(){ active_=previous_; }
    NewtonAuditScope(const NewtonAuditScope&)=delete;
    NewtonAuditScope& operator=(const NewtonAuditScope&)=delete;
    template<class Copy> static void record(Copy&& copy) noexcept {
        if(!active_ || active_->failed || active_->unsupported) return;
        try { copy(*active_); } catch(...) { active_->failed=true; }
    }
    static void begin(bool supported) noexcept {
        record([&](auto& a){a.unsupported=!supported || a.started; a.started=true;});
    }
    static void system(int iteration,const VectorXd& x,const VectorXd& r,
                       const SparseMatrixd& j,const VectorXd& weights) noexcept {
        record([&](auto& a){
            if(a.iterations.size()>=a.limit){a.truncated=true;return;}
            NewtonAuditIteration e;e.iteration=iteration;e.x=x;e.residual=r;
            e.jacobian=j;e.rowWeights=weights;a.iterations.push_back(std::move(e));
        });
    }
    static void direction(int iteration,const VectorXd& step) noexcept {
        record([&](auto& a){
            if(a.truncated || a.iterations.empty())return;
            auto& e=a.iterations.back();
            if(e.iteration!=iteration){a.failed=true;return;}
            e.linearDirection=step;e.hasLinearDirection=true;
        });
    }
    static void search(int iteration,const VectorXd& step,const LineSearchResult& search) noexcept {
        record([&](auto& a){
            if(a.truncated || a.iterations.empty())return;
            auto& e=a.iterations.back();
            if(e.iteration!=iteration){a.failed=true;return;}
            e.protectedDirection=step;e.searched=search;e.hasSearch=true;
        });
    }
private:
    NewtonAudit* previous_;
    static inline thread_local NewtonAudit* active_=nullptr;
};
} // namespace vela::experimental
