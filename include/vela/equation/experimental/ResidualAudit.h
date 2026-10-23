#pragma once
#include "vela/core/Types.h"
#include <vector>
#include <cstddef>

namespace vela::experimental {
struct ResidualOperation { int row; char operation; Real operand; };
struct ResidualAudit {
    std::size_t target=1, calls=0, limit=1000000;
    bool recording=false, failed=false, truncated=false, finished=false;
    int rows=0;
    VectorXd final;
    std::vector<ResidualOperation> operations;
};
// Passive observations only. No replacement residual, population, or cache.
class ResidualAuditScope {
public:
    explicit ResidualAuditScope(ResidualAudit* sink) noexcept : previous_(active_) { active_=sink; }
    ~ResidualAuditScope() { active_=previous_; }
    ResidualAuditScope(const ResidualAuditScope&)=delete;
    ResidualAuditScope& operator=(const ResidualAuditScope&)=delete;
    static void begin(int rows) noexcept {
        if(!active_ || active_->failed)return;
        active_->recording=++active_->calls==active_->target;
        if(active_->recording)active_->rows=rows;
    }
    template<class F> static void operation(int row,char op,F&& operand) noexcept {
        if(!active_ || !active_->recording || active_->failed)return;
        try {
            if(active_->operations.size()>=active_->limit) {
                active_->truncated=true;active_->recording=false;return;
            }
            active_->operations.push_back({row,op,operand()});
        } catch(...) { active_->failed=true;active_->recording=false; }
    }
    static void end(const VectorXd& value) noexcept {
        if(!active_ || !active_->recording || active_->failed)return;
        try {active_->final=value;active_->finished=true;active_->recording=false;}
        catch(...) {active_->failed=true;active_->recording=false;}
    }
private:
    ResidualAudit* previous_;
    static inline thread_local ResidualAudit* active_=nullptr;
};
}
