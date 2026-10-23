#pragma once
#include "vela/core/Types.h"
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace vela::experimental {
struct CacheAuditEvent {
    int iteration=0;
    std::string name;
    std::uint64_t amount=0;
    bool previousValid=false;
    VectorXd key;
};
struct CacheAudit {
    std::size_t limit=100000;
    int iteration=0;
    bool failed=false,truncated=false;
    std::vector<CacheAuditEvent> events;
};
// Observes existing counters and exact continuity-cache keys. No profiler is
// enabled and no cache lookup, refresh or physics evaluation is added.
class CacheAuditScope {
public:
    explicit CacheAuditScope(CacheAudit* sink) noexcept:previous_(active_){active_=sink;}
    ~CacheAuditScope(){active_=previous_;}
    CacheAuditScope(const CacheAuditScope&)=delete;
    CacheAuditScope& operator=(const CacheAuditScope&)=delete;
    static void iteration(int value) noexcept {if(active_)active_->iteration=value;}
    static void event(std::string_view name,std::uint64_t amount=1,
                      const VectorXd* key=nullptr,bool valid=false) noexcept {
        if(!active_||active_->failed||active_->truncated)return;
        try {
            if(active_->events.size()>=active_->limit){active_->truncated=true;return;}
            CacheAuditEvent e;e.iteration=active_->iteration;e.name=name;
            e.amount=amount;e.previousValid=valid;if(key)e.key=*key;
            active_->events.push_back(std::move(e));
        }catch(...){active_->failed=true;}
    }
private:
    CacheAudit* previous_;
    static inline thread_local CacheAudit* active_=nullptr;
};
}
