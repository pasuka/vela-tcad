#pragma once
#include "vela/core/Types.h"
#include <cstddef>
#include <vector>

namespace vela::experimental {
struct AssemblyContribution {
    Index edge,cell,node;
    int row,col;
    Real value;
};
struct AssemblyAudit {
    std::size_t target=1,calls=0,limit=1000000;
    bool recording=false,failed=false,truncated=false,finished=false;
    VectorXd x,n,p;
    SparseMatrixd initial,final;
    std::vector<unsigned char> constrained;
    std::vector<AssemblyContribution> contributions;
};
// Passive copies of actual scatter operands. One selected assembly per scope.
class AssemblyAuditScope {
public:
    explicit AssemblyAuditScope(AssemblyAudit* sink) noexcept:previous_(active_){active_=sink;}
    ~AssemblyAuditScope(){active_=previous_;}
    AssemblyAuditScope(const AssemblyAuditScope&)=delete;
    AssemblyAuditScope& operator=(const AssemblyAuditScope&)=delete;
    template<class F>static void record(F&& f) noexcept {
        if(!active_||active_->failed)return;
        try{f(*active_);}catch(...){active_->failed=true;active_->recording=false;}
    }
    static void begin(const VectorXd& x,const VectorXd& n,const VectorXd& p,
                      const SparseMatrixd& j,const std::vector<bool>& mask) noexcept {
        record([&](auto& a){
            a.recording=++a.calls==a.target;
            if(!a.recording)return;
            a.x=x;a.n=n;a.p=p;a.initial=j;
            a.constrained.assign(mask.begin(),mask.end());
        });
    }
    static void add(Index edge,Index cell,Index node,int row,int col,Real value) noexcept {
        if(!active_||!active_->recording)return;
        record([&](auto& a){
            if(a.contributions.size()>=a.limit){a.truncated=true;a.recording=false;return;}
            a.contributions.push_back({edge,cell,node,row,col,value});
        });
    }
    static void end(const SparseMatrixd& j) noexcept {
        if(!active_||!active_->recording)return;
        record([&](auto& a){a.final=j;a.finished=true;a.recording=false;});
    }
private:
    AssemblyAudit* previous_;
    static inline thread_local AssemblyAudit* active_=nullptr;
};
}
