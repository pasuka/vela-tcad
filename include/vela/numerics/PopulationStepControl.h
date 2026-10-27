#pragma once
#include "vela/core/Types.h"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace vela::population_step {
// Selected from 0ee79c7dd4b06a5c31676e701be571bccc6c3a85,
// scripts/diagnostics/BvMethodsPopulationNewton.h. Shared without changing
// the diagnostic API; no density floor or acceptance predicate is introduced.
struct FractionToBoundary {
    double alpha=1.;
    int limitingNode=-1, limitingCarrier=-1;
};
inline FractionToBoundary populationFractionToBoundary(
    const VectorXd& n,const VectorXd& p,const VectorXd& dn,const VectorXd& dp,
    const std::vector<bool>& active,double tau) {
    if(!(tau>0.&&tau<1.)||n.size()!=p.size()||n.size()!=dn.size()||n.size()!=dp.size()||
       active.size()!=static_cast<std::size_t>(n.size()))throw std::invalid_argument("invalid fraction-to-boundary inputs");
    FractionToBoundary out;
    for(int i=0;i<n.size();++i)if(active[i])for(int b=0;b<2;++b) {
        const double c=b?p[i]:n[i],dc=b?dp[i]:dn[i];
        if(!(c>0.)||!std::isfinite(c)||!std::isfinite(dc))throw std::domain_error("invalid positive population direction");
        if(dc<0.) {const double bound=tau*(c/(-dc));if(bound<out.alpha){out.alpha=bound;out.limitingNode=i;out.limitingCarrier=b;}}
    }
    if(!(out.alpha>0.)||!std::isfinite(out.alpha))throw std::domain_error("no representable feasible population step");
    return out;
}
// The same c/(-dc) bound evaluated as 1/(-dc/c), so forming a dimensional
// population increment cannot underflow before it limits the common step.
inline FractionToBoundary relativeFractionToBoundary(
    const VectorXd& rn,const VectorXd& rp,const std::vector<bool>& active,double tau) {
    const VectorXd ones=VectorXd::Ones(rn.size());
    return populationFractionToBoundary(ones,ones,rn,rp,active,tau);
}
inline double log1pRemainder(double x) {
    if(!std::isfinite(x)||!(x>-1.))throw std::domain_error("invalid finite population increment");
    if(std::abs(x)<1e-4)return x*x*(-.5+x*(1./3.+x*(-.25+x*(.2-x/6.))));
    return std::log1p(x)-x;
}
}
