#pragma once

#include <array>
#include <cmath>
#include <stdexcept>

namespace vela::experimental {
template<class Scalar> struct StatisticsJet {
    Scalar value, first, second;
};
template<class Scalar> struct SgDifferential {
    Scalar value;
    // Physical potentials: psi0, psi1, qf0, qf1. Mobility is held fixed.
    std::array<Scalar,4> partial;
    bool midpoint = false;
};

// Classical, isothermal Fermi SG control. The statistics callback supplies a
// mutually consistent value and two eta derivatives. No quantum correction,
// mobility chain, clipping, source, or acceptance policy is supplied here.
template<class S, class Statistics>
SgDifferential<S> fermiSgDifferential(
    const S& psi0, const S& psi1, const S& qf0, const S& qf1,
    const S& vt, const S& coefficient, const S& log0, const S& log1,
    const S& dos0, const S& dos1, const S& driftOffset, bool electron,
    Statistics&& statistics) {
    using std::abs; using std::exp; using std::expm1; using std::log;
    if (vt <= 0 || dos0 <= 0 || dos1 <= 0)
        throw std::invalid_argument("Invalid Fermi SG material input");
    const S sign = electron ? 1 : -1;
    const S eta0 = sign*(psi0-qf0)/vt+log0;
    const S eta1 = sign*(psi1-qf1)/vt+log1;
    const auto f0 = statistics(eta0), f1 = statistics(eta1);
    if (f0.value <= 0 || f1.value <= 0)
        throw std::invalid_argument("Nonpositive Fermi SG population");
    const S n0=dos0*f0.value, n1=dos1*f1.value;
    const S l=log(S(n1/n0));
    S g=1, g0=0, g1=0;
    bool midpoint=true;
    if (abs(l)>S(1e-10)) {
        g=(eta1-eta0)/l;
        if (g>0) {
            g0=(-1+g*f0.first/f0.value)/l;
            g1=(1-g*f1.first/f1.value)/l;
            midpoint=false;
        }
    }
    if (midpoint) {
        const auto fm=statistics(S((eta0+eta1)/2));
        g=1;g0=0;g1=0;
        if(fm.first>0 && fm.value/fm.first>1) {
            g=fm.value/fm.first;
            g0=g1=(1-fm.value*fm.second/(fm.first*fm.first))/2;
        }
    }
    const S a=sign*(psi1-psi0+driftOffset)/(vt*g);
    const S t=sign*(qf1-qf0)/(vt*g);
    S b,db;
    if(abs(a)<S(1e-12)) {
        b=1-a/2+a*a/12-a*a*a*a/720;
        db=-S(1)/2+a/6-a*a*a/180;
    } else {
        const S em=expm1(a);
        b=a/em;db=(em-a*exp(a))/(em*em);
    }
    const S em=expm1(t), et=exp(t);
    SgDifferential<S> result{coefficient*g*b*n1*em,{},midpoint};
    for(int k=0;k<4;++k) {
        const S e0=k==0?S(sign/vt):k==2?S(-sign/vt):S(0);
        const S e1=k==1?S(sign/vt):k==3?S(-sign/vt):S(0);
        const S dg=g0*e0+g1*e1;
        const S da=sign*S(k==1?1:k==0?-1:0)/(vt*g)-a*dg/g;
        const S dt=sign*S(k==3?1:k==2?-1:0)/(vt*g)-t*dg/g;
        const S dn=dos1*f1.first*e1;
        result.partial[k]=coefficient*(dg*b*n1*em+g*db*da*n1*em+
            g*b*dn*em+g*b*n1*et*dt);
    }
    return result;
}
} // namespace vela::experimental
