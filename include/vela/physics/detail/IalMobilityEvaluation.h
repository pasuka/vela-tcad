#pragma once

#include "vela/physics/IalMobility.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace vela::ial_detail {

// Six simultaneous forward directions, in SI input units. This implementation
// is private to the IALMob kernel; the equation layer owns potential DOFs.
struct Dual {
    Real value = 0.;
    std::array<Real, 6> derivative{};
    Dual() = default;
    Dual(Real v) : value(v) {}
};
inline Real value(Real x) { return x; }
inline Real value(const Dual& x) { return x.value; }
inline Dual operator+(const Dual& a, const Dual& b) {
    Dual r(a.value+b.value);
    for (int i=0;i<6;++i) r.derivative[i]=a.derivative[i]+b.derivative[i];
    return r;
}
inline Dual operator-(const Dual& a, const Dual& b) {
    Dual r(a.value-b.value);
    for (int i=0;i<6;++i) r.derivative[i]=a.derivative[i]-b.derivative[i];
    return r;
}
inline Dual operator-(const Dual& a) { return Dual(0.)-a; }
inline Dual operator*(const Dual& a, const Dual& b) {
    Dual r(a.value*b.value);
    for (int i=0;i<6;++i) r.derivative[i]=a.derivative[i]*b.value+a.value*b.derivative[i];
    return r;
}
inline Dual operator/(const Dual& a, const Dual& b) {
    Dual r(a.value/b.value);
    for (int i=0;i<6;++i) r.derivative[i]=(a.derivative[i]-r.value*b.derivative[i])/b.value;
    return r;
}
inline Real power(Real a, Real b) { return std::pow(a,b); }
inline Dual power(const Dual& a, Real b) {
    Dual r(std::pow(a.value,b));
    if (a.value>0.) {
        const Real slope=b*std::pow(a.value,b-1.);
        for (int i=0;i<6;++i) r.derivative[i]=slope*a.derivative[i];
    }
    return r;
}
inline Dual power(const Dual& a, const Dual& b) {
    Dual r=power(a,b.value);
    if (a.value>0.)
        for (int i=0;i<6;++i) r.derivative[i]+=r.value*std::log(a.value)*b.derivative[i];
    return r;
}
inline Real exponential(Real x) { return std::exp(x); }
inline Dual exponential(const Dual& x) {
    Dual r(std::exp(x.value));
    for (int i=0;i<6;++i) r.derivative[i]=r.value*x.derivative[i];
    return r;
}
template<class S> S inverse(const S& x) {
    return std::isinf(value(x)) ? S(0.) : S(1.)/x;
}
template<class S> S hypot(const S& a, const S& b) {
    const Real scale=std::max(std::abs(value(a)),std::abs(value(b)));
    if (scale==0.) return S(0.);
    return S(scale)*power((a/S(scale))*(a/S(scale))+(b/S(scale))*(b/S(scale)),.5);
}
template<class S> S clustered(const S& doping, Real reference, Real coefficient) {
    const S ratio=doping/S(reference);
    return doping*(S(1.)+ratio*ratio/(S(1.)+S(coefficient)*ratio*ratio));
}
template<class S> S coulomb2d(const S& doping, const S& carrier, Real d1, Real d2,
                             Real nu0, Real nu1, Real nu2, const IalMobilityParameters& p) {
    if (value(doping)==0.) return S(std::numeric_limits<Real>::infinity());
    const S normalized=doping/S(p.nDopRef);
    return hypot(S(d1)*power(carrier/S(p.nScRef),nu0)/power(normalized,nu1),
                 S(d2)/power(normalized,nu2));
}

// Same value path for scalar diagnostics and differentiated evaluations.
template<class S> std::array<S,6> evaluate(const std::array<S,6>& state,
    const IalMobilityParameters& p, bool electron, Real pMin) {
    const S inf(std::numeric_limits<Real>::infinity());
    const S nd=state[0]/S(1e6), na=state[1]/S(1e6);
    const S n=state[2]/S(1e6), h=state[3]/S(1e6);
    const S field=state[4]/S(100.), distance=state[5]*S(100.);
    const S c=electron?n:h, other=electron?h:n;
    const S inv=electron?na:nd, acc=electron?nd:na;
    const S ndStar=clustered(nd,p.nRefD,p.cRefD), naStar=clustered(na,p.nRefA,p.cRefA);
    const S nsc=ndStar+naStar+other;
    S mu3=inf, g(1.);
    if (value(nsc)>0.) {
        const S screening=S(1.)/(S(2.459/3.97e13)*power(nsc,2./3.)+
                                      S(3.828/(p.mass*1.36e20))*(c+other));
        const S clamped=value(screening)>pMin?screening:S(pMin);
        g=S(1.)-S(.89233)/power(S(.41372)+clamped*S(std::pow(1./p.mass,.28227)),.19778)
            +S(.005978)/power(clamped*S(std::pow(p.mass,.72169)),1.80618);
        const S pPower=power(screening,.6478);
        const S f=(S(.7643)*pPower+S(2.2999+6.5502*p.mass/p.otherMass))/
                  (pPower+S(2.3670-.8552*p.mass/p.otherMass));
        const S effective=(electron?ndStar:naStar)+g*(electron?naStar:ndStar)+other/f;
        const Real muN=p.muMax*p.muMax/(p.muMax-p.muMin);
        const Real muC=p.muMax*p.muMin/(p.muMax-p.muMin);
        mu3=S(muN)*nsc/effective*power(S(p.nRef)/nsc,p.alpha)+S(muC)*(c+other)/effective;
    }
    const S muInv=coulomb2d(inv,c,p.d1Inv,p.d2Inv,p.nu0Inv,p.nu1Inv,p.nu2Inv,p);
    S muAcc=coulomb2d(acc,c,p.d1Acc,p.d2Acc,p.nu0Acc,p.nu1Acc,p.nu2Acc,p);
    if (std::isfinite(value(muAcc))) muAcc=muAcc*g;
    const S inverse2=inverse(muInv)+inverse(muAcc);
    const S mu2=value(inverse2)==0.?inf:inverse(inverse2);
    const S transition=S(p.S/300.)*power(field,2./3.)-S(p.transitionP);
    const S f=value(transition)>700.?S(0.):S(1.)/(S(1.)+exponential(transition));
    const S weight=exponential(-distance/S(p.lCritC))*(S(1.)-f);
    // Infinite component mobilities mean absent scattering, not infinite
    // derivatives. Do not propagate infinity*zero through the AD arithmetic.
    const S muC=value(weight)==0.?mu3:value(weight)==1.?mu2:
        (!std::isfinite(value(mu3))||!std::isfinite(value(mu2)))?inf:
        (S(1.)-weight)*mu3+weight*mu2;
    const S damping=exponential(-distance/S(p.lCrit));
    S muPh(p.muMax), muSr=inf;
    if (value(field)>0. && value(damping)>0.) {
        const S muPh2=S(p.B)/field+S(p.C)*power(na+nd+S(p.N2),p.lambda)/power(field,1./3.);
        muPh=S(1.)/(damping/muPh2+S(1./p.muMax));
        const S exponent=S(p.A)+S(p.alphaSr)*(n+h)/power(na+nd+S(p.N1),p.nu);
        const S inverseSr=power(field,exponent)/S(p.delta)+field*field*field/S(p.eta);
        muSr=power(na+nd+S(p.N2),p.lambdaSr)/inverseSr;
    }
    const S mobility=S(1.)/(inverse(muC)+inverse(muPh)+damping*inverse(muSr));
    return {mobility*S(1e-4),mu3*S(1e-4),mu2*S(1e-4),muC*S(1e-4),muPh*S(1e-4),muSr*S(1e-4)};
}
} // namespace vela::ial_detail
