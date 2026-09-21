#pragma once

#include "vela/physics/IalMobility.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace vela::ial_detail {

// Seven simultaneous forward directions, in SI input units. This implementation
// is private to the IALMob kernel; the equation layer owns potential DOFs.
struct Dual {
    Real value = 0.;
    std::array<Real, 7> derivative{};
    Dual() = default;
    Dual(Real v) : value(v) {}
};
inline Real value(Real x) { return x; }
inline Real value(const Dual& x) { return x.value; }
inline Dual operator+(const Dual& a, const Dual& b) {
    Dual r(a.value+b.value);
    for (int i=0;i<7;++i) r.derivative[i]=a.derivative[i]+b.derivative[i];
    return r;
}
inline Dual operator-(const Dual& a, const Dual& b) {
    Dual r(a.value-b.value);
    for (int i=0;i<7;++i) r.derivative[i]=a.derivative[i]-b.derivative[i];
    return r;
}
inline Dual operator-(const Dual& a) { return Dual(0.)-a; }
inline Dual operator*(const Dual& a, const Dual& b) {
    Dual r(a.value*b.value);
    for (int i=0;i<7;++i) r.derivative[i]=a.derivative[i]*b.value+a.value*b.derivative[i];
    return r;
}
inline Dual operator/(const Dual& a, const Dual& b) {
    Dual r(a.value/b.value);
    for (int i=0;i<7;++i) r.derivative[i]=(a.derivative[i]-r.value*b.derivative[i])/b.value;
    return r;
}
inline Real power(Real a, Real b) { return std::pow(a,b); }
inline Dual power(const Dual& a, Real b) {
    Dual r(std::pow(a.value,b));
    if (a.value>0.) {
        const Real slope=b*std::pow(a.value,b-1.);
        for (int i=0;i<7;++i) r.derivative[i]=slope*a.derivative[i];
    }
    return r;
}
inline Dual power(const Dual& a, const Dual& b) {
    Dual r=power(a,b.value);
    if (a.value>0.)
        for (int i=0;i<7;++i) r.derivative[i]+=r.value*std::log(a.value)*b.derivative[i];
    return r;
}
inline Real exponential(Real x) { return std::exp(x); }
inline Dual exponential(const Dual& x) {
    Dual r(std::exp(x.value));
    for (int i=0;i<7;++i) r.derivative[i]=r.value*x.derivative[i];
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
                             Real nu0, Real nu1, Real nu2, Real a1, Real a2, const S& t, const IalMobilityParameters& p) {
    if (value(doping)==0.) return S(std::numeric_limits<Real>::infinity());
    const S normalized=doping/S(p.nDopRef);
    return hypot(S(d1)*power(t,a1)*power(carrier/S(p.nScRef),nu0)/power(normalized,nu1),
                 S(d2)*power(t,a2)/power(normalized,nu2));
}

// Terms independent of the cell's normal field. All derivatives retain their
// original seven SI directions, including local temperature and distance.
template<class S> struct Preparation {
    S mu3,mu2,transitionScale,coulombDamping,damping,phonon3d;
    S phononNumerator,phononTemperaturePower,roughnessExponent,roughnessNumerator;
};
template<class S> Preparation<S> prepare(const std::array<S,7>& state,
    const IalMobilityParameters& p, bool electron, Real pMin) {
    const S inf(std::numeric_limits<Real>::infinity());
    const S nd=state[0]/S(1e6), na=state[1]/S(1e6);
    const S n=state[2]/S(1e6), h=state[3]/S(1e6);
    const S distance=state[5]*S(100.);
    const S temperature=state[6], t=temperature/S(300.);
    const S c=electron?n:h, other=electron?h:n;
    const S inv=electron?na:nd, acc=electron?nd:na;
    const S ndStar=clustered(nd,p.nRefD,p.cRefD), naStar=clustered(na,p.nRefA,p.cRefA);
    const S nsc=ndStar+naStar+other;
    S mu3=inf, g(1.);
    if (value(nsc)>0.) {
        const S screening=t*t/(S(2.459/3.97e13)*power(nsc,2./3.)+
                                      S(3.828/(p.mass*1.36e20))*(c+other));
        const S clamped=value(screening)>pMin?screening:S(pMin);
        g=S(1.)-S(.89233)/power(S(.41372)+clamped*power(t/S(p.mass),.28227),.19778)
            +S(.005978)/power(clamped*power(S(p.mass)/t,.72169),1.80618);
        const S pPower=power(screening,.6478);
        const S f=(S(.7643)*pPower+S(2.2999+6.5502*p.mass/p.otherMass))/
                  (pPower+S(2.3670-.8552*p.mass/p.otherMass));
        const S effective=(electron?ndStar:naStar)+g*(electron?naStar:ndStar)+other/f;
        const Real muN=p.muMax*p.muMax/(p.muMax-p.muMin);
        const Real muC=p.muMax*p.muMin/(p.muMax-p.muMin);
        mu3=S(muN)*power(t,3.*p.alpha-1.5)*nsc/effective*power(S(p.nRef)/nsc,p.alpha)
            +S(muC)*power(t,-.5)*(c+other)/effective;
    }
    const S muInv=coulomb2d(inv,c,p.d1Inv,p.d2Inv,p.nu0Inv,p.nu1Inv,p.nu2Inv,p.alpha1Inv,p.alpha2Inv,t,p);
    S muAcc=coulomb2d(acc,c,p.d1Acc,p.d2Acc,p.nu0Acc,p.nu1Acc,p.nu2Acc,p.alpha1Acc,p.alpha2Acc,t,p);
    if (std::isfinite(value(muAcc))) muAcc=muAcc*g;
    const S inverse2=inverse(muInv)+inverse(muAcc);
    const S mu2=value(inverse2)==0.?inf:inverse(inverse2);
    return {mu3,mu2,S(p.S)/temperature,exponential(-distance/S(p.lCritC)),
        exponential(-distance/S(p.lCrit)),S(p.muMax)*power(t,-p.theta),
        S(p.C)*power(na+nd+S(p.N2),p.lambda),power(t,p.k),
        S(p.A)+S(p.alphaSr)*(n+h)/power(na+nd+S(p.N1),p.nu),
        power(na+nd+S(p.N2),p.lambdaSr)};
}

template<class S> std::array<S,6> evaluatePrepared(const Preparation<S>& q,
    const S& field_SI,const IalMobilityParameters& p) {
    const S inf(std::numeric_limits<Real>::infinity());
    const S field=field_SI/S(100.);
    const S& mu3=q.mu3;const S& mu2=q.mu2;
    const S transition=q.transitionScale*power(field,2./3.)-S(p.transitionP);
    const S f=value(transition)>700.?S(0.):S(1.)/(S(1.)+exponential(transition));
    const S weight=q.coulombDamping*(S(1.)-f);
    // Infinite component mobilities mean absent scattering, not infinite
    // derivatives. Do not propagate infinity*zero through the AD arithmetic.
    const S muC=value(weight)==0.?mu3:value(weight)==1.?mu2:
        (!std::isfinite(value(mu3))||!std::isfinite(value(mu2)))?inf:
        (S(1.)-weight)*mu3+weight*mu2;
    const S& damping=q.damping;
    const S& muPh3=q.phonon3d;
    S muPh=muPh3, muSr=inf;
    if (value(field)>0. && value(damping)>0.) {
        const S muPh2=S(p.B)/field+q.phononNumerator/(power(field,1./3.)*q.phononTemperaturePower);
        muPh=S(1.)/(damping/muPh2+S(1.)/muPh3);
        const S inverseSr=power(field,q.roughnessExponent)/S(p.delta)+field*field*field/S(p.eta);
        muSr=q.roughnessNumerator/inverseSr;
    }
    const S mobility=S(1.)/(inverse(muC)+inverse(muPh)+damping*inverse(muSr));
    return {mobility*S(1e-4),mu3*S(1e-4),mu2*S(1e-4),muC*S(1e-4),muPh*S(1e-4),muSr*S(1e-4)};
}
// Same arithmetic for uncached scalar diagnostics and differentiated callers.
template<class S> std::array<S,6> evaluate(const std::array<S,7>& state,
    const IalMobilityParameters& p,bool electron,Real pMin) {
    return evaluatePrepared(prepare(state,p,electron,pMin),state[4],p);
}
} // namespace vela::ial_detail
