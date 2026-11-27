#pragma once

#include "vela/physics/IalMobility.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace vela::ial_detail {

// Seven simultaneous forward directions, in SI input units. This implementation
// is private to the IALMob kernel; the equation layer owns potential DOFs.
template<class Scalar>
struct BasicDual {
    Scalar value = 0;
    std::array<Scalar, 7> derivative{};
    BasicDual() = default;
    BasicDual(Scalar v) : value(v) {}
};
using Dual = BasicDual<Real>;

template<class Scalar> requires std::is_floating_point_v<Scalar>
Scalar value(Scalar x) { return x; }
template<class Scalar>
Scalar value(const BasicDual<Scalar>& x) { return x.value; }

template<class T> BasicDual<T> operator+(const BasicDual<T>& a, const BasicDual<T>& b) {
    BasicDual<T> r(a.value+b.value);
    for (int i=0;i<7;++i) r.derivative[i]=a.derivative[i]+b.derivative[i];
    return r;
}
template<class T> BasicDual<T> operator-(const BasicDual<T>& a, const BasicDual<T>& b) {
    BasicDual<T> r(a.value-b.value);
    for (int i=0;i<7;++i) r.derivative[i]=a.derivative[i]-b.derivative[i];
    return r;
}
template<class T> BasicDual<T> operator-(const BasicDual<T>& a) { return BasicDual<T>(0)-a; }
template<class T> BasicDual<T> operator*(const BasicDual<T>& a, const BasicDual<T>& b) {
    BasicDual<T> r(a.value*b.value);
    for (int i=0;i<7;++i) r.derivative[i]=a.derivative[i]*b.value+a.value*b.derivative[i];
    return r;
}
template<class T> BasicDual<T> operator/(const BasicDual<T>& a, const BasicDual<T>& b) {
    BasicDual<T> r(a.value/b.value);
    for (int i=0;i<7;++i) r.derivative[i]=(a.derivative[i]-r.value*b.derivative[i])/b.value;
    return r;
}
template<class T> requires std::is_floating_point_v<T>
T power(T a, std::type_identity_t<T> b) { return std::pow(a,b); }
template<class T> BasicDual<T> power(const BasicDual<T>& a, std::type_identity_t<T> b) {
    BasicDual<T> r(std::pow(a.value,b));
    if (a.value>0) {
        const T slope=b*std::pow(a.value,b-T(1));
        for (int i=0;i<7;++i) r.derivative[i]=slope*a.derivative[i];
    }
    return r;
}
template<class T> BasicDual<T> power(const BasicDual<T>& a, const BasicDual<T>& b) {
    auto r=power(a,b.value);
    if (a.value>0)
        for (int i=0;i<7;++i) r.derivative[i]+=r.value*std::log(a.value)*b.derivative[i];
    return r;
}

// Version-bound numerical convention observed in U-2022.12 binary64. Keep the
// recorded auxiliary constants and single update; converging this iteration
// instead would select the ordinary stationary-minimum policy.
template<class S, class T> S u202212ScreeningCutoff(T mass,const S& temperature) {
    const S u=temperature/S(T(300.L)*mass), y=power(u,T(.28227L));
    const S x=S(T(.15L))*(S(T(1.L))+u), c=S(T(7.324855641L))*power(u,-T(.727845002L));
    const S base=S(T(.41372L))+x*y, f=power(base,T(.853672965L));
    const S slope=S(T(.853672965L))*y*f/base;
    return x+(f-c*x*x)/(S(T(2.L))*c*x-slope);
}
template<class T> requires std::is_floating_point_v<T>
T exponential(T x) { return std::exp(x); }
template<class T> BasicDual<T> exponential(const BasicDual<T>& x) {
    BasicDual<T> r(std::exp(x.value));
    for (int i=0;i<7;++i) r.derivative[i]=r.value*x.derivative[i];
    return r;
}
template<class S> S inverse(const S& x) {
    return std::isinf(value(x)) ? S(0.) : S(1.)/x;
}
template<class S> S hypot(const S& a, const S& b) {
    const auto scale=std::max(std::abs(value(a)),std::abs(value(b)));
    if (scale==0.) return S(0.);
    return S(scale)*power((a/S(scale))*(a/S(scale))+(b/S(scale))*(b/S(scale)),.5);
}
template<class S, class T> S clustered(const S& doping, T reference, T coefficient) {
    const S ratio=doping/S(reference);
    return doping*(S(1.)+ratio*ratio/(S(1.)+S(coefficient)*ratio*ratio));
}
template<class S, class T> S coulomb2d(const S& doping, const S& carrier, T d1, T d2,
                             T nu0, T nu1, T nu2, T a1, T a2, const S& t,
                             const BasicIalMobilityParameters<T>& p) {
    if (value(doping)==0.) return S(std::numeric_limits<T>::infinity());
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
template<class S, class T> Preparation<S> prepare(const std::array<S,7>& state,
    const BasicIalMobilityParameters<T>& p, bool electron, std::type_identity_t<T> pMin,
    std::type_identity_t<T> minimumTemperatureDerivative=0) {
    const S inf(std::numeric_limits<T>::infinity());
    const S nd=state[0]/S(T(1e6L)), na=state[1]/S(T(1e6L));
    const S n=state[2]/S(T(1e6L)), h=state[3]/S(T(1e6L));
    const S distance=state[5]*S(T(100.L));
    const S temperature=state[6], t=temperature/S(T(300.L));
    const S c=electron?n:h, other=electron?h:n;
    const S inv=electron?na:nd, acc=electron?nd:na;
    const S ndStar=clustered(nd,p.nRefD,p.cRefD), naStar=clustered(na,p.nRefA,p.cRefA);
    const S nsc=ndStar+naStar+other;
    S mu3=inf, g(T(1.L));
    if (value(nsc)>T(0.L)) {
        const S screening=t*t/(S(T(2.459L)/T(3.97e13L))*power(nsc,T(2.L)/T(3.L))+
                                      S(T(3.828L)/(p.mass*T(1.36e20L)))*(c+other));
        S cutoff(pMin);
        if constexpr(std::is_same_v<S,BasicDual<T>>)cutoff.derivative[6]=minimumTemperatureDerivative;
        const S clamped=value(screening)>pMin?screening:cutoff;
        g=S(T(1.L))-S(T(.89233L))/power(S(T(.41372L))+clamped*power(t/S(p.mass),T(.28227L)),T(.19778L))
            +S(T(.005978L))/power(clamped*power(S(p.mass)/t,T(.72169L)),T(1.80618L));
        const S pPower=power(screening,T(.6478L));
        const S f=(S(T(.7643L))*pPower+S(T(2.2999L)+T(6.5502L)*p.mass/p.otherMass))/
                  (pPower+S(T(2.3670L)-T(.8552L)*p.mass/p.otherMass));
        const S effective=(electron?ndStar:naStar)+g*(electron?naStar:ndStar)+other/f;
        const T muN=p.muMax*p.muMax/(p.muMax-p.muMin);
        const T muC=p.muMax*p.muMin/(p.muMax-p.muMin);
        mu3=S(muN)*power(t,T(3.L)*p.alpha-T(1.5L))*nsc/effective*power(S(p.nRef)/nsc,p.alpha)
            +S(muC)*power(t,-T(.5L))*(c+other)/effective;
    }
    const S muInv=coulomb2d(inv,c,p.d1Inv,p.d2Inv,p.nu0Inv,p.nu1Inv,p.nu2Inv,p.alpha1Inv,p.alpha2Inv,t,p);
    S muAcc=coulomb2d(acc,c,p.d1Acc,p.d2Acc,p.nu0Acc,p.nu1Acc,p.nu2Acc,p.alpha1Acc,p.alpha2Acc,t,p);
    if (std::isfinite(value(muAcc))) muAcc=muAcc*g;
    const S inverse2=inverse(muInv)+inverse(muAcc);
    const S mu2=value(inverse2)==T(0.L)?inf:inverse(inverse2);
    const S impurityPower=power(na+nd+S(p.N2),p.lambda);
    return {mu3,mu2,S(p.S)/temperature,exponential(-distance/S(p.lCritC)),
        exponential(-distance/S(p.lCrit)),S(p.muMax)*power(t,-p.theta),
        S(p.C)*impurityPower,power(t,p.k),
        S(p.A)+S(p.alphaSr)*(n+h)/power(na+nd+S(p.N1),p.nu),
        p.lambdaSr==p.lambda ? impurityPower : power(na+nd+S(p.N2),p.lambdaSr)};
}

template<class S, class T> std::array<S,6> evaluatePrepared(const Preparation<S>& q,
    const S& field_SI,const BasicIalMobilityParameters<T>& p) {
    const S inf(std::numeric_limits<T>::infinity());
    const S field=field_SI/S(T(100.L));
    const S& mu3=q.mu3;const S& mu2=q.mu2;
    const S transition=q.transitionScale*power(field,T(2.L)/T(3.L))-S(p.transitionP);
    const S f=value(transition)>T(700.L)?S(T(0.L)):S(T(1.L))/(S(T(1.L))+exponential(transition));
    const S weight=q.coulombDamping*(S(T(1.L))-f);
    // Infinite component mobilities mean absent scattering, not infinite
    // derivatives. Do not propagate infinity*zero through the AD arithmetic.
    const S muC=value(weight)==T(0.L)?mu3:value(weight)==T(1.L)?mu2:
        (!std::isfinite(value(mu3))||!std::isfinite(value(mu2)))?inf:
        (S(T(1.L))-weight)*mu3+weight*mu2;
    const S& damping=q.damping;
    const S& muPh3=q.phonon3d;
    S muPh=muPh3, muSr=inf;
    if (value(field)>T(0.L) && value(damping)>T(0.L)) {
        const S muPh2=S(p.B)/field+q.phononNumerator/(power(field,T(1.L)/T(3.L))*q.phononTemperaturePower);
        muPh=S(T(1.L))/(damping/muPh2+S(T(1.L))/muPh3);
        const S inverseSr=power(field,q.roughnessExponent)/S(p.delta)+field*field*field/S(p.eta);
        muSr=q.roughnessNumerator/inverseSr;
    }
    const S mobility=S(T(1.L))/(inverse(muC)+inverse(muPh)+damping*inverse(muSr));
    return {mobility*S(T(1e-4L)),mu3*S(T(1e-4L)),mu2*S(T(1e-4L)),muC*S(T(1e-4L)),muPh*S(T(1e-4L)),muSr*S(T(1e-4L))};
}
// Same arithmetic for uncached scalar diagnostics and differentiated callers.
template<class S, class T> std::array<S,6> evaluate(const std::array<S,7>& state,
    const BasicIalMobilityParameters<T>& p,bool electron,std::type_identity_t<T> pMin) {
    return evaluatePrepared(prepare(state,p,electron,pMin),state[4],p);
}
} // namespace vela::ial_detail
