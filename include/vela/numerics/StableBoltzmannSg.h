#pragma once
#include <cmath>
namespace vela::numerics {
// Derivative of C B(eta) population expm1(q), for strictly constant C.
// Form relative directions before multiplying large equilibrium terms.
template<class T> T stableBoltzmannSgDirection(const T& cp,const T& eta,const T& q,
    const T& dLogPopulation,const T& dEta,const T& dQ) {
    using std::abs;using std::exp;using std::expm1;
    const T b=eta==0?T(1):T(eta/expm1(eta));
    T db;
    if(abs(eta)<T(1e-8)) {
        const T e2=eta*eta;
        db=-T(1)/2+eta/6-eta*e2/180+eta*e2*e2/5040;
    } else {
        const T em=expm1(eta);db=(em-eta*exp(eta))/(em*em);
    }
    return cp*(b*expm1(q)*dLogPopulation+db*expm1(q)*dEta+b*exp(q)*dQ);
}
}
