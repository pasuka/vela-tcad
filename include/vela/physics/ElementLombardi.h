#pragma once
#include "vela/physics/ElementHighField.h"
#include <cmath>
#include <stdexcept>

namespace vela {
inline bool usesElementDistanceLombardi(const MobilityModelConfig& c) {
    return (c.model=="phumob_lombardi" || usesElementHighField(c)) && c.edgeAveraging=="element_box_phumob" &&
        c.surface.discretization=="element_distance_gradient";
}

// The explicit element profile currently qualifies the carrier-independent
// Enhanced Lombardi exponent (alpha=0). Carrier dependence still enters through
// the full live PhuMob model. SI inputs and inverse SI mobility output.
template<class T>
T elementLombardiInverse(const T& field,const T& distance,const T& impurity,
                        double temperature,const LombardiParameters& p) {
    using std::exp;using std::pow;
    if(p.alpha!=0 || p.B<=0 || p.C<0 || p.N0<=0 || p.N2<0 || p.delta<=0 ||
       p.eta<=0 || p.criticalLength<=0 || p.A<1 || p.acousticFactor<0 ||
       p.roughnessFactor<0 || temperature<=0)
        throw std::invalid_argument("Element Lombardi requires valid alpha=0 parameters");
    if(field<=0)return T(0);
    T damp=p.criticalLength>1.?T(1):T(exp(-distance/T(p.criticalLength)));
    T C=T(p.C)*pow(T(temperature/300.),T(-p.k))*pow((impurity+T(p.N2))/T(p.N0),T(p.lambda));
    T ac=field/(T(p.B)+C*pow(field,T(2)/3));
    T sr=pow(field/T(100),T(p.A))/T(p.delta)+field*field*field/T(p.eta);
    return T(damp*(T(p.acousticFactor)*ac+T(p.roughnessFactor)*sr));
}

struct ElementLombardiResult { Real mobility,bulkDerivative,fieldDerivative; };
inline ElementLombardiResult evaluateElementLombardi(Real bulk,Real field,
    Real distance,Real impurity,Real temperature,const LombardiParameters& p,
    const MobilityModelConfig& config) {
    if(bulk<=0)return {0,0,0};
    if(!std::isfinite(distance) || !std::isfinite(field))return {bulk,1,0};
    const long double F=std::abs(field)*config.internalFieldToVPerM;
    const long double D=std::max(0.,distance)*config.internalLengthToM;
    const long double N=impurity*config.internalConcentrationToM3;
    const long double M=bulk*config.internalMobilityToM2PerVS;
    const long double inv=elementLombardiInverse(F,D,N,temperature,p);
    const long double mu=1/(1/M+inv),factor=(mu/M)*(mu/M);
    long double derivative=0;
    if(F>0) {
        const long double C=p.C*std::pow((long double)temperature/300,-p.k)*std::pow((N+p.N2)/p.N0,p.lambda);
        const long double z=C*std::pow(F,2.L/3),den=p.B+z;
        const long double ac=(p.B+z/3)/(den*den);
        const long double sr=p.A*std::pow(F/100,p.A)/(p.delta*F)+3*F*F/p.eta;
        const long double damp=p.criticalLength>1?1:std::exp(-D/p.criticalLength);
        derivative=-mu*mu*damp*(p.acousticFactor*ac+p.roughnessFactor*sr)*
            config.internalFieldToVPerM/config.internalMobilityToM2PerVS;
    }
    return {static_cast<Real>(mu/config.internalMobilityToM2PerVS),
        static_cast<Real>(factor),static_cast<Real>(derivative)};
}
}
