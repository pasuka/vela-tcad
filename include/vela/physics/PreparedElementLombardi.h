#pragma once
#include "vela/physics/ElementLombardi.h"

namespace vela {
// Immutable coefficients for one material/node/configuration. Both constant
// expressions are retained because the original value and tangent use slightly
// different temperature arithmetic away from 300 K.
struct PreparedElementLombardi {
    long double valueCoefficient=0.,tangentCoefficient=0.,damping=1.;
    bool finiteDistance=true,validParameters=true;
};
inline PreparedElementLombardi prepareElementLombardi(
    Real distance,Real impurity,Real temperature,const LombardiParameters& p,
    const MobilityModelConfig& config) {
    PreparedElementLombardi out;out.finiteDistance=std::isfinite(distance);
    out.validParameters=!(p.alpha!=0 || p.B<=0 || p.C<0 || p.N0<=0 || p.N2<0 ||
        p.delta<=0 || p.eta<=0 || p.criticalLength<=0 || p.A<1 ||
        p.acousticFactor<0 || p.roughnessFactor<0 || temperature<=0);
    if(!out.finiteDistance || !out.validParameters)return out;
    const long double D=std::max(0.,distance)*config.internalLengthToM;
    const long double N=impurity*config.internalConcentrationToM3;
    out.damping=p.criticalLength>1.?1.L:std::exp(-D/(long double)p.criticalLength);
    out.valueCoefficient=(long double)p.C*std::pow((long double)(temperature/300.),(long double)(-p.k))*
        std::pow((N+(long double)p.N2)/(long double)p.N0,(long double)p.lambda);
    out.tangentCoefficient=p.C*std::pow((long double)temperature/300,-p.k)*
        std::pow((N+p.N2)/p.N0,p.lambda);
    return out;
}
// The residual transport visits three vertices for each carrier on the same
// cell. Only the powers depend on the common field and the carrier exponent.
struct LombardiCellPowers {
    struct Entry {
        long double base = 0.;
        long double exponent = 0.;
        long double value = 0.;

        long double evaluate(long double nextBase, long double nextExponent) {
            if (!(nextBase > 0.) || !std::isfinite(nextBase) || !std::isfinite(nextExponent))
                return std::pow(nextBase, nextExponent);
            if (base != nextBase || exponent != nextExponent) {
                value = std::pow(nextBase, nextExponent);
                base = nextBase;
                exponent = nextExponent;
            }
            return value;
        }
    };

    Entry acoustic;
    Entry electronRoughness;
    Entry holeRoughness;
};

inline ElementLombardiResult evaluatePreparedElementLombardi(Real bulk,Real field,
    const PreparedElementLombardi& prepared,const LombardiParameters& p,
    const MobilityModelConfig& config, LombardiCellPowers* powers = nullptr,
    bool electron = false) {
    if(bulk<=0)return {0,0,0};
    if(!prepared.finiteDistance || !std::isfinite(field))return {bulk,1,0};
    if(!prepared.validParameters)
        throw std::invalid_argument("Element Lombardi requires valid alpha=0 parameters");
    const long double F=std::abs(field)*config.internalFieldToVPerM;
    const long double M=bulk*config.internalMobilityToM2PerVS;
    long double inv=0.,power=0.,roughness=0.;
    if(!(F<=0.)) {
        if (powers) {
            power = powers->acoustic.evaluate(F, (long double)2/3);
            auto& entry = electron ? powers->electronRoughness : powers->holeRoughness;
            roughness = entry.evaluate(F/(long double)100, (long double)p.A);
        } else {
            power=std::pow(F,(long double)2/3);
            roughness=std::pow(F/(long double)100,(long double)p.A);
        }
        const long double ac=F/((long double)p.B+prepared.valueCoefficient*power);
        const long double sr=roughness/(long double)p.delta+F*F*F/(long double)p.eta;
        inv=prepared.damping*((long double)p.acousticFactor*ac+(long double)p.roughnessFactor*sr);
    }
    const long double mu=1/(1/M+inv),factor=(mu/M)*(mu/M);
    long double derivative=0.;
    if(F>0.) {
        const long double z=prepared.tangentCoefficient*power,den=p.B+z;
        const long double ac=(p.B+z/3)/(den*den);
        const long double sr=p.A*roughness/(p.delta*F)+3*F*F/p.eta;
        derivative=-mu*mu*prepared.damping*(p.acousticFactor*ac+p.roughnessFactor*sr)*
            config.internalFieldToVPerM/config.internalMobilityToM2PerVS;
    }
    return {static_cast<Real>(mu/config.internalMobilityToM2PerVS),
        static_cast<Real>(factor),static_cast<Real>(derivative)};
}
} // namespace vela
