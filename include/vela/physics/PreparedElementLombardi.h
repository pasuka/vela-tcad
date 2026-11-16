#pragma once
#include "vela/core/ArithmeticPrecision.h"
#include "vela/physics/ElementLombardi.h"

namespace vela {
// Immutable coefficients for one material/node/configuration. Both constant
// expressions are retained because the original value and tangent use slightly
// different temperature arithmetic away from 300 K.
struct PreparedElementLombardi {
    vela::ArithmeticReal valueCoefficient=0.,tangentCoefficient=0.,damping=1.;
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
    const vela::ArithmeticReal D=std::max(0.,distance)*config.internalLengthToM;
    const vela::ArithmeticReal N=impurity*config.internalConcentrationToM3;
    out.damping=p.criticalLength>1.?VELA_ARITHMETIC_LITERAL(1.):std::exp(-D/(vela::ArithmeticReal)p.criticalLength);
    out.valueCoefficient=(vela::ArithmeticReal)p.C*std::pow((vela::ArithmeticReal)(temperature/300.),(vela::ArithmeticReal)(-p.k))*
        std::pow((N+(vela::ArithmeticReal)p.N2)/(vela::ArithmeticReal)p.N0,(vela::ArithmeticReal)p.lambda);
    out.tangentCoefficient=p.C*std::pow((vela::ArithmeticReal)temperature/300,-p.k)*
        std::pow((N+p.N2)/p.N0,p.lambda);
    return out;
}
// The transport assembly visits three vertices for each carrier on the same
// cell. Only the powers depend on the common field and the carrier exponent.
struct LombardiCellPowers {
    struct Entry {
        vela::ArithmeticReal base = 0.;
        vela::ArithmeticReal exponent = 0.;
        vela::ArithmeticReal value = 0.;

        vela::ArithmeticReal evaluate(vela::ArithmeticReal nextBase, vela::ArithmeticReal nextExponent) {
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
    const vela::ArithmeticReal F=std::abs(field)*config.internalFieldToVPerM;
    const vela::ArithmeticReal M=bulk*config.internalMobilityToM2PerVS;
    vela::ArithmeticReal inv=0.,power=0.,roughness=0.;
    if(!(F<=0.)) {
        if (powers) {
            power = powers->acoustic.evaluate(F, (vela::ArithmeticReal)2/3);
            auto& entry = electron ? powers->electronRoughness : powers->holeRoughness;
            roughness = entry.evaluate(F/(vela::ArithmeticReal)100, (vela::ArithmeticReal)p.A);
        } else {
            power=std::pow(F,(vela::ArithmeticReal)2/3);
            roughness=std::pow(F/(vela::ArithmeticReal)100,(vela::ArithmeticReal)p.A);
        }
        const vela::ArithmeticReal ac=F/((vela::ArithmeticReal)p.B+prepared.valueCoefficient*power);
        const vela::ArithmeticReal sr=roughness/(vela::ArithmeticReal)p.delta+F*F*F/(vela::ArithmeticReal)p.eta;
        inv=prepared.damping*((vela::ArithmeticReal)p.acousticFactor*ac+(vela::ArithmeticReal)p.roughnessFactor*sr);
    }
    const vela::ArithmeticReal mu=1/(1/M+inv),factor=(mu/M)*(mu/M);
    vela::ArithmeticReal derivative=0.;
    if(F>0.) {
        const vela::ArithmeticReal z=prepared.tangentCoefficient*power,den=p.B+z;
        const vela::ArithmeticReal ac=(p.B+z/3)/(den*den);
        const vela::ArithmeticReal sr=p.A*roughness/(p.delta*F)+3*F*F/p.eta;
        derivative=-mu*mu*prepared.damping*(p.acousticFactor*ac+p.roughnessFactor*sr)*
            config.internalFieldToVPerM/config.internalMobilityToM2PerVS;
    }
    return {static_cast<Real>(mu/config.internalMobilityToM2PerVS),
        static_cast<Real>(factor),static_cast<Real>(derivative)};
}
} // namespace vela
