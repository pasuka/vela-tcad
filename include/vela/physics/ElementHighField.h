#pragma once
#include "vela/core/ArithmeticPrecision.h"
#include "vela/physics/MobilityModel.h"
#include <cmath>
#include <stdexcept>

namespace vela {
inline bool usesMasettiElementHighField(const MobilityModelConfig& config)
{
    return config.model == "masetti_field" && config.edgeAveraging == "element_box" &&
        config.highFieldGradientDiscretization == "element_box_vertex";
}

inline bool usesElementHighField(const MobilityModelConfig& c) {
    return c.model=="phumob_field_lombardi" &&
        c.edgeAveraging=="element_box_phumob" &&
        c.surface.discretization=="element_distance_gradient" &&
        c.highFieldGradientDiscretization=="element_vertex_partial_layer";
}
// F and velocity use the active unit system. The calibrated 1 V/cm
// branch is deliberately discontinuous; derivatives below are branch-local.
template<class T> T elementCanali(const T& low,const T& field,
    const FieldMobilityParameters& p,Real cutoff) {
    using std::pow;
    if(!(p.saturationVelocity>0.) || !(p.beta>1.) || !(cutoff>0.))
        throw std::invalid_argument("Element Canali requires positive vsat, beta > 1 and cutoff");
    if(field<T(cutoff) || low<=0)return low;
    T z=low*field/T(p.saturationVelocity);
    return low/pow(T(1)+pow(z,T(p.beta)),T(1)/T(p.beta));
}
struct ElementCanaliValue {Real mobility,lowDerivative,fieldDerivative;};
inline ElementCanaliValue evaluateElementCanali(Real low,Real field,
    const FieldMobilityParameters& p,Real cutoff) {
    const vela::ArithmeticReal mu=elementCanali((vela::ArithmeticReal)low,(vela::ArithmeticReal)field,p,cutoff);
    if(field<cutoff || low<=0)return {static_cast<Real>(mu),1.,0.};
    const vela::ArithmeticReal z=(vela::ArithmeticReal)low*field/p.saturationVelocity;
    const vela::ArithmeticReal dm=std::pow(1+std::pow(z,(vela::ArithmeticReal)p.beta),-1-1/(vela::ArithmeticReal)p.beta);
    return {static_cast<Real>(mu),static_cast<Real>(dm),
        static_cast<Real>(-(vela::ArithmeticReal)low*low/p.saturationVelocity*std::pow(z,(vela::ArithmeticReal)p.beta-1)*dm)};
}
inline ElementCanaliValue evaluateElementCanaliSharedPower(Real low,Real field,
    const FieldMobilityParameters& p,Real cutoff) {
    if(!(p.saturationVelocity>0.) || !(p.beta>1.) || !(cutoff>0.))
        throw std::invalid_argument("Element Canali requires positive vsat, beta > 1 and cutoff");
    if(field<cutoff || low<=0)return {low,1.,0.};
    const vela::ArithmeticReal beta=p.beta;
    const vela::ArithmeticReal z=static_cast<vela::ArithmeticReal>(low)*field/p.saturationVelocity;
    const vela::ArithmeticReal base=1+std::pow(z,beta);
    const vela::ArithmeticReal mu=static_cast<vela::ArithmeticReal>(low)/std::pow(base,1/beta);
    const vela::ArithmeticReal dm=std::pow(base,-1-1/beta);
    return {static_cast<Real>(mu),static_cast<Real>(dm),
        static_cast<Real>(-static_cast<vela::ArithmeticReal>(low)*low/p.saturationVelocity*std::pow(z,beta-1)*dm)};
}

}
