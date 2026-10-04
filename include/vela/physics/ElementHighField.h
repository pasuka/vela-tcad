#pragma once
#include "vela/physics/MobilityModel.h"
#include <cmath>
#include <stdexcept>

namespace vela {
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
    const long double mu=elementCanali((long double)low,(long double)field,p,cutoff);
    if(field<cutoff || low<=0)return {static_cast<Real>(mu),1.,0.};
    const long double z=(long double)low*field/p.saturationVelocity;
    const long double dm=std::pow(1+std::pow(z,(long double)p.beta),-1-1/(long double)p.beta);
    return {static_cast<Real>(mu),static_cast<Real>(dm),
        static_cast<Real>(-(long double)low*low/p.saturationVelocity*std::pow(z,(long double)p.beta-1)*dm)};
}
}
