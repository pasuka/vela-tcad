#include "vela/physics/IalHighFieldMobility.h"
#include "vela/core/IalKernelProfiling.h"
#include "vela/physics/detail/IalMobilityEvaluation.h"
#include <cmath>
#include <stdexcept>

namespace vela {
namespace {
void validate(Real low, Real field, Real temperature,const IalHighFieldParameters& p) {
    for(Real v:{low,p.saturationVelocity300_m_per_s,p.beta300})
        if(!std::isfinite(v)||v<=0.)throw std::invalid_argument("IALMob HFS requires positive scales");
    if(!std::isfinite(p.saturationVelocityTemperatureExponent)||!std::isfinite(p.betaTemperatureExponent)
       ||!std::isfinite(field)||field<0.||!std::isfinite(temperature)||temperature<50.)
        throw std::invalid_argument("Invalid explicit IALMob HFS input");
}
template<bool derivatives> IalHighFieldResult explicitHighField(Real low,Real field,Real temperature,
    const IalHighFieldParameters& p) {
    ++ialKernelProfile.highFieldEvaluations;validate(low,field,temperature,p);
    if(field==0.)return {low,1.,0.,0.}; // Preserve the reference cusp convention.
    const Real t=temperature/300.;
    const Real vsat=p.saturationVelocity300_m_per_s*std::pow(t,-p.saturationVelocityTemperatureExponent);
    const Real beta=p.beta300*std::pow(t,p.betaTemperatureExponent);
    const Real x=low*field/vsat,s=std::pow(x,beta),den=std::pow(1.+s,1./beta);
    const Real mu=low/den;
    IalHighFieldResult out{mu,0.,0.,0.};
    if constexpr(derivatives) {
        const Real r=s/(1.+s);
        out.lowFieldDerivative=1./(den*(1.+s));
        out.drivingFieldDerivative_m3_per_V2s=-mu*r/field;
        // d(log(mu))/dT at fixed low and field. Retain both beta(T) and vsat(T).
        const Real entropy=std::log1p(s)-(s==0.?0.:r*std::log(s));
        out.temperatureDerivative_m2_per_Vs_K=mu/temperature*
            (p.betaTemperatureExponent*entropy/beta-p.saturationVelocityTemperatureExponent*r);
    }
    for(Real v:{out.mobility_m2_per_Vs,out.lowFieldDerivative,
                out.drivingFieldDerivative_m3_per_V2s,out.temperatureDerivative_m2_per_Vs_K})
        if(!std::isfinite(v))throw std::runtime_error("Explicit IALMob HFS produced nonfinite result");
    if(mu<=0.)throw std::runtime_error("Explicit IALMob HFS produced nonpositive mobility");
    return out;
}
}
IalHighFieldResult evaluateIalHighFieldMobilityExplicit(Real m,Real e,Real t,const IalHighFieldParameters& p) {
    return explicitHighField<true>(m,e,t,p);
}
Real evaluateIalHighFieldMobilityValue(Real m,Real e,Real t,const IalHighFieldParameters& p) {
    return explicitHighField<false>(m,e,t,p).mobility_m2_per_Vs;
}
IalHighFieldResult evaluateIalHighFieldMobility(Real low, Real field, Real temperature,
    const IalHighFieldParameters& p)
{
    ++ialKernelProfile.highFieldEvaluations;
    for (Real v : {low,p.saturationVelocity300_m_per_s,p.beta300})
        if (!std::isfinite(v)||v<=0.) throw std::invalid_argument("IALMob HFS requires positive scales");
    for (Real v : {p.saturationVelocityTemperatureExponent,p.betaTemperatureExponent})
        if (!std::isfinite(v)) throw std::invalid_argument("IALMob HFS requires finite exponents");
    if (!std::isfinite(field)||field<0.||!std::isfinite(temperature)||temperature<50.)
        throw std::invalid_argument("IALMob HFS requires nonnegative field and temperature >= 50 K");
    if (field==0.) return {low,1.,0.,0.};
    using namespace ial_detail;
    Dual mu(low), drive(field), t(temperature/300.);
    mu.derivative[0]=1.;drive.derivative[1]=1.;t.derivative[2]=1./300.;
    const auto vsat=Dual(p.saturationVelocity300_m_per_s)*power(t,-p.saturationVelocityTemperatureExponent);
    const auto beta=Dual(p.beta300)*power(t,p.betaTemperatureExponent);
    const auto result=mu/power(Dual(1.)+power(mu*drive/vsat,beta),Dual(1.)/beta);
    for (Real v : {result.value,result.derivative[0],result.derivative[1],result.derivative[2]})
        if (!std::isfinite(v)) throw std::runtime_error("IALMob HFS produced nonfinite result");
    if (result.value<=0.) throw std::runtime_error("IALMob HFS produced nonpositive mobility");
    return {result.value,result.derivative[0],result.derivative[1],result.derivative[2]};
}
} // namespace vela
