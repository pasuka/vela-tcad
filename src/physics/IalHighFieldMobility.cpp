#include "vela/physics/IalHighFieldMobility.h"
#include "vela/physics/detail/IalMobilityEvaluation.h"
#include <cmath>
#include <stdexcept>

namespace vela {
IalHighFieldResult evaluateIalHighFieldMobility(Real low, Real field, Real temperature,
    const IalHighFieldParameters& p)
{
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
