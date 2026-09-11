#include "vela/physics/IalMobility.h"
#include "vela/physics/detail/IalMobilityEvaluation.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vela {
namespace {
Real screeningGDerivative(Real p, Real mass)
{
    const Real a = std::pow(1. / mass, .28227);
    const Real b = std::pow(mass, .72169);
    return .89233 * .19778 * a / std::pow(.41372 + a * p, 1.19778)
        - .005978 * 1.80618 * b / std::pow(b * p, 2.80618);
}


}

IalMobilityParameters IalMobility::siliconDefaults(bool electron)
{
    IalMobilityParameters p;
    if (!electron) {
        p.muMax = 470.5; p.muMin = 44.9; p.theta = 2.247;
        p.alpha = .719; p.nRef = 2.23e17;
        p.mass = 1.258; p.otherMass = 1.;
    }
    return p;
}

IalMobility::IalMobility(IalMobilityParameters p, bool electron)
    : params_(p), electron_(electron)
{
    for (Real value : {p.muMax,p.muMin,p.alpha,p.nRef,p.mass,p.otherMass,
                      p.nRefD,p.nRefA,p.cRefD,p.cRefA,p.nDopRef,p.nScRef,
                      p.B,p.C,p.delta,p.eta,p.N1,p.lCrit,p.lCritC,
                      p.d1Inv,p.d2Inv,p.d1Acc,p.d2Acc}) {
        if (!std::isfinite(value) || value <= 0.)
            throw std::invalid_argument("IALMob requires finite positive scale parameters");
    }
    for (Real value : {p.S,p.transitionP,p.lambda,p.lambdaSr,p.A,p.alphaSr,p.nu,p.N2,
                      p.nu0Inv,p.nu1Inv,p.nu2Inv,p.nu0Acc,p.nu1Acc,p.nu2Acc}) {
        if (!std::isfinite(value) || value < 0.)
            throw std::invalid_argument("IALMob requires finite nonnegative exponents and offsets");
    }
    if (p.muMax <= p.muMin)
        throw std::invalid_argument("IALMob requires muMax > muMin");
    // Locate the minimum once, not during each mobility evaluation.
    Real lower = 1e-12, upper = 1e12;
    for (int i = 0; i < 100; ++i) {
        const Real midpoint = std::sqrt(lower * upper);
        if (screeningGDerivative(midpoint, p.mass) < 0.) lower = midpoint;
        else upper = midpoint;
    }
    pMin_ = std::sqrt(lower * upper);
}

IalMobilityResult IalMobility::evaluate(const IalMobilityState& state) const
{
    for (Real value : {state.donors_m3,state.acceptors_m3,state.electrons_m3,
                      state.holes_m3,state.normalField_V_per_m,state.interfaceDistance_m}) {
        if (!std::isfinite(value) || value < 0.)
            throw std::invalid_argument("IALMob state must be finite and nonnegative");
    }
    const auto values = ial_detail::evaluate<Real>({state.donors_m3,state.acceptors_m3,
        state.electrons_m3,state.holes_m3,state.normalField_V_per_m,state.interfaceDistance_m},
        params_,electron_,pMin_);
    if (!std::isfinite(values[0]) || values[0] <= 0.)
        throw std::runtime_error("IALMob produced invalid mobility");
    return {values[0],values[1],values[2],values[3],values[4],values[5]};
}

IalMobilityDifferential IalMobility::evaluateWithDerivatives(const IalMobilityState& state) const
{
    IalMobilityDifferential result{evaluate(state),{}};
    const std::array<Real,6> inputs{state.donors_m3,state.acceptors_m3,
        state.electrons_m3,state.holes_m3,state.normalField_V_per_m,state.interfaceDistance_m};
    std::array<ial_detail::Dual,6> variables;
    for (std::size_t i=0;i<6;++i) {
        variables[i]=ial_detail::Dual(inputs[i]);
        variables[i].derivative[i]=1.;
    }
    result.derivative_SI=ial_detail::evaluate(variables,params_,electron_,pMin_)[0].derivative;
    for (Real derivative:result.derivative_SI)
        if (!std::isfinite(derivative))
            throw std::runtime_error("IALMob produced a nonfinite derivative");
    return result;
}

int IalMobility::orientationFamily(const std::array<Real,3>& normal)
{
    std::array<Real,3> values;
    for (std::size_t i=0;i<3;++i) {
        if (!std::isfinite(normal[i])) throw std::invalid_argument("Invalid crystal normal");
        values[i]=std::abs(normal[i]);
    }
    std::sort(values.rbegin(),values.rend());
    if (values[0]==0.) throw std::invalid_argument("Zero crystal normal");
    // Rescale first to avoid overflow for large, non-unit normals.
    const Real b=values[1]/values[0],c=values[2]/values[0];
    const std::array<Real,3> score{1.,(1.+b)/std::sqrt(2.),(1.+b+c)/std::sqrt(3.)};
    const auto selected=std::max_element(score.begin(),score.end())-score.begin();
    return std::array<int,3>{100,110,111}[selected];
}
} // namespace vela
