#include "vela/physics/IalMobility.h"
#include "vela/physics/detail/IalMobilityEvaluation.h"
#include "vela/core/PhysicsCallCounters.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vela {
namespace {
Real screeningGDerivative(Real p, Real a, Real b)
{
    return .89233 * .19778 * a / std::pow(.41372 + a * p, 1.19778)
        - .005978 * 1.80618 * b / std::pow(b * p, 2.80618);
}
Real screeningMinimum(Real mass, Real temperature)
{
    ++physicsCallCounters.ialScreeningMinimumSolves;
    const Real a = std::pow(temperature / (300. * mass), .28227);
    const Real b = std::pow(mass * 300. / temperature, .72169);
    Real predictedLower = 0., predictedUpper = 0.;
    if (mass >= .01 && mass <= 100. && temperature <= 1e5) {
        // g'(p)=0 is equivalent to p^v/(c+a*p)^u=K. Its logarithmic
        // derivative v-u*a*p/(c+a*p) is positive (v>u), so there is one
        // sign change. Predict a bracket, but retain the original bisection
        // sequence and evaluate its original derivative near that root.
        constexpr Real u = 1.19778, v = 2.80618, c = .41372;
        const Real logK = std::log(.005978 * 1.80618 / (.89233 * .19778))
            + (1. - v) * std::log(b) - std::log(a);
        Real x = 0.;
        for (int i = 0; i < 6; ++i) {
            const Real ap = a * std::exp(x);
            x -= (v * x - u * std::log(c + ap) - logK)
                / (v - u * ap / (c + ap));
        }
        const Real predicted = std::exp(x);
        const Real lo = predicted * (1. - 1e-12), hi = predicted * (1. + 1e-12);
        // The margin is well outside cancellation at the floating root.
        // Nonfinite/unsuccessful predictions use the original full search.
        if (lo > 1e-12 && hi < 1e12 && screeningGDerivative(lo,a,b) < 0.
            && screeningGDerivative(hi,a,b) > 0.) {
            predictedLower = lo;
            predictedUpper = hi;
        }
    }
    Real lower = 1e-12, upper = 1e12;
    for (int i = 0; i < 100; ++i) {
        const Real midpoint = std::sqrt(lower * upper);
        // Once geometric bisection reaches an endpoint, another identical
        // iteration cannot refine the representable root.
        if (midpoint == lower || midpoint == upper) return midpoint;
        const bool negative = predictedLower > 0. && midpoint < predictedLower ? true
            : predictedUpper > 0. && midpoint > predictedUpper ? false
            : screeningGDerivative(midpoint, a, b) < 0.;
        if (negative) lower = midpoint;
        else upper = midpoint;
    }
    return std::sqrt(lower * upper);
}

void validateState(const IalMobilityState& state) {
    for(Real value:{state.donors_m3,state.acceptors_m3,state.electrons_m3,
                   state.holes_m3,state.normalField_V_per_m,state.interfaceDistance_m})
        if(!std::isfinite(value)||value<0.)throw std::invalid_argument("IALMob state must be finite and nonnegative");
    if(!std::isfinite(state.temperature_K)||state.temperature_K<50.)
        throw std::invalid_argument("IALMob requires finite temperature >= 50 K");
}



}

struct IalMobilityPreparationCache::Impl {
    using Key=std::pair<const IalMobility*,std::array<Real,6>>;
    struct Less {
        bool operator()(const Key& a,const Key& b) const {
            return a.first!=b.first?std::less<const IalMobility*>{}(a.first,b.first):a.second<b.second;
        }
    };
    std::map<Key,ial_detail::Preparation<Real>,Less> scalar;
    std::map<Key,ial_detail::Preparation<ial_detail::Dual>,Less> differentiated;
    std::size_t hits=0;
};
IalMobilityPreparationCache::IalMobilityPreparationCache():impl_(std::make_unique<Impl>()) {}
IalMobilityPreparationCache::~IalMobilityPreparationCache()=default;
std::size_t IalMobilityPreparationCache::hits() const {return impl_->hits;}
std::size_t IalMobilityPreparationCache::size() const {return impl_->scalar.size()+impl_->differentiated.size();}

Real IalScreeningCache::minimum(Real mass,Real temperature) {
    if(!std::isfinite(mass)||mass<=0.||!std::isfinite(temperature)||temperature<50.)
        throw std::invalid_argument("Invalid IALMob screening mass/temperature");
    const auto key=std::make_pair(mass,temperature);
    if(const auto found=roots_.find(key);found!=roots_.end())return found->second;
    const Real result=screeningMinimum(mass,temperature);roots_.emplace(key,result);return result;
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
    for (Real value : {p.theta,p.k,p.alpha1Inv,p.alpha2Inv,p.alpha1Acc,p.alpha2Acc})
        if (!std::isfinite(value))
            throw std::invalid_argument("IALMob requires finite temperature exponents");
    if (p.muMax <= p.muMin)
        throw std::invalid_argument("IALMob requires muMax > muMin");
    // Keep the original fast path for the isothermal coupled solver.
    pMin_ = screeningMinimum(p.mass, 300.);
}

IalMobilityResult IalMobility::evaluate(const IalMobilityState& state,IalScreeningCache* cache,
    IalMobilityPreparationCache* preparation) const
{
    validateState(state);
    const auto minimum=[&]{return state.temperature_K==300.?pMin_:
        cache?cache->minimum(params_.mass,state.temperature_K):screeningMinimum(params_.mass,state.temperature_K);};
    if(!preparation)return evaluatePrepared(state,minimum());
    const IalMobilityPreparationCache::Impl::Key key{this,{state.donors_m3,state.acceptors_m3,
        state.electrons_m3,state.holes_m3,state.interfaceDistance_m,state.temperature_K}};
    auto& entries=preparation->impl_->scalar;auto found=entries.find(key);
    if(found==entries.end())found=entries.emplace(key,ial_detail::prepare<Real>({state.donors_m3,
        state.acceptors_m3,state.electrons_m3,state.holes_m3,state.normalField_V_per_m,
        state.interfaceDistance_m,state.temperature_K},params_,electron_,minimum())).first;
    else ++preparation->impl_->hits;
    const auto result=ial_detail::evaluatePrepared(found->second,state.normalField_V_per_m,params_);
    if(!std::isfinite(result[0])||result[0]<=0.)throw std::runtime_error("IALMob produced invalid mobility");
    return {result[0],result[1],result[2],result[3],result[4],result[5]};
}

IalMobilityResult IalMobility::evaluatePrepared(const IalMobilityState& state,Real minimum) const
{
    const auto values = ial_detail::evaluate<Real>({state.donors_m3,state.acceptors_m3,
        state.electrons_m3,state.holes_m3,state.normalField_V_per_m,state.interfaceDistance_m,state.temperature_K},
        params_,electron_,minimum);
    if (!std::isfinite(values[0]) || values[0] <= 0.)
        throw std::runtime_error("IALMob produced invalid mobility");
    return {values[0],values[1],values[2],values[3],values[4],values[5]};
}

IalMobilityDifferential IalMobility::evaluateWithDerivatives(const IalMobilityState& state,IalScreeningCache* cache,
    IalMobilityPreparationCache* preparation) const
{
    validateState(state);
    IalMobilityDifferential result;
    const std::array<Real,7> inputs{state.donors_m3,state.acceptors_m3,
        state.electrons_m3,state.holes_m3,state.normalField_V_per_m,state.interfaceDistance_m,state.temperature_K};
    std::array<ial_detail::Dual,7> variables;
    for (std::size_t i=0;i<7;++i) {
        variables[i]=ial_detail::Dual(inputs[i]);
        variables[i].derivative[i]=1.;
    }
    // At the clamped minimum, dG/dP = 0. The envelope derivative is the
    // partial dG/dT; differentiating the numerical minimizer is unnecessary.
    const auto prepare=[&]{
        const Real minimum=state.temperature_K==300.?pMin_:
            cache?cache->minimum(params_.mass,state.temperature_K):screeningMinimum(params_.mass,state.temperature_K);
        return ial_detail::prepare(variables,params_,electron_,minimum);
    };
    const auto differentiated=[&]{
        if(!preparation)return ial_detail::evaluatePrepared(prepare(),variables[4],params_);
        const IalMobilityPreparationCache::Impl::Key key{this,{state.donors_m3,state.acceptors_m3,
            state.electrons_m3,state.holes_m3,state.interfaceDistance_m,state.temperature_K}};
        auto& entries=preparation->impl_->differentiated;auto found=entries.find(key);
        if(found==entries.end())found=entries.emplace(key,prepare()).first;
        else ++preparation->impl_->hits;
        return ial_detail::evaluatePrepared(found->second,variables[4],params_);
    }();
    result.result={differentiated[0].value,differentiated[1].value,differentiated[2].value,
        differentiated[3].value,differentiated[4].value,differentiated[5].value};
    if(!std::isfinite(result.result.mobility_m2_per_Vs)||result.result.mobility_m2_per_Vs<=0.)
        throw std::runtime_error("IALMob produced invalid mobility");
    const auto& derivatives=differentiated[0].derivative;
    std::copy_n(derivatives.begin(),6,result.derivative_SI.begin());
    result.temperatureDerivative_m2_per_Vs_K=derivatives[6];
    if (!std::isfinite(derivatives[6]))
        throw std::runtime_error("IALMob produced a nonfinite temperature derivative");
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
