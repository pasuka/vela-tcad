#include "vela/physics/IalMobility.h"
#include "vela/core/IalKernelProfiling.h"
#include "vela/core/ElectrothermalCostProfile.h"
#include "vela/physics/detail/IalMobilityEvaluation.h"
#include "vela/physics/detail/IalMobilityGenerated.h"
#include "vela/core/PhysicsCallCounters.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <boost/math/tools/toms748_solve.hpp>

namespace vela {
namespace {
template<class T>
T screeningGDerivative(T p, T a, T b)
{
    return T(.89233L) * T(.19778L) * a / std::pow(T(.41372L) + a * p, T(1.19778L))
        - T(.005978L) * T(1.80618L) * b / std::pow(b * p, T(2.80618L));
}
template<class T>
T legacyScreeningMinimum(T mass, T temperature)
{
    const T a = std::pow(temperature / (T(300.L) * mass), T(.28227L));
    const T b = std::pow(mass * T(300.L) / temperature, T(.72169L));
    T predictedLower = T(0.L), predictedUpper = T(0.L);
    if (mass >= T(.01L) && mass <= T(100.L) && temperature <= T(1e5L)) {
        // g'(p)=0 is equivalent to p^v/(c+a*p)^u=K. Its logarithmic
        // derivative v-u*a*p/(c+a*p) is positive (v>u), so there is one
        // sign change. Predict a bracket, but retain the original bisection
        // sequence and evaluate its original derivative near that root.
        constexpr T u = T(1.19778L), v = T(2.80618L), c = T(.41372L);
        const T logK = std::log(T(.005978L) * T(1.80618L) / (T(.89233L) * T(.19778L)))
            + (T(1.L) - v) * std::log(b) - std::log(a);
        T x = T(0.L);
        for (int i = 0; i < 6; ++i) {
            const T ap = a * std::exp(x);
            x -= (v * x - u * std::log(c + ap) - logK)
                / (v - u * ap / (c + ap));
        }
        const T predicted = std::exp(x);
        const T lo = predicted * (T(1.L) - T(1e-12L)), hi = predicted * (T(1.L) + T(1e-12L));
        // The margin is well outside cancellation at the floating root.
        // Nonfinite/unsuccessful predictions use the original full search.
        if (lo > T(1e-12L) && hi < T(1e12L) && screeningGDerivative(lo,a,b) < T(0.L)
            && screeningGDerivative(hi,a,b) > T(0.L)) {
            predictedLower = lo;
            predictedUpper = hi;
        }
    }
    T lower = T(1e-12L), upper = T(1e12L);
    for (int i = 0; i < 100; ++i) {
        const T midpoint = std::sqrt(lower * upper);
        // Once geometric bisection reaches an endpoint, another identical
        // iteration cannot refine the representable root.
        if (midpoint == lower || midpoint == upper) return midpoint;
        const bool negative = predictedLower > T(0.L) && midpoint < predictedLower ? true
            : predictedUpper > T(0.L) && midpoint > predictedUpper ? false
            : screeningGDerivative(midpoint, a, b) < T(0.L);
        if (negative) lower = midpoint;
        else upper = midpoint;
    }
    return std::sqrt(lower * upper);
}

template<class T>
T candidateScreeningMinimum(T mass,T temperature,IalScreeningMethod method) {
    ++ialKernelProfile.screeningCandidateCalls;
    const auto fallback=[&]{++ialKernelProfile.screeningFallbacks;return legacyScreeningMinimum(mass,temperature);};
    if(mass<T(.01L) || mass>T(100.L) || temperature>T(1e5L))return fallback();
    constexpr T u=T(1.19778L),v=T(2.80618L),c=T(.41372L),alpha=T(.28227L);
    constexpr T eps=std::numeric_limits<T>::epsilon();
    const T logr=std::log(temperature/(T(300.L)*mass));
    const T d=std::log(T(.005978L)*T(1.80618L)/(T(.89233L)*T(.19778L)))
        +(v-T(1.L))*(alpha+T(.72169L))*logr-(v-u)*std::log(c);
    const auto inverse=[](T z){return z<T(0.L)?z/v:z/(v-u);};
    T lo=inverse(d),hi=inverse(d+u*std::log(T(2.L)));
    // Outward slack plus actual sign checks account for endpoint rounding.
    lo-=T(16.L)*eps*(T(1.L)+std::abs(lo));hi+=T(16.L)*eps*(T(1.L)+std::abs(hi));
    const auto evaluate=[&](T y){
        ++ialKernelProfile.screeningFunctionEvaluations;
        const T e=std::exp(-std::abs(y));
        const T soft=std::max(y,T(0.L))+std::log1p(e);
        const T sigmoid=y>=T(0.L)?T(1.L)/(T(1.L)+e):e/(T(1.L)+e);
        return std::array<T,4>{v*y-u*soft-d,v-u*sigmoid,
            -u*sigmoid*(T(1.L)-sigmoid),T(4.L)*eps*(T(1.L)+std::abs(d)+v*std::abs(y)+u*soft)};
    };
    const T flo=evaluate(lo)[0],fhi=evaluate(hi)[0];
    if(!std::isfinite(flo)||!std::isfinite(fhi)||flo>T(0.L)||fhi<T(0.L))return fallback();
    T y=lo+(hi-lo)/T(2.L);
    if(method==IalScreeningMethod::Toms748){
        std::uintmax_t budget=64;
        const auto interval=boost::math::tools::toms748_solve(
            [&](T x){return evaluate(x)[0];},lo,hi,flo,fhi,
            [](T a,T b){return b-a<=T(8.L)*eps*(T(1.L)+std::abs(a)+std::abs(b));},budget);
        y=interval.first+(interval.second-interval.first)/T(2.L);
    }else{
        for(int iteration=0;iteration<32;++iteration){
            const auto q=evaluate(y);
            if(std::abs(q[0])<=q[3])break;
            if(q[0]<T(0.L))lo=y;else hi=y;
            T next=y-q[0]/q[1];
            if(method==IalScreeningMethod::Halley){
                const T denominator=T(2.L)*q[1]*q[1]-q[0]*q[2];
                if(std::isfinite(denominator)&&denominator>T(0.L))next=y-T(2.L)*q[0]*q[1]/denominator;
            }
            if(!std::isfinite(next)||next<=lo||next>=hi)next=lo+(hi-lo)/T(2.L);
            if(next==y)break;
            y=next;
        }
    }
    const auto final=evaluate(y);
    const T result=std::exp(std::log(c)+y-alpha*logr);
    if(!std::isfinite(result)||result<=T(0.L)||std::abs(final[0])>final[3]*T(4.L))return fallback();
    return result;
}

template<class Scalar>
void validateParameters(const BasicIalMobilityParameters<Scalar>& p) {
    for (Scalar value : {p.muMax,p.muMin,p.alpha,p.nRef,p.mass,p.otherMass,
                      p.nRefD,p.nRefA,p.cRefD,p.cRefA,p.nDopRef,p.nScRef,
                      p.B,p.C,p.delta,p.eta,p.N1,p.lCrit,p.lCritC,
                      p.d1Inv,p.d2Inv,p.d1Acc,p.d2Acc}) {
        if (!std::isfinite(value) || value <= 0.)
            throw std::invalid_argument("IALMob requires finite positive scale parameters");
    }
    for (Scalar value : {p.S,p.transitionP,p.lambda,p.lambdaSr,p.A,p.alphaSr,p.nu,p.N2,
                      p.nu0Inv,p.nu1Inv,p.nu2Inv,p.nu0Acc,p.nu1Acc,p.nu2Acc}) {
        if (!std::isfinite(value) || value < 0.)
            throw std::invalid_argument("IALMob requires finite nonnegative exponents and offsets");
    }
    for (Scalar value : {p.theta,p.k,p.alpha1Inv,p.alpha2Inv,p.alpha1Acc,p.alpha2Acc})
        if (!std::isfinite(value))
            throw std::invalid_argument("IALMob requires finite temperature exponents");
    if (p.muMax <= p.muMin)
        throw std::invalid_argument("IALMob requires muMax > muMin");
}

template<class Scalar>
void validateState(const BasicIalMobilityState<Scalar>& state) {
    for(Scalar value:{state.donors_m3,state.acceptors_m3,state.electrons_m3,
                   state.holes_m3,state.normalField_V_per_m,state.interfaceDistance_m})
        if(!std::isfinite(value)||value<0.)throw std::invalid_argument("IALMob state must be finite and nonnegative");
    if(!std::isfinite(state.temperature_K)||state.temperature_K<50.)
        throw std::invalid_argument("IALMob requires finite temperature >= 50 K");
}



}

IalScreeningMethod ialScreeningMethod(const std::string& name) {
    if(name=="legacy")return IalScreeningMethod::Legacy;
    if(name=="newton")return IalScreeningMethod::Newton;
    if(name=="halley")return IalScreeningMethod::Halley;
    if(name=="toms748")return IalScreeningMethod::Toms748;
    if(name=="u2022_12_single_newton")return IalScreeningMethod::U202212SingleNewton;
    throw std::invalid_argument("Unknown IALMob screening root method");
}
const char* ialScreeningMethodName(IalScreeningMethod method) {
    switch(method) {
    case IalScreeningMethod::Legacy: return "legacy";
    case IalScreeningMethod::Newton: return "newton";
    case IalScreeningMethod::Halley: return "halley";
    case IalScreeningMethod::Toms748: return "toms748";
    case IalScreeningMethod::U202212SingleNewton: return "u2022_12_single_newton";
    }
    throw std::invalid_argument("Unknown IALMob screening root method");
}
Real ialScreeningMinimum(Real mass,Real temperature,IalScreeningMethod method) {
    if(!std::isfinite(mass)||mass<=0.||!std::isfinite(temperature)||temperature<50.)
        throw std::invalid_argument("Invalid IALMob screening mass/temperature");
    ElectrothermalCostTimer cost(ElectrothermalCostProfile::ScreeningRoot,
        electrothermalCostProfile.mode==ElectrothermalCostProfile::Screening);
    ++physicsCallCounters.ialScreeningMinimumSolves;
    if(method==IalScreeningMethod::U202212SingleNewton){
        // Runtime identity was established for Si at 300 K. This bounded
        // thermal interval is for opt-in consistency tests, not a claim of
        // native validation at other temperatures or arbitrary material mass.
        if((mass!=1. && mass!=1.258)||temperature<250.||temperature>500.)
            throw std::invalid_argument("u2022_12_single_newton requires Si mass 1 or 1.258 and 250 <= T <= 500 K");
        const Real cutoff=ial_detail::u202212ScreeningCutoff(mass,temperature);
        if(!std::isfinite(cutoff)||cutoff<=0.)
            throw std::runtime_error("Invalid U-2022.12 screening cutoff");
        return cutoff;
    }
    return method==IalScreeningMethod::Legacy?legacyScreeningMinimum(mass,temperature):candidateScreeningMinimum(mass,temperature,method);
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
    std::map<Key,ial_detail::Preparation<Real>,Less> generatedScalar;
    std::map<Key,ial_detail::Preparation<ial_detail::Dual>,Less> generatedDifferentiated;
    std::size_t hits=0;
};
IalMobilityPreparationCache::IalMobilityPreparationCache():impl_(std::make_unique<Impl>()) {}
IalMobilityPreparationCache::~IalMobilityPreparationCache()=default;
std::size_t IalMobilityPreparationCache::hits() const {return impl_->hits;}
std::size_t IalMobilityPreparationCache::size() const {return impl_->scalar.size()+impl_->differentiated.size()+impl_->generatedScalar.size()+impl_->generatedDifferentiated.size();}

Real IalScreeningCache::minimum(Real mass,Real temperature,IalScreeningMethod method) {
    ++ialKernelProfile.screeningRequests;
    if(!std::isfinite(mass)||mass<=0.||!std::isfinite(temperature)||temperature<50.)
        throw std::invalid_argument("Invalid IALMob screening mass/temperature");
    const auto key=std::make_tuple(mass,temperature,method);
    if(const auto found=roots_.find(key);found!=roots_.end()){++ialKernelProfile.screeningHits;return found->second;}
    const Real result=ialScreeningMinimum(mass,temperature,method);roots_.emplace(key,result);return result;
}

IalMobilityParameters IalMobility::siliconDefaults(bool electron)
{
    return ialSiliconDefaults<Real>(electron);
}

BasicIalMobilityDifferential<long double> evaluateExtendedSlotIalMobility(
    const BasicIalMobilityState<long double>& state,
    const BasicIalMobilityParameters<long double>& parameters, bool electron,
    IalScreeningMethod method)
{
    validateState(state);
    validateParameters(parameters);
    const auto compatibleMass = static_cast<Real>(parameters.mass);
    if (method == IalScreeningMethod::U202212SingleNewton &&
        ((compatibleMass != 1. && compatibleMass != 1.258) ||
         state.temperature_K < 250.L || state.temperature_K > 500.L))
        throw std::invalid_argument("extended U-2022.12 IALMob requires silicon and 250 <= T <= 500 K");

    using Dual80 = ial_detail::BasicDual<long double>;
    const std::array inputs{state.donors_m3,state.acceptors_m3,state.electrons_m3,
        state.holes_m3,state.normalField_V_per_m,state.interfaceDistance_m,state.temperature_K};
    std::array<Dual80,7> variables;
    for (std::size_t i=0; i<variables.size(); ++i) {
        variables[i] = Dual80(inputs[i]);
        variables[i].derivative[i] = 1.L;
    }
    const auto cutoff = [&] {
        if (method == IalScreeningMethod::U202212SingleNewton)
            return ial_detail::u202212ScreeningCutoff(parameters.mass,variables[6]);
        const auto minimum = method == IalScreeningMethod::Legacy
            ? legacyScreeningMinimum(parameters.mass,state.temperature_K)
            : candidateScreeningMinimum(parameters.mass,state.temperature_K,method);
        return Dual80(minimum);
    }();
    if (!std::isfinite(cutoff.value) || cutoff.value <= 0.L)
        throw std::runtime_error("Invalid extended U-2022.12 IALMob screening cutoff");
    const auto prepared = ial_detail::prepare(
        variables,parameters,electron,cutoff.value,cutoff.derivative[6]);
    const auto values = ial_detail::evaluatePrepared(prepared,variables[4],parameters);
    if (!std::isfinite(values[0].value) || values[0].value <= 0.L ||
        !std::ranges::all_of(values[0].derivative, [](long double x) { return std::isfinite(x); }))
        throw std::runtime_error("Extended IALMob produced invalid mobility or derivative");

    BasicIalMobilityDifferential<long double> result;
    result.result = {values[0].value,values[1].value,values[2].value,values[3].value,
                     values[4].value,values[5].value};
    std::copy_n(values[0].derivative.begin(),result.derivative_SI.size(),result.derivative_SI.begin());
    result.temperatureDerivative_m2_per_Vs_K = values[0].derivative[6];
    return result;
}

IalMobility::IalMobility(IalMobilityParameters p, bool electron,IalScreeningMethod method)
    : params_(p), electron_(electron),screeningMethod_(method)
{
    validateParameters(p);
    // Keep the original fast path for the isothermal coupled solver.
    pMin_ = ialScreeningMinimum(p.mass, 300.,screeningMethod_);
}

IalMobilityResult IalMobility::evaluate(const IalMobilityState& state,IalScreeningCache* cache,
    IalMobilityPreparationCache* preparation,bool generated) const
{
    validateState(state);
    if(generated&&screeningMethod_==IalScreeningMethod::U202212SingleNewton)
        throw std::invalid_argument("U-2022.12 screening compatibility requires the reference IALMob kernel");
    const auto minimum=[&]{return state.temperature_K==300.?pMin_:
        cache?cache->minimum(params_.mass,state.temperature_K,screeningMethod_):ialScreeningMinimum(params_.mass,state.temperature_K,screeningMethod_);};
    if(generated) {
        const std::array<Real,7> input{state.donors_m3,state.acceptors_m3,state.electrons_m3,
            state.holes_m3,state.normalField_V_per_m,state.interfaceDistance_m,state.temperature_K};
        const auto make=[&]{return ial_detail::generatedPrepareValue(input,params_,electron_,minimum());};
        const auto result=[&]{
            if(!preparation)return ial_detail::generatedEvaluateValue(make(),state.normalField_V_per_m,params_);
            const IalMobilityPreparationCache::Impl::Key key{this,{state.donors_m3,state.acceptors_m3,
                state.electrons_m3,state.holes_m3,state.interfaceDistance_m,state.temperature_K}};
            auto& entries=preparation->impl_->generatedScalar;auto found=entries.find(key);
            if(found==entries.end())found=entries.emplace(key,make()).first;else ++preparation->impl_->hits;
            return ial_detail::generatedEvaluateValue(found->second,state.normalField_V_per_m,params_);
        }();
        if(!std::isfinite(result[0])||result[0]<=0.)throw std::runtime_error("Generated IALMob produced invalid mobility");
        return {result[0],result[1],result[2],result[3],result[4],result[5]};
    }
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
    IalMobilityPreparationCache* preparation,bool generated) const
{
    validateState(state);
    if(generated&&screeningMethod_==IalScreeningMethod::U202212SingleNewton)
        throw std::invalid_argument("U-2022.12 screening compatibility requires the reference IALMob kernel");
    IalMobilityDifferential result;
    const std::array<Real,7> inputs{state.donors_m3,state.acceptors_m3,
        state.electrons_m3,state.holes_m3,state.normalField_V_per_m,state.interfaceDistance_m,state.temperature_K};
    std::array<ial_detail::Dual,7> variables;
    for (std::size_t i=0;i<7;++i) {
        variables[i]=ial_detail::Dual(inputs[i]);
        variables[i].derivative[i]=1.;
    }
    // Ordinary methods use the stationary envelope derivative. The explicit
    // one-step compatibility cutoff is not stationary and must carry dPcut/dT.
    const auto prepare=[&]{
        const Real minimum=state.temperature_K==300.?pMin_:
            cache?cache->minimum(params_.mass,state.temperature_K,screeningMethod_):ialScreeningMinimum(params_.mass,state.temperature_K,screeningMethod_);
        const Real slope=screeningMethod_==IalScreeningMethod::U202212SingleNewton?
            ial_detail::u202212ScreeningCutoff(params_.mass,variables[6]).derivative[6]:0.;
        return ial_detail::prepare(variables,params_,electron_,minimum,slope);
    };
    const auto differentiated=[&]{
        if(generated){
            const auto make=[&]{
                const Real minimum=state.temperature_K==300.?pMin_:
                    cache?cache->minimum(params_.mass,state.temperature_K,screeningMethod_):ialScreeningMinimum(params_.mass,state.temperature_K,screeningMethod_);
                return ial_detail::generatedPreparePartials(inputs,params_,electron_,minimum);
            };
            if(!preparation)return ial_detail::generatedEvaluatePartials(make(),state.normalField_V_per_m,params_);
            const IalMobilityPreparationCache::Impl::Key key{this,{state.donors_m3,state.acceptors_m3,
                state.electrons_m3,state.holes_m3,state.interfaceDistance_m,state.temperature_K}};
            auto& entries=preparation->impl_->generatedDifferentiated;auto found=entries.find(key);
            if(found==entries.end())found=entries.emplace(key,make()).first;else ++preparation->impl_->hits;
            return ial_detail::generatedEvaluatePartials(found->second,state.normalField_V_per_m,params_);
        }
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
