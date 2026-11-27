#pragma once

#include "vela/core/Types.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace vela {

// The nonlinear state contains physical internal concentrations, not their logs.
// Linear column scaling is independent of this storage contract.
struct PopulationUpdateConfig {
    std::string protection = "local_nonpositive";
    Real localRetainedFraction = 0.01;
    Real fractionTau = 0.99;
    Real minimumDensity = 0.; // Internal concentration units; opt-in, not a reference error.
};

inline void validatePopulationUpdate(const PopulationUpdateConfig& c) {
    if (!std::isfinite(c.minimumDensity) || c.minimumDensity<0.)
        throw std::invalid_argument("population minimum density must be finite and nonnegative");
    if (c.protection != "local_nonpositive" && c.protection != "fraction_to_boundary")
        throw std::invalid_argument("population protection must be local_nonpositive or fraction_to_boundary");
    if (!(c.localRetainedFraction > 0. && c.localRetainedFraction < 1.) ||
        !(c.fractionTau > 0. && c.fractionTau < 1.))
        throw std::invalid_argument("population protection ratios must be finite and between zero and one");
}

// Predictor coordinates are independent of the Newton unknowns. Compute logs
// separately to avoid overflow of current/previous. The optional cutoff
// suppresses prediction; it does not impose a floor on the state.
template<class Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> predictPopulation(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& previous,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& current,
    std::type_identity_t<Scalar> ratio, std::type_identity_t<Scalar> lowDensityLimit = 0) {
    if(previous.size()!=current.size() || current.size()==0 || current.size()%3 ||
       !previous.allFinite() || !current.allFinite() || !std::isfinite(ratio) ||
       !std::isfinite(lowDensityLimit) || lowDensityLimit<0.)
        throw std::invalid_argument("invalid population predictor history");
    const int n=current.size()/3;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> next = current + ratio * (current - previous);
    for(int k=n;k<current.size();++k) {
        if(previous[k]==0. && current[k]==0.){next[k]=0.;continue;}
        if(!(previous[k]>0.) || !(current[k]>0.))throw std::domain_error("invalid population predictor domain");
        if(lowDensityLimit>0. && current[k]<=lowDensityLimit){next[k]=current[k];continue;}
        next[k]=std::exp(std::log(current[k])+ratio*(std::log(current[k])-std::log(previous[k])));
        if(!(next[k]>0.) || !std::isfinite(next[k]))throw std::domain_error("unrepresentable population prediction");
    }
    if(!next.allFinite())throw std::domain_error("nonfinite population prediction");
    return next;
}

template<class Scalar>
struct BasicPopulationTrialLimits {
    Scalar localRetainedFraction, fractionTau, minimumDensity;
};

template<class Scalar>
struct BasicPopulationTrial {
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> state;
    std::vector<int> protectedElectrons, protectedHoles;
    std::vector<int> flooredElectrons, flooredHoles;
    Scalar alpha = 1;
};
using PopulationTrial = BasicPopulationTrial<Real>;

// Re-evaluate the original candidate for every alpha. Contact and oxide slots
// are not projected: prescribed equations and their exact updates own them.
template<class Scalar>
BasicPopulationTrial<Scalar> makePopulationTrial(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& before,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& raw, std::type_identity_t<Scalar> alpha,
    const std::vector<bool>& electronFree, const std::vector<bool>& holeFree,
    const PopulationUpdateConfig& c,
    const std::type_identity_t<BasicPopulationTrialLimits<Scalar>>* exact = nullptr) {
    validatePopulationUpdate(c);
    const auto limits = exact ? *exact : BasicPopulationTrialLimits<Scalar>{
        Scalar(c.localRetainedFraction), Scalar(c.fractionTau), Scalar(c.minimumDensity)};
    if (!(limits.localRetainedFraction > 0 && limits.localRetainedFraction < 1) ||
        !(limits.fractionTau > 0 && limits.fractionTau < 1) ||
        !(limits.minimumDensity >= 0) || !std::isfinite(limits.minimumDensity))
        throw std::invalid_argument("invalid arithmetic population trial limits");
    const int n = static_cast<int>(electronFree.size());
    if (!n || holeFree.size() != electronFree.size() || before.size() != 3*n ||
        raw.size() != before.size() || !before.allFinite() || !raw.allFinite() ||
        !(alpha > 0. && alpha <= 1.))
        throw std::invalid_argument("invalid population trial");
    BasicPopulationTrial<Scalar> result;
    result.alpha = alpha;
    for (int b=1; b<=2; ++b) for (int i=0; i<n; ++i) {
        if (!(b==1 ? electronFree[i] : holeFree[i])) continue;
        const int k=b*n+i;
        if (!(before[k]>0.)) throw std::domain_error("population trial requires positive free carriers");
        if (c.protection=="fraction_to_boundary" && raw[k]<0.)
            result.alpha=std::min(result.alpha,limits.fractionTau * (before[k]/(-raw[k])));
    }
    if (!(result.alpha>0.)) throw std::domain_error("unrepresentable population step");
    result.state=before+result.alpha*raw;
    for (int b=1; b<=2; ++b) for (int i=0; i<n; ++i) {
        if (!(b==1 ? electronFree[i] : holeFree[i])) continue;
        const int k=b*n+i;
        if(!std::isfinite(result.state[k]))throw std::domain_error("nonfinite raw population candidate");
        if (c.protection=="local_nonpositive" && result.state[k]<=0.) {
            result.state[k]=limits.localRetainedFraction * before[k];
            (b==1 ? result.protectedElectrons : result.protectedHoles).push_back(i);
        }
        if(result.state[k]<limits.minimumDensity) {
            result.state[k]=limits.minimumDensity;
            (b==1 ? result.flooredElectrons : result.flooredHoles).push_back(i);
        }
        if (!(result.state[k]>0.) || !std::isfinite(result.state[k]))
            throw std::domain_error("unrepresentable positive population candidate");
    }
    if (!result.state.allFinite()) throw std::domain_error("nonfinite population state");
    return result;
}
// Preserve the expression conversions offered by the binary64 entry points.
inline VectorXd predictPopulation(const VectorXd& previous, const VectorXd& current,
                                  Real ratio, Real lowDensityLimit = 0.) {
    return predictPopulation<Real>(previous, current, ratio, lowDensityLimit);
}

inline PopulationTrial makePopulationTrial(
    const VectorXd& before, const VectorXd& raw, Real alpha,
    const std::vector<bool>& electronFree, const std::vector<bool>& holeFree,
    const PopulationUpdateConfig& config) {
    return makePopulationTrial<Real>(before, raw, alpha, electronFree, holeFree, config);
}
} // namespace vela
