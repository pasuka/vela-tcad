#pragma once

#include "vela/core/Types.h"
#include "vela/physics/IalScreeningRoot.h"
#include <array>
#include <map>
#include <memory>
#include <tuple>

namespace vela {

/// Independent low-field IALMob kernel, T-2022.03 equations 290--311.
/// Parameters use the manual's cm-based units; state/result boundaries use SI.
/// FullPhuMob, PhononCombination=1, no stress/thin-layer corrections.
/// Impurity weighting factors/exponents in the phonon/roughness terms stay at
/// their manual defaults of one; ClusteringEverywhere is disabled. Temperature
/// dependencies follow equations 293, 295--299, 305--309. T >= 50 K only;
/// the optional low-temperature phonon exponent correction is not enabled.
/// This kernel does not enable IALMob in the coupled solver configuration.
template<class Scalar>
struct BasicIalMobilityParameters {
    Scalar muMax = 1417.L, muMin = 52.2L, theta = 2.285L;
    Scalar alpha = .68L, nRef = 9.68e16L;
    Scalar mass = 1.L, otherMass = 1.258L;
    Scalar nRefD = 4e20L, nRefA = 7.2e20L, cRefD = .21L, cRefA = .5L;
    Scalar nDopRef = 1e18L, nScRef = 1e18L;
    Scalar S = .3042L, transitionP = 4.L;
    Scalar B = 9e5L, C = 4400.L, lambda = .057L, k = 1.L;
    Scalar delta = 3.97e13L, eta = 1e50L, lambdaSr = .057L;
    Scalar A = 2.L, alphaSr = 0.L, nu = 0.L, N1 = 1.L, N2 = 1.L;
    Scalar lCrit = 1e3L, lCritC = 1e3L;
    Scalar d1Inv = 135.L, d2Inv = 40.L, nu0Inv = 1.5L, nu1Inv = 2.L, nu2Inv = .5L;
    Scalar alpha1Inv = 0.L, alpha2Inv = 0.L, alpha1Acc = 0.L, alpha2Acc = 0.L;
    Scalar d1Acc = 135.L, d2Acc = 40.L, nu0Acc = 1.5L, nu1Acc = 2.L, nu2Acc = .5L;
};

template<class Scalar>
struct BasicIalMobilityState {
    Scalar donors_m3 = 0.L, acceptors_m3 = 0.L;
    Scalar electrons_m3 = 0.L, holes_m3 = 0.L;
    Scalar normalField_V_per_m = 0.L, interfaceDistance_m = 0.L;
    Scalar temperature_K = 300.L;
};

template<class Scalar>
struct BasicIalMobilityResult {
    Scalar mobility_m2_per_Vs;
    Scalar coulomb3d_m2_per_Vs, coulomb2d_m2_per_Vs, coulomb_m2_per_Vs;
    Scalar phonon_m2_per_Vs, roughness_m2_per_Vs;
};

/// Derivatives of the final low-field mobility with respect to SI state fields
/// in this order: Nd, Na, n, p, |Enormal|, interface distance. At a zero input
/// with a fractional-power cusp the selected branch slope is zero, not a
/// claim of a classical derivative. Geometry and orientation are fixed.
template<class Scalar>
struct BasicIalMobilityDifferential {
    BasicIalMobilityResult<Scalar> result;
    std::array<Scalar, 6> derivative_SI{};
    /// Partial at fixed Nd, Na, n, p, field and distance (m2 / (V s K)).
    Scalar temperatureDerivative_m2_per_Vs_K = 0.L;
};

using IalMobilityParameters = BasicIalMobilityParameters<Real>;
using IalMobilityState = BasicIalMobilityState<Real>;
using IalMobilityResult = BasicIalMobilityResult<Real>;
using IalMobilityDifferential = BasicIalMobilityDifferential<Real>;

template<class Scalar>
BasicIalMobilityParameters<Scalar> ialSiliconDefaults(bool electron) {
    BasicIalMobilityParameters<Scalar> p;
    if (!electron) {
        p.muMax = 470.5L; p.muMin = 44.9L; p.theta = 2.247L;
        p.alpha = .719L; p.nRef = 2.23e17L;
        p.mass = 1.258L; p.otherMass = 1.L;
    }
    return p;
}

// Evaluate the SLOT reference kernel from wide state and PAR parameters using
// the configured screening branch, independently of the Real caches.
BasicIalMobilityDifferential<long double> evaluateExtendedSlotIalMobility(
    const BasicIalMobilityState<long double>& state,
    const BasicIalMobilityParameters<long double>& parameters, bool electron,
    IalScreeningMethod method);

/// Scoped to one transport-state preparation. The screening minimizer depends
/// only on mass and temperature; exact keys never reuse another temperature.
class IalScreeningCache {
public:
    Real minimum(Real mass,Real temperature,IalScreeningMethod method=defaultIalScreeningMethod);
    std::size_t size() const {return roots_.size();}
private:
    std::map<std::tuple<Real,Real,IalScreeningMethod>,Real> roots_;
};

/// Local to one transport-state preparation. Model objects must remain alive
/// and immutable; exact keys include doping, densities, distance and temperature.
class IalMobilityPreparationCache {
public:
    IalMobilityPreparationCache();
    ~IalMobilityPreparationCache();
    IalMobilityPreparationCache(const IalMobilityPreparationCache&)=delete;
    IalMobilityPreparationCache& operator=(const IalMobilityPreparationCache&)=delete;
    std::size_t hits() const;
    std::size_t size() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    friend class IalMobility;
};

class IalMobility {
public:
    /// Local temperature is explicit; coupled electrical callers still supply 300 K.
    explicit IalMobility(IalMobilityParameters parameters, bool electron,
                         IalScreeningMethod method=defaultIalScreeningMethod);
    IalScreeningMethod screeningMethod() const { return screeningMethod_; }
    IalMobility withScreeningMethod(IalScreeningMethod method) const {
        return IalMobility(params_,electron_,method);
    }
    IalMobilityResult evaluate(const IalMobilityState& state,IalScreeningCache* cache=nullptr,
                              IalMobilityPreparationCache* preparation=nullptr,bool generated=false) const;
    IalMobilityDifferential evaluateWithDerivatives(const IalMobilityState& state,
                                                  IalScreeningCache* cache=nullptr,
                                                  IalMobilityPreparationCache* preparation=nullptr,bool generated=false) const;
    static IalMobilityParameters siliconDefaults(bool electron);
    /// Immutable model data for evidence fingerprints; never carrier state.
    const IalMobilityParameters& parameters() const { return params_; }
    /// Closest cubic plane family, including sign/permutation symmetry.
    /// Input is a nonzero normal already transformed into crystal coordinates.
    static int orientationFamily(const std::array<Real, 3>& crystalNormal);
private:
    IalMobilityResult evaluatePrepared(const IalMobilityState& state,Real minimum) const;
    IalMobilityParameters params_;
    bool electron_;
    IalScreeningMethod screeningMethod_;
    Real pMin_;
};

} // namespace vela
