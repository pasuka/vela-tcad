#pragma once

#include "vela/core/Types.h"
#include <array>
#include <map>
#include <memory>

namespace vela {

/// Independent low-field IALMob kernel, T-2022.03 equations 290--311.
/// Parameters use the manual's cm-based units; state/result boundaries use SI.
/// FullPhuMob, PhononCombination=1, no stress/thin-layer corrections.
/// Impurity weighting factors/exponents in the phonon/roughness terms stay at
/// their manual defaults of one; ClusteringEverywhere is disabled. Temperature
/// dependencies follow equations 293, 295--299, 305--309. T >= 50 K only;
/// the optional low-temperature phonon exponent correction is not enabled.
/// This kernel does not enable IALMob in the coupled solver configuration.
struct IalMobilityParameters {
    Real muMax = 1417., muMin = 52.2, theta = 2.285;
    Real alpha = .68, nRef = 9.68e16;
    Real mass = 1., otherMass = 1.258;
    Real nRefD = 4e20, nRefA = 7.2e20, cRefD = .21, cRefA = .5;
    Real nDopRef = 1e18, nScRef = 1e18;
    Real S = .3042, transitionP = 4.;
    Real B = 9e5, C = 4400., lambda = .057, k = 1.;
    Real delta = 3.97e13, eta = 1e50, lambdaSr = .057;
    Real A = 2., alphaSr = 0., nu = 0., N1 = 1., N2 = 1.;
    Real lCrit = 1e3, lCritC = 1e3;
    Real d1Inv = 135., d2Inv = 40., nu0Inv = 1.5, nu1Inv = 2., nu2Inv = .5;
    Real alpha1Inv = 0., alpha2Inv = 0., alpha1Acc = 0., alpha2Acc = 0.;
    Real d1Acc = 135., d2Acc = 40., nu0Acc = 1.5, nu1Acc = 2., nu2Acc = .5;
};

struct IalMobilityState {
    Real donors_m3 = 0., acceptors_m3 = 0.;
    Real electrons_m3 = 0., holes_m3 = 0.;
    Real normalField_V_per_m = 0., interfaceDistance_m = 0.;
    Real temperature_K = 300.;
};

struct IalMobilityResult {
    Real mobility_m2_per_Vs;
    Real coulomb3d_m2_per_Vs, coulomb2d_m2_per_Vs, coulomb_m2_per_Vs;
    Real phonon_m2_per_Vs, roughness_m2_per_Vs;
};

/// Derivatives of the final low-field mobility with respect to SI state fields
/// in this order: Nd, Na, n, p, |Enormal|, interface distance. At a zero input
/// with a fractional-power cusp the selected branch slope is zero, not a
/// claim of a classical derivative. Geometry and orientation are fixed.
struct IalMobilityDifferential {
    IalMobilityResult result;
    std::array<Real, 6> derivative_SI{};
    /// Partial at fixed Nd, Na, n, p, field and distance (m2 / (V s K)).
    Real temperatureDerivative_m2_per_Vs_K = 0.;
};

/// Scoped to one transport-state preparation. The screening minimizer depends
/// only on mass and temperature; exact keys never reuse another temperature.
class IalScreeningCache {
public:
    Real minimum(Real mass,Real temperature);
    std::size_t size() const {return roots_.size();}
private:
    std::map<std::pair<Real,Real>,Real> roots_;
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
    explicit IalMobility(IalMobilityParameters parameters, bool electron);
    IalMobilityResult evaluate(const IalMobilityState& state,IalScreeningCache* cache=nullptr,
                              IalMobilityPreparationCache* preparation=nullptr) const;
    IalMobilityDifferential evaluateWithDerivatives(const IalMobilityState& state,
                                                  IalScreeningCache* cache=nullptr,
                                                  IalMobilityPreparationCache* preparation=nullptr) const;
    static IalMobilityParameters siliconDefaults(bool electron);
    /// Closest cubic plane family, including sign/permutation symmetry.
    /// Input is a nonzero normal already transformed into crystal coordinates.
    static int orientationFamily(const std::array<Real, 3>& crystalNormal);
private:
    IalMobilityResult evaluatePrepared(const IalMobilityState& state,Real minimum) const;
    IalMobilityParameters params_;
    bool electron_;
    Real pMin_;
};

} // namespace vela
