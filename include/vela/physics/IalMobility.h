#pragma once

#include "vela/core/Types.h"
#include <array>

namespace vela {

/// Independent low-field IALMob kernel, T-2022.03 equations 290--311.
/// Parameters use the manual's cm-based units; state/result boundaries use SI.
/// FullPhuMob, PhononCombination=1, no stress/thin-layer corrections.
/// Impurity weighting factors/exponents in the phonon/roughness terms stay at
/// their manual defaults of one; ClusteringEverywhere is disabled. Temperature
/// exponents theta/k are retained as metadata and have no effect at 300 K.
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
    Real d1Acc = 135., d2Acc = 40., nu0Acc = 1.5, nu1Acc = 2., nu2Acc = .5;
};

struct IalMobilityState {
    Real donors_m3 = 0., acceptors_m3 = 0.;
    Real electrons_m3 = 0., holes_m3 = 0.;
    Real normalField_V_per_m = 0., interfaceDistance_m = 0.;
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
};

class IalMobility {
public:
    /// The first implementation contract is explicitly isothermal at 300 K.
    explicit IalMobility(IalMobilityParameters parameters, bool electron);
    IalMobilityResult evaluate(const IalMobilityState& state) const;
    IalMobilityDifferential evaluateWithDerivatives(const IalMobilityState& state) const;
    static IalMobilityParameters siliconDefaults(bool electron);
    /// Closest cubic plane family, including sign/permutation symmetry.
    /// Input is a nonzero normal already transformed into crystal coordinates.
    static int orientationFamily(const std::array<Real, 3>& crystalNormal);
private:
    IalMobilityParameters params_;
    bool electron_;
    Real pMin_;
};

} // namespace vela
