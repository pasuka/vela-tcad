#pragma once
#include "vela/physics/NonlocalPathSourceAssembly.h"

namespace vela::experimental {
// Explicit private numerical policy. The reference is a scalar energy-window
// integral [m^-2 s^-1], not a physical generation floor or an accuracy bound.
struct TunnelingCellIntegrationPolicy {
    Real minimumWidth_kT = .2, maximumWidth_kT = 1, variationTolerance = .1;
    Real variationReference_per_m2_s = 0;
    Index maximumDepth = 6;
};
struct TunnelingFractionIntegrals {
    Real valenceAction = 0, conductionAction = 0;
    Real valenceInverse_m2 = 0, conductionInverse_m2 = 0;
};
// P1 in dimensionless band fraction; intentionally distinct from the N3
// piecewise-linear barrier-energy kernel. Original segment scales survive clips.
TunnelingFractionIntegrals tunnelingFractionIntegrals(
    const std::vector<TunnelingBandKnot>&, Real energy_J,
    Real conductionMass_kg, Real valenceMass_kg, Real hbar_J_s);
std::vector<Real> tunnelingFirstCrossingWeights(const std::vector<Real>& energies_J,
    Real lower_J, Real upper_J);
TunnelingSourceAssembly assembleNonlocalPathCellSources(const TunnelingGeometryMesh&,
    const TunnelingBandState&, const NonlocalPathBandToBandConfig&,
    const TunnelingCellIntegrationPolicy&, const TunnelingGeometryOptions& = {});
} // namespace vela::experimental
