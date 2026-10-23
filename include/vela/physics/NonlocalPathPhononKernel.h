#pragma once
#include "vela/physics/NonlocalPathGeometry.h"

namespace vela::experimental {
struct TunnelingBandKnot {
    Real distance_m, valence_J, conduction_J, temperature_K;
};
struct TunnelingBandProfile {
    TunnelingStateIdentity identity;
    Real initialValenceGradient_J_per_m = 0;
    std::vector<TunnelingBandKnot> knots;
};
struct TunnelingReservoirs {
    Real electronFermi_J, holeFermi_J;
    // Explicit available energies beyond the turning endpoints, NOT extrema
    // of the barrier segment (which would incorrectly make both windows zero).
    Real valenceAvailableEnergy_J, conductionAvailableEnergy_J;
};
struct PhononTunnelingResult {
    Real crossover_m = 0, action = 0; // action = int(kv)+int(kc), dimensionless.
    Real valenceReciprocalIntegral_m2 = 0, conductionReciprocalIntegral_m2 = 0;
    Real valenceWindowFactor = 0, conductionWindowFactor = 0;
    Real logPrefactor = 0, logAbsoluteGeneration = 0;
    Real generation_per_m3_s = 0; // Positive generation, negative recombination.
    int sign = 0;
    bool underflow = false; // Log remains finite; this is NOT a physical zero.
    Real estimatedIntegrationRoundoff = 0;
};
// Only complete, state-matched N2 geometry is accepted. No source scatter.
TunnelingBandProfile tunnelingProfileFromGeometry(const TunnelingGeometryResult&,
    const TunnelingBandState&, const NonlocalPathBandToBandConfig&);
// Homogeneous Si, isothermal, P1 bands with a unique Keldysh crossover.
// Invalid/unresolved geometry and unsupported barrier pockets throw; never E2.
PhononTunnelingResult evaluatePhononTunneling(const TunnelingBandProfile&,
    const TunnelingReservoirs&, const NonlocalPathBandToBandConfig&);
} // namespace vela::experimental
