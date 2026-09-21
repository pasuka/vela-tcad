#pragma once

namespace vela {

/// Control-volume contract for semiconductor mobile and ionized-dopant charge
/// in the Poisson equation.  The legacy global policy uses the full geometric
/// node volume, including adjacent dielectric cells at a conforming interface.
enum class PoissonChargeVolumePolicy {
    Global,
    MaterialLocalBarycentric,
};

inline const char* poissonChargeVolumePolicyName(
    PoissonChargeVolumePolicy policy)
{
    switch (policy) {
        case PoissonChargeVolumePolicy::Global:
            return "global";
        case PoissonChargeVolumePolicy::MaterialLocalBarycentric:
            return "material_local";
    }
    return "unknown";
}

} // namespace vela
