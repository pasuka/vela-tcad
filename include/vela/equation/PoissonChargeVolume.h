#pragma once

namespace vela {

/// Control-volume contract for semiconductor mobile and ionized-dopant charge
/// in the Poisson equation.  The legacy global policy uses the full geometric
/// node volume, including adjacent dielectric cells at a conforming interface.
enum class PoissonChargeVolumePolicy {
    Global,
    MaterialLocalBarycentric,
    MaterialLocalMixedVoronoi,
    MaterialLocalCircumcentric,
    MaterialLocalBox,
};

inline const char* poissonChargeVolumePolicyName(
    PoissonChargeVolumePolicy policy)
{
    switch (policy) {
        case PoissonChargeVolumePolicy::Global:
            return "global";
        case PoissonChargeVolumePolicy::MaterialLocalBarycentric:
            return "material_local";
    case PoissonChargeVolumePolicy::MaterialLocalMixedVoronoi:
        return "material_local_mixed_voronoi";
    case PoissonChargeVolumePolicy::MaterialLocalBox:
        return "material_local_box";
    case PoissonChargeVolumePolicy::MaterialLocalCircumcentric:
        return "material_local_circumcentric";
    }
    return "unknown";
}

} // namespace vela
