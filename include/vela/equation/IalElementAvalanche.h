#pragma once

#include "vela/equation/IalElementMobility.h"

namespace vela {

struct ElementIonizationCurve {
    Real lowA_per_m, lowB_V_per_m, highA_per_m, highB_V_per_m;
    Real switchField_V_per_m = 4e7;
};

/// Research kernel for the identified 300 K AvalDensGradQF construction.
/// This API does not select a solver model or infer native mesh conventions.
struct IalElementAvalancheOptions {
    /// Caller supplies the geometry-selected vertex, independently of state.
    std::size_t mobilityVertex = 3; ///< Unset; evaluation rejects until supplied.
    /// Explicit source shares, distinct from transport/mobility box measures.
    std::array<Real,3> sourceMeasure_m2{};
    ElementIonizationCurve electron{7.03e7,1.231e8,7.03e7,1.231e8};
    ElementIonizationCurve hole{1.582e8,2.036e8,6.71e7,1.693e8};
    Real minimumField_V_per_m = 100.;
    bool highFieldMobility = true;
};

struct IalElementAvalancheCarrier {
    detail::Tri3LocalForwardDual density_m3, mobility_m2_per_Vs;
    detail::Tri3LocalForwardDual currentField_V_per_m, ionizationField_V_per_m;
    detail::Tri3LocalForwardDual alpha_per_m, generation_m3_per_s;
    std::array<detail::Tri3LocalForwardDual,3> source_per_m_s{};
};

struct IalElementAvalancheResult {
    IalElementAvalancheCarrier electron, hole;
};

/// All nine (psi, phin, phip) columns include supplied carrier-statistics
/// responses, local IALMob/HFS, QF/normal fields, and cosine alpha mixing.
/// Geometry, doping, orientation and temperature are fixed. Equal extrema use
/// the first local index on each selected min/max branch; zero norms and the
/// cutoff select zero slope. At such nonsmooth points this is a branch Jacobian,
/// not a claim that a unique classical derivative exists.
IalElementAvalancheResult evaluateIalElementAvalanche(
    const IalElementGeometry& geometry,
    const std::array<IalElementVertexState,3>& state,
    const std::array<const IalMobility*,3>& electronModels,
    const std::array<const IalMobility*,3>& holeModels,
    const IalElementAvalancheOptions& options);

} // namespace vela
