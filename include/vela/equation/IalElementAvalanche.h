#pragma once

#include "vela/equation/IalElementMobility.h"

namespace vela {

struct ElementIonizationCurve {
    Real lowA_per_m, lowB_V_per_m, highA_per_m, highB_V_per_m;
    Real switchField_V_per_m = 4e7;
};

struct ArithmeticElementIonizationCurve {
    ArithmeticReal lowA_per_m, lowB_V_per_m, highA_per_m, highB_V_per_m, switchField_V_per_m;
};

struct ArithmeticIalIonizationCoefficients {
    ArithmeticElementIonizationCurve electron, hole;
    ArithmeticReal minimumField_V_per_m;
};

/// Research kernel for the identified 300 K AvalDensGradQF construction.
/// This API does not select a solver model or infer native mesh conventions.
struct IalElementAvalancheOptions {
    GradQfContactField contactField = GradQfContactField::ElectricField;
    GradQfContactField mobilityContactField = GradQfContactField::ElectricField;
    /// Caller supplies the geometry-selected vertex, independently of state.
    std::size_t mobilityVertex = 3; ///< Unset; evaluation rejects until supplied.
    /// Explicit source shares, distinct from transport/mobility box measures.
    std::array<Real,3> sourceMeasure_m2{};
    ElementIonizationCurve electron{7.03e7,1.231e8,7.03e7,1.231e8};
    ElementIonizationCurve hole{1.582e8,2.036e8,6.71e7,1.693e8};
    Real minimumField_V_per_m = 100.;
    bool highFieldMobility = true;
    Real populationThermalVoltage_V = 0.;
};

template<class Scalar>
struct BasicIalElementAvalancheCarrier {
    using Dual = detail::BasicTri3LocalForwardDual<Scalar>;
    Dual density_m3, mobility_m2_per_Vs;
    Dual currentField_V_per_m, ionizationField_V_per_m;
    Dual alpha_per_m, generation_m3_per_s;
    std::array<Dual,3> source_per_m_s{};
};

template<class Scalar>
struct BasicIalElementAvalancheResult {
    BasicIalElementAvalancheCarrier<Scalar> electron, hole;
};
using IalElementAvalancheCarrier = BasicIalElementAvalancheCarrier<Real>;
using IalElementAvalancheResult = BasicIalElementAvalancheResult<Real>;

using IalAvalanchePopulationVertex = IalPopulationVertex;

/// Evaluate the same source with prepared mobility and authoritative population.
/// The potential is in the solver coordinate, with a common voltage scale.
/// Mobility inputs may carry derivatives in the same [psi, relative n, relative p]
/// basis. A constant dual freezes mobility to isolate the source response.
BasicIalElementAvalancheResult<ArithmeticReal> evaluateIalPopulationAvalancheFromMobility(
    const BasicIalElementGeometry<ArithmeticReal>& geometry,
    const std::array<IalAvalanchePopulationVertex, 3>& state,
    ArithmeticReal potentialScale_V, ArithmeticReal thermalVoltage_V,
    const std::array<detail::BasicTri3LocalForwardDual<ArithmeticReal>, 2>& mobility,
    const IalElementAvalancheOptions& options,
    const std::array<long double, 3>* arithmeticSourceMeasures = nullptr,
    const ArithmeticIalIonizationCoefficients* arithmeticCoefficients = nullptr);

/// All nine (psi, phin, phip) columns include supplied carrier-statistics
/// responses, local IALMob/HFS, QF/normal fields, and cosine alpha mixing.
/// Contact elements retain the full QF current gradient, including at interface
/// overlaps; non-contact interface elements use its tangential component.
/// The contact alpha drive uses the electric field, separately from the current.
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
