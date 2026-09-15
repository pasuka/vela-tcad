#pragma once

#include "vela/equation/Tri3LocalForwardAD.h"
#include "vela/physics/IalMobility.h"
#include "vela/physics/MobilityModel.h"

namespace vela {

/// Frozen geometry in SI. The distance gradient is deliberately not normalized.
/// Orientation is selected upstream, independently of the electrical state.
struct IalElementGeometry {
    std::array<std::array<Real,2>,3> coordinates_m{};
    std::array<Real,3> interfaceDistance_m{};
    std::array<Real,3> vertexMeasure_m2{};
    bool partialBoundaryLayer = false;
    std::array<Real,2> boundaryTangent{}; ///< Unit tangent, if projection is enabled.
    bool touchesEffectiveElectrode = false;
};

struct IalElementVertexState {
    Real potential_V = 0., electronQf_V = 0., holeQf_V = 0.;
    Real donors_m3 = 0., acceptors_m3 = 0.;
    Real electrons_m3 = 0., holes_m3 = 0.;
    /// Positive derivatives with respect to psi-phin and phip-psi, respectively.
    /// The caller supplies derivatives from its actual carrier statistics,
    /// including Fermi-Dirac statistics when selected. No Boltzmann assumption.
    Real electronResponse_m3_per_V = 0., holeResponse_m3_per_V = 0.;
    Real temperature_K = 300.;
    /// Carrier temperature partials at fixed local potentials; may be negative.
    Real electronTemperatureResponse_m3_per_K = 0., holeTemperatureResponse_m3_per_K = 0.;
};

struct IalElementMobilityOptions {
    bool highField = true;
    Real referenceDensity_m3 = 1e18; ///< D4 RefDens = 1e12 cm^-3.
    FieldMobilityParameters electronField{1.07e5,1.109};
    FieldMobilityParameters holeField{8.37e4,1.213};
    bool temperatureDependentHighField = false;
    bool temperatureDerivatives = false;
    /// False requests values only and requires temperatureDerivatives=false.
    bool spatialDerivatives = true;
    /// Opt-in transport preparation policies; no physical parameter changes.
    bool reuseLocalPreparation = false;
    bool residualValuesOnly = false;
    /// Spatial and temperature passes have identical local physical inputs.
    bool reuseThermalLocalDifferentials = true;
    Real electronVelocityTemperatureExponent = .87, holeVelocityTemperatureExponent = .52;
    Real electronBetaTemperatureExponent = .66, holeBetaTemperatureExponent = .17;
    /// Non-owning optional preparation cache, valid for this evaluation only.
    IalScreeningCache* screeningCache=nullptr;
    IalMobilityPreparationCache* preparationCache=nullptr;
};

struct IalElementMobilityResult {
    /// Derivative ordering: (psi, phin, phip) at vertex 0, then 1, then 2.
    /// Values in m^2/(V s), derivatives per physical volt.
    detail::Tri3LocalForwardDual electronLowField, holeLowField;
    detail::Tri3LocalForwardDual electron, hole;
    /// Three additional T columns, one per vertex, when requested in options.
    std::array<Real,3> electronTemperatureDerivative{},holeTemperatureDerivative{};
    std::array<Real,3> electronLowTemperatureDerivative{},holeLowTemperatureDerivative{};
};

/// Element-vertex mobility and its coupled potential Jacobian, prior to
/// transport-flux assembly. Models are already selected by each vertex's
/// independent crystal orientation. This API does not enable a solver model.
IalElementMobilityResult evaluateIalElementMobility(
    const IalElementGeometry& geometry,
    const std::array<IalElementVertexState,3>& state,
    const std::array<const IalMobility*,3>& electronModels,
    const std::array<const IalMobility*,3>& holeModels,
    const IalElementMobilityOptions& options = {});

} // namespace vela
