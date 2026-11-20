#pragma once
#include <optional>
#include <stdexcept>

#include "vela/equation/Tri3LocalForwardAD.h"
#include "vela/physics/IalMobility.h"
#include "vela/physics/MobilityModel.h"

namespace vela {

namespace detail {
// The population basis uses unit relative increments in n and p.
template<class Scalar>
void seedBoltzmannPopulationVertex(std::size_t i, std::type_identity_t<Scalar> thermalVoltage,
    BasicTri3LocalForwardDual<Scalar>& potential, BasicTri3LocalForwardDual<Scalar>& electronQf,
    BasicTri3LocalForwardDual<Scalar>& holeQf, BasicTri3LocalForwardDual<Scalar>& electrons,
    BasicTri3LocalForwardDual<Scalar>& holes) {
    electronQf.derivative.fill(0.);
    holeQf.derivative.fill(0.);
    electrons.derivative.fill(0.);
    holes.derivative.fill(0.);
    if (electrons.value > 0.) {
        electronQf.derivative[3*i]=potential.derivative[3*i];
        electronQf.derivative[3*i+1]=-thermalVoltage;
        electrons.derivative[3*i+1]=electrons.value;
    }
    if (holes.value > 0.) {
        holeQf.derivative[3*i]=potential.derivative[3*i];
        holeQf.derivative[3*i+2]=thermalVoltage;
        holes.derivative[3*i+2]=holes.value;
    }
}
} // namespace detail

/// Frozen geometry in SI. The distance gradient is deliberately not normalized.
/// Orientation is selected upstream, independently of the electrical state.
template<class Scalar>
struct BasicIalElementGeometry {
    std::array<std::array<Scalar,2>,3> coordinates_m{};
    std::array<Scalar,3> interfaceDistance_m{};
    std::array<Scalar,3> vertexMeasure_m2{};
    bool partialBoundaryLayer = false;
    std::array<Scalar,2> boundaryTangent{}; ///< Unit tangent, if projection is enabled.
    bool touchesEffectiveElectrode = false;
};
using IalElementGeometry = BasicIalElementGeometry<Real>;

struct IalPopulationVertex {
    ArithmeticReal potential;
    ArithmeticReal electrons_m3, holes_m3;
    ArithmeticReal electronIntrinsic_m3, holeIntrinsic_m3;
};

namespace detail {
using IalPopulationDual = BasicTri3LocalForwardDual<ArithmeticReal>;
struct IalPopulationFields {
    std::array<IalPopulationDual,3> psi, electronQf, holeQf, electrons, holes;
};
IalPopulationFields buildIalPopulationFields(
    const std::array<IalPopulationVertex,3>& state,
    ArithmeticReal potentialScale_V, ArithmeticReal thermalVoltage_V);

template<class Scalar>
class IalTri3Gradient {
public:
    explicit IalTri3Gradient(const std::array<std::array<Scalar,2>,3>& xy)
        : ax_(xy[1][0]-xy[0][0]), ay_(xy[1][1]-xy[0][1]),
          bx_(xy[2][0]-xy[0][0]), by_(xy[2][1]-xy[0][1]),
          det_(ax_*by_-ay_*bx_) {
        if (!std::isfinite(det_) || det_ == 0)
            throw std::invalid_argument("IALMob requires a nondegenerate triangle");
    }

    template<class D>
    std::array<D,2> operator()(const std::array<D,3>& v) const {
        const D a=v[1]-v[0], b=v[2]-v[0];
        return {(a*D(by_)-b*D(ay_))/D(det_), (b*D(ax_)-a*D(bx_))/D(det_)};
    }

private:
    Scalar ax_, ay_, bx_, by_, det_;
};
} // namespace detail

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

struct IalHighFieldDecimal {
    ArithmeticReal velocity_m_s;
    ArithmeticReal beta;
};

struct IalElementMobilityOptions {
    GradQfContactField contactField = GradQfContactField::ElectricField;
    bool highField = true;
    Real referenceDensity_m3 = 1e18; ///< D4 RefDens = 1e12 cm^-3.
    FieldMobilityParameters electronField{1.07e5,1.109};
    FieldMobilityParameters holeField{8.37e4,1.213};
    // The two entries are electron and hole, in SI, for the 300 K SLOT path.
    std::optional<std::array<IalHighFieldDecimal, 2>> decimalHighField;
    bool temperatureDependentHighField = false;
    bool temperatureDerivatives = false;
    /// False requests values only and requires temperatureDerivatives=false.
    bool spatialDerivatives = true;
    // Positive enables direct (psi, relative n, relative p) spatial seeds.
    Real populationThermalVoltage_V = 0.;
    /// Opt-in transport preparation policies; no physical parameter changes.
    bool reuseLocalPreparation = false;
    bool residualValuesOnly = false;
    /// Spatial and temperature passes have identical local physical inputs.
    bool reuseThermalLocalDifferentials = true;
    /// Opt-in reuse of identical high-field values/partials between the spatial
    /// and temperature passes of one element. Never reused across states.
    bool reuseThermalHighField = false;
    bool explicitHighFieldPartials = false;
    bool generatedLowFieldPartials = false;
    Real electronVelocityTemperatureExponent = .87, holeVelocityTemperatureExponent = .52;
    Real electronBetaTemperatureExponent = .66, holeBetaTemperatureExponent = .17;
    /// Non-owning optional preparation cache, valid for this evaluation only.
    IalScreeningCache* screeningCache=nullptr;
    IalMobilityPreparationCache* preparationCache=nullptr;
};

struct IalElementMobilityResult {
    /// Derivative ordering: (psi, phin, phip) at vertex 0, then 1, then 2.
    /// Values in m^2/(V s). When population seeds are selected, density
    /// columns are per unit relative increment and psi remains per volt.
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

/// Isothermal SLOT mobility in SI, differentiated in the population basis.
/// The supplied model at each vertex is the already-selected crystal family.
struct IalPopulationMobility {
    std::array<detail::IalPopulationDual,3> electron, hole;
    std::array<detail::IalPopulationDual,3> electronLowField, holeLowField;
};
IalPopulationMobility evaluateIalPopulationMobility(
    const BasicIalElementGeometry<ArithmeticReal>& geometry,
    const detail::IalPopulationFields& fields,
    const std::array<std::array<ArithmeticReal,2>,3>& doping_m3,
    const std::array<const BasicIalMobilityParameters<long double>*,3>& electrons,
    const std::array<const BasicIalMobilityParameters<long double>*,3>& holes,
    ArithmeticReal temperature_K, ArithmeticReal referenceDensity_m3,
    const std::optional<std::array<IalHighFieldDecimal,2>>& highField,
    GradQfContactField contactField, IalScreeningMethod method);

} // namespace vela
