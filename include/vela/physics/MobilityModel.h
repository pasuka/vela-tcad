#pragma once

#include "vela/physics/GradQfContactField.h"
#include <optional>

#include "vela/core/UnitScaling.h"
#include "vela/core/Types.h"
#include "vela/material/Material.h"
#include "vela/physics/IALMobilityParameters.h"
#include <nlohmann/json_fwd.hpp>
#include <memory>
#include <limits>
#include <string>
#include <vector>
#include <array>

namespace vela {

struct IalTransportOptions;
struct IalTransportGeometry;
struct IalTransportState;

enum class CarrierType {
    Electron,
    Hole,
};

/** Separate ionized dopants in active internal concentration units.
 * Required by IALMob; a single net or total concentration loses compensation.
 */
struct MobilityDopants {
    Real donors = std::numeric_limits<Real>::quiet_NaN();
    Real acceptors = std::numeric_limits<Real>::quiet_NaN();
};
enum class PhuMobDonorSpecies {
    Arsenic,
    Phosphorus,
};

/** Carrier-specific parameters from the Sentaurus T-2022.03 PhuMob model. */
template<class Scalar>
struct BasicPhuMobCarrierParameters {
    Scalar muMax = 0; ///< Lattice mobility at 300 K [active mobility units].
    Scalar muMin = 0; ///< Carrier-scattering mobility floor [active mobility units].
    Scalar theta = 0; ///< Lattice temperature exponent [-].
    Scalar nRef = 1;  ///< Reference concentration [active concentration units].
    Scalar alpha = 1; ///< Impurity-scattering exponent [-].
};

/**
 * Sentaurus T-2022.03 Silicon Philips unified mobility parameters.
 *
 * The defaults are expressed in the legacy SI unit system.  The mobility
 * configuration parser converts them to the active unit system together with
 * the other mobility parameters.
 */
template<class Scalar>
struct BasicPhuMobParameters {
    BasicPhuMobCarrierParameters<Scalar> electronArsenic{
        Scalar(0.14170L), Scalar(0.00522L), Scalar(2.285L), Scalar(9.68e22L), Scalar(0.68L)};
    BasicPhuMobCarrierParameters<Scalar> electronPhosphorus{
        Scalar(0.14140L), Scalar(0.00685L), Scalar(2.285L), Scalar(9.20e22L), Scalar(0.711L)};
    BasicPhuMobCarrierParameters<Scalar> holeBoron{
        Scalar(0.04705L), Scalar(0.00449L), Scalar(2.247L), Scalar(2.23e23L), Scalar(0.719L)};

    Scalar donorClusterReference = Scalar(4.0e26L);   ///< nref_D [active concentration units].
    Scalar acceptorClusterReference = Scalar(7.2e26L); ///< nref_A [active concentration units].
    Scalar donorClusterCoefficient = Scalar(0.21L);   ///< cref_D [-].
    Scalar acceptorClusterCoefficient = Scalar(0.5L); ///< cref_A [-].
    Scalar electronMassRatio = 1;          ///< m*_e / m0 [-].
    Scalar holeMassRatio = Scalar(1.258L);            ///< m*_h / m0 [-].
    Scalar conwellWeisskopfFactor = Scalar(2.459L);   ///< f_CW [-].
    Scalar brooksHerringFactor = Scalar(3.828L);      ///< f_BH [-].
    Scalar electronHoleScatteringFactor = 1; ///< f_e [-].
    Scalar holeElectronScatteringFactor = 1; ///< f_h [-].

    // Built-in coefficients for Sentaurus Device User Guide equation 282.
    Scalar gA = Scalar(0.89233L);
    Scalar gB = Scalar(0.41372L);
    Scalar gC = Scalar(0.005978L);
    Scalar gAlpha = Scalar(0.28227L);
    Scalar gAlphaPrime = Scalar(0.72169L);
    Scalar gBeta = Scalar(0.19778L);
    Scalar gGamma = Scalar(1.80618L);

    PhuMobDonorSpecies donorSpecies = PhuMobDonorSpecies::Arsenic;
    Scalar internalConcentrationToCm3 = Scalar(1.0e-6L);
    Scalar internalMobilityToCm2PerVS = Scalar(1.0e4L);
};

template<class Scalar>
struct BasicPhuMobScalarState {
    Scalar donors = 0;
    Scalar acceptors = 0;
    Scalar electrons = 0;
    Scalar holes = 0;
    Scalar temperature_K = 300;
};

/** Component-resolved result of one scalar PhuMob carrier evaluation. */
template<class Scalar>
struct BasicPhuMobScalarResult {
    Scalar mobility = 0;
    Scalar latticeMobility = 0;
    Scalar scatteringMobility = 0;
    Scalar scatteringConcentration = 0;          ///< N_i,sc [active concentration units].
    Scalar effectiveScatteringConcentration = 0; ///< N_i,sc,eff [active concentration units].
    Scalar screeningParameter = 0;               ///< P_i [-].
    Scalar screeningF = 0;                       ///< F(P_i) [-].
    Scalar screeningG = 0;                       ///< G(P_i) after the documented clamp [-].
    Scalar screeningGDerivative = 0;             ///< dG/dP; zero on the floor branch.
};

/** Carrier chain terms at fixed impurities and temperature, in active mobility units. */
template<class Scalar>
struct BasicPhuMobLogDensityDerivatives {
    Scalar electrons = 0; ///< d(mu)/d(log n), retaining small n contributions.
    Scalar holes = 0;     ///< d(mu)/d(log p), retaining small p contributions.
};

using PhuMobCarrierParameters = BasicPhuMobCarrierParameters<Real>;
using PhuMobParameters = BasicPhuMobParameters<Real>;
using PhuMobScalarState = BasicPhuMobScalarState<Real>;
using PhuMobScalarResult = BasicPhuMobScalarResult<Real>;
using PhuMobLogDensityDerivatives = BasicPhuMobLogDensityDerivatives<Real>;

using ExtendedPhuMobParameters = BasicPhuMobParameters<long double>;
using ExtendedPhuMobScalarState = BasicPhuMobScalarState<long double>;
using ExtendedPhuMobScalarResult = BasicPhuMobScalarResult<long double>;
using ExtendedPhuMobLogDensityDerivatives = BasicPhuMobLogDensityDerivatives<long double>;

struct CaugheyThomasParameters {
    Real muMin = 0.0; ///< Low-field mobility floor [m^2/V/s]
    Real nRef  = 1.0; ///< Reference doping concentration [m^-3]
    Real alpha = 1.0; ///< Empirical roll-off exponent [-]
};

struct MasettiParameters {
    Real muConst = 0.0; ///< ConstantMobility reference mobility [m^2/V/s]
    Real muMin1 = 0.0;  ///< Low-doping exponential floor term [m^2/V/s]
    Real muMin2 = 0.0;  ///< Reference floor in the concentration roll-off term [m^2/V/s]
    Real mu1 = 0.0;     ///< High-doping correction term [m^2/V/s]
    Real pc = 0.0;      ///< Exponential activation concentration [m^-3]
    Real cr = 1.0;      ///< Reference doping concentration [m^-3]
    Real cs = 1.0;      ///< High-doping correction concentration [m^-3]
    Real alpha = 1.0;   ///< Concentration roll-off exponent [-]
    Real beta = 1.0;    ///< High-doping correction exponent [-]
};

struct FieldMobilityParameters {
    Real saturationVelocity = 1.0e5; ///< Saturation velocity [m/s]
    Real beta = 2.0;                 ///< High-field roll-off exponent [-]
};

struct SurfaceMobilityParameters {
    std::string discretization = "legacy_cell_centroid";
    Real thetaElectron = 0.0; ///< Electron vertical-field degradation coefficient [m/V]
    Real thetaHole = 0.0;     ///< Hole vertical-field degradation coefficient [m/V]
    Real beta = 1.0;          ///< Vertical-field roll-off exponent [-]
    Real referenceField = 0.0; ///< Field offset before degradation starts [V/m]
    Real minFactor = 0.0;     ///< Optional lower clamp for mu_surface / mu_bulk [-]
    Real maxFactor = 1.0;     ///< Optional upper clamp for mu_surface / mu_bulk [-]
    std::string surfaceRegion; ///< Optional semiconductor region where degradation is active.
    std::vector<std::string> surfaceInterface; ///< Optional two-region interface selector.
    Real coordinateFieldFactor = 1.0; ///< Internal coordinate-gradient to internal electric-field units.
    /// Per-cell nearest-interface field and distance, populated by assemblers
    /// for interface-distance-aware models such as Enhanced Lombardi.
    std::vector<Real> cellNormalFields;
    std::vector<Real> cellDistances;
    std::vector<Real> cellNormalX;
    std::vector<Real> cellNormalY;
    // True vertex distances; the element profile differentiates their linear
    // interpolant without normalizing its gradient.
    std::vector<Real> nodeDistances;
};

/** Sentaurus T-2022.03 Enhanced Lombardi coefficients, stored in SI units. */
struct LombardiParameters {
    Real B = 0.0;       ///< Acoustic-phonon coefficient [m/s].
    Real C = 0.0;       ///< Acoustic-phonon coefficient [m^(5/3)/(V^(2/3)s)].
    Real N0 = 1.0e6;   ///< Reference concentration [m^-3].
    Real N2 = 1.0e6;   ///< Concentration offset [m^-3].
    Real lambda = 0.0;
    Real k = 1.0;
    Real delta = 0.0;  ///< Surface-roughness coefficient in the SI-reduced formula.
    Real A = 2.0;
    Real alpha = 0.0;  ///< Carrier-dependent exponent coefficient [m^3].
    Real aOther = 0.0;
    Real N1 = 1.0e6;   ///< Carrier exponent concentration offset [m^-3].
    Real nu = 1.0;
    Real eta = 0.0;    ///< Cubic-field coefficient [V^2/(m s)].
    Real criticalLength = 1.0e-8; ///< Interface damping length [m].
    Real acousticFactor = 1.0;
    Real roughnessFactor = 1.0;
};

struct MobilityModelConfig {
    std::shared_ptr<const IalTransportOptions> ialmob;
    std::shared_ptr<const IalTransportGeometry> ialmobGeometry;
    std::shared_ptr<const IalTransportState> ialmobState;
    /// ``constant_field`` uses the material low-field mobility and applies
    /// only the configured high-field saturation limiter.
    std::string model = "constant";
    // Automatic sampling evaluates nodal Masetti before averaging.
    // Explicit element-box policies also select the cell and edge quadrature.
    std::string edgeAveraging = "auto";
    bool prepareBoxMasettiLombardi = false;
    /// Explicit element-box material override; skip PhuMob, surface and HFS.
    std::vector<std::string> constantMaterials;
    std::string highFieldDrivingForce = "electric_field";
    /// Spatial discretization for a quasi-Fermi-gradient high-field drive.
    /// ``edge_projection`` preserves the historical edge-aligned difference;
    /// ``transport_cell_vector`` recovers the full vector gradient from
    /// adjacent transport cells before evaluating the edge mobility.
    std::string highFieldGradientDiscretization = "edge_projection";
    /// Experimental Sentaurus-compatible HFS support: retain the configured
    /// quasi-Fermi-gradient drive in the device interior, but use the P1
    /// electrostatic-field magnitude in transport cells touching a contact.
    /// This is deliberately independent of the impact-ionization contact
    /// fallback and remains disabled unless a template opts in explicitly.
    bool contactElectricFieldFallback = false;
    // JSON resolves this explicitly. Direct C++ callers can retain their kernel default.
    std::optional<GradQfContactField> gradQfContactField;

    bool useElectricFieldAtGradQfContacts(bool kernelDefault = false) const {
        return highFieldDrivingForce == "quasi_fermi_gradient" &&
            (gradQfContactField ? *gradQfContactField == GradQfContactField::ElectricField
                               : kernelDefault || contactElectricFieldFallback);
    }
    std::string contactElectricFieldFallbackScope = "contact_node_cell";
    std::string contactElectricFieldFallbackMode = "cell_gradient_magnitude";
    /// Spatial support used by the carrier continuity operator and terminal
    /// current integration. ``scharfetter_gummel_edge`` preserves the
    /// historical edge flux. ``element_qf_gradient`` uses a Tri3 P1
    /// quasi-Fermi-gradient current reconstructed inside each transport cell.
    std::string carrierCurrentDiscretization = "scharfetter_gummel_edge";
    std::string dopingConcentrationBasis = "net_doping";
    bool jacobianFieldDerivatives = true;
    // Omit only transport feedback; preserve values and source AD caches.
    bool ialmobTransportFeedbackDerivatives = true;
    // Geometry only: two independent gradient columns, in active field/V.
    struct HighFieldCell {bool contact=false; std::array<Real,4> stencil{};};
    std::vector<HighFieldCell> highFieldCells;

    // 300 K silicon defaults converted from common Caughey-Thomas parameter
    // sets expressed in cm^2/(V s) and cm^-3.
    CaugheyThomasParameters electronCT{0.00522, 9.68e22, 0.68};
    CaugheyThomasParameters holeCT{0.00449, 2.23e23, 0.70};
    // Sentaurus 2018 Silicon DopingDependence Formula 1 (Masetti) defaults,
    // converted from sdevice -P:Silicon output.
    MasettiParameters electronMasetti{
        0.14170, 0.00522, 0.00522, 0.00434, 0.0, 9.68e22, 3.43e26, 0.68, 2.0};
    MasettiParameters holeMasetti{
        0.04705, 0.00449, 0.0, 0.00290, 9.23e22, 2.23e23, 6.10e26, 0.719, 2.0};
    PhuMobParameters phuMob{};
    // Sentaurus 2018 Silicon HighFieldDependence defaults from sdevice -P:
    // vsat0 = 1.07e7, 8.37e6 cm/s and beta0 = 1.109, 1.213 at 300 K.
    FieldMobilityParameters electronField{1.07e5, 1.109};
    FieldMobilityParameters holeField{8.37e4, 1.213};
    SurfaceMobilityParameters surface{};
    // Legacy SLOT endpoint model; separate from the element IAL transport above.
    IALMobilityParameters slotIalMobility{};
    // Sentaurus T-2022.03 Silicon EnormalDependence defaults.  C is
    // converted from cm^(5/3)/(V^(2/3)s), delta from cm^2/(V s), and eta
    // from V^2/(cm s).
    LombardiParameters electronLombardi{
        4.7500e5, 5.8000e2 * 4.641588833612778e-4,
        1.0e6, 1.0e6, 0.125, 1.0, 5.8200e10, 2.0,
        0.0, 0.0, 1.0e6, 1.0, 5.8200e32, 1.0e-8, 1.0, 1.0};
    LombardiParameters holeLombardi{
        9.9250e4, 2.9470e3 * 4.641588833612778e-4,
        1.0e6, 1.0e6, 0.0317, 1.0, 2.0546e10, 2.0,
        0.0, 0.0, 1.0e6, 1.0, 2.0546e32, 1.0e-8, 1.0, 1.0};
    Real internalFieldToVPerM = 1.0;
    Real internalConcentrationToM3 = 1.0;
    Real internalMobilityToM2PerVS = 1.0;
    Real internalLengthToM = 1.0;
};

class MobilityModel {
public:
    virtual ~MobilityModel() = default;

    virtual Real electronMobility(const Material& material,
                                  Real netDoping,
                                  Real n,
                                  Real p,
                                  Real electricField = 0.0,
                                  Real surfaceNormalField = 0.0,
                                  Real surfaceDistance = 0.0,
                                  MobilityDopants dopants = {}) const = 0;

    virtual Real holeMobility(const Material& material,
                              Real netDoping,
                              Real n,
                              Real p,
                              Real electricField = 0.0,
                              Real surfaceNormalField = 0.0,
                              Real surfaceDistance = 0.0,
                                  MobilityDopants dopants = {}) const = 0;
    /**
     * Mobility evaluation retaining separate ionized donor and acceptor
     * concentrations.  Models that only depend on net doping use the default
     * forwarding implementation; PhuMob overrides it because compensation and
     * carrier-carrier scattering are part of its state.
     */
    virtual Real electronMobilityWithIonizedImpurities(
        const Material& material,
        Real donors,
        Real acceptors,
        Real n,
        Real p,
        Real electricField = 0.0,
        Real surfaceNormalField = 0.0,
        Real surfaceDistance = 0.0) const;

    virtual Real holeMobilityWithIonizedImpurities(
        const Material& material,
        Real donors,
        Real acceptors,
        Real n,
        Real p,
        Real electricField = 0.0,
        Real surfaceNormalField = 0.0,
        Real surfaceDistance = 0.0) const;
};

class ConstantMobility final : public MobilityModel {
public:
    Real electronMobility(const Material& material,
                          Real netDoping,
                          Real n,
                          Real p,
                          Real electricField = 0.0,
                          Real surfaceNormalField = 0.0,
                          Real surfaceDistance = 0.0,
                                  MobilityDopants dopants = {}) const override;

    Real holeMobility(const Material& material,
                      Real netDoping,
                      Real n,
                      Real p,
                      Real electricField = 0.0,
                      Real surfaceNormalField = 0.0,
                      Real surfaceDistance = 0.0,
                                  MobilityDopants dopants = {}) const override;
};

class DopingDependentMobility final : public MobilityModel {
public:
    explicit DopingDependentMobility(MobilityModelConfig config = {});

    // One isothermal Masetti/Lombardi vertex sample with zero carrier inputs.
    // The enclosing assembler owns the material, doping, and model identity.
    struct PreparedMasettiLombardi {
        Real bulkMobility = 0.0;
        Real acousticCoefficient = 0.0;
        Real roughnessExponent = 0.0;
        bool validSurfaceParameters = false;
    };
    PreparedMasettiLombardi prepareMasettiLombardi(
        const Material& material, Real impurity, CarrierType carrier) const;
    Real evaluatePreparedMasettiLombardi(
        const PreparedMasettiLombardi& prepared, CarrierType carrier,
        Real electricField, Real surfaceNormalField, Real surfaceDistance) const;

    Real electronMobility(const Material& material,
                          Real netDoping,
                          Real n,
                          Real p,
                          Real electricField = 0.0,
                          Real surfaceNormalField = 0.0,
                          Real surfaceDistance = 0.0,
                                  MobilityDopants dopants = {}) const override;

    Real holeMobility(const Material& material,
                      Real netDoping,
                      Real n,
                      Real p,
                      Real electricField = 0.0,
                      Real surfaceNormalField = 0.0,
                      Real surfaceDistance = 0.0,
                                  MobilityDopants dopants = {}) const override;

    Real electronMobilityWithIonizedImpurities(
        const Material& material,
        Real donors,
        Real acceptors,
        Real n,
        Real p,
        Real electricField = 0.0,
        Real surfaceNormalField = 0.0,
        Real surfaceDistance = 0.0) const override;

    Real holeMobilityWithIonizedImpurities(
        const Material& material,
        Real donors,
        Real acceptors,
        Real n,
        Real p,
        Real electricField = 0.0,
        Real surfaceNormalField = 0.0,
        Real surfaceDistance = 0.0) const override;

private:
    static Real caugheyThomas(Real muMax,
                              Real netDoping,
                              const CaugheyThomasParameters& params);
    static Real masetti(Real netDoping,
                        const MasettiParameters& params);
    static Real fieldLimit(Real lowFieldMobility,
                           Real electricField,
                           const FieldMobilityParameters& params);
    static Real surfaceLimit(Real bulkMobility,
                             Real surfaceNormalField,
                             Real theta,
                             const SurfaceMobilityParameters& params);
    Real lombardiLimit(Real bulkMobility,
                       Real netDoping,
                       Real n,
                       Real p,
                       Real surfaceNormalField,
                       Real surfaceDistance,
                       Real latticeTemperature_K,
                       CarrierType carrier,
                       const LombardiParameters& params) const;
    Real lombardiMobilityFromField(
        Real bulkMobility, Real field, Real damping, Real acousticCoefficient,
        Real exponent, const LombardiParameters& params) const;

    enum class BulkLaw { Constant, CaugheyThomas, Masetti, PhuMob };
    enum class SurfaceLaw { None, Generic, Lombardi };

    MobilityModelConfig config_;
    BulkLaw bulkLaw_;
    SurfaceLaw surfaceLaw_;
    bool fieldMobility_;
};

MobilityModelConfig mobilityModelConfig(std::string modelName);
MobilityModelConfig mobilityModelConfigFromJson(
    const nlohmann::json& value,
    UnitScalingConfig scaling = {});
bool usesVertexSampledMasetti(const MobilityModelConfig& config);
bool isSurfaceMobilityModel(const MobilityModelConfig& config);
// Identifies only the legacy SLOT endpoint models, not native model="ialmob".
bool isIALMobilityModel(const MobilityModelConfig& config);
bool mobilityDependsOnCarrierDensity(const MobilityModelConfig& config);
bool isPhuMobModel(const MobilityModelConfig& config);
bool surfaceMobilityAppliesToRegionPair(const MobilityModelConfig& config,
                                        const std::string& regionName,
                                        const std::vector<std::string>& adjacentRegionNames);
std::unique_ptr<MobilityModel> makeMobilityModel(const MobilityModelConfig& config);

/**
 * Evaluate the T-2022.03 Philips unified low-field mobility equations
 * (User Guide equations 270-283) for one carrier and one local scalar state.
 * M6 wires this scalar kernel into DD assembly through the separate-ionized-
 * impurity mobility interface above.
 */
PhuMobScalarResult evaluatePhuMobScalar(
    CarrierType carrier,
    const PhuMobScalarState& state,
    const PhuMobParameters& params = {});

PhuMobLogDensityDerivatives evaluatePhuMobLogDensityDerivatives(
    CarrierType carrier,
    const PhuMobScalarState& state,
    const PhuMobParameters& params = {});

ExtendedPhuMobScalarResult evaluateExtendedPhuMobScalar(
    CarrierType carrier, const ExtendedPhuMobScalarState& state,
    const ExtendedPhuMobParameters& params);

ExtendedPhuMobLogDensityDerivatives evaluateExtendedPhuMobLogDensityDerivatives(
    CarrierType carrier, const ExtendedPhuMobScalarState& state,
    const ExtendedPhuMobParameters& params);

} // namespace vela
