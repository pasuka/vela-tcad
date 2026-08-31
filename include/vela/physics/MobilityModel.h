#pragma once

#include "vela/core/UnitScaling.h"
#include "vela/core/Types.h"
#include "vela/material/Material.h"
#include <nlohmann/json_fwd.hpp>
#include <memory>
#include <string>
#include <vector>

namespace vela {

enum class CarrierType {
    Electron,
    Hole,
};

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
};

/** Sentaurus O-2018.06 Enhanced Lombardi coefficients, stored in SI units. */
struct LombardiParameters {
    Real B = 0.0;       ///< Acoustic-phonon coefficient [m/s].
    Real C = 0.0;       ///< Acoustic-phonon coefficient [m^(5/3)/(V^(2/3)s)].
    Real N0 = 1.0e6;   ///< Reference concentration [m^-3].
    Real N2 = 1.0e6;   ///< Concentration offset [m^-3].
    Real lambda = 0.0;
    Real k = 1.0;
    Real delta = 0.0;  ///< Surface-roughness mobility coefficient [m^2/(V s)].
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
    /// ``constant_field`` uses the material low-field mobility and applies
    /// only the configured high-field saturation limiter.
    std::string model = "constant";
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
    std::string contactElectricFieldFallbackScope = "contact_node_cell";
    std::string contactElectricFieldFallbackMode = "cell_gradient_magnitude";
    /// Spatial support used by the carrier continuity operator and terminal
    /// current integration. ``scharfetter_gummel_edge`` preserves the
    /// historical edge flux. ``element_qf_gradient`` uses a Tri3 P1
    /// quasi-Fermi-gradient current reconstructed inside each transport cell.
    std::string carrierCurrentDiscretization = "scharfetter_gummel_edge";
    std::string dopingConcentrationBasis = "net_doping";
    bool jacobianFieldDerivatives = true;

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
    // Sentaurus 2018 Silicon HighFieldDependence defaults from sdevice -P:
    // vsat0 = 1.07e7, 8.37e6 cm/s and beta0 = 1.109, 1.213 at 300 K.
    FieldMobilityParameters electronField{1.07e5, 1.109};
    FieldMobilityParameters holeField{8.37e4, 1.213};
    SurfaceMobilityParameters surface{};
    // Sentaurus O-2018.06 Silicon EnormalDependence defaults.  C is
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
                                  Real surfaceDistance = 0.0) const = 0;

    virtual Real holeMobility(const Material& material,
                              Real netDoping,
                              Real n,
                              Real p,
                              Real electricField = 0.0,
                              Real surfaceNormalField = 0.0,
                              Real surfaceDistance = 0.0) const = 0;
};

class ConstantMobility final : public MobilityModel {
public:
    Real electronMobility(const Material& material,
                          Real netDoping,
                          Real n,
                          Real p,
                          Real electricField = 0.0,
                          Real surfaceNormalField = 0.0,
                          Real surfaceDistance = 0.0) const override;

    Real holeMobility(const Material& material,
                      Real netDoping,
                      Real n,
                      Real p,
                      Real electricField = 0.0,
                      Real surfaceNormalField = 0.0,
                      Real surfaceDistance = 0.0) const override;
};

class DopingDependentMobility final : public MobilityModel {
public:
    explicit DopingDependentMobility(MobilityModelConfig config = {});

    Real electronMobility(const Material& material,
                          Real netDoping,
                          Real n,
                          Real p,
                          Real electricField = 0.0,
                          Real surfaceNormalField = 0.0,
                          Real surfaceDistance = 0.0) const override;

    Real holeMobility(const Material& material,
                      Real netDoping,
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
                       CarrierType carrier,
                       const LombardiParameters& params) const;

    MobilityModelConfig config_;
};

MobilityModelConfig mobilityModelConfig(std::string modelName);
MobilityModelConfig mobilityModelConfigFromJson(
    const nlohmann::json& value,
    UnitScalingConfig scaling = {});
bool isSurfaceMobilityModel(const MobilityModelConfig& config);
bool surfaceMobilityAppliesToRegionPair(const MobilityModelConfig& config,
                                        const std::string& regionName,
                                        const std::vector<std::string>& adjacentRegionNames);
std::unique_ptr<MobilityModel> makeMobilityModel(const MobilityModelConfig& config);

} // namespace vela
