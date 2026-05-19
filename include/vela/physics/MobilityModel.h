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
};

struct MobilityModelConfig {
    std::string model = "constant";

    // 300 K silicon defaults converted from common Caughey-Thomas parameter
    // sets expressed in cm^2/(V s) and cm^-3.
    CaugheyThomasParameters electronCT{0.00522, 9.68e22, 0.68};
    CaugheyThomasParameters holeCT{0.00449, 2.23e23, 0.70};
    FieldMobilityParameters electronField{};
    FieldMobilityParameters holeField{};
    SurfaceMobilityParameters surface{};
};

class MobilityModel {
public:
    virtual ~MobilityModel() = default;

    virtual Real electronMobility(const Material& material,
                                  Real netDoping,
                                  Real n,
                                  Real p,
                                  Real electricField = 0.0,
                                  Real surfaceNormalField = 0.0) const = 0;

    virtual Real holeMobility(const Material& material,
                              Real netDoping,
                              Real n,
                              Real p,
                              Real electricField = 0.0,
                              Real surfaceNormalField = 0.0) const = 0;
};

class ConstantMobility final : public MobilityModel {
public:
    Real electronMobility(const Material& material,
                          Real netDoping,
                          Real n,
                          Real p,
                          Real electricField = 0.0,
                          Real surfaceNormalField = 0.0) const override;

    Real holeMobility(const Material& material,
                      Real netDoping,
                      Real n,
                      Real p,
                      Real electricField = 0.0,
                      Real surfaceNormalField = 0.0) const override;
};

class DopingDependentMobility final : public MobilityModel {
public:
    explicit DopingDependentMobility(MobilityModelConfig config = {});

    Real electronMobility(const Material& material,
                          Real netDoping,
                          Real n,
                          Real p,
                          Real electricField = 0.0,
                          Real surfaceNormalField = 0.0) const override;

    Real holeMobility(const Material& material,
                      Real netDoping,
                      Real n,
                      Real p,
                      Real electricField = 0.0,
                      Real surfaceNormalField = 0.0) const override;

private:
    static Real caugheyThomas(Real muMax,
                              Real netDoping,
                              const CaugheyThomasParameters& params);
    static Real fieldLimit(Real lowFieldMobility,
                           Real electricField,
                           const FieldMobilityParameters& params);
    static Real surfaceLimit(Real bulkMobility,
                             Real surfaceNormalField,
                             Real theta,
                             const SurfaceMobilityParameters& params);

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
