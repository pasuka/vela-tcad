#pragma once

#include "vela/core/Types.h"
#include "vela/material/Material.h"
#include <memory>
#include <string>

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

struct MobilityModelConfig {
    std::string model = "constant";

    // 300 K silicon defaults converted from common Caughey-Thomas parameter
    // sets expressed in cm^2/(V s) and cm^-3.
    CaugheyThomasParameters electronCT{0.00522, 9.68e22, 0.68};
    CaugheyThomasParameters holeCT{0.00449, 2.23e23, 0.70};
    FieldMobilityParameters electronField{};
    FieldMobilityParameters holeField{};
};

class MobilityModel {
public:
    virtual ~MobilityModel() = default;

    virtual Real electronMobility(const Material& material,
                                  Real netDoping,
                                  Real n,
                                  Real p,
                                  Real electricField = 0.0) const = 0;

    virtual Real holeMobility(const Material& material,
                              Real netDoping,
                              Real n,
                              Real p,
                              Real electricField = 0.0) const = 0;
};

class ConstantMobility final : public MobilityModel {
public:
    Real electronMobility(const Material& material,
                          Real netDoping,
                          Real n,
                          Real p,
                          Real electricField = 0.0) const override;

    Real holeMobility(const Material& material,
                      Real netDoping,
                      Real n,
                      Real p,
                      Real electricField = 0.0) const override;
};

class DopingDependentMobility final : public MobilityModel {
public:
    explicit DopingDependentMobility(MobilityModelConfig config = {});

    Real electronMobility(const Material& material,
                          Real netDoping,
                          Real n,
                          Real p,
                          Real electricField = 0.0) const override;

    Real holeMobility(const Material& material,
                      Real netDoping,
                      Real n,
                      Real p,
                      Real electricField = 0.0) const override;

private:
    static Real caugheyThomas(Real muMax,
                              Real netDoping,
                              const CaugheyThomasParameters& params);
    static Real fieldLimit(Real lowFieldMobility,
                           Real electricField,
                           const FieldMobilityParameters& params);

    MobilityModelConfig config_;
};

MobilityModelConfig mobilityModelConfig(std::string modelName);
std::unique_ptr<MobilityModel> makeMobilityModel(const MobilityModelConfig& config);

} // namespace vela
