#pragma once

#include "vela/physics/MobilityModel.h"

namespace vela {

/** Diagnostic component mobilities are in cm^2/(V s). */
struct IALMobilityComponents {
    Real screening = 0.0;
    Real F = 0.0;
    Real G = 0.0;
    Real coulomb3D = 0.0;
    Real coulomb2D = 0.0;
    Real coulomb = 0.0;
    Real phonon3D = 0.0;
    Real phonon2D = 0.0;
    Real phonon = 0.0;
    Real roughness = 0.0;
    Real interfaceDamping = 0.0;
    Real coulombDamping = 0.0;
    Real bulkWeight = 1.0;
    Real lowField = 0.0;
};

/** Isotropic IALMob (T-2022.03 equations 290--311), without thin-layer,
 * stress, auto-orientation, or low-temperature phonon corrections.
 * Field saturation is applied to the complete low-field IAL mobility.
 */
class IALMobility final : public MobilityModel {
public:
    explicit IALMobility(MobilityModelConfig config);
    Real electronMobility(const Material&, Real, Real n, Real p,
                          Real electricField = 0.0, Real surfaceNormalField = 0.0,
                          Real surfaceDistance = 0.0, MobilityDopants dopants = {}) const override;
    Real holeMobility(const Material&, Real, Real n, Real p,
                      Real electricField = 0.0, Real surfaceNormalField = 0.0,
                      Real surfaceDistance = 0.0, MobilityDopants dopants = {}) const override;
    IALMobilityComponents components(const Material&, CarrierType, MobilityDopants,
                                    Real n, Real p, Real surfaceNormalField,
                                    Real surfaceDistance) const;

private:
    Real evaluate(const Material&, CarrierType, MobilityDopants, Real n, Real p,
                  Real field, Real normal, Real distance) const;
    MobilityModelConfig config_;
    Real electronScreeningMinimum_;
    Real holeScreeningMinimum_;
};

void parseIALMobilityParameters(const nlohmann::json&, IALMobilityParameters&);

} // namespace vela
