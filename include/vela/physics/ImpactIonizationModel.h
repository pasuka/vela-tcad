#pragma once

#include "vela/core/Types.h"
#include <memory>
#include <string>

namespace vela {

struct ImpactIonizationModelConfig {
    std::string model = "none";
    Real electronA = 7.03e7; ///< Selberherr electron prefactor [1/m]
    Real electronB = 1.231e8; ///< Selberherr electron critical field [V/m]
    Real holeA = 1.582e8; ///< Selberherr hole prefactor [1/m]
    Real holeB = 2.036e8; ///< Selberherr hole critical field [V/m]
    Real carrierVelocity = 1.0e5; ///< Effective saturated carrier speed [m/s]
};


class ImpactIonizationModel {
public:
    virtual ~ImpactIonizationModel() = default;
    virtual Real electronCoefficient(Real electricField) const = 0;
    virtual Real holeCoefficient(Real electricField) const = 0;
    virtual Real generationRate(Real electricField, Real n, Real p) const = 0;
};

class NoImpactIonization final : public ImpactIonizationModel {
public:
    Real electronCoefficient(Real electricField) const override;
    Real holeCoefficient(Real electricField) const override;
    Real generationRate(Real electricField, Real n, Real p) const override;
};

class SelberherrImpactIonization final : public ImpactIonizationModel {
public:
    explicit SelberherrImpactIonization(ImpactIonizationModelConfig config = {});
    Real electronCoefficient(Real electricField) const override;
    Real holeCoefficient(Real electricField) const override;
    Real generationRate(Real electricField, Real n, Real p) const override;

private:
    static Real coefficient(Real electricField, Real prefactor, Real criticalField);
    ImpactIonizationModelConfig config_;
};

ImpactIonizationModelConfig impactIonizationModelConfig(std::string modelName);
std::unique_ptr<ImpactIonizationModel> makeImpactIonizationModel(
    const ImpactIonizationModelConfig& config);

} // namespace vela
