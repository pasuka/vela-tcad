#pragma once

#include "vela/core/Types.h"
#include <memory>
#include <string>

namespace vela {

struct ImpactIonizationModelConfig {
    std::string model = "none";
    std::string drivingForce = "electric_field";
    std::string generation = "carrier_density";
    Real electronA = 7.03e7; ///< Selberherr electron prefactor [1/m]
    Real electronB = 1.231e8; ///< Selberherr electron critical field [V/m]
    Real holeA = 1.582e8; ///< Selberherr hole prefactor [1/m]
    Real holeB = 2.036e8; ///< Selberherr hole critical field [V/m]
    Real carrierVelocity = 1.0e5; ///< Effective saturated carrier speed [m/s]
    Real electronALow = 7.03e7; ///< Van Overstraeten electron low-field prefactor [1/m]
    Real electronAHigh = 7.03e7; ///< Van Overstraeten electron high-field prefactor [1/m]
    Real electronBLow = 1.231e8; ///< Van Overstraeten electron low-field critical field [V/m]
    Real electronBHigh = 1.231e8; ///< Van Overstraeten electron high-field critical field [V/m]
    Real holeALow = 1.582e8; ///< Van Overstraeten hole low-field prefactor [1/m]
    Real holeAHigh = 6.71e7; ///< Van Overstraeten hole high-field prefactor [1/m]
    Real holeBLow = 2.036e8; ///< Van Overstraeten hole low-field critical field [V/m]
    Real holeBHigh = 1.693e8; ///< Van Overstraeten hole high-field critical field [V/m]
    Real switchField = 4.0e7; ///< Van Overstraeten low/high switch field [V/m]
    Real phononEnergy = 0.063; ///< Optical phonon energy for temperature factor [eV]
    Real referenceTemperature_K = 300.0; ///< Reference temperature for gamma factor [K]
    Real temperature_K = 300.0; ///< Lattice temperature for gamma factor [K]
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

class VanOverstraetenImpactIonization final : public ImpactIonizationModel {
public:
    explicit VanOverstraetenImpactIonization(ImpactIonizationModelConfig config = {});
    Real electronCoefficient(Real electricField) const override;
    Real holeCoefficient(Real electricField) const override;
    Real generationRate(Real electricField, Real n, Real p) const override;

private:
    static Real coefficient(Real electricField,
                            Real switchField,
                            Real lowPrefactor,
                            Real highPrefactor,
                            Real lowCriticalField,
                            Real highCriticalField,
                            Real gamma);
    Real gamma() const;

    ImpactIonizationModelConfig config_;
};

ImpactIonizationModelConfig impactIonizationModelConfig(std::string modelName);
std::unique_ptr<ImpactIonizationModel> makeImpactIonizationModel(
    const ImpactIonizationModelConfig& config);

} // namespace vela
