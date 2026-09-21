#pragma once

#include "vela/core/Types.h"
#include "vela/core/UnitScaling.h"
#include "vela/physics/BandToBandTunnelingModel.h"
#include <nlohmann/json_fwd.hpp>
#include <string>
#include <vector>

namespace vela {

struct RecombinationRateDerivatives {
    Real dRateDn = 0.0;      ///< Partial derivative at fixed excess product [s^-1]
    Real dRateDp = 0.0;      ///< Partial derivative at fixed excess product [s^-1]
    Real dRateDExcess = 0.0; ///< Partial derivative wrt n*p-ni^2 [m^3/s]
};

struct GeneralizedRecombinationRateDerivatives {
    Real dRateDn = 0.0; ///< Partial derivative at fixed generalized excess product [s^-1]
    Real dRateDp = 0.0; ///< Partial derivative at fixed generalized excess product [s^-1]
    Real dRateDExcess = 0.0; ///< Partial derivative wrt generalized excess product [m^3/s]
    Real dRateDElectronDegeneracy = 0.0; ///< Partial derivative wrt gamma_n [m^-3 s^-1]
    Real dRateDHoleDegeneracy = 0.0; ///< Partial derivative wrt gamma_p [m^-3 s^-1]
};

struct RecombinationLinearization {
    Real diagonal = 0.0; ///< Coefficient multiplying the solved carrier [s^-1]
    Real rhs = 0.0;      ///< Source contribution moved to the RHS [m^-3 s^-1]
};

struct SRHLifetimeParameters {
    Real tauMin = 0.0;  ///< High-doping lifetime limit [s]
    Real tauMax = 1.0e-7; ///< Low-doping lifetime limit [s]
    Real referenceDoping = 1.0e22; ///< Reference concentration [m^-3]
    Real gamma = 1.0; ///< Doping roll-off exponent [-]
};

struct SRHDopingDependenceConfig {
    bool enabled = false;
    /// ``total_impurity`` uses Nd+Na; ``net_doping`` uses |Nd-Na|.
    std::string concentrationBasis = "total_impurity";
    SRHLifetimeParameters electron{};
    SRHLifetimeParameters hole{};
    bool temperatureDependence = false;
    Real temperature_K = 300.0;
    Real referenceTemperature_K = 300.0;
    Real electronTemperatureExponent = 0.0;
    Real holeTemperatureExponent = 0.0;
    /// Electron density used by SRH under density-gradient quantum correction:
    /// ``quantum`` preserves the legacy Vela behavior; ``sentaurus_default``
    /// evaluates the SRH denominator and generalized Fermi factors with the
    /// classical density while transport and Poisson retain the quantum density.
    std::string densityCoupling = "quantum";
};

struct AugerDensityDependenceConfig {
    bool enabled = false;
    Real electronEnhancement = 0.0;
    Real holeEnhancement = 0.0;
    Real electronReferenceDensity = 1.0e24; ///< Internal concentration units; SI by default.
    Real holeReferenceDensity = 1.0e24;
};

struct RecombinationModelConfig {
    std::vector<std::string> mechanisms = {"srh"};
    Real taun = 1.0e-5; ///< Electron SRH lifetime [s]
    Real taup = 3.0e-6; ///< Hole SRH lifetime [s]

    // Sentaurus 2018 silicon Auger defaults at 300 K [m^6/s], converted from
    // Cn = A + B + C in cm^6/s and Cp = A + B + C in cm^6/s.
    Real augerCn = 2.90e-43;
    Real augerCp = 1.028e-43;
    /// Excess-product model used by Auger under Fermi statistics:
    /// ``generalized_fermi`` preserves the historical Vela behavior;
    /// ``classical_np`` uses n*p-ni_eff^2, matching the conventional model.
    std::string augerExcessProduct = "generalized_fermi";
    bool augerWithGeneration = true; ///< Compatibility default; false clips Auger generation.
    AugerDensityDependenceConfig augerDensityDependence{};
    SRHDopingDependenceConfig srhDopingDependence{};
    BandToBandTunnelingConfig bandToBand{};
};

class RecombinationModel {
public:
    explicit RecombinationModel(RecombinationModelConfig config = {});

    bool srhEnabled() const { return srhEnabled_; }
    bool augerEnabled() const { return augerEnabled_; }
    bool usesClassicalAugerExcessProduct() const
    {
        return config_.augerExcessProduct == "classical_np";
    }
    bool bandToBandEnabled() const { return bandToBand_.enabled(); }
    const BandToBandTunnelingModel& bandToBand() const { return bandToBand_; }

    Real electronLifetime(Real dopingConcentration) const;
    Real holeLifetime(Real dopingConcentration) const;
    Real srhDopingConcentration(Real donors, Real acceptors) const;
    Real srhRate(Real n, Real p, Real ni,
                 Real dopingConcentration = 0.0) const;
    Real srhRateFromExcessProduct(Real excessProduct,
                                  Real n,
                                  Real p,
                                  Real ni,
                                  Real dopingConcentration = 0.0) const;
    Real srhRateGeneralizedFromExcessProduct(
        Real excessProduct,
        Real n,
        Real p,
        Real n1,
        Real p1,
        Real electronDegeneracy,
        Real holeDegeneracy,
        Real dopingConcentration = 0.0) const;
    Real augerRate(Real n, Real p, Real ni) const;
    Real augerRateFromExcessProduct(Real excessProduct,
                                    Real n,
                                    Real p) const;
    Real totalRate(Real n, Real p, Real ni,
                   Real dopingConcentration = 0.0) const;
    Real totalRateFromExcessProduct(Real excessProduct,
                                    Real n,
                                    Real p,
                                    Real ni,
                                    Real dopingConcentration = 0.0) const;
    Real totalRateGeneralizedFromExcessProduct(
        Real excessProduct,
        Real n,
        Real p,
        Real n1,
        Real p1,
        Real electronDegeneracy,
        Real holeDegeneracy,
        Real dopingConcentration = 0.0) const;
    RecombinationRateDerivatives totalRateDerivativesFromExcessProduct(
        Real excessProduct,
        Real n,
        Real p,
        Real ni,
        Real dopingConcentration = 0.0) const;
    GeneralizedRecombinationRateDerivatives
    srhRateGeneralizedDerivativesFromExcessProduct(
        Real excessProduct,
        Real n,
        Real p,
        Real n1,
        Real p1,
        Real electronDegeneracy,
        Real holeDegeneracy,
        Real dopingConcentration = 0.0) const;
    RecombinationRateDerivatives augerRateDerivativesFromExcessProduct(
        Real excessProduct,
        Real n,
        Real p) const;

    RecombinationLinearization electronLinearization(
        Real n, Real p, Real ni, Real dopingConcentration = 0.0) const;
    RecombinationLinearization holeLinearization(
        Real n, Real p, Real ni, Real dopingConcentration = 0.0) const;

private:
    Real srhDenominator(Real n, Real p, Real ni,
                        Real dopingConcentration) const;
    Real srhGeneralizedDenominator(
        Real n, Real p, Real n1, Real p1,
        Real electronDegeneracy, Real holeDegeneracy,
        Real dopingConcentration) const;

    RecombinationModelConfig config_;
    BandToBandTunnelingModel bandToBand_;
    bool srhEnabled_ = false;
    bool augerEnabled_ = false;
};

RecombinationModelConfig recombinationModelConfig(
    std::vector<std::string> mechanisms,
    Real taun = 1.0e-5,
    Real taup = 3.0e-6,
    SRHDopingDependenceConfig srhDopingDependence = {});

SRHDopingDependenceConfig srhDopingDependenceConfigFromJson(
    const nlohmann::json& value,
    UnitScalingConfig scaling = {});

AugerDensityDependenceConfig augerDensityDependenceConfigFromJson(
    const nlohmann::json& value, UnitScalingConfig scaling = {});

} // namespace vela
