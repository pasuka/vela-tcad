#pragma once

#include "vela/core/Types.h"
#include <string>
#include <vector>

namespace vela {

struct RecombinationRateDerivatives {
    Real dRateDn = 0.0;      ///< Partial derivative at fixed excess product [s^-1]
    Real dRateDp = 0.0;      ///< Partial derivative at fixed excess product [s^-1]
    Real dRateDExcess = 0.0; ///< Partial derivative wrt n*p-ni^2 [m^3/s]
};

struct RecombinationLinearization {
    Real diagonal = 0.0; ///< Coefficient multiplying the solved carrier [s^-1]
    Real rhs = 0.0;      ///< Source contribution moved to the RHS [m^-3 s^-1]
};

struct RecombinationModelConfig {
    std::vector<std::string> mechanisms = {"srh"};
    Real taun = 1.0e-5; ///< Electron SRH lifetime [s]
    Real taup = 3.0e-6; ///< Hole SRH lifetime [s]

    // Sentaurus 2018 silicon Auger defaults at 300 K [m^6/s], converted from
    // Cn = A + B + C in cm^6/s and Cp = A + B + C in cm^6/s.
    Real augerCn = 2.90e-43;
    Real augerCp = 1.028e-43;
};

class RecombinationModel {
public:
    explicit RecombinationModel(RecombinationModelConfig config = {});

    bool srhEnabled() const { return srhEnabled_; }
    bool augerEnabled() const { return augerEnabled_; }

    Real srhRate(Real n, Real p, Real ni) const;
    Real srhRateFromExcessProduct(Real excessProduct,
                                  Real n,
                                  Real p,
                                  Real ni) const;
    Real augerRate(Real n, Real p, Real ni) const;
    Real augerRateFromExcessProduct(Real excessProduct,
                                    Real n,
                                    Real p) const;
    Real totalRate(Real n, Real p, Real ni) const;
    Real totalRateFromExcessProduct(Real excessProduct,
                                    Real n,
                                    Real p,
                                    Real ni) const;
    RecombinationRateDerivatives totalRateDerivativesFromExcessProduct(
        Real excessProduct,
        Real n,
        Real p,
        Real ni) const;

    RecombinationLinearization electronLinearization(Real n, Real p, Real ni) const;
    RecombinationLinearization holeLinearization(Real n, Real p, Real ni) const;

private:
    Real srhDenominator(Real n, Real p, Real ni) const;

    RecombinationModelConfig config_;
    bool srhEnabled_ = false;
    bool augerEnabled_ = false;
};

RecombinationModelConfig recombinationModelConfig(
    std::vector<std::string> mechanisms,
    Real taun = 1.0e-5,
    Real taup = 3.0e-6);

} // namespace vela
