#pragma once

#include "vela/core/Types.h"
#include <string>
#include <vector>

namespace vela {

struct RecombinationLinearization {
    Real diagonal = 0.0; ///< Coefficient multiplying the solved carrier [s^-1]
    Real rhs = 0.0;      ///< Source contribution moved to the RHS [m^-3 s^-1]
};

struct RecombinationModelConfig {
    std::vector<std::string> mechanisms = {"srh"};
    Real taun = 1.0e-7; ///< Electron SRH lifetime [s]
    Real taup = 1.0e-7; ///< Hole SRH lifetime [s]

    // Silicon Auger defaults [m^6/s].  Values are intentionally modest so the
    // default examples remain stable while high-injection tests see growth.
    Real augerCn = 2.8e-43;
    Real augerCp = 9.9e-44;
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
    Real taun = 1.0e-7,
    Real taup = 1.0e-7);

} // namespace vela
