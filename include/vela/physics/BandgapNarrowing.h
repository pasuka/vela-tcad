#pragma once

#include "vela/core/Types.h"
#include <memory>
#include <string>

namespace vela {

struct BandgapNarrowingConfig {
    std::string model = "none"; ///< "none", "slotboom", or "old_slotboom"
    Real referenceDoping = 1.0e23; ///< Slotboom reference concentration [m^-3]
    Real coefficient = 9.0e-3; ///< Slotboom narrowing coefficient [eV]
    Real smoothing = 0.5; ///< Dimensionless Slotboom smoothing term
    Real offset = 0.0; ///< Optional additive narrowing offset [eV]
};

class BandgapNarrowing {
public:
    virtual ~BandgapNarrowing() = default;

    /// Return the effective bandgap narrowing DeltaEg [eV] at a node.
    virtual Real deltaEg(Real impurityConcentration, Real n, Real p) const;
};

class NoBandgapNarrowing final : public BandgapNarrowing {
public:
    Real deltaEg(Real impurityConcentration, Real n, Real p) const override;
};

class SlotboomBandgapNarrowing final : public BandgapNarrowing {
public:
    explicit SlotboomBandgapNarrowing(BandgapNarrowingConfig config = {});

    Real deltaEg(Real impurityConcentration, Real n, Real p) const override;

private:
    BandgapNarrowingConfig config_;
};

/// Return ni_eff = ni * exp(DeltaEg / (2 Vt)) for a narrowing in eV.
Real effectiveIntrinsicDensity(Real ni, Real thermalVoltage, Real deltaEg);

BandgapNarrowingConfig bandgapNarrowingConfig(std::string modelName);
std::unique_ptr<BandgapNarrowing> makeBandgapNarrowingModel(
    const BandgapNarrowingConfig& config);

} // namespace vela
