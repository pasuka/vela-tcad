#pragma once

#include "vela/core/Types.h"
#include "vela/core/UnitScaling.h"
#include <memory>
#include <string>

namespace vela {

struct BandgapNarrowingConfig {
    std::string model = "none"; ///< "none", "slotboom", or "old_slotboom"
    /// Numerical evaluation used on Boltzmann, equal-ni transport edges when
    /// BGN is disabled.  The legacy path subtracts two quasi-Fermi
    /// exponentials; the compensated path uses the algebraically equivalent
    /// log/expm1 VariableNi kernel.  The legacy default preserves existing
    /// solver results unless a diagnostic deck opts in explicitly.
    std::string equalNiFluxEvaluation = "legacy_factor_difference";
    Real referenceDoping = 1.0e23; ///< Slotboom reference concentration [m^-3]
    Real coefficient = 9.0e-3; ///< Slotboom narrowing coefficient [eV]
    Real smoothing = 0.5; ///< Dimensionless Slotboom smoothing term
    Real offset = 0.0; ///< Optional additive narrowing offset [eV]
    /// Apply the Sentaurus correction used when EffectiveIntrinsicDensity is
    /// combined with Fermi statistics.  This remains opt-in so existing Vela
    /// calibration decks retain their previous semantics.
    bool fermiStatisticsCorrection = false;
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

/// Sentaurus-compatible apparent-BGN correction for Fermi statistics [eV].
///
/// OldSlotboom parameters were extracted using Maxwell-Boltzmann statistics.
/// Sentaurus adds the difference between Fermi-Dirac and Boltzmann majority
/// carrier energies, evaluated at 300 K, unless NoFermi is requested.
Real fermiStatisticsBandgapCorrection(Real donors, Real acceptors,
                                      Real Nc, Real Nv,
                                      Real thermalVoltage);

BandgapNarrowingConfig bandgapNarrowingConfig(
    std::string modelName,
    UnitScalingConfig scaling = UnitScalingConfig{});
std::unique_ptr<BandgapNarrowing> makeBandgapNarrowingModel(
    const BandgapNarrowingConfig& config);

} // namespace vela
