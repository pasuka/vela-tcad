#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/CarrierStatistics.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <utility>

namespace vela {

Real BandgapNarrowing::deltaEg(Real, Real, Real) const
{
    return 0.0;
}

Real NoBandgapNarrowing::deltaEg(Real, Real, Real) const
{
    return 0.0;
}

SlotboomBandgapNarrowing::SlotboomBandgapNarrowing(BandgapNarrowingConfig config)
    : config_(std::move(config))
{
    if (config_.referenceDoping <= 0.0)
        throw std::invalid_argument(
            "SlotboomBandgapNarrowing: referenceDoping must be positive.");
    if (config_.coefficient < 0.0)
        throw std::invalid_argument(
            "SlotboomBandgapNarrowing: coefficient cannot be negative.");
    if (config_.smoothing < 0.0)
        throw std::invalid_argument(
            "SlotboomBandgapNarrowing: smoothing cannot be negative.");
}

Real SlotboomBandgapNarrowing::deltaEg(Real impurityConcentration, Real n, Real p) const
{
    const Real effectiveConcentration = std::max({std::abs(impurityConcentration), n, p});
    if (effectiveConcentration <= 0.0 || config_.coefficient == 0.0)
        return 0.0;

    const Real x = std::log(effectiveConcentration / config_.referenceDoping);
    const Real delta =
        config_.offset +
        config_.coefficient * (x + std::sqrt(x * x + config_.smoothing));
    return std::max(delta, 0.0);
}

Real effectiveIntrinsicDensity(Real ni, Real thermalVoltage, Real deltaEg)
{
    if (ni <= 0.0 || deltaEg <= 0.0)
        return ni;
    if (thermalVoltage <= 0.0)
        throw std::invalid_argument("effectiveIntrinsicDensity: thermalVoltage must be positive.");

    const Real exponent = deltaEg / (2.0 * thermalVoltage);
    const Real maxExponent = std::log(std::numeric_limits<Real>::max() / ni);
    if (exponent >= maxExponent)
        return std::numeric_limits<Real>::max();
    return ni * std::exp(exponent);
}

Real fermiStatisticsBandgapCorrection(Real donors, Real acceptors,
                                      Real Nc, Real Nv,
                                      Real thermalVoltage)
{
    if (!(thermalVoltage > 0.0) || !std::isfinite(thermalVoltage))
        throw std::invalid_argument(
            "fermiStatisticsBandgapCorrection: thermal voltage must be positive and finite.");
    if (!(Nc > 0.0) || !(Nv > 0.0) || !std::isfinite(Nc) || !std::isfinite(Nv))
        throw std::invalid_argument(
            "fermiStatisticsBandgapCorrection: density of states must be positive and finite.");

    auto speciesCorrection = [](Real concentration, Real densityOfStates) {
        if (!(concentration > 0.0) || !std::isfinite(concentration))
            return Real{0.0};
        const Real ratio = concentration / densityOfStates;
        if (!(ratio > 0.0) || !std::isfinite(ratio))
            return Real{0.0};
        return std::max<Real>(
            0.0, inverseFermiDiracHalf(ratio) - std::log(ratio));
    };

    return thermalVoltage * (
        speciesCorrection(donors, Nc) + speciesCorrection(acceptors, Nv));
}

BandgapNarrowingConfig bandgapNarrowingConfig(
    std::string modelName,
    UnitScalingConfig scaling)
{
    BandgapNarrowingConfig config;
    config.model = std::move(modelName);
    // Slotboom and OldSlotboom share the same Sentaurus reference
    // concentration, so both must express the compiled SI default in the
    // active internal concentration unit.
    if (config.model == "slotboom" || config.model == "old_slotboom") {
        config.referenceDoping =
            scaling.unitSystem().m3ToInternalConcentration(1.0e23);
        config.coefficient = 9.0e-3;
        config.smoothing = 0.5;
        config.offset = 0.0;
    }
    return config;
}

std::unique_ptr<BandgapNarrowing> makeBandgapNarrowingModel(
    const BandgapNarrowingConfig& config)
{
    if (config.model == "none")
        return std::make_unique<NoBandgapNarrowing>();
    if (config.model == "slotboom" || config.model == "old_slotboom")
        return std::make_unique<SlotboomBandgapNarrowing>(config);

    throw std::invalid_argument(
        "makeBandgapNarrowingModel: unknown bandgap narrowing model '" + config.model + "'.");
}

} // namespace vela
