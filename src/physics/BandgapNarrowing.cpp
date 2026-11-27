#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/BandgapNarrowingArithmetic.h"
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
    return detail::slotboomNarrowing(effectiveConcentration, config_.referenceDoping,
                                   config_.coefficient, config_.smoothing, config_.offset);
}

/// Sentaurus Bennett-Wilson: Eref * log(Ntotal/Nref)^2 above Nref, zero below.
/// Ntotal is donors + acceptors in the active concentration unit, not net
/// doping or injected carriers. Material dEg0 belongs in the base ni, not here.
BennettWilsonBandgapNarrowing::BennettWilsonBandgapNarrowing(BandgapNarrowingConfig config)
    : config_(std::move(config))
{
    if (!(config_.referenceDoping > 0.0) || !std::isfinite(config_.referenceDoping))
        throw std::invalid_argument("BennettWilson: referenceDoping must be positive and finite.");
    if (!(config_.coefficient >= 0.0) || !std::isfinite(config_.coefficient))
        throw std::invalid_argument("BennettWilson: coefficient must be nonnegative and finite.");
    if (config_.offset != 0.0)
        throw std::invalid_argument("BennettWilson: offset must be zero; material dEg0 belongs in base ni.");
    if (config_.smoothing != 0.0)
        throw std::invalid_argument("BennettWilson: smoothing is not part of this model; use zero.");
}

Real BennettWilsonBandgapNarrowing::deltaEg(Real impurityConcentration, Real, Real) const
{
    if (!(impurityConcentration >= 0.0) || !std::isfinite(impurityConcentration))
        throw std::invalid_argument("BennettWilson: total impurity must be nonnegative and finite.");
    return detail::bennettWilsonNarrowing(
        impurityConcentration, config_.referenceDoping, config_.coefficient);
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
                                      Real thermalVoltage, CarrierStatisticsModel model)
{
    if (!(thermalVoltage > 0.0) || !std::isfinite(thermalVoltage))
        throw std::invalid_argument(
            "fermiStatisticsBandgapCorrection: thermal voltage must be positive and finite.");
    if (!(Nc > 0.0) || !(Nv > 0.0) || !std::isfinite(Nc) || !std::isfinite(Nv))
        throw std::invalid_argument(
            "fermiStatisticsBandgapCorrection: density of states must be positive and finite.");

    auto speciesCorrection = [model](Real concentration, Real densityOfStates) {
        if (!(concentration > 0.0) || !std::isfinite(concentration))
            return Real{0.0};
        const Real ratio = concentration / densityOfStates;
        if (!(ratio > 0.0) || !std::isfinite(ratio))
            return Real{0.0};
        return std::max<Real>(
            0.0, inverseFermiDiracHalf(ratio, model) - std::log(ratio));
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
    if (config.model == "bennett_wilson") {
        // Sentaurus Device User Guide Eq. 165 / Bennett parameter set.
        config.referenceDoping =
            scaling.unitSystem().m3ToInternalConcentration(3.162e24);
        config.coefficient = 6.84e-3;
        config.smoothing = 0.0;
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
    if (config.model == "bennett_wilson")
        return std::make_unique<BennettWilsonBandgapNarrowing>(config);

    throw std::invalid_argument(
        "makeBandgapNarrowingModel: unknown bandgap narrowing model '" + config.model + "'.");
}

} // namespace vela
