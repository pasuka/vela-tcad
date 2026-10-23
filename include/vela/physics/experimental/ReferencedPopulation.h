#pragma once

#include <cmath>
#include <stdexcept>

namespace vela::experimental {

// Frozen material/reference inputs. No fitted offset or derivative substitute.
// The Fermi callback must evaluate the value of the selected statistics model
// in Scalar precision. Passing a double-only callback defeats this contract.
struct PopulationReference {
    double intrinsic = 0.0;
    double densityOfStates = 0.0;
    double thermalVoltage = 0.0;
    double potentialScale = 1.0;
    double quasiFermiReference = 0.0;
    double quantumPotential = 0.0;
    bool electron = true;
    bool fermi = false;
};

template<class Scalar, class FermiValue>
Scalar referencedPopulation(const PopulationReference& r,
                            const Scalar& packedPsi, const Scalar& packedQf,
                            FermiValue&& fermiValue) {
    using std::exp;
    using std::isfinite;
    if (!std::isfinite(r.intrinsic) || r.intrinsic < 0.0 ||
        !std::isfinite(r.thermalVoltage) || r.thermalVoltage <= 0.0 ||
        !std::isfinite(r.potentialScale) || r.potentialScale <= 0.0 ||
        !std::isfinite(r.quasiFermiReference) || !std::isfinite(r.quantumPotential) ||
        !isfinite(packedPsi) || !isfinite(packedQf))
        throw std::invalid_argument("Invalid referenced population inputs");
    if (r.intrinsic == 0.0) return Scalar(0);
    // Do not form either physical potential, their difference, or eta in double.
    Scalar argument = (packedPsi * Scalar(r.potentialScale) -
        Scalar(r.quasiFermiReference) - packedQf * Scalar(r.potentialScale));
    if (r.electron) argument -= Scalar(r.quantumPotential);
    else argument = -argument;
    argument /= Scalar(r.thermalVoltage);
    Scalar value;
    if (r.fermi) {
        if (!std::isfinite(r.densityOfStates) || r.densityOfStates <= 0.0)
            throw std::invalid_argument("Invalid Fermi density of states");
        // Keep the binary64 material log used by the ordinary evaluator fixed;
        // this control changes state arithmetic, not material normalization.
        const double logRatio = std::log(r.intrinsic / r.densityOfStates);
        if (!std::isfinite(logRatio))
            throw std::invalid_argument("Unrepresentable material log ratio");
        value = Scalar(r.densityOfStates) * fermiValue(Scalar(argument + Scalar(logRatio)));
    } else {
        if (argument < -500) argument = -500;
        if (argument > 500) argument = 500;
        value = Scalar(r.intrinsic) * exp(argument);
    }
    if (!isfinite(value) || value < 0)
        throw std::runtime_error("Unrepresentable referenced population");
    return value;
}

template<class Scalar, class FermiValue>
Scalar referencedPopulationIncrement(const PopulationReference& r,
    const Scalar& beforePsi, const Scalar& beforeQf,
    const Scalar& afterPsi, const Scalar& afterQf, FermiValue&& fermiValue) {
    const Scalar before = referencedPopulation(r, beforePsi, beforeQf, fermiValue);
    const Scalar after = referencedPopulation(r, afterPsi, afterQf, fermiValue);
    return after - before;
}
} // namespace vela::experimental
