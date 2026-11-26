#include "vela/equation/ExtendedBoltzmannMaterials.h"

#include "vela/physics/BandgapNarrowingArithmetic.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace vela {

ExtendedEquilibriumMaterials buildExtendedEquilibriumMaterials(
    const DeviceMesh& mesh, const MaterialDatabase& database, const DopingModel& doping,
    const BandgapNarrowingConfig& narrowing, long double temperature_K, long double thermalVoltage_V)
{
    return buildExtendedEquilibriumMaterials(mesh, database, doping, narrowing, temperature_K, thermalVoltage_V,
        {narrowing.referenceDoping, narrowing.coefficient, narrowing.smoothing, narrowing.offset});
}

ExtendedEquilibriumMaterials buildExtendedEquilibriumMaterials(
    const DeviceMesh& mesh, const MaterialDatabase& database, const DopingModel& doping,
    const BandgapNarrowingConfig& narrowing, long double temperature_K, long double thermalVoltage_V,
    const ExtendedBandgapCoefficients& coefficients, CarrierStatisticsModel statistics)
{
    if constexpr (compiledArithmeticPrecision != ArithmeticPrecision::Extended80)
        throw std::logic_error("Extended Boltzmann materials require an extended80 build");
    const bool slotboom = narrowing.model == "slotboom" || narrowing.model == "old_slotboom";
    const bool bennett = narrowing.model == "bennett_wilson";
    const bool fermi = usesFermiDirac(statistics);
    if (doping.numNodes() != mesh.numNodes() ||
        !(temperature_K > 0 && thermalVoltage_V > 0) ||
        !std::isfinite(temperature_K) || !std::isfinite(thermalVoltage_V) ||
        (!fermi && statistics != CarrierStatisticsModel::Boltzmann) ||
        (narrowing.fermiStatisticsCorrection && (!fermi || temperature_K != 300.L)) ||
        (!slotboom && !bennett && narrowing.model != "none"))
        throw std::invalid_argument("Invalid extended Boltzmann material inputs");
    if ((slotboom || bennett) &&
        (!(coefficients.referenceDoping > 0 && coefficients.coefficient >= 0 && coefficients.smoothing >= 0) ||
         !std::isfinite(coefficients.referenceDoping) || !std::isfinite(coefficients.coefficient) ||
         !std::isfinite(coefficients.smoothing) || !std::isfinite(coefficients.offset) ||
         (bennett && (coefficients.smoothing != 0 || coefficients.offset != 0))))
        throw std::invalid_argument("Invalid extended bandgap narrowing parameters");

    const auto count = static_cast<Eigen::Index>(mesh.numNodes());
    ExtendedEquilibriumMaterials result;
    result.statistics = statistics;
    for (auto* field : {&result.intrinsic, &result.electronIntrinsic, &result.holeIntrinsic,
                        &result.intrinsicPotential_V, &result.electronDOS, &result.holeDOS})
        field->setZero(count);
    std::vector<const MaterialDatabase::DecimalSemiconductor*> owners(mesh.numNodes(), nullptr);
    for (const auto& cell : mesh.cells()) {
        const auto& name = mesh.getRegion(cell.region_id).material;
        const auto& material = database.getMaterial(name);
        if (!(material.ni > 0 || material.mun > 0 || material.mup > 0))
            continue;
        const auto* bound = database.decimalSemiconductor(name, static_cast<ArithmeticReal>(temperature_K));
        if (!bound || !(bound->intrinsicDensity_internal > 0))
            throw std::invalid_argument("Missing extended Boltzmann material at operating temperature: " + name);
        if (!std::isfinite(bound->boundReferencePotential_V) ||
            (bound->boundReferencePotential_V != 0 && std::abs(temperature_K - 300.L) > 1e-8L))
            throw std::invalid_argument("Invalid isothermal intrinsic reference for material: " + name);
        for (const auto node : cell.node_ids) {
            const auto* current = owners.at(node);
            if (current && current->boundReferencePotential_V != bound->boundReferencePotential_V)
                throw std::invalid_argument("Mixed semiconductor intrinsic references are unsupported");
            if (!current)
                owners[node] = bound;
        }
    }
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const auto* bound = owners[node];
        if (!bound)
            continue;
        const auto impurities = doping.dopingOperandsForAssembly(node);
        const long double totalDoping = static_cast<long double>(impurities[0]) + impurities[1];
        long double delta = 0;
        if (slotboom)
            delta = detail::slotboomNarrowing(totalDoping, coefficients.referenceDoping,
                coefficients.coefficient, coefficients.smoothing, coefficients.offset);
        if (bennett)
            delta = detail::bennettWilsonNarrowing(
                totalDoping, coefficients.referenceDoping, coefficients.coefficient);
        if (narrowing.fermiStatisticsCorrection) {
            delta += thermalVoltage_V * (
                extendedFermiDegeneracy(impurities[0] / bound->electronsDOS_internal, statistics).logCorrection +
                extendedFermiDegeneracy(impurities[1] / bound->holesDOS_internal, statistics).logCorrection);
        }
        const auto intrinsic = bound->intrinsicDensity_internal * std::exp(delta / (2 * thermalVoltage_V));
        const auto reference = bound->boundReferencePotential_V / thermalVoltage_V;
        result.intrinsic[node] = intrinsic;
        result.electronIntrinsic[node] = intrinsic * std::exp(-reference);
        result.holeIntrinsic[node] = intrinsic * std::exp(reference);
        result.intrinsicPotential_V[node] = bound->boundReferencePotential_V;
        result.electronDOS[node] = bound->electronsDOS_internal;
        result.holeDOS[node] = bound->holesDOS_internal;
        if (!(result.electronIntrinsic[node] > 0 && result.holeIntrinsic[node] > 0) ||
            !std::isfinite(result.electronIntrinsic[node]) || !std::isfinite(result.holeIntrinsic[node]))
            throw std::runtime_error("Unrepresentable extended Boltzmann intrinsic populations");
    }
    return result;
}

ExtendedOhmicPopulation extendedBoltzmannOhmicTarget(
    long double netDoping, long double intrinsic, long double intrinsicPotential_V,
    long double thermalVoltage_V, long double appliedBias_V,
    long double potentialScale_V, long double potentialOrigin_V)
{
    for (const auto value : {netDoping, intrinsic, intrinsicPotential_V, thermalVoltage_V,
                             appliedBias_V, potentialScale_V, potentialOrigin_V})
        if (!std::isfinite(value))
            throw std::invalid_argument("Nonfinite extended Boltzmann Ohmic input");
    if (!(intrinsic > 0 && thermalVoltage_V > 0 && potentialScale_V > 0))
        throw std::invalid_argument("Invalid extended Boltzmann Ohmic scale");

    const auto half = netDoping / 2;
    const auto major = std::hypot(half, intrinsic) + std::abs(half);
    const auto minor = intrinsic * (intrinsic / major);
    const auto ratio = half / intrinsic;
    const auto reduced = std::isfinite(ratio) ? std::asinh(ratio)
        : std::copysign(std::log(2.L) + std::log(std::abs(half)) - std::log(intrinsic), half);
    const auto potential = (thermalVoltage_V * reduced + intrinsicPotential_V +
                            (appliedBias_V - potentialOrigin_V)) / potentialScale_V;
    if (!std::isfinite(potential) || !std::isfinite(major) || !(minor > 0))
        throw std::runtime_error("Unrepresentable extended Boltzmann Ohmic target");
    return netDoping >= 0 ? ExtendedOhmicPopulation{potential, major, minor}
                         : ExtendedOhmicPopulation{potential, minor, major};
}

ExtendedOhmicPopulation extendedEquilibriumOhmicTarget(
    long double netDoping, const ExtendedEquilibriumMaterials& materials, Index node,
    long double thermalVoltage_V)
{
    const auto estimate = extendedBoltzmannOhmicTarget(
        netDoping, materials.intrinsic[node], materials.intrinsicPotential_V[node],
        thermalVoltage_V, 0, thermalVoltage_V, 0);
    if (materials.statistics == CarrierStatisticsModel::Boltzmann) return estimate;

    const auto populations = [&](long double reduced) {
        return std::pair{
            extendedFermiPopulation(materials.electronDOS[node], materials.electronIntrinsic[node], reduced,
                                    materials.statistics),
            extendedFermiPopulation(materials.holeDOS[node], materials.holeIntrinsic[node], -reduced,
                                    materials.statistics)};
    };
    auto reduced = estimate.potential;
    auto lower = reduced - 1;
    auto upper = reduced + 1;
    for (long double distance = 2; distance <= 8192; distance *= 2) {
        const auto [leftN, leftP] = populations(lower);
        const auto [rightN, rightP] = populations(upper);
        if (leftN.first - leftP.first <= netDoping && rightN.first - rightP.first >= netDoping)
            break;
        if (distance == 8192) throw std::runtime_error("Cannot bracket extended Fermi Ohmic target");
        lower = estimate.potential - distance;
        upper = estimate.potential + distance;
    }
    constexpr auto epsilon = std::numeric_limits<long double>::epsilon();
    for (int iteration = 0; iteration < 100; ++iteration) {
        const auto [electron, hole] = populations(reduced);
        const auto residual = electron.first - hole.first - netDoping;
        if (std::abs(residual) <= 8 * epsilon * std::max(electron.first, hole.first))
            return {reduced, electron.first, hole.first};
        if (residual > 0) upper = reduced;
        else lower = reduced;
        const auto midpoint = lower + (upper - lower) / 2;
        const auto next = reduced - residual / (electron.second + hole.second);
        if (next == reduced || midpoint == lower || midpoint == upper)
            return {reduced, electron.first, hole.first};
        reduced = next > lower && next < upper ? next : midpoint;
    }
    throw std::runtime_error("Extended Fermi Ohmic target failed to converge");
}

} // namespace vela
