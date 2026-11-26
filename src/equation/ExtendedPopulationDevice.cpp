#include "vela/equation/ExtendedPopulationDevice.h"

#include <cmath>
#include <stdexcept>

namespace vela {

ExtendedPopulationDeviceSystem assembleExtendedPopulationDevice(
    const BasicPoissonAssembler<long double>& poisson, const ExtendedPopulationState& state,
    const ExtendedEquilibriumMaterials& material, const ExtendedOhmicTargets& contacts,
    std::span<const Index> drivenContactNodes, const ExtendedIalPopulationContinuity& continuity,
    const ExtendedAveragePoissonGeometry* geometry,
    const ExtendedPopulationNormalization* normalization)
{
    const auto count = poisson.matrix().rows();
    if (count <= 0 || state.coordinates.size() != 3 * count || !state.coordinates.allFinite() ||
        !std::isfinite(state.deviceBias_V) || state.potentialOrigin_V != 0 ||
        state.potentialScale_V != poisson.potentialScale_V() ||
        material.electronIntrinsic.size() != count || material.holeIntrinsic.size() != count ||
        std::ssize(contacts.potential) != count || contacts.electrons.size() != count ||
        contacts.holes.size() != count || continuity.residual.size() != 2 * count ||
        continuity.jacobian.rows() != 2 * count || continuity.jacobian.cols() != 3 * count)
        throw std::invalid_argument("Inconsistent extended population device assembly inputs");

    std::vector<bool> driven(count, false);
    for (const auto node : drivenContactNodes) {
        if (node >= contacts.potential.size() || !contacts.potential[node] || driven[node])
            throw std::invalid_argument("Invalid driven extended Ohmic node");
        driven[node] = true;
    }
    ExtendedPopulationDeviceSystem result;
    result.residual.resize(3 * count);
    result.residual.head(count) = -poisson.rhs();
    result.residual.tail(2 * count) = continuity.residual;
    result.voltageDerivative.setZero(3 * count);
    result.electronFree.resize(count);
    result.holeFree.resize(count);
    std::vector<Eigen::Triplet<long double>> terms;
    terms.reserve(poisson.matrix().nonZeros() + continuity.jacobian.nonZeros() + 5 * count);
    for (Eigen::Index col = 0; col < count; ++col) {
        for (ExtendedSparseMatrix::InnerIterator term(poisson.matrix(), col); term; ++term) {
            if (contacts.potential[term.row()])
                continue;
            terms.emplace_back(term.row(), col, term.value());
            if (term.row() != col)
                result.residual[term.row()] -= term.value() *
                    (state.coordinates[term.row()] - state.coordinates[col]);
        }
    }
    for (Eigen::Index node = 0; node < count; ++node) {
        if (const auto& target = contacts.potential[node]) {
            terms.emplace_back(node, node, 1);
            if (driven[node]) {
                const auto biasedTarget = *target + state.deviceBias_V / state.potentialScale_V;
                result.residual[node] = state.coordinates[node] - biasedTarget;
                result.voltageDerivative[node] = -1 / state.potentialScale_V;
            } else {
                result.residual[node] = state.coordinates[node] - *target;
                result.retainedIndependentPotentialRows.push_back(static_cast<int>(node));
            }
        } else {
            const auto weight = geometry ? poisson.mobileChargeWeight(node, *geometry)
                                         : poisson.mobileChargeWeight(node);
            if (weight != 0) {
                const auto electron = material.electronIntrinsic[node] > 0;
                const auto hole = material.holeIntrinsic[node] > 0;
                const auto n = electron ? state.coordinates[count + node] : 0.L;
                const auto p = hole ? state.coordinates[2 * count + node] : 0.L;
                result.residual[node] += weight * (n - p);
                if (electron)
                    terms.emplace_back(node, count + node, weight);
                if (hole)
                    terms.emplace_back(node, 2 * count + node, -weight);
            }
        }
        for (const bool electron : {true, false}) {
            const auto offset = electron ? count : 2 * count;
            const auto& intrinsic = electron ? material.electronIntrinsic : material.holeIntrinsic;
            const auto& target = electron ? contacts.electrons : contacts.holes;
            auto& free = electron ? result.electronFree : result.holeFree;
            free[node] = !contacts.potential[node] && intrinsic[node] > 0;
            if (!free[node]) {
                const auto row = offset + node;
                result.residual[row] = state.coordinates[row] - target[node];
                terms.emplace_back(row, row, 1);
                result.prescribedPopulationRows.push_back(static_cast<int>(row));
            }
        }
    }
    for (Eigen::Index col = 0; col < continuity.jacobian.outerSize(); ++col) {
        for (ExtendedSparseMatrix::InnerIterator term(continuity.jacobian, col); term; ++term) {
            const auto row = term.row();
            if (row < count ? result.electronFree[row] : result.holeFree[row - count])
                terms.emplace_back(count + row, col, term.value());
        }
    }
    result.jacobian.resize(3 * count, 3 * count);
    result.jacobian.setFromTriplets(terms.begin(), terms.end());
    if (normalization) {
        const auto [concentration, poissonFactor, continuityFactor] = *normalization;
        for (auto scale : {concentration, poissonFactor, continuityFactor})
            if (!(scale > 0) || !std::isfinite(scale))
                throw std::invalid_argument("Invalid extended population equation scale");

        ExtendedVector rows = ExtendedVector::Constant(3 * count, continuityFactor);
        for (Eigen::Index i = 0; i < count; ++i) {
            rows[i] = contacts.potential[i] ? 1.L : poissonFactor;
            if (!result.electronFree[i])
                rows[count + i] = 1.L / concentration;
            if (!result.holeFree[i])
                rows[2 * count + i] = 1.L / concentration;
        }
        result.residual.array() *= rows.array();
        result.voltageDerivative.array() *= rows.array();
        for (Eigen::Index col = 0; col < result.jacobian.outerSize(); ++col) {
            for (ExtendedSparseMatrix::InnerIterator entry(result.jacobian, col); entry; ++entry) {
                const auto row = entry.row();
                const bool populationIdentity = row >= count && (row < 2 * count
                    ? !result.electronFree[row - count] : !result.holeFree[row - 2 * count]);
                if (!populationIdentity)
                    entry.valueRef() *= rows[row] * (col < count ? 1.L : concentration);
            }
        }
    }
    if (!result.residual.allFinite())
        throw std::runtime_error("Nonfinite extended population device residual");
    return result;
}

} // namespace vela
