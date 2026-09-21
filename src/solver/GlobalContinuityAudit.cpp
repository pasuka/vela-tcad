#include "vela/solver/NewtonSolver.h"
#include <nlohmann/json.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <algorithm>
#include <bit>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>

namespace vela {
namespace {
Real pairwise(const std::vector<Real>& values, std::size_t first, std::size_t last)
{
    if (last == first) return 0;
    if (last - first == 1) return values[first];
    const auto mid = first + (last - first) / 2;
    return pairwise(values, first, mid) + pairwise(values, mid, last);
}
Real neumaier(const std::vector<Real>& values)
{
    Real sum = 0, correction = 0;
    for (Real value : values) {
        const Real next = sum + value;
        correction += std::abs(sum) >= std::abs(value)
            ? (sum - next) + value : (value - next) + sum;
        sum = next;
    }
    return sum + correction;
}
}

nlohmann::json globalContinuityReductionAudit(const std::vector<Real>& values)
{
    using High = boost::multiprecision::cpp_bin_float_quad;
    High reference = 0;
    long double extended = 0;
    Real kahan = 0, correction = 0;
    for (Real value : values) {
        reference += High(value);
        extended += static_cast<long double>(value);
        const Real adjusted = value - correction;
        const Real next = kahan + adjusted;
        correction = (next - kahan) - adjusted;
        kahan = next;
    }
    auto sum = [](const auto& v) { return std::accumulate(v.begin(), v.end(), Real{0}); };
    auto ordered = values;
    std::stable_sort(ordered.begin(), ordered.end(), [](Real a, Real b) { return std::abs(a) < std::abs(b); });
    nlohmann::json results = {{"double", sum(values)}, {"pairwise", pairwise(values, 0, values.size())},
        {"kahan", kahan}, {"neumaier", neumaier(values)}, {"long_double", static_cast<Real>(extended)},
        {"ascending", sum(ordered)}, {"quad", static_cast<Real>(reference)}};
    std::reverse(ordered.begin(), ordered.end()); results["descending"] = sum(ordered);
    std::mt19937 generator(20260905); std::shuffle(ordered.begin(), ordered.end(), generator);
    results["permuted"] = sum(ordered);
    Real maxError = 0;
    for (const auto& item : results.items())
        maxError = std::max(maxError, static_cast<Real>(abs(High(item.value().get<Real>()) - reference)));
    return {{"sums", results}, {"max_error_vs_quad", maxError},
        {"quad_reference", reference.str(36, std::ios_base::scientific)},
        {"long_double_digits", std::numeric_limits<long double>::digits}, {"count", values.size()}};
}

nlohmann::json auditGlobalContinuityRoundoff(const CoupledDDAssembler& assembler,
    const VectorXd& state, const CoupledDDBoundaryConditions& bcs,
    const NewtonGlobalContinuityClosureConfig& config, Index edgeCount)
{
    const int N = static_cast<int>(assembler.numNodes());
    std::vector<Index> ec, hc;
    for (const auto& entry : bcs.phin) ec.push_back(entry.first);
    for (const auto& entry : bcs.phip) hc.push_back(entry.first);
    auto freeValues = [&](const VectorXd& r, int carrier) {
        std::vector<Real> values;
        for (int i = 0; i < N; ++i)
            if (!(carrier == 0 ? bcs.phin : bcs.phip).contains(i)) values.push_back(r((carrier + 1) * N + i));
        return values;
    };
    const VectorXd baseline = assembler.residual(state, bcs);
    nlohmann::json result = {{"event", "roundoff_audit"},
        {"fixed_residual_reduction", {globalContinuityReductionAudit(freeValues(baseline, 0)),
            globalContinuityReductionAudit(freeValues(baseline, 1))}}, {"variants", nlohmann::json::array()}};
    std::vector<Index> reverse(edgeCount), permuted;
    std::iota(reverse.begin(), reverse.end(), Index{0}); permuted = reverse;
    std::reverse(reverse.begin(), reverse.end());
    std::mt19937 generator(20260905); std::shuffle(permuted.begin(), permuted.end(), generator);
    for (const std::string order : {"forward", "reverse", "permuted"}) {
        for (bool compensated : {false, true}) {
            CoupledDDAssemblyAuditOptions options;
            options.compensated = compensated;
            if (order == "reverse") options.edgeOrder = reverse;
            if (order == "permuted") options.edgeOrder = permuted;
            const VectorXd residual = assembler.auditAssemblyResidual(state, bcs, options);
            const VectorXd repeat = assembler.auditAssemblyResidual(state, bcs, options);
            int changed = 0, repeatChanged = 0;
            for (int i = 0; i < residual.size(); ++i) {
                changed += std::bit_cast<std::uint64_t>(residual(i)) != std::bit_cast<std::uint64_t>(baseline(i));
                repeatChanged += std::bit_cast<std::uint64_t>(residual(i)) != std::bit_cast<std::uint64_t>(repeat(i));
            }
            const auto terms = assembler.auditAssemblyCarrierTerms(state, {}, options);
            const auto closure = evaluateGlobalContinuityClosure(terms, ec, hc, config);
            std::array<long double, 2> gauge{}, physical{}, inactiveCount{};
            for (int i = 0; i < N; ++i) {
                if (!bcs.phin.contains(i)) {
                    gauge[0] += terms[i].electronGauge;
                    physical[0] += static_cast<long double>(terms[i].electronFlux) + terms[i].electronRecombination + terms[i].electronImpact;
                    inactiveCount[0] += !terms[i].electronContinuityActive;
                }
                if (!bcs.phip.contains(i)) {
                    gauge[1] += terms[i].holeGauge;
                    physical[1] += static_cast<long double>(terms[i].holeFlux) + terms[i].holeRecombination + terms[i].holeImpact;
                    inactiveCount[1] += !terms[i].holeContinuityActive;
                }
            }
            const auto er = globalContinuityReductionAudit(freeValues(residual, 0));
            const auto hr = globalContinuityReductionAudit(freeValues(residual, 1));
            result["variants"].push_back({{"order", order}, {"compensated", compensated},
                {"residual_linf_delta", (residual - baseline).lpNorm<Eigen::Infinity>()},
                {"changed_entries", changed}, {"repeat_bitwise_changed_entries", repeatChanged},
                {"electron_reduction", er}, {"hole_reduction", hr},
                {"physical_mismatch", {closure.electron.mismatch, closure.hole.mismatch}},
                {"gauge_free_sum", gauge}, {"inactive_free_rows", inactiveCount},
                {"physical_term_free_sum", physical}});
        }
    }
    return result;
}
} // namespace vela
