#pragma once
#include "vela/core/Types.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace vela::shadow_convergence {
using Json = nlohmann::json;
// P1 observer only. These values never feed a solver predicate.
inline Json unavailable(const std::string& reason) {
    return {{"status", "not_evaluated"}, {"reason", reason}};
}
inline Json uncalibrated(double value) {
    if (!std::isfinite(value)) return unavailable("nonfinite diagnostic value");
    return {{"status", "not_calibrated"}, {"value", value}};
}
inline std::string combine(const std::vector<std::string>& states) {
    if (states.empty()) return "not_evaluated";
    for (const auto* s : {"fail", "not_calibrated", "not_evaluated"})
        if (std::find(states.begin(), states.end(), s) != states.end()) return s;
    for (const auto& s : states)
        if (s != "pass") return "not_evaluated";
    return "pass";
}
struct Sum {
    long double sum = 0, correction = 0;
    void add(long double x) {
        const auto t = sum + x;
        correction += std::abs(sum) >= std::abs(x) ? (sum - t) + x : (x - t) + sum;
        sum = t;
    }
    long double value() const { return sum + correction; }
};
inline Json vector(const VectorXd& x) {
    Json a = Json::array();
    for (Eigen::Index i = 0; i < x.size(); ++i) a.push_back(x[i]);
    return a;
}
// The caller must identify J: a pre-closure matrix is not the factorized matrix.
inline Json linearQuality(const SparseMatrixd& J, const VectorXd& r, const VectorXd& d) {
    if (J.rows() != r.size() || J.cols() != d.size())
        return unavailable("matrix/vector dimensions differ");
    if (!r.allFinite() || !d.allFinite()) return unavailable("nonfinite operand");
    std::vector<Sum> defect(r.size()), denominator(r.size());
    for (Eigen::Index i = 0; i < r.size(); ++i) {
        defect[i].add(r[i]); denominator[i].add(std::abs(r[i]));
    }
    for (int k = 0; k < J.outerSize(); ++k)
        for (SparseMatrixd::InnerIterator it(J, k); it; ++it) {
            if (!std::isfinite(it.value())) return unavailable("nonfinite matrix");
            const long double product = static_cast<long double>(it.value()) * d[it.col()];
            defect[it.row()].add(product); denominator[it.row()].add(std::abs(product));
        }
    long double maxRatio = 0, defectInf = 0, denominatorInf = 0;
    Eigen::Index worst = -1;
    std::size_t exactZeroRows = 0;
    for (Eigen::Index i = 0; i < r.size(); ++i) {
        const auto a = std::abs(defect[i].value()), b = denominator[i].value();
        if (b == 0) { ++exactZeroRows; continue; }
        if (a / b > maxRatio) { maxRatio = a / b; worst = i; }
        defectInf = std::max(defectInf, a); denominatorInf = std::max(denominatorInf, b);
    }
    return {{"beta_comp", uncalibrated(static_cast<double>(maxRatio))},
        {"beta_norm", denominatorInf > 0
            ? uncalibrated(static_cast<double>(defectInf / denominatorInf))
            : unavailable("zero norm denominator")},
        {"norm_definition", "max|Jd+r| / max(|J||d|+|r|)"},
        {"worst_row", worst}, {"exact_zero_denominator_rows", exactZeroRows},
        {"arithmetic", "long double products and Neumaier compensated row sums"},
        {"forward_error_guaranteed", false}};
}
} // namespace vela::shadow_convergence
