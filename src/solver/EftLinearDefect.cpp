
#include "vela/solver/LinearRefinement.h"

#include <bit>
#include <cfenv>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <optional>
#include <string_view>
#if defined(__SSE__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 1)
#include <xmmintrin.h>
#define VELA_EFT_MXCSR
#endif

namespace vela::linear_refinement {
namespace {
#if defined(VELA_STRICT_EFT_FP) && defined(VELA_EFT_MXCSR)
constexpr double splitter = 134217729.;

class FloatingEnvironment {
public:
    FloatingEnvironment() : held_(std::feholdexcept(&saved_) == 0) {}
    ~FloatingEnvironment() { if (held_) std::fesetenv(&saved_); }
    bool supported() const {
        constexpr unsigned flushOrDenormalsAreZero = (1u << 15) | (1u << 6);
        return held_ && std::numeric_limits<double>::is_iec559 &&
               std::fegetround() == FE_TONEAREST && !(_mm_getcsr() & flushOrDenormalsAreZero);
    }

private:
    std::fenv_t saved_{};
    bool held_;
};

void twoSum(double a, double b, double& high, double& low) {
    high = a + b;
    const double shifted = high - a;
    low = (a - (high - shifted)) + (b - shifted);
}

void twoProduct(double a, double b, double& high, double& low) {
    const double splitA = a * splitter;
    const double splitB = b * splitter;
    const double ah = splitA - (splitA - a), al = a - ah;
    const double bh = splitB - (splitB - b), bl = b - bh;
    high = a * b;
    low = ((ah * bh - high) + ah * bl + al * bh) + al * bl;
}

int leastSetBit(double v) {
    const auto u = std::bit_cast<std::uint64_t>(v) & 0x7fffffffffffffffULL;
    const int exponent = static_cast<int>(u >> 52);
    const auto mantissa = (u & 0xfffffffffffffULL) | (exponent ? 0x10000000000000ULL : 0);
    return (exponent ? exponent - 1023 - 52 : -1074) + std::countr_zero(mantissa);
}

struct Accumulator {
    double sum = 0, first = 0, second = 0;
    bool discarded = false;

    void add(double value, double tail) {
        double q1, q2, q3;
        twoSum(sum, value, sum, q1);
        twoSum(first, q1, first, q2);
        twoSum(first, tail, first, q3);
        const double combined = q2 + q3;
        const double updated = second + combined;
        // Detect the cancellation witness where the third level loses a
        // nonzero term and later cancels to zero. Keep the K3 sum order.
        discarded |= (combined == q2 && q3 != 0.) ||
                     (combined == q3 && q2 != 0.) ||
                     (updated == second && combined != 0.) ||
                     (updated == combined && second != 0.);
        second = updated;
    }

    bool finite() const {
        return std::isfinite(sum) && std::isfinite(first) && std::isfinite(second);
    }

    bool finish(MP& value) const {
        double sumHigh, sumLow, high, low;
        twoSum(sum, first, sumHigh, sumLow);
        twoSum(sumHigh, sumLow + second, high, low);
        if (!std::isfinite(high) || !std::isfinite(low) ||
            (discarded && high == 0. && low == 0.))
            return false;
        value = MP(high) + MP(low);
        return true;
    }
};

enum class FallbackReason { none, range, cancellation };
struct Row {
    Accumulator residual, denominator;
    FallbackReason fallback = FallbackReason::none;

    explicit Row(double f) : residual{f}, denominator{std::abs(f)} {}

    void add(double a, double b) {
        if (fallback != FallbackReason::none) return;
        if (!std::isfinite(a * splitter) || !std::isfinite(b * splitter) || !std::isfinite(a * b) ||
            (a != 0. && b != 0. &&
             (std::fpclassify(a) == FP_SUBNORMAL || std::fpclassify(b) == FP_SUBNORMAL ||
              leastSetBit(a) + leastSetBit(b) < -1074))) {
            fallback = FallbackReason::range;
            return;
        }
        double high, low;
        twoProduct(a, b, high, low);
        residual.add(high, low);
        if (high < 0.) { high = -high; low = -low; }
        denominator.add(high, low);
        if (!residual.finite() || !denominator.finite())
            fallback = FallbackReason::range;
    }
};

double convert(const MP& value) {
    // Match cpp_dec_float's binary64 range policy before formatting.
    if (!boost::multiprecision::isfinite(value)) return static_cast<double>(value);
    static const MP minimum((std::numeric_limits<double>::min)());
    static const MP maximum((std::numeric_limits<double>::max)());
    const MP magnitude = abs(value);
    if (value == 0 || magnitude < minimum) return 0.;
    if (magnitude > maximum) return static_cast<double>(value);

    const std::string text = value.str(std::numeric_limits<double>::max_digits10, std::ios_base::scientific);
    double parsed = 0.;
    const auto end = text.data() + text.size();
    const auto result = std::from_chars(text.data(), end, parsed, std::chars_format::scientific);
    return result.ec == std::errc{} && result.ptr == end ? parsed : static_cast<double>(value);
}
#endif
} // namespace

bool useBvMethodsEftDefect() {
    const char* value = std::getenv("VELA_BVMETHODS_MP100_DEFECT");
    return !value || std::string_view(value) != "1";
}

Defect eftDefect(const SparseMatrixd& J, const VectorXd& F, const VectorXd& step) {
#if defined(VELA_STRICT_EFT_FP) && defined(VELA_EFT_MXCSR)
    if (J.rows() != F.size() || J.cols() != step.size() || !F.allFinite() || !step.allFinite())
        throw std::invalid_argument("Invalid refinement defect dimensions or input");
    FloatingEnvironment environment;
    if (environment.supported()) {
        ScopedPerformanceTimer timer("linear.mp_defect.total");
        incrementPerformanceCounter("linear.mp_defect.calls");
        incrementPerformanceCounter("linear.eft_defect.calls");
        std::vector<Row> rows;
        {
            ScopedPerformanceTimer stage("linear.eft_defect.accumulate");
            rows.reserve(F.size());
            for (Eigen::Index i = 0; i < F.size(); ++i) rows.emplace_back(F(i));
            for (Eigen::Index col = 0; col < J.outerSize(); ++col)
                for (SparseMatrixd::InnerIterator it(J, col); it; ++it)
                    rows[it.row()].add(it.value(), step(col));
        }

        Defect result{VectorXd(F.size()), 0., VectorXd::Zero(F.size())};
        struct FallbackRow {
            Eigen::Index row;
            MP residual;
            MP denominator;
        };
        std::vector<FallbackRow> fallback;
        std::vector<Eigen::Index> index(F.size(), -1);
        std::uint64_t rangeRows = 0, cancellationRows = 0;
        auto writeRow = [&](Eigen::Index row, const MP& residual, const MP& denominator) {
            result.values(row) = convert(residual);
            if (denominator != 0) {
                const double error = convert(MP(abs(residual) / denominator));
                result.componentwiseBackwardErrors(row) = error;
                result.backwardMax = std::max(result.backwardMax, error);
            }
        };
        {
            ScopedPerformanceTimer stage("linear.eft_defect.finish");
            for (Eigen::Index i = 0; i < F.size(); ++i) {
                const auto& row = rows[i];
                if (row.fallback == FallbackReason::none) {
                    MP residual, denominator;
                    if (row.residual.finish(residual) && row.denominator.finish(denominator) && denominator >= 0) {
                        writeRow(i, residual, denominator);
                        continue;
                    }
                    ++cancellationRows;
                } else {
                    ++rangeRows;
                }
                index[i] = static_cast<Eigen::Index>(fallback.size());
                const MP value(F(i));
                fallback.push_back({i, value, abs(value)});
            }
        }
        if (!fallback.empty()) {
            ScopedPerformanceTimer stage("linear.eft_defect.fallback");
            for (Eigen::Index col = 0; col < J.outerSize(); ++col) {
                std::optional<MP> direction;
                for (SparseMatrixd::InnerIterator it(J, col); it; ++it) {
                    const auto row = index[it.row()];
                    if (row < 0) continue;
                    if (!direction) direction.emplace(step(col));
                    const MP product = MP(it.value()) * *direction;
                    fallback[row].residual += product;
                    fallback[row].denominator += abs(product);
                }
            }
            for (const auto& row : fallback) writeRow(row.row, row.residual, row.denominator);
        }
        if (!result.values.allFinite() || !result.componentwiseBackwardErrors.allFinite())
            throw std::runtime_error("Nonfinite linear refinement defect");
        incrementPerformanceCounter("linear.eft_defect.rows", rows.size() - fallback.size());
        incrementPerformanceCounter("linear.eft_defect.range_fallback_rows", rangeRows);
        incrementPerformanceCounter("linear.eft_defect.cancellation_fallback_rows", cancellationRows);
        return result;
    }
#endif
    incrementPerformanceCounter("linear.eft_defect.whole_call_fallbacks");
    return highPrecisionDefect(J, F, step);
}
} // namespace vela::linear_refinement
