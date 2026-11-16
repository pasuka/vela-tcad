
#include "vela/solver/LinearRefinement.h"
#include "vela/core/ArithmeticPrecision.h"
#include "vela/numerics/ExactBinary64.h"

#include <bit>
#include <cfenv>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <optional>
#include <string_view>
#include <type_traits>
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

    template<class Number>
    bool finish(Number& value) const {
        double sumHigh, sumLow, high, low;
        twoSum(sum, first, sumHigh, sumLow);
        twoSum(sumHigh, sumLow + second, high, low);
        if (!std::isfinite(high) || !std::isfinite(low) ||
            (discarded && high == 0. && low == 0.))
            return false;
        value = Number(high) + Number(low);
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
using Exact = exact_binary64::Rational;

double convert(ArithmeticReal value) {
    // Preserve the existing result policy for binary64 subnormals.
    const ArithmeticReal magnitude = std::abs(value);
    if (magnitude < std::numeric_limits<double>::min()) return 0.;
    if (magnitude > std::numeric_limits<double>::max()) {
        const double infinity = std::numeric_limits<double>::infinity();
        return value < 0. ? -infinity : infinity;
    }
    return static_cast<double>(value);
}

double convert(const Exact& value) {
    return exact_binary64::toNormalDouble(value);
}

template<class Number>
Number exactInput(double value) {
    if constexpr (std::is_same_v<Number, Exact>)
        return exact_binary64::fromDouble(value);
    else
        return Number(value);
}

template<class Number>
struct FallbackRow {
    Eigen::Index row;
    Number residual;
    Number denominator;
};

template<class Number, class WriteRow>
void accumulateFallback(const SparseMatrixd& J, const VectorXd& step,
                        std::vector<FallbackRow<Number>>& fallback, const WriteRow& writeRow) {
    if (fallback.empty()) return;
    using std::abs;
    std::vector<Eigen::Index> index(J.rows(), -1);
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(fallback.size()); ++i)
        index[fallback[i].row] = i;
    for (Eigen::Index col = 0; col < J.outerSize(); ++col) {
        std::optional<Number> direction;
        for (SparseMatrixd::InnerIterator it(J, col); it; ++it) {
            const auto row = index[it.row()];
            if (row < 0) continue;
            if (!direction) direction.emplace(exactInput<Number>(step(col)));
            const Number product = exactInput<Number>(it.value()) * *direction;
            fallback[row].residual += product;
            fallback[row].denominator += abs(product);
        }
    }
    for (const auto& row : fallback) writeRow(row.row, row.residual, row.denominator);
}

template<class Number>
Defect eftDefectImpl(const SparseMatrixd& J, const VectorXd& F, const VectorXd& step) {
    constexpr bool nativeFinish = std::is_same_v<Number, ArithmeticReal>;
    constexpr bool extended = compiledArithmeticPrecision == ArithmeticPrecision::Extended80;
    using std::abs;
    ScopedPerformanceTimer timer("linear.mp_defect.total");
    incrementPerformanceCounter("linear.mp_defect.calls");
    incrementPerformanceCounter("linear.eft_defect.calls");
    if constexpr (nativeFinish)
        incrementPerformanceCounter(extended ? "linear.eft_defect.long_double_calls"
                                             : "linear.eft_defect.double_calls");
    else
        incrementPerformanceCounter("linear.eft_defect.mp100_finish_calls");
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
    std::vector<FallbackRow<Number>> fallback;
    std::vector<FallbackRow<Exact>> exactFallback;
    std::uint64_t rangeRows = 0, cancellationRows = 0;
    auto writeRow = [&]<class Value>(Eigen::Index row, const Value& residual, const Value& denominator) {
        using std::abs;
        result.values(row) = convert(residual);
        if (denominator != 0) {
            const double error = convert(Value(abs(residual) / denominator));
            result.componentwiseBackwardErrors(row) = error;
            result.backwardMax = std::max(result.backwardMax, error);
        }
    };
    {
        ScopedPerformanceTimer stage("linear.eft_defect.finish");
        for (Eigen::Index i = 0; i < F.size(); ++i) {
            const auto& row = rows[i];
            if (row.fallback == FallbackReason::none) {
                Number residual, denominator;
                if (row.residual.finish(residual) && row.denominator.finish(denominator) && denominator >= 0) {
                    writeRow(i, residual, denominator);
                    continue;
                }
                ++cancellationRows;
            } else {
                ++rangeRows;
            }
            if constexpr (nativeFinish) {
                // Recompute discarded tails exactly. Binary64 also needs
                // exact products for rows that exceed its exponent range.
                if (!extended || row.fallback == FallbackReason::none) {
                    const Exact value = exact_binary64::fromDouble(F(i));
                    exactFallback.push_back({i, value, abs(value)});
                    continue;
                }
            }
            const Number value(F(i));
            fallback.push_back({i, value, abs(value)});
        }
    }
    if (!fallback.empty() || !exactFallback.empty()) {
        ScopedPerformanceTimer stage("linear.eft_defect.fallback");
        accumulateFallback(J, step, fallback, writeRow);
        if constexpr (nativeFinish) accumulateFallback(J, step, exactFallback, writeRow);
    }
    if (!result.values.allFinite() || !result.componentwiseBackwardErrors.allFinite())
        throw std::runtime_error("Nonfinite linear refinement defect");
    incrementPerformanceCounter("linear.eft_defect.rows", rows.size() - fallback.size() - exactFallback.size());
    incrementPerformanceCounter("linear.eft_defect.range_fallback_rows", rangeRows);
    incrementPerformanceCounter("linear.eft_defect.cancellation_fallback_rows", cancellationRows);
    return result;
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
        const char* forceMp100 = std::getenv("VELA_BVMETHODS_EFT_MP100_FINISH");
        if (!forceMp100 || std::string_view(forceMp100) != "1")
            return eftDefectImpl<ArithmeticReal>(J, F, step);
        return eftDefectImpl<MP>(J, F, step);
    }
#endif
    incrementPerformanceCounter("linear.eft_defect.whole_call_fallbacks");
    return highPrecisionDefect(J, F, step);
}
} // namespace vela::linear_refinement
