
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

template<class Scalar>
void twoSum(Scalar a, Scalar b, Scalar& high, Scalar& low) {
    high = a + b;
    const Scalar shifted = high - a;
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

template<class Scalar>
struct Accumulator {
    Scalar sum = 0, first = 0, second = 0;
    bool discarded = false;

    void add(Scalar value, Scalar tail) {
        Scalar q1, q2, q3;
        twoSum(sum, value, sum, q1);
        twoSum(first, q1, first, q2);
        twoSum(first, tail, first, q3);
        const Scalar combined = q2 + q3;
        const Scalar updated = second + combined;
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
        Scalar sumHigh, sumLow, high, low;
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
    Accumulator<double> residual, denominator;
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
#if defined(VELA_ARITHMETIC_EXTENDED80)
using Wide = long double;
using WideMatrix = Eigen::SparseMatrix<Wide>;
using WideVector = Eigen::VectorX<Wide>;
using exact_binary64::Integer;
using exact_binary64::Rational;
constexpr int wideDigits = std::numeric_limits<Wide>::digits;
static_assert(wideDigits == 64 && std::numeric_limits<Wide>::radix == 2);

Rational fromWide(Wide value) {
    if (!std::isfinite(value))
        throw std::invalid_argument("Nonfinite extended refinement operand");
    if (value == 0)
        return 0;
    int exponent;
    const auto fraction = std::frexp(std::abs(value), &exponent);
    Integer mantissa = static_cast<std::uint64_t>(std::ldexp(fraction, wideDigits));
    if (std::signbit(value))
        mantissa = -mantissa;
    const int shift = exponent - wideDigits;
    if (shift >= 0)
        return Rational(mantissa << shift);
    return Rational(mantissa, Integer(1) << -shift);
}

// Keep the subnormal range and round once to nearest, ties to even.
Wide toWide(const Rational& value) {
    Integer numerator = boost::multiprecision::numerator(value);
    if (numerator == 0)
        return 0;
    const bool negative = numerator < 0;
    if (negative)
        numerator = -numerator;
    Integer denominator = boost::multiprecision::denominator(value);
    int exponent = static_cast<int>(msb(numerator)) - static_cast<int>(msb(denominator));
    if (exponent >= 0 ? numerator < (denominator << exponent) : (numerator << -exponent) < denominator)
        --exponent;
    if (exponent >= std::numeric_limits<Wide>::max_exponent)
        return negative ? -std::numeric_limits<Wide>::infinity() : std::numeric_limits<Wide>::infinity();
    int shift = wideDigits - 1 - std::max(exponent, std::numeric_limits<Wide>::min_exponent - 1);
    if (shift >= 0)
        numerator <<= shift;
    else
        denominator <<= -shift;
    Integer quotient, remainder;
    divide_qr(numerator, denominator, quotient, remainder);
    remainder <<= 1;
    if (remainder > denominator || (remainder == denominator && bit_test(quotient, 0)))
        ++quotient;
    if (quotient != 0 && msb(quotient) >= wideDigits) {
        quotient >>= 1;
        --shift;
    }
    const auto rounded = std::scalbn(static_cast<Wide>(quotient.convert_to<std::uint64_t>()), -shift);
    return negative ? -rounded : rounded;
}

bool finishWide(const Accumulator<Wide>& accumulator, Wide& value) {
    Wide high, low;
    twoSum(accumulator.sum, accumulator.first, high, low);
    value = high + (low + accumulator.second);
    return std::isfinite(value) && !(accumulator.discarded && value == 0);
}

Wide wideProductTail(Wide a, Wide b, Wide product) {
    // Split the 64-bit significands in a range where their partial products
    // stay normal and the splitting multiplication cannot overflow.
    if (std::abs(a) >= 0x1p-16300L && std::abs(a) <= 0x1p16340L &&
        std::abs(b) >= 0x1p-16300L && std::abs(b) <= 0x1p16340L &&
        std::abs(product) >= 0x1p-16250L && std::abs(product) <= 0x1p16380L) {
        constexpr Wide splitter = 0x1.00000001p32L;
        const Wide sa = splitter * a;
        const Wide sb = splitter * b;
        const Wide ah = sa - (sa - a);
        const Wide bh = sb - (sb - b);
        const Wide al = a - ah;
        const Wide bl = b - bh;
        const Wide e1 = product - ah * bh;
        const Wide e2 = e1 - al * bh;
        const Wide e3 = e2 - ah * bl;
        const Wide tail = al * bl - e3;
        return tail == 0.L ? 0.L : tail;
    }
    return std::fma(a, b, -product);
}

BasicDefect<Wide> extendedEftDefectImpl(const WideMatrix& matrix, const WideVector& residual,
                                       const WideVector& direction) {
    ScopedPerformanceTimer timer("linear.extended_eft_defect.total");
    incrementPerformanceCounter("linear.extended_eft_defect.calls");
    struct WideRow { Accumulator<Wide> residual, denominator; bool range = false; };
    std::vector<WideRow> rows(residual.size());
    for (Eigen::Index i = 0; i < residual.size(); ++i) {
        rows[i].residual.sum = residual[i];
        rows[i].denominator.sum = std::abs(residual[i]);
    }
    // A normal leading product is insufficient near underflow because the
    // fused remainder may need another two mantissas of exponent range.
    constexpr auto exactProductFloor = std::numeric_limits<Wide>::min() * 0x1p128L;
    for (Eigen::Index col = 0; col < matrix.outerSize(); ++col) {
        for (WideMatrix::InnerIterator entry(matrix, col); entry; ++entry) {
            auto& row = rows[entry.row()];
            if (row.range)
                continue;
            const auto a = entry.value(), b = direction[col];
            const auto product = a * b;
            if (!std::isfinite(product) || (a != 0 && b != 0 && std::abs(product) < exactProductFloor)) {
                row.range = true;
                continue;
            }
            const auto tail = wideProductTail(a, b, product);
            row.residual.add(product, tail);
            row.denominator.add(std::abs(product), std::signbit(product) ? -tail : tail);
        }
    }

    BasicDefect<Wide> result{WideVector(residual.size()), 0, WideVector(residual.size())};
    const auto writeRow = [&](Eigen::Index i, Wide value, Wide error) {
        if (!std::isfinite(value) || !std::isfinite(error))
            throw std::runtime_error("Nonfinite extended linear refinement defect");
        result.values[i] = value;
        result.componentwiseBackwardErrors[i] = error;
        result.backwardMax = std::max(result.backwardMax, error);
    };
    struct ExactRow { Eigen::Index index; Rational residual, denominator; };
    std::vector<Eigen::Index> exactIndex(residual.size(), -1);
    std::vector<ExactRow> exactRows;
    std::uint64_t rangeFallbacks = 0;
    for (Eigen::Index i = 0; i < residual.size(); ++i) {
        const auto& row = rows[i];
        Wide sum, scale;
        if (!row.range && finishWide(row.residual, sum) && finishWide(row.denominator, scale) && scale >= 0) {
            writeRow(i, sum, scale == 0 ? 0 : std::abs(sum) / scale);
            continue;
        }
        rangeFallbacks += row.range;
        exactIndex[i] = static_cast<Eigen::Index>(exactRows.size());
        const auto exact = fromWide(residual[i]);
        exactRows.push_back({i, exact, abs(exact)});
    }
    if (!exactRows.empty()) {
        for (Eigen::Index col = 0; col < matrix.outerSize(); ++col) {
            std::optional<Rational> exactDirection;
            for (WideMatrix::InnerIterator entry(matrix, col); entry; ++entry) {
                const auto index = exactIndex[entry.row()];
                if (index < 0)
                    continue;
                if (!exactDirection)
                    exactDirection = fromWide(direction[col]);
                const Rational term = fromWide(entry.value()) * *exactDirection;
                auto& row = exactRows[index];
                row.residual += term;
                row.denominator += abs(term);
            }
        }
        for (const auto& row : exactRows) {
            const auto error = row.denominator == 0 ? 0.L : toWide(Rational(abs(row.residual) / row.denominator));
            writeRow(row.index, toWide(row.residual), error);
        }
    }
    incrementPerformanceCounter("linear.extended_eft_defect.range_fallback_rows", rangeFallbacks);
    incrementPerformanceCounter("linear.extended_eft_defect.cancellation_fallback_rows",
                                exactRows.size() - rangeFallbacks);
    return result;
}
#endif

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

BasicDefect<long double> extendedEftDefect(const Eigen::SparseMatrix<long double>& matrix,
    const Eigen::VectorX<long double>& residual, const Eigen::VectorX<long double>& direction) {
#if defined(VELA_ARITHMETIC_EXTENDED80) && defined(VELA_STRICT_EFT_FP) && defined(VELA_EFT_MXCSR)
    if (matrix.rows() != residual.size() || matrix.cols() != direction.size() ||
        !residual.allFinite() || !direction.allFinite())
        throw std::invalid_argument("Invalid extended refinement defect dimensions or input");
    {
        FloatingEnvironment environment;
        if (environment.supported())
            return extendedEftDefectImpl(matrix, residual, direction);
    }
#endif
    incrementPerformanceCounter("linear.extended_eft_defect.whole_call_fallbacks");
    return computeHighPrecisionDefect(matrix, residual, direction);
}
} // namespace vela::linear_refinement
