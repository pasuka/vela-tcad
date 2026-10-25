#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vela::experimental {

struct JoyceDixonJet {
    double value;
    double first;
    double second;
};

/// Explicit Joyce-Dixon approximation used by fermi_dirac_joyce_dixon.
///
/// The low-density inverse is the four-term Joyce-Dixon expression (1977,
/// doi:10.1063/1.89697). A shifted Sommerfeld expression is joined at the
/// common slope, making the inverse C1. Both eta derivatives and the forward
/// function are derived from this same inverse. The second derivative is
/// piecewise: it must not be differenced across the join.
///
/// This construction is a candidate for a measured native mapping; it does
/// not change Vela's defining-integral Fermi-Dirac default. Apparent FD-BGN
/// must call this inverse too when the candidate is used in a private replay.
class JoyceDixonStatistics {
public:
    struct Join {
        double ratio;
        double eta;
        double shift;
    };

    static const Join& join()
    {
        static const Join value = [] {
            double lo = 8.0, hi = 9.0;
            for (int i = 0; i < 64; ++i) {
                const double mid = lo + (hi - lo) * 0.5;
                if (lowSlope(mid) > highSlope(mid)) lo = mid;
                else hi = mid;
            }
            const double ratio = lo + (hi - lo) * 0.5;
            return Join{ratio, lowInverse(ratio),
                        highInverse(ratio) - lowInverse(ratio)};
        }();
        return value;
    }

    static double inverse(double ratio)
    {
        if (ratio == std::numeric_limits<double>::infinity()) return ratio;
        if (!(ratio > 0.0) || !std::isfinite(ratio))
            throw std::invalid_argument("Joyce-Dixon inverse requires a positive ratio");
        return ratio < join().ratio ? lowInverse(ratio)
                                    : highInverse(ratio) - join().shift;
    }

    static JoyceDixonJet evaluate(double eta)
    {
        if (std::isnan(eta)) return {eta, eta, eta};
        if (eta == -std::numeric_limits<double>::infinity()) return {0., 0., 0.};
        if (eta == std::numeric_limits<double>::infinity()) return {eta, eta, 0.};
        if (eta >= join().eta) {
            const double t = eta + join().shift;
            const double radius = std::hypot(t, std::sqrt(c));
            const double root = std::sqrt(radius);
            const double value = norm * radius * root;
            const double first = 1.5 * norm * t / root;
            const double second = 1.5 * norm / root
                * (0.5 + 0.5 * (c / radius) / radius);
            return {value, first, second};
        }
        if (eta < -36.0) {
            const double value = std::exp(eta);
            return {value, value, value};
        }
        double lo = eta - a1 * join().ratio;
        double hi = std::min(eta, std::log(join().ratio));
        double z = std::clamp(eta - a1 * std::exp(std::min(eta, 0.)), lo, hi);
        for (int i = 0; i < 40; ++i) {
            const double ratio = std::exp(z);
            const double f = z + correction(ratio) - eta;
            if (std::abs(f) <= 4. * std::numeric_limits<double>::epsilon()
                    * std::max(1., std::abs(eta))) break;
            if (f > 0.) hi = z;
            else lo = z;
            const double trial = z - f / logSlope(ratio);
            const double next = trial > lo && trial < hi
                ? trial : lo + (hi - lo) * 0.5;
            if (next == z) break;
            z = next;
        }
        const double ratio = std::exp(z);
        const double d = logSlope(ratio);
        const double dprime = a1 + ratio * (4*a2 + ratio*(9*a3 + ratio*16*a4));
        return {ratio, ratio / d, ratio * (d - ratio*dprime) / (d*d*d)};
    }

private:
    static constexpr double a1 = 0.35355339059327379;
    static constexpr double a2 = -0.0049500897298752622;
    static constexpr double a3 = 1.4838577128872821e-4;
    static constexpr double a4 = -4.4256301190009895e-6;
    static constexpr double norm = 0.75225277806367504926; // 4/(3 sqrt(pi))
    static constexpr double c = 1.64493406684822643647; // pi^2/6

    static double correction(double ratio)
    {
        return ratio*(a1 + ratio*(a2 + ratio*(a3 + ratio*a4)));
    }
    static double logSlope(double ratio)
    {
        return 1. + ratio*(a1 + ratio*(2*a2 + ratio*(3*a3 + ratio*4*a4)));
    }
    static double lowInverse(double ratio)
    {
        return std::log(ratio) + correction(ratio);
    }
    static double lowSlope(double ratio) { return logSlope(ratio) / ratio; }
    static double highInverse(double ratio)
    {
        const double u = std::pow(ratio, 2./3.) * std::pow(1./norm, 2./3.);
        return u * std::sqrt(1. - (c / u) / u);
    }
    static double highSlope(double ratio)
    {
        const double u = std::pow(ratio, 2./3.) * std::pow(1./norm, 2./3.);
        return (2./3.) * (u / ratio) * (u / highInverse(ratio));
    }
};

} // namespace vela::experimental
