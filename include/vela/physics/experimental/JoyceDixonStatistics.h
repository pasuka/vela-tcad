#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace vela::experimental {

template<class Scalar>
struct BasicJoyceDixonJet {
    Scalar value;
    Scalar first;
    Scalar second;
};

/// Explicit Joyce-Dixon approximation used by fermi_dirac_joyce_dixon.
///
/// The low-density inverse is the four-term Joyce-Dixon expression (1977,
/// doi:10.1063/1.89697). A shifted Sommerfeld expression is joined at the
/// common slope, making the inverse C1. Both eta derivatives and the forward
/// function are derived from this same inverse. The second derivative is
/// piecewise: it must not be differenced across the join.
///
/// Population, transport, and apparent bandgap narrowing can share this law.
template<class Scalar>
class BasicJoyceDixonStatistics {
public:
    struct Join {
        Scalar ratio;
        Scalar eta;
        Scalar shift;
    };

    static const Join& join()
    {
        static const Join value = [] {
            Scalar lo = 8.0, hi = 9.0;
            for (int i = 0; i < 64; ++i) {
                const Scalar mid = lo + (hi - lo) * 0.5;
                if (lowSlope(mid) > highSlope(mid)) lo = mid;
                else hi = mid;
            }
            const Scalar ratio = lo + (hi - lo) * 0.5;
            return Join{ratio, lowInverse(ratio),
                        highInverse(ratio) - lowInverse(ratio)};
        }();
        return value;
    }

    static Scalar inverse(Scalar ratio)
    {
        if (ratio == std::numeric_limits<Scalar>::infinity()) return ratio;
        if (!(ratio > 0.0) || !std::isfinite(ratio))
            throw std::invalid_argument("Joyce-Dixon inverse requires a positive ratio");
        return ratio < join().ratio ? lowInverse(ratio)
                                    : highInverse(ratio) - join().shift;
    }

    // Explicit stable primitive contract; existing inverse/forward dispatch is unchanged.
    // g=log(r)-eta(r) and dg/dr have finite nondegenerate limits at r=0.
    static Scalar primitiveLogGamma(Scalar ratio) {
        if (!(ratio>=0.) || !std::isfinite(ratio))
            throw std::invalid_argument("JD gamma requires a nonnegative finite ratio");
        return ratio<join().ratio ? -correction(ratio) : std::log(ratio)-inverse(ratio);
    }
    static Scalar primitiveLogGammaSlope(Scalar ratio) {
        if (!(ratio>=0.) || !std::isfinite(ratio))
            throw std::invalid_argument("JD gamma slope requires a nonnegative finite ratio");
        if(ratio<join().ratio)return -(a1+ratio*(2*a2+ratio*(3*a3+ratio*4*a4)));
        return 1./ratio-highSlope(ratio);
    }

    // Scale before the nondegenerate exponential reaches the subnormal
    // range. DOS*exp(eta) can be normal although exp(eta) has lost digits
    // or underflowed. These are the same low-density JD law, without a floor.
    static Scalar density(Scalar eta,Scalar dos)
    {
        if(eta<std::log(std::numeric_limits<Scalar>::min()))
            return std::exp(std::log(dos)+eta);
        return dos*evaluate(eta).value;
    }
    static Scalar densityDerivativeEta(Scalar eta,Scalar dos)
    {
        if(eta<std::log(std::numeric_limits<Scalar>::min()))return density(eta,dos);
        return dos*evaluate(eta).first;
    }
    static Scalar inverseDensity(Scalar population,Scalar dos)
    {
        if(!(population>0.) || !(dos>0.) || !std::isfinite(population) || !std::isfinite(dos))
            throw std::invalid_argument("Joyce-Dixon density inverse requires positive finite population and DOS");
        const Scalar ratio=population/dos;
        if(ratio<std::numeric_limits<Scalar>::min())return std::log(population)-std::log(dos);
        return inverse(ratio);
    }
    static Scalar densitySlopeFromPopulation(Scalar population,Scalar dos)
    {
        if(population/dos<std::numeric_limits<Scalar>::min())return population;
        return dos*evaluate(inverseDensity(population,dos)).first;
    }

    static BasicJoyceDixonJet<Scalar> evaluate(Scalar eta)
    {
        if (std::isnan(eta)) return {eta, eta, eta};
        if (eta == -std::numeric_limits<Scalar>::infinity()) return {0., 0., 0.};
        if (eta == std::numeric_limits<Scalar>::infinity()) return {eta, eta, 0.};
        if (eta >= join().eta) {
            const Scalar t = eta + join().shift;
            const Scalar radius = std::hypot(t, std::sqrt(c));
            const Scalar root = std::sqrt(radius);
            const Scalar value = norm * radius * root;
            const Scalar first = 1.5 * norm * t / root;
            const Scalar second = 1.5 * norm / root
                * (0.5 + 0.5 * (c / radius) / radius);
            return {value, first, second};
        }
        if (eta < (std::is_same_v<Scalar, double> ? -36 : -46)) {
            const Scalar value = std::exp(eta);
            return {value, value, value};
        }
        Scalar lo = eta - a1 * join().ratio;
        Scalar hi = std::min(eta, std::log(join().ratio));
        Scalar z = std::clamp(eta - a1 * std::exp(std::min(eta, Scalar{0})), lo, hi);
        for (int i = 0; i < 40; ++i) {
            const Scalar ratio = std::exp(z);
            const Scalar f = z + correction(ratio) - eta;
            if (std::abs(f) <= 4. * std::numeric_limits<Scalar>::epsilon()
                    * std::max(Scalar{1}, std::abs(eta))) break;
            if (f > 0.) hi = z;
            else lo = z;
            const Scalar trial = z - f / logSlope(ratio);
            const Scalar next = trial > lo && trial < hi
                ? trial : lo + (hi - lo) * 0.5;
            if (next == z) break;
            z = next;
        }
        const Scalar ratio = std::exp(z);
        const Scalar d = logSlope(ratio);
        const Scalar dprime = a1 + ratio * (4*a2 + ratio*(9*a3 + ratio*16*a4));
        return {ratio, ratio / d, ratio * (d - ratio*dprime) / (d*d*d)};
    }

private:
    // Correctly rounded from high-precision reversion of the F_{1/2} fugacity
    // series. Computing high-order coefficients by cancelling doubles changes
    // the dense-branch C1 shift by about 2e-14.
    static constexpr Scalar a1 = std::is_same_v<Scalar, double>
        ? static_cast<Scalar>(0.3535533905932738) : static_cast<Scalar>(0.35355339059327376220042218105242452L);
    static constexpr Scalar a2 = std::is_same_v<Scalar, double>
        ? static_cast<Scalar>(-0.004950089729875255) : static_cast<Scalar>(-0.004950089729875254836382926833985819L);
    static constexpr Scalar a3 = std::is_same_v<Scalar, double>
        ? static_cast<Scalar>(0.00014838577128872333) : static_cast<Scalar>(0.0001483857712887233393665671378556174L);
    static constexpr Scalar a4 = std::is_same_v<Scalar, double>
        ? static_cast<Scalar>(-4.425630118996707e-06) : static_cast<Scalar>(-0.000004425630118996706728478812616541828L);
    static constexpr Scalar norm = std::is_same_v<Scalar, double>
        ? static_cast<Scalar>(0.75225277806367504926) : static_cast<Scalar>(0.75225277806367504926410593541436345L); // 4/(3 sqrt(pi))
    static constexpr Scalar c = std::is_same_v<Scalar, double>
        ? static_cast<Scalar>(1.64493406684822643647) : static_cast<Scalar>(1.644934066848226436472415166646025189L); // pi^2/6

    static Scalar correction(Scalar ratio)
    {
        return ratio*(a1 + ratio*(a2 + ratio*(a3 + ratio*a4)));
    }
    static Scalar logSlope(Scalar ratio)
    {
        return 1. + ratio*(a1 + ratio*(2*a2 + ratio*(3*a3 + ratio*4*a4)));
    }
    static Scalar lowInverse(Scalar ratio)
    {
        return std::log(ratio) + correction(ratio);
    }
    static Scalar lowSlope(Scalar ratio) { return logSlope(ratio) / ratio; }
    static Scalar highInverse(Scalar ratio)
    {
        const Scalar u = std::pow(ratio, Scalar{2}/3) * std::pow(Scalar{1}/norm, Scalar{2}/3);
        return u * std::sqrt(1. - (c / u) / u);
    }
    static Scalar highSlope(Scalar ratio)
    {
        const Scalar u = std::pow(ratio, Scalar{2}/3) * std::pow(Scalar{1}/norm, Scalar{2}/3);
        return (Scalar{2}/3) * (u / ratio) * (u / highInverse(ratio));
    }
};

using JoyceDixonJet = BasicJoyceDixonJet<double>;
using JoyceDixonStatistics = BasicJoyceDixonStatistics<double>;

} // namespace vela::experimental
