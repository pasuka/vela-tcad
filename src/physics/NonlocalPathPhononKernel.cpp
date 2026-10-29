#include "vela/physics/NonlocalPathPhononKernel.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace vela::experimental {
namespace {
using Wide = long double;
constexpr Wide representationTolerance = 128 * std::numeric_limits<Real>::epsilon();
void require(bool ok, const char* message) { if (!ok) throw std::invalid_argument(message); }
bool positive(Real x) { return std::isfinite(x) && x > 0; }
Real narrow(Wide x) {
    const Real result = static_cast<Real>(x);
    require(std::isfinite(result), "unrepresentable phonon kernel diagnostic");
    return result;
}
// Analytic P1 primitives; no endpoint epsilon cutoff or singular sampling.
Wide sqrtIntegral(Wide a, Wide b, Wide length) {
    const Wide sa = std::sqrt(a), sb = std::sqrt(b);
    require(sa + sb > 0, "zero barrier interval in selected Keldysh leg");
    return (2 * length / 3) * (a + sa * sb + b) / (sa + sb);
}
Wide inverseSqrtIntegral(Wide a, Wide b, Wide length) {
    const Wide sum = std::sqrt(a) + std::sqrt(b);
    require(sum > 0, "classically allowed plateau requires a different branch");
    return 2 * length / sum;
}
Wide gapIntegral(Wide a, Wide b, Wide length) {
    if (a == b) return length * std::pow(a, -1.25L);
    const Wide relative = (b - a) / a;
    const Wide logarithm = std::abs(relative) < .5L ? std::log1p(relative) : std::log(b) - std::log(a);
    return length * std::pow(a, -1.25L) * (-4 * std::expm1(-.25L * logarithm) / std::expm1(logarithm));
}
Wide softplus(Wide x) { return std::max(x, Wide(0)) + std::log1p(std::exp(-std::abs(x))); }
Wide logWindow(Wide exponent) {
    require(exponent >= 0 && std::isfinite(exponent), "invalid transverse window");
    return exponent == 0 ? -std::numeric_limits<Wide>::infinity() : std::log(-std::expm1(-exponent));
}
struct BarrierKnot { Wide x, a, b; };
}

TunnelingBandProfile tunnelingProfileFromGeometry(const TunnelingGeometryResult& geometry,
    const TunnelingBandState& state, const NonlocalPathBandToBandConfig& config) {
    state.validate(); config.validate();
    require(geometry.path.termination == TunnelingTermination::Complete,
        "phonon kernel requires resolved complete geometry");
    require(geometry.path.identity == state.identity, "stale phonon geometry/state identity");
    geometry.path.validate(state, config.maximumPathLength_m);
    require(!geometry.segments.empty() && geometry.path.points.size() == geometry.segments.size() + 1,
        "incomplete geometry segment support");
    TunnelingBandProfile result;
    result.identity = state.identity;
    result.initialValenceGradient_J_per_m = geometry.initialValenceGradient_J_per_m;
    auto sample = [&](Real x, const TunnelingEndpoint& support) {
        support.validate(state);
        Wide ev = 0, ec = 0, temperature = 0;
        for (const auto& [node, weight] : support.weights) {
            ev += Wide(weight) * state.nodes[node].valence_J;
            ec += Wide(weight) * state.nodes[node].conduction_J;
            temperature += Wide(weight) * state.nodes[node].temperature_K;
        }
        return TunnelingBandKnot{x, narrow(ev), narrow(ec), narrow(temperature)};
    };
    result.knots.push_back(sample(0, geometry.segments[0].start));
    const Index material = geometry.segments[0].start.material;
    Real previous = 0;
    for (const auto& segment : geometry.segments) {
        require(segment.begin_m == previous && segment.end_m > previous &&
            segment.start.cell == segment.end.cell && segment.start.material == material &&
            segment.end.material == material, "invalid or heterogeneous phonon segment support");
        const auto start = sample(segment.begin_m, segment.start);
        const auto& last = result.knots.back();
        const Wide continuityScale = std::abs(Wide(start.valence_J)) + std::abs(Wide(last.valence_J)) +
            std::abs(Wide(start.conduction_J)) + std::abs(Wide(last.conduction_J));
        require(std::abs(Wide(start.valence_J) - last.valence_J) <= representationTolerance *
            continuityScale &&
            std::abs(Wide(start.conduction_J) - last.conduction_J) <= representationTolerance *
            continuityScale, "discontinuous band support");
        result.knots.push_back(sample(segment.end_m, segment.end));
        previous = segment.end_m;
    }
    return result;
}

PhononTunnelingResult evaluatePhononTunneling(const TunnelingBandProfile& profile,
    const TunnelingReservoirs& reservoirs, const NonlocalPathBandToBandConfig& config) {
    config.validate();
    require(profile.knots.size() >= 2 && profile.knots.front().distance_m == 0 &&
        positive(profile.initialValenceGradient_J_per_m), "invalid phonon band profile or initial gradient");
    require(std::isfinite(reservoirs.electronFermi_J) && std::isfinite(reservoirs.holeFermi_J) &&
        std::isfinite(reservoirs.valenceAvailableEnergy_J) && reservoirs.valenceAvailableEnergy_J >= 0 &&
        std::isfinite(reservoirs.conductionAvailableEnergy_J) && reservoirs.conductionAvailableEnergy_J >= 0,
        "finite Fermi energies and nonnegative explicit reservoir windows required");
    require(profile.knots.size() <= config.maximumIntegrationIntervals,
        "phonon integration interval resource limit exceeded");
    const auto constants = tunnelingPhysicalConstants(config.constants);
    const auto masses = resolveIndirectParameters(config);
    const Wide h = constants.planck_J_s, hbar = h / (2 * std::numbers::pi_v<Wide>);
    const Wide mc = masses.conductionMass_kg, mv = masses.valenceMass_kg;
    const Wide mr = mc / (1 + mc / mv);
    const Wide energy = Wide(profile.knots[0].valence_J) - config.valenceOffset_J;
    const Wide temperature = profile.knots[0].temperature_K;
    require(temperature > 0 && std::isfinite(temperature), "positive physical temperature required");
    const Real t = profile.knots[0].temperature_K;
    const Wide temperatureAllowance = 4 * (std::nextafter(t, std::numeric_limits<Real>::infinity()) - t);
    std::vector<BarrierKnot> knots;
    for (Index i = 0; i < profile.knots.size(); ++i) {
        const auto& p = profile.knots[i];
        require(std::isfinite(p.distance_m) && std::isfinite(p.valence_J) && std::isfinite(p.conduction_J) &&
            positive(p.temperature_K) && std::abs(Wide(p.temperature_K) - temperature) <= temperatureAllowance,
            "nonfinite or nonisothermal phonon profile");
        require(i == 0 || p.distance_m > profile.knots[i - 1].distance_m, "unordered phonon band knots");
        Wide a = energy - p.valence_J + config.valenceOffset_J;
        Wide b = Wide(p.conduction_J) + config.conductionOffset_J - energy;
        const Wide scale = std::abs(energy) + std::abs(Wide(p.valence_J)) + std::abs(Wide(p.conduction_J)) +
            std::abs(Wide(config.conductionOffset_J)) + std::abs(Wide(config.valenceOffset_J));
        const Wide allowance = representationTolerance * scale;
        require(a >= -allowance && b >= -allowance, "classically allowed barrier pocket is not a qualified branch");
        if (i == 0) { require(std::abs(a) <= allowance, "invalid valence turning endpoint"); a = 0; }
        if (i + 1 == profile.knots.size()) { require(std::abs(b) <= allowance, "invalid conduction turning endpoint"); b = 0; }
        a = std::max(a, Wide(0)); b = std::max(b, Wide(0));
        require(a + b > 0, "nonpositive effective tunneling gap");
        knots.push_back({p.distance_m, a, b});
    }
    const Wide length = knots.back().x;
    require(length <= Wide(config.maximumPathLength_m) * (1 + representationTolerance), "phonon path exceeds length limit");
    const Wide massScale = std::max(mc, mv);
    std::vector<Wide> roots;
    for (Index i = 1; i < knots.size(); ++i) {
        const auto& left = knots[i - 1]; const auto& right = knots[i];
        const Wide z0 = mv / massScale * left.a - mc / massScale * left.b;
        const Wide z1 = mv / massScale * right.a - mc / massScale * right.b;
        require(z0 != 0 || z1 != 0, "degenerate Keldysh crossover interval");
        if ((z0 <= 0 && z1 >= 0) || (z0 >= 0 && z1 <= 0)) {
            const Wide root = left.x + (right.x - left.x) * (-z0 / (z1 - z0));
            if (roots.empty() || std::abs(root - roots.back()) > representationTolerance * length) roots.push_back(root);
        }
    }
    require(roots.size() == 1 && roots[0] > 0 && roots[0] < length, "phonon profile requires a unique interior crossover");
    const Wide x0 = roots[0];
    Wide actionV = 0, actionC = 0, reciprocalV = 0, reciprocalC = 0, gapSum = 0;
    for (Index i = 1; i < knots.size(); ++i) {
        const auto& left = knots[i - 1]; const auto& right = knots[i];
        const Wide span = right.x - left.x;
        gapSum += gapIntegral(left.a + left.b, right.a + right.b, span);
        auto integrate = [&](Wide lo, Wide hi, bool valence) {
            if (hi <= lo) return;
            const Wide p0 = valence ? left.a : left.b, p1 = valence ? right.a : right.b;
            const Wide a = std::max(Wide(0), p0 + (p1 - p0) * ((lo - left.x) / span));
            const Wide b = std::max(Wide(0), p0 + (p1 - p0) * ((hi - left.x) / span));
            const Wide k = std::sqrt(2 * (valence ? mv : mc)) / hbar;
            if (valence) {
                actionV += k * sqrtIntegral(a, b, hi - lo);
                reciprocalV += inverseSqrtIntegral(a, b, hi - lo) / k;
            } else {
                actionC += k * sqrtIntegral(a, b, hi - lo);
                reciprocalC += inverseSqrtIntegral(a, b, hi - lo) / k;
            }
        };
        integrate(left.x, std::min(right.x, x0), true);
        integrate(std::max(left.x, x0), right.x, false);
    }
    require(actionV >= 0 && actionC >= 0 && reciprocalV > 0 && reciprocalC > 0 && gapSum > 0,
        "invalid phonon integral");
    const Wide estimated = 128 * std::numeric_limits<Wide>::epsilon() * (knots.size() + 1);
    require(estimated <= config.relativeIntegrationTolerance, "requested phonon integration precision is unattainable");
    const Wide kbT = Wide(constants.boltzmann_J_per_K) * temperature;
    const Wide phononArgument = Wide(config.phononEnergy_J) / kbT;
    const Wide logThermal = phononArgument > 50 ? std::log1p(2 * std::exp(-phononArgument) /
        (-std::expm1(-phononArgument))) : std::log1p(2 / std::expm1(phononArgument));
    const Wide logBase = std::log(Wide(masses.coupling_J2_m_per_kg)) + logThermal -
        std::log(64 * std::numbers::pi_v<Wide> * std::numbers::pi_v<Wide> * Wide(config.phononEnergy_J)) +
        .5L * (std::log(mv) + std::log(mc) - std::log(h) - std::log(length) - .5L * std::log(2 * mr)) +
        std::log(gapSum) - std::log(reciprocalV) - std::log(reciprocalC);
    const Wide logV = logWindow(2 * mv * reservoirs.valenceAvailableEnergy_J / (hbar * hbar) * reciprocalV);
    const Wide logC = logWindow(2 * mc * reservoirs.conductionAvailableEnergy_J / (hbar * hbar) * reciprocalC);
    const Wide delta = (Wide(reservoirs.holeFermi_J) - reservoirs.electronFermi_J) / kbT;
    PhononTunnelingResult result;
    result.crossover_m = narrow(x0); result.action = narrow(actionV + actionC);
    result.valenceReciprocalIntegral_m2 = narrow(reciprocalV);
    result.conductionReciprocalIntegral_m2 = narrow(reciprocalC);
    result.valenceWindowFactor = narrow(std::exp(logV)); result.conductionWindowFactor = narrow(std::exp(logC));
    result.logPrefactor = reservoirs.valenceAvailableEnergy_J == 0 || reservoirs.conductionAvailableEnergy_J == 0 ?
        -std::numeric_limits<Real>::infinity() : narrow(logBase + logV + logC);
    result.estimatedIntegrationRoundoff = narrow(estimated);
    if (delta == 0 || reservoirs.valenceAvailableEnergy_J == 0 || reservoirs.conductionAvailableEnergy_J == 0) {
        result.logAbsoluteGeneration = -std::numeric_limits<Real>::infinity(); return result;
    }
    const Wide high = std::max(Wide(reservoirs.holeFermi_J), Wide(reservoirs.electronFermi_J));
    const Wide low = std::min(Wide(reservoirs.holeFermi_J), Wide(reservoirs.electronFermi_J));
    const Wide logOccupation = std::log(-std::expm1(-std::abs(delta))) -
        softplus((energy - high) / kbT) - softplus((low - energy) / kbT);
    const Wide logRate = std::log(Wide(profile.initialValenceGradient_J_per_m)) + logBase + logV + logC -
        2 * (actionV + actionC) + logOccupation;
    result.logAbsoluteGeneration = narrow(logRate);
    result.sign = delta > 0 ? 1 : -1;
    result.generation_per_m3_s = narrow(result.sign * std::exp(logRate));
    result.underflow = result.generation_per_m3_s == 0;
    return result;
}
} // namespace vela::experimental
