#include "vela/physics/NonlocalPathBandToBand.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>
#include <tuple>

namespace vela::experimental {
namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
bool positive(Real x) { return std::isfinite(x) && x > 0; }
bool nonnegative(Real x) { return std::isfinite(x) && x >= 0; }
void keys(const nlohmann::json& j, std::initializer_list<const char*> allowed) {
    check(j.is_object(), "NonLocalPath private config must be an object");
    for (auto it = j.begin(); it != j.end(); ++it) {
        check(std::find(allowed.begin(), allowed.end(), it.key()) != allowed.end(),
              "unknown NonLocalPath private config key or unit");
    }
}
Index unsignedIndex(const nlohmann::json& j) {
    check(j.is_number_unsigned(), "NonLocalPath indices must be unsigned integers");
    return j.get<Index>();
}
struct PhysicalConstants { Real q, h, kb, m0; };
PhysicalConstants physical(TunnelingConstants c) {
    if (c == TunnelingConstants::Codata2018)
        return {1.602176634e-19, 6.62607015e-34, 1.380649e-23, 9.1093837015e-31};
    check(c == TunnelingConstants::SentaurusT202203, "unknown constant profile");
    return {1.602192e-19, 6.626176e-34, 1.380662e-23, 9.109534e-31};
}
long double sum(const std::vector<Real>& x) {
    long double total = 0;
    for (Real v : x) { check(std::isfinite(v), "nonfinite source integral"); total += v; }
    return total;
}
}
Real tunnelingElementaryCharge(TunnelingConstants c) { return physical(c).q; }
TunnelingPhysicalConstants tunnelingPhysicalConstants(TunnelingConstants c) {
    const auto k = physical(c);
    return {k.q, k.h, k.kb, k.m0};
}

void NonlocalPathBandToBandConfig::validate(TunnelingSolveMode mode) const {
    check(mode == TunnelingSolveMode::CoupledDC, "NonLocalPath MVP requires coupled DC");
    check(material == "Si", "NonLocalPath N1 parameter contract supports indirect Si only");
    check(!regions.empty(), "explicit NonLocalPath regions required");
    check(std::set<Index>(regions.begin(), regions.end()).size() == regions.size(), "duplicate region");
    check(ab.has_value() != microscopic.has_value(), "exactly one AB or microscopic mode required");
    check(positive(phononEnergy_J) && positive(maximumPathLength_m), "positive phonon and length required");
    check(std::isfinite(conductionOffset_J) && std::isfinite(valenceOffset_J), "nonfinite band offset");
    check(positive(relativeIntegrationTolerance) && relativeIntegrationTolerance < 1 &&
          maximumIntegrationIntervals > 0, "invalid integration control");
    physical(constants);
    check(derivatives == TunnelingDerivativeMode::NativeDefault ||
          derivatives == TunnelingDerivativeMode::NativeExplicit ||
          derivatives == TunnelingDerivativeMode::FullMathematical, "unknown derivative contract");
    if (ab) check(positive(ab->A_per_m3_s) && positive(ab->B_V_per_m) &&
                  positive(ab->referenceGap_J) && nonnegative(ab->massRatio), "invalid AB parameters");
    if (microscopic) check(positive(microscopic->conductionMass_kg) &&
                           positive(microscopic->valenceMass_kg) &&
                           positive(microscopic->coupling_J2_m_per_kg), "invalid microscopic parameters");
}

IndirectMicroscopicParameters resolveIndirectParameters(const NonlocalPathBandToBandConfig& c) {
    c.validate();
    if (c.microscopic) return *c.microscopic;
    const auto k = physical(c.constants);
    const auto& a = *c.ab;
    const Real mr = std::pow(3 * k.q * k.h * a.B_V_per_m /
        (std::pow(2., 3.5) * std::numbers::pi * std::pow(a.referenceGap_J, 1.5)), 2);
    Real mc, mv;
    if (a.massRatio == 0) {
        check(2 * mr < k.m0, "automatic valence mass is nonpositive");
        mc = 1 / (1 / (2 * mr) + 1 / k.m0);
        mv = 1 / (1 / (2 * mr) - 1 / k.m0);
    } else {
        mc = mr * (1 + 1 / a.massRatio);
        mv = mr * (1 + a.massRatio);
    }
    const Real phonons = 1 / std::expm1(c.phononEnergy_J / (k.kb * 300));
    const Real coupling = a.A_per_m3_s * std::pow(2., 5.25) * std::pow(k.h, 2.5) *
        std::pow(mr, 1.25) * c.phononEnergy_J * std::pow(a.referenceGap_J, 1.75) /
        (std::pow(mc * mv, 1.5) * (1 + 2 * phonons) * std::pow(k.q * 100, 2.5));
    check(positive(mc) && positive(mv) && positive(coupling), "unrepresentable indirect parameters");
    return {mc, mv, coupling};
}

NonlocalPathBandToBandConfig nonlocalPathConfigFromJson(const nlohmann::json& j) {
    keys(j, {"material", "regions", "ab", "microscopic", "phonon_energy_J",
         "conduction_offset_J", "valence_offset_J", "maximum_path_length_m",
         "relative_integration_tolerance", "maximum_integration_intervals", "constants", "derivatives"});
    NonlocalPathBandToBandConfig c;
    c.material = j.at("material").get<std::string>();
    check(j.at("regions").is_array(), "regions must be an array");
    for (const auto& region : j.at("regions")) c.regions.push_back(unsignedIndex(region));
    if (j.contains("ab")) {
        const auto& a = j.at("ab");
        keys(a, {"A_per_m3_s", "B_V_per_m", "reference_gap_J", "mass_ratio"});
        c.ab = IndirectABParameters{a.at("A_per_m3_s").get<Real>(), a.at("B_V_per_m").get<Real>(),
             a.at("reference_gap_J").get<Real>(), a.value("mass_ratio", 0.)};
    }
    if (j.contains("microscopic")) {
        const auto& m = j.at("microscopic");
        keys(m, {"conduction_mass_kg", "valence_mass_kg", "coupling_J2_m_per_kg"});
        c.microscopic = IndirectMicroscopicParameters{m.at("conduction_mass_kg").get<Real>(),
            m.at("valence_mass_kg").get<Real>(), m.at("coupling_J2_m_per_kg").get<Real>()};
    }
    c.phononEnergy_J = j.at("phonon_energy_J").get<Real>();
    c.maximumPathLength_m = j.at("maximum_path_length_m").get<Real>();
    c.conductionOffset_J = j.value("conduction_offset_J", 0.);
    c.valenceOffset_J = j.value("valence_offset_J", 0.);
    c.relativeIntegrationTolerance = j.value("relative_integration_tolerance", 1e-8);
    if (j.contains("maximum_integration_intervals"))
        c.maximumIntegrationIntervals = unsignedIndex(j.at("maximum_integration_intervals"));
    const auto constants = j.at("constants").get<std::string>();
    check(constants == "codata2018" || constants == "sentaurus_t2022.03", "unknown constants");
    c.constants = constants == "codata2018" ? TunnelingConstants::Codata2018 : TunnelingConstants::SentaurusT202203;
    const auto d = j.at("derivatives").get<std::string>();
    check(d == "native_default" || d == "native_explicit" || d == "full_mathematical", "unknown derivatives");
    c.derivatives = d == "native_default" ? TunnelingDerivativeMode::NativeDefault :
        d == "native_explicit" ? TunnelingDerivativeMode::NativeExplicit : TunnelingDerivativeMode::FullMathematical;
    c.validate();
    return c;
}
nlohmann::json nonlocalPathConfigToJson(const NonlocalPathBandToBandConfig& c) {
    c.validate();
    nlohmann::json j = {{"material", c.material}, {"regions", c.regions},
        {"phonon_energy_J", c.phononEnergy_J}, {"maximum_path_length_m", c.maximumPathLength_m},
        {"conduction_offset_J", c.conductionOffset_J}, {"valence_offset_J", c.valenceOffset_J},
        {"relative_integration_tolerance", c.relativeIntegrationTolerance},
        {"maximum_integration_intervals", c.maximumIntegrationIntervals},
        {"constants", c.constants == TunnelingConstants::Codata2018 ? "codata2018" : "sentaurus_t2022.03"},
        {"derivatives", c.derivatives == TunnelingDerivativeMode::NativeDefault ? "native_default" :
            c.derivatives == TunnelingDerivativeMode::NativeExplicit ? "native_explicit" : "full_mathematical"}};
    if (c.ab) j["ab"] = {{"A_per_m3_s", c.ab->A_per_m3_s}, {"B_V_per_m", c.ab->B_V_per_m},
        {"reference_gap_J", c.ab->referenceGap_J}, {"mass_ratio", c.ab->massRatio}};
    if (c.microscopic) j["microscopic"] = {{"conduction_mass_kg", c.microscopic->conductionMass_kg},
        {"valence_mass_kg", c.microscopic->valenceMass_kg}, {"coupling_J2_m_per_kg", c.microscopic->coupling_J2_m_per_kg}};
    return j;
}

void TunnelingBandState::validate() const {
    check(!nodes.empty(), "empty tunneling state");
    for (const auto& n : nodes) {
        std::set<Index> regions{n.region};
        for (Index r : n.additionalRegions) check(regions.insert(r).second, "duplicate node region membership");
        check(std::isfinite(n.conduction_J) && std::isfinite(n.valence_J) &&
            std::isfinite(n.electronFermi_J) && std::isfinite(n.holeFermi_J) && std::isfinite(n.narrowing_J) &&
            n.conductionGradient_J_per_m.allFinite() && n.valenceGradient_J_per_m.allFinite() &&
            positive(n.temperature_K) && nonnegative(n.electrons_per_m3) && nonnegative(n.holes_per_m3) &&
            nonnegative(n.conductionDOS_per_m3) && nonnegative(n.valenceDOS_per_m3), "invalid band state");
        if (n.semiconductor) check(n.conduction_J > n.valence_J && positive(n.conductionDOS_per_m3) &&
            positive(n.valenceDOS_per_m3), "semiconductor bands and DOS required");
        std::set<std::pair<Index, int>> seen;
        for (const auto& d : n.derivatives) {
            check(d.column.node < nodes.size() && static_cast<int>(d.column.variable) >= 0 &&
                static_cast<int>(d.column.variable) <= static_cast<int>(TunnelingVariable::Temperature_K) &&
                seen.emplace(d.column.node, static_cast<int>(d.column.variable)).second &&
                std::isfinite(d.conduction_J) && std::isfinite(d.valence_J) &&
                std::isfinite(d.electronFermi_J) && std::isfinite(d.holeFermi_J) && std::isfinite(d.narrowing_J) &&
                std::isfinite(d.electrons_per_m3) && std::isfinite(d.holes_per_m3) &&
                std::isfinite(d.conductionDOS_per_m3) && std::isfinite(d.valenceDOS_per_m3) &&
                std::isfinite(d.temperature_K) &&
                d.conductionGradient_J_per_m.allFinite() && d.valenceGradient_J_per_m.allFinite(), "invalid band derivative");
        }
    }
}
void TunnelingEndpoint::validate(const TunnelingBandState& s) const {
    check(!weights.empty(), "empty endpoint support");
    long double total = 0;
    std::set<Index> seen;
    for (const auto& [node, weight] : weights) {
        check(node < s.nodes.size() && nonnegative(weight) && seen.insert(node).second, "invalid endpoint weight");
        const auto& n = s.nodes[node];
        check(n.semiconductor && n.material == material && (n.region == region ||
            std::find(n.additionalRegions.begin(), n.additionalRegions.end(), region) != n.additionalRegions.end()),
            "endpoint outside selected semiconductor");
        total += weight;
    }
    check(std::abs(total - 1) <= 32 * std::numeric_limits<Real>::epsilon(), "endpoint weights must sum to one");
}
void TunnelingPath::validate(const TunnelingBandState& s, Real maximumLength_m) const {
    check(identity == s.identity, "stale tunneling path");
    check(positive(maximumLength_m) && !points.empty() &&
        (termination != TunnelingTermination::Complete || points.size() >= 2) &&
        std::isfinite(energyLower_J) && std::isfinite(energyUpper_J) && energyLower_J <= energyUpper_J,
        "invalid path bounds");
    check(static_cast<int>(termination) >= 0 && static_cast<int>(termination) <=
        static_cast<int>(TunnelingTermination::NoBarrier), "invalid path termination");
    Real length = 0;
    for (Index i = 0; i < points.size(); ++i) {
        check(points[i].position_m.allFinite(), "nonfinite path point");
        points[i].support.validate(s);
        if (i) {
            const Real step = (points[i].position_m - points[i - 1].position_m).norm();
            check(positive(step), "degenerate path segment");
            length += step;
        }
    }
    check(std::isfinite(length) && length <= maximumLength_m * (1 + 32 * std::numeric_limits<Real>::epsilon()),
        "path exceeds maximum length");
    Real previous = 0;
    for (const auto& [lo, hi] : integrationIntervals_m) {
        check(nonnegative(lo) && positive(hi) && lo >= previous && hi > lo && hi <= length,
            "invalid path integration interval");
        previous = hi;
    }
}

NonlocalPairSources distributePairEvents(const TunnelingBandState& s, const std::vector<NonlocalPairEvent>& events) {
    s.validate();
    NonlocalPairSources result;
    result.identity = s.identity;
    result.electronIntegrals_per_m_s.assign(s.nodes.size(), 0);
    result.holeIntegrals_per_m_s.assign(s.nodes.size(), 0);
    std::set<std::uint64_t> ids;
    for (const auto& event : events) {
        check(ids.insert(event.pathId).second, "duplicate pair event path ID");
        check(std::isfinite(event.signedRate_per_m_s), "nonfinite pair event");
        event.holeStart.validate(s);
        event.electronEnd.validate(s);
        for (const auto& [node, weight] : event.holeStart.weights)
            result.holeIntegrals_per_m_s[node] += event.signedRate_per_m_s * weight;
        for (const auto& [node, weight] : event.electronEnd.weights)
            result.electronIntegrals_per_m_s[node] += event.signedRate_per_m_s * weight;
    }
    result.ledger = events;
    check(result.balanceError() <= 1e-12, "pair transfer conservation failure");
    return result;
}
void NonlocalPairSources::validate(const TunnelingBandState& s) const {
    check(identity == s.identity, "stale nonlocal pair sources");
    check(electronIntegrals_per_m_s.size() == s.nodes.size() &&
        holeIntegrals_per_m_s.size() == s.nodes.size(), "nonlocal source size mismatch");
    check(balanceError() <= 1e-12, "pair transfer conservation failure");
    const auto expected = distributePairEvents(s, ledger);
    long double absolute = 0, discrepancy = 0;
    for (const auto& event : ledger) absolute += std::abs(static_cast<long double>(event.signedRate_per_m_s));
    for (Index i = 0; i < s.nodes.size(); ++i) {
        discrepancy += std::abs(static_cast<long double>(electronIntegrals_per_m_s[i]) - expected.electronIntegrals_per_m_s[i]);
        discrepancy += std::abs(static_cast<long double>(holeIntegrals_per_m_s[i]) - expected.holeIntegrals_per_m_s[i]);
    }
    check(absolute == 0 ? discrepancy == 0 : discrepancy / absolute <= 1e-12,
        "source distribution differs from event ledger");
    check(evaluatedDerivativeMode.has_value() || derivatives.empty(), "unlabeled source derivatives");
    if (evaluatedDerivativeMode) check(*evaluatedDerivativeMode == TunnelingDerivativeMode::NativeDefault ||
        *evaluatedDerivativeMode == TunnelingDerivativeMode::NativeExplicit ||
        *evaluatedDerivativeMode == TunnelingDerivativeMode::FullMathematical, "unknown evaluated derivative mode");
    std::set<std::tuple<bool, Index, Index, int>> seen;
    for (const auto& d : derivatives) {
        check(d.rowNode < s.nodes.size() && d.column.node < s.nodes.size() &&
            s.nodes[d.rowNode].semiconductor && std::isfinite(d.integralDerivative) &&
            static_cast<int>(d.column.variable) >= 0 &&
            static_cast<int>(d.column.variable) <= static_cast<int>(TunnelingVariable::Temperature_K) &&
            seen.emplace(d.electron, d.rowNode, d.column.node, static_cast<int>(d.column.variable)).second,
            "invalid source derivative support");
    }
}
Real NonlocalPairSources::balanceError() const {
    long double absolute = 0;
    for (const auto& e : ledger) {
        check(std::isfinite(e.signedRate_per_m_s), "nonfinite pair ledger");
        absolute += std::abs(static_cast<long double>(e.signedRate_per_m_s));
    }
    const auto e = sum(electronIntegrals_per_m_s), h = sum(holeIntegrals_per_m_s);
    if (absolute == 0) {
        check(std::all_of(electronIntegrals_per_m_s.begin(), electronIntegrals_per_m_s.end(), [](Real x) { return x == 0; }) &&
            std::all_of(holeIntegrals_per_m_s.begin(), holeIntegrals_per_m_s.end(), [](Real x) { return x == 0; }),
            "nonzero source without pair events");
        return 0;
    }
    return static_cast<Real>(std::abs(e - h) / absolute);
}
std::vector<Real> NonlocalPairSources::volumeRates(bool electron, const std::vector<Real>& areas) const {
    const auto& integrals = electron ? electronIntegrals_per_m_s : holeIntegrals_per_m_s;
    check(areas.size() == integrals.size(), "source area size mismatch");
    std::vector<Real> rates(areas.size());
    for (Index i = 0; i < areas.size(); ++i) {
        check(nonnegative(areas[i]) && std::isfinite(integrals[i]), "nonnegative node area required");
        check(areas[i] != 0 || integrals[i] == 0, "nonzero source on zero-area node");
        rates[i] = areas[i] == 0 ? 0 : integrals[i] / areas[i];
        check(std::isfinite(rates[i]), "unrepresentable volume rate");
    }
    return rates;
}
Real NonlocalPairSources::currentPerUm(bool electron, Real factor, TunnelingConstants c) const {
    check(positive(factor), "positive AreaFactor required");
    const Real result = static_cast<Real>(sum(electron ? electronIntegrals_per_m_s : holeIntegrals_per_m_s) *
        tunnelingElementaryCharge(c) * 1e-6L * factor);
    check(std::isfinite(result), "unrepresentable source current");
    return result;
}
bool NonlocalJacobianSupport::replace(Index rows, Index columns, std::vector<Entry> entries) {
    check(rows > 0 && columns > 0, "empty Jacobian dimensions");
    for (const auto& [r, c] : entries) check(r < rows && c < columns, "Jacobian support out of range");
    std::sort(entries.begin(), entries.end());
    entries.erase(std::unique(entries.begin(), entries.end()), entries.end());
    if (rows == rows_ && columns == columns_ && entries == entries_) return false;
    rows_ = rows; columns_ = columns; entries_ = std::move(entries); ++revision_;
    return true;
}
void NonlocalJacobianSupport::require(Index row, Index column) const {
    if (!std::binary_search(entries_.begin(), entries_.end(), Entry{row, column}))
        throw std::logic_error("NonLocalPath sparse support missing: rebuild symbolic structure");
}
} // namespace vela::experimental
