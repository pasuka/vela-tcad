#include "vela/physics/NonlocalPathBandToBand.h"
#include "vela/physics/BandToBandTunnelingModel.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>
#include <cmath>
#include <limits>

using namespace vela;
using namespace vela::experimental;
using Catch::Approx;
namespace {
NonlocalPathBandToBandConfig config() {
    NonlocalPathBandToBandConfig c;
    c.constants = TunnelingConstants::SentaurusT202203;
    const Real q = tunnelingElementaryCharge(c.constants);
    c.regions = {1};
    c.ab = IndirectABParameters{4e20, 1.9e9, 1.1082092307692308 * q};
    c.phononEnergy_J = .037 * q;
    c.maximumPathLength_m = 1e-5;
    return c;
}
TunnelingBandState state() {
    TunnelingBandState s;
    s.identity.trial = 7;
    s.nodes.resize(4);
    for (auto& n : s.nodes) {
        n.material = 0; n.region = 1;
        n.conduction_J = 1.6e-19;
        n.conductionDOS_per_m3 = 2.8e25;
        n.valenceDOS_per_m3 = 1.04e25;
    }
    s.nodes[3].contact = true;
    return s;
}
TunnelingEndpoint endpoint(Index node) { return {0, 0, 1, {{node, 1}}}; }
}

TEST_CASE("NLP N1 strict SI configuration round trips both input modes", "[nonlocalpath][contract]") {
    auto c = config();
    auto j = nonlocalPathConfigToJson(c);
    REQUIRE(nonlocalPathConfigToJson(nonlocalPathConfigFromJson(j)) == j);
    c.microscopic = resolveIndirectParameters(c);
    c.ab.reset();
    c.conductionOffset_J = -1e-21;
    j = nonlocalPathConfigToJson(c);
    REQUIRE(nonlocalPathConfigToJson(nonlocalPathConfigFromJson(j)) == j);
    REQUIRE(resolveIndirectParameters(c).conductionMass_kg == c.microscopic->conductionMass_kg);
    for (const auto* mode : {"native_default", "native_explicit", "full_mathematical"}) {
        j["derivatives"] = mode;
        REQUIRE(nonlocalPathConfigToJson(nonlocalPathConfigFromJson(j)) == j);
    }
}
TEST_CASE("NLP N1 rejects ambiguous units invalid parameters and unsupported solvers", "[nonlocalpath][contract]") {
    const auto c = config();
    for (const auto* key : {"Apath", "maximum_path_length_cm", "provenance", "strategy", "temperature_K"}) {
        auto j = nonlocalPathConfigToJson(c); j[key] = 1;
        REQUIRE_THROWS(nonlocalPathConfigFromJson(j));
    }
    auto j = nonlocalPathConfigToJson(c);
    j["ab"]["A_cm3_s"] = 4e14;
    REQUIRE_THROWS(nonlocalPathConfigFromJson(j));
    for (auto mode : {TunnelingSolveMode::Gummel, TunnelingSolveMode::Transient,
                      TunnelingSolveMode::AC, TunnelingSolveMode::Noise}) REQUIRE_THROWS(c.validate(mode));
    for (Real v : {0., -1., std::numeric_limits<Real>::infinity(), std::nan("")}) {
        auto bad = c; bad.phononEnergy_J = v; REQUIRE_THROWS(bad.validate());
        bad = c; bad.ab->referenceGap_J = v; REQUIRE_THROWS(bad.validate());
    }
    auto bad = c; bad.microscopic = resolveIndirectParameters(c); REQUIRE_THROWS(bad.validate());
    bad = c; bad.ab.reset(); REQUIRE_THROWS(bad.validate());
    bad = c; bad.regions = {1, 1}; REQUIRE_THROWS(bad.validate());
    bad = c; bad.material = "Ge"; REQUIRE_THROWS(bad.validate());
    bad = c; bad.ab->B_V_per_m = 1e12; REQUIRE_THROWS(resolveIndirectParameters(bad));
    j = nonlocalPathConfigToJson(c); j["regions"] = {-1}; REQUIRE_THROWS(nonlocalPathConfigFromJson(j));
    j = nonlocalPathConfigToJson(c); j["maximum_integration_intervals"] = 2.5;
    REQUIRE_THROWS(nonlocalPathConfigFromJson(j));
}
TEST_CASE("NLP N1 parameter conversion matches independently qualified N0 SI values", "[nonlocalpath][contract]") {
    auto c = config();
    const auto m = resolveIndirectParameters(c);
    // N0 independently verified immutable native AB/microscopic source controls.
    REQUIRE(m.conductionMass_kg / 9.109534e-31 == Approx(.10208356687911158).epsilon(2e-13));
    REQUIRE(m.valenceMass_kg / 9.109534e-31 == Approx(.12827262005550075).epsilon(2e-13));
    REQUIRE(m.coupling_J2_m_per_kg / .01 == Approx(6.456560487646502e-17).epsilon(2e-13));
    c.ab->A_per_m3_s *= 2;
    const auto twice = resolveIndirectParameters(c);
    REQUIRE(twice.conductionMass_kg == m.conductionMass_kg);
    REQUIRE(twice.coupling_J2_m_per_kg == 2 * m.coupling_J2_m_per_kg);
    c.ab->massRatio = 2;
    const auto ratio = resolveIndirectParameters(c);
    REQUIRE(ratio.valenceMass_kg / ratio.conductionMass_kg == Approx(2));
    c.ab->B_V_per_m *= 2;
    REQUIRE(resolveIndirectParameters(c).conductionMass_kg / ratio.conductionMass_kg == Approx(4));
    c.constants = TunnelingConstants::Codata2018;
    REQUIRE(tunnelingElementaryCharge(c.constants) == 1.602176634e-19);
    REQUIRE(resolveIndirectParameters(c).conductionMass_kg != twice.conductionMass_kg);
}
TEST_CASE("NLP N1 actual band state preserves gauge gradients and dependency identity", "[nonlocalpath][contract]") {
    auto s = state();
    s.nodes[0].derivatives.push_back({{3, TunnelingVariable::HoleQF_J}, 0, 0, 0, 1, 0});
    REQUIRE_NOTHROW(s.validate());
    for (auto& n : s.nodes) {
        n.conduction_J += 1e-19; n.valence_J += 1e-19;
        n.electronFermi_J += 1e-19; n.holeFermi_J += 1e-19;
    }
    REQUIRE_NOTHROW(s.validate());
    auto bad = s; bad.nodes[0].temperature_K = 0; REQUIRE_THROWS(bad.validate());
    bad = s; bad.nodes[0].derivatives[0].column.node = 4; REQUIRE_THROWS(bad.validate());
    bad = s; bad.nodes[0].derivatives.push_back(bad.nodes[0].derivatives[0]); REQUIRE_THROWS(bad.validate());
    bad = s; bad.nodes[0].valenceGradient_J_per_m[0] = std::nan(""); REQUIRE_THROWS(bad.validate());
    bad = s; bad.nodes[0].derivatives[0].conductionDOS_per_m3 = std::nan(""); REQUIRE_THROWS(bad.validate());
    auto changed = s.identity;
    ++changed.trial; REQUIRE(changed != s.identity);
    changed = s.identity; ++changed.contacts; REQUIRE(changed != s.identity);
    changed = s.identity; ++changed.population; REQUIRE(changed != s.identity);
    changed = s.identity; ++changed.temperature; REQUIRE(changed != s.identity);
}
TEST_CASE("NLP N1 endpoint transfer separates spatial sources and retains contact events", "[nonlocalpath][contract]") {
    auto s = state();
    TunnelingEndpoint end{1, 0, 1, {{2, .25}, {3, .75}}};
    auto out = distributePairEvents(s, {{1, 8, endpoint(0), end}, {2, -2, endpoint(1), endpoint(3)}});
    REQUIRE(out.electronIntegrals_per_m_s == std::vector<Real>{0, 0, 2, 4});
    REQUIRE(out.holeIntegrals_per_m_s == std::vector<Real>{8, -2, 0, 0});
    REQUIRE(out.balanceError() == 0);
    REQUIRE(out.electronIntegrals_per_m_s[3] == 4); // Dirichlet ledger is retained.
    REQUIRE(out.volumeRates(true, {1, 2, 4, 8}) == std::vector<Real>{0, 0, .5, .5});
    REQUIRE(out.currentPerUm(true, 1, TunnelingConstants::Codata2018) / 1.602176634e-19 == Approx(6e-6));
    REQUIRE(out.currentPerUm(false, 3, TunnelingConstants::Codata2018) ==
        Approx(3 * out.currentPerUm(true, 1, TunnelingConstants::Codata2018)).epsilon(1e-14).margin(0));
    REQUIRE_THROWS(out.volumeRates(true, {1, 1, 0, 1}));
    REQUIRE_THROWS(out.currentPerUm(true, -1, TunnelingConstants::Codata2018));
    REQUIRE_THROWS(distributePairEvents(s, {{1, 8, endpoint(0), end}, {1, 1, endpoint(0), end}}));
    end.weights[0].second = .3; REQUIRE_THROWS(distributePairEvents(s, {{1, 8, endpoint(0), end}}));
    s.nodes[3].semiconductor = false;
    REQUIRE_THROWS(distributePairEvents(s, {{1, 1, endpoint(0), endpoint(3)}}));
}
TEST_CASE("NLP N1 signed tiny events use an absolute event budget without a unit floor", "[nonlocalpath][contract]") {
    const auto s = state();
    auto out = distributePairEvents(s, {{1, 1e-100, endpoint(0), endpoint(3)}});
    REQUIRE(out.balanceError() == 0);
    out.electronIntegrals_per_m_s[3] = 0;
    REQUIRE(out.balanceError() == 1);
    out = distributePairEvents(s, {{1, -1e-100, endpoint(0), endpoint(3)}});
    REQUIRE(out.currentPerUm(true, 1, TunnelingConstants::Codata2018) < 0);
    out = distributePairEvents(s, {});
    REQUIRE(out.balanceError() == 0);
    out.electronIntegrals_per_m_s[0] = 1e-100;
    REQUIRE_THROWS(out.balanceError());
}
TEST_CASE("NLP N1 insulator zero boxes stay empty and contact deletion cannot pass ledger validation", "[nonlocalpath][contract]") {
    auto s = state();
    s.nodes[1].semiconductor = false;
    s.nodes[1].conductionDOS_per_m3 = s.nodes[1].valenceDOS_per_m3 = 0;
    s.nodes[1].conduction_J = s.nodes[1].valence_J = 0;
    REQUIRE_NOTHROW(s.validate());
    auto out = distributePairEvents(s, {{1, 2, endpoint(0), endpoint(3)}});
    REQUIRE(out.volumeRates(true, {1, 0, 1, 1})[1] == 0);
    REQUIRE_NOTHROW(out.validate(s));
    out.electronIntegrals_per_m_s[3] = 0;
    out.holeIntegrals_per_m_s[0] = 0;
    REQUIRE(out.balanceError() == 0); // Global e/h equality alone is insufficient.
    REQUIRE_THROWS(out.validate(s));
    REQUIRE_THROWS(distributePairEvents(s, {{1, 1, endpoint(0), endpoint(1)}}));
}
TEST_CASE("NLP N1 path bounds and state revisions prevent stale trial reuse", "[nonlocalpath][contract]") {
    auto s = state();
    TunnelingPath p;
    p.identity = s.identity;
    p.points = {{Point2(0, 0), endpoint(0)}, {Point2(1e-6, 0), endpoint(3)}};
    p.integrationIntervals_m = {{0, 1e-6}};
    p.termination = TunnelingTermination::Complete;
    REQUIRE_NOTHROW(p.validate(s, 1e-6));
    auto bad = p; bad.identity.trial++; REQUIRE_THROWS(bad.validate(s, 1e-6));
    REQUIRE_THROWS(p.validate(s, .9e-6));
    bad = p; bad.integrationIntervals_m = {{0, .7e-6}, {.6e-6, 1e-6}};
    REQUIRE_THROWS(bad.validate(s, 1e-6));
    bad = p; bad.points[1].position_m = Point2::Zero(); REQUIRE_THROWS(bad.validate(s, 1e-6));
}
TEST_CASE("NLP N1 remote sparse support changes invalidate symbolic caches", "[nonlocalpath][contract]") {
    NonlocalJacobianSupport graph;
    REQUIRE(graph.replace(12, 12, {{0, 11}, {1, 7}, {0, 11}}));
    const auto revision = graph.revision();
    REQUIRE_NOTHROW(graph.require(0, 11));
    REQUIRE_FALSE(graph.replace(12, 12, {{1, 7}, {0, 11}}));
    REQUIRE(graph.revision() == revision);
    REQUIRE_THROWS(graph.require(0, 10));
    REQUIRE_THROWS(graph.replace(12, 12, {{0, 12}}));
    REQUIRE(graph.revision() == revision); // Failed update is atomic.
    REQUIRE(graph.replace(12, 12, {{1, 7}}));
    REQUIRE(graph.revision() == revision + 1);
    REQUIRE_THROWS(graph.require(0, 11));
}
TEST_CASE("NLP N1 source derivative columns carry physical units and an explicit mode", "[nonlocalpath][contract]") {
    const auto s = state();
    auto out = distributePairEvents(s, {{1, 2, endpoint(0), endpoint(3)}});
    REQUIRE_FALSE(out.evaluatedDerivativeMode.has_value());
    REQUIRE_NOTHROW(out.validate(s));
    out.derivatives = {{true, 3, {0, TunnelingVariable::ElectronQF_J}, 1e20}};
    REQUIRE_THROWS(out.validate(s));
    out.evaluatedDerivativeMode = TunnelingDerivativeMode::FullMathematical;
    REQUIRE_NOTHROW(out.validate(s));
    out.derivatives[0].column.node = 4; REQUIRE_THROWS(out.validate(s));
    out.derivatives[0].column.node = 0; out.identity.trial++; REQUIRE_THROWS(out.validate(s));
}
TEST_CASE("NLP N1 does not silently map unsupported production models to E2", "[nonlocalpath][contract]") {
    BandToBandTunnelingConfig c;
    c.model = "NonLocalPath";
    REQUIRE_THROWS(BandToBandTunnelingModel(c));
    c.model = "e2";
    const BandToBandTunnelingModel local(c);
    const Real field = 1e8;
    REQUIRE(local.generationRate(field) / (c.prefactorA_SI * field * field *
        std::exp(-c.exponentialB_V_per_m / field)) == Approx(1));
}
