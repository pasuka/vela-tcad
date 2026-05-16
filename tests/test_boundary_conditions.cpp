#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <nlohmann/json.hpp>

#include "vela/boundary/BoundaryCondition.h"

#include <stdexcept>
#include <string>

using namespace vela;
using Catch::Approx;

namespace {

nlohmann::json makeContactsJson(const nlohmann::json& contacts)
{
    return nlohmann::json{{"contacts", contacts}};
}

} // namespace

TEST_CASE("contactTypeFromString accepts canonical names", "[boundary]")
{
    REQUIRE(contactTypeFromString("ohmic")     == ContactType::Ohmic);
    REQUIRE(contactTypeFromString("dirichlet") == ContactType::Dirichlet);
    REQUIRE(contactTypeFromString("schottky")  == ContactType::Schottky);
    REQUIRE(contactTypeFromString("floating")  == ContactType::Floating);
}

TEST_CASE("contactTypeFromString normalises hyphen/underscore/case", "[boundary]")
{
    REQUIRE(contactTypeFromString("metal_gate") == ContactType::MetalGate);
    REQUIRE(contactTypeFromString("metal-gate") == ContactType::MetalGate);
    REQUIRE(contactTypeFromString("Metal-Gate") == ContactType::MetalGate);
    REQUIRE(contactTypeFromString("METALGATE")  == ContactType::MetalGate);
    REQUIRE(contactTypeFromString("Ohmic")      == ContactType::Ohmic);
    REQUIRE(contactTypeFromString("DIRICHLET")  == ContactType::Dirichlet);
}

TEST_CASE("contactTypeFromString rejects unknown values", "[boundary]")
{
    REQUIRE_THROWS_AS(contactTypeFromString("mystery"), std::invalid_argument);
    REQUIRE_THROWS_AS(contactTypeFromString(""), std::invalid_argument);
}

TEST_CASE("boundaryTypeFromString normalises and validates", "[boundary]")
{
    REQUIRE(boundaryTypeFromString("dirichlet")   == BoundaryType::Dirichlet);
    REQUIRE(boundaryTypeFromString("Neumann")     == BoundaryType::Neumann);
    REQUIRE(boundaryTypeFromString("insulating")  == BoundaryType::Insulating);
    REQUIRE(boundaryTypeFromString("Symmetry")    == BoundaryType::Symmetry);
    REQUIRE_THROWS_AS(boundaryTypeFromString("robin"), std::invalid_argument);
}

TEST_CASE("parseContactBoundarySpecs treats untyped contacts as Ohmic", "[boundary]")
{
    const nlohmann::json cfg = makeContactsJson({
        {{"name", "cathode"}, {"bias", 0.0}},
        {{"name", "anode"},   {"bias", 0.5}},
    });

    const auto specs = parseContactBoundarySpecs(cfg);
    REQUIRE(specs.size() == 2);
    REQUIRE(specs[0].name == "cathode");
    REQUIRE(specs[0].type == ContactType::Ohmic);
    REQUIRE(specs[0].rawType.empty());
    REQUIRE(specs[0].bias == Approx(0.0));
    REQUIRE_FALSE(specs[0].flatbandVoltage.has_value());
    REQUIRE_FALSE(specs[0].workFunction_eV.has_value());

    REQUIRE(specs[1].name == "anode");
    REQUIRE(specs[1].bias == Approx(0.5));
    REQUIRE(specs[1].type == ContactType::Ohmic);
}

TEST_CASE("parseContactBoundarySpecs honours explicit type with normalised names", "[boundary]")
{
    const nlohmann::json cfg = makeContactsJson({
        {{"name", "gate"},   {"bias", 1.0}, {"type", "Metal-Gate"}},
        {{"name", "anode"},  {"bias", 0.0}, {"type", "DIRICHLET"}},
        {{"name", "field"},  {"bias", 0.0}, {"type", "schottky"},
         {"barrier_eV", 0.6}},
    });

    const auto specs = parseContactBoundarySpecs(cfg);
    REQUIRE(specs.size() == 3);
    REQUIRE(specs[0].type == ContactType::MetalGate);
    REQUIRE(specs[0].rawType == "Metal-Gate");
    REQUIRE(specs[1].type == ContactType::Dirichlet);
    REQUIRE(specs[2].type == ContactType::Schottky);
    REQUIRE(specs[2].barrier_eV.has_value());
    REQUIRE(*specs[2].barrier_eV == Approx(0.6));
}

TEST_CASE("parseContactBoundarySpecs captures flatband and work-function offsets", "[boundary]")
{
    const nlohmann::json cfg = makeContactsJson({
        {{"name", "gate_fb"}, {"bias", 1.0}, {"flatband_voltage", 0.3}},
        {{"name", "gate_wf"}, {"bias", 1.0}, {"work_function_eV", 4.1}},
    });

    const auto specs = parseContactBoundarySpecs(cfg);
    REQUIRE(specs.size() == 2);
    REQUIRE(specs[0].flatbandVoltage.has_value());
    REQUIRE(*specs[0].flatbandVoltage == Approx(0.3));
    REQUIRE_FALSE(specs[0].workFunction_eV.has_value());
    REQUIRE(specs[1].workFunction_eV.has_value());
    REQUIRE(*specs[1].workFunction_eV == Approx(4.1));
    REQUIRE_FALSE(specs[1].flatbandVoltage.has_value());
}

TEST_CASE("parseContactBoundarySpecs rejects mutually exclusive offsets", "[boundary]")
{
    const nlohmann::json cfg = makeContactsJson({
        {{"name", "gate"},
         {"bias", 1.0},
         {"flatband_voltage", 0.2},
         {"work_function_eV", 4.1}},
    });
    REQUIRE_THROWS_AS(parseContactBoundarySpecs(cfg), std::runtime_error);
}

TEST_CASE("parseContactBoundarySpecs rejects unknown types early", "[boundary]")
{
    const nlohmann::json cfg = makeContactsJson({
        {{"name", "weird"}, {"bias", 0.0}, {"type", "mystery"}},
    });
    REQUIRE_THROWS_AS(parseContactBoundarySpecs(cfg), std::invalid_argument);
}

TEST_CASE("effectivePoissonDirichletPotential matches legacy formula", "[boundary]")
{
    SECTION("bias only")
    {
        ContactBoundarySpec s;
        s.name = "anode";
        s.bias = 0.7;
        REQUIRE(effectivePoissonDirichletPotential(s) == Approx(0.7));
    }
    SECTION("flatband shift")
    {
        ContactBoundarySpec s;
        s.name = "gate";
        s.bias = 1.0;
        s.flatbandVoltage = 0.3;
        REQUIRE(effectivePoissonDirichletPotential(s) == Approx(0.7));
    }
    SECTION("work-function shift in V")
    {
        ContactBoundarySpec s;
        s.name = "gate";
        s.bias = 1.0;
        s.workFunction_eV = 4.1; // 1 eV/q == 1 V
        REQUIRE(effectivePoissonDirichletPotential(s) == Approx(1.0 - 4.1));
    }
    SECTION("both offsets rejected")
    {
        ContactBoundarySpec s;
        s.name = "gate";
        s.bias = 1.0;
        s.flatbandVoltage = 0.3;
        s.workFunction_eV = 4.1;
        REQUIRE_THROWS_AS(effectivePoissonDirichletPotential(s), std::runtime_error);
    }
}

TEST_CASE("parseContactBoundarySpecs handles missing or empty contacts", "[boundary]")
{
    SECTION("missing contacts field returns empty vector")
    {
        nlohmann::json cfg = nlohmann::json::object();
        REQUIRE(parseContactBoundarySpecs(cfg).empty());
    }
    SECTION("empty contacts array returns empty vector")
    {
        nlohmann::json cfg = makeContactsJson(nlohmann::json::array());
        REQUIRE(parseContactBoundarySpecs(cfg).empty());
    }
    SECTION("non-array contacts raises invalid_argument")
    {
        nlohmann::json cfg;
        cfg["contacts"] = "oops";
        REQUIRE_THROWS_AS(parseContactBoundarySpecs(cfg), std::invalid_argument);
    }
}

TEST_CASE("toString round-trips contact and boundary enums", "[boundary]")
{
    REQUIRE(toString(ContactType::Ohmic)     == "ohmic");
    REQUIRE(toString(ContactType::Dirichlet) == "dirichlet");
    REQUIRE(toString(ContactType::MetalGate) == "metal_gate");
    REQUIRE(toString(ContactType::Schottky)  == "schottky");
    REQUIRE(toString(ContactType::Floating)  == "floating");

    REQUIRE(toString(BoundaryType::Dirichlet)  == "dirichlet");
    REQUIRE(toString(BoundaryType::Neumann)    == "neumann");
    REQUIRE(toString(BoundaryType::Insulating) == "insulating");
    REQUIRE(toString(BoundaryType::Symmetry)   == "symmetry");
}

// ---------------------------------------------------------------------------
// Boundary segment parsing tests
// ---------------------------------------------------------------------------

TEST_CASE("parseBoundarySegmentSpecs handles missing boundaries field", "[boundary]")
{
    nlohmann::json cfg = nlohmann::json::object();
    REQUIRE(parseBoundarySegmentSpecs(cfg).empty());
}

TEST_CASE("parseBoundarySegmentSpecs parses Neumann boundary", "[boundary]")
{
    nlohmann::json cfg;
    cfg["boundaries"] = nlohmann::json::array({
        {{"name", "top"}, {"type", "neumann"}, {"node_ids", {0, 1, 2}},
         {"normal_displacement_C_per_m2", 1.5e-8}}
    });

    const auto specs = parseBoundarySegmentSpecs(cfg);
    REQUIRE(specs.size() == 1);
    REQUIRE(specs[0].name == "top");
    REQUIRE(specs[0].type == BoundaryType::Neumann);
    REQUIRE(specs[0].node_ids.size() == 3);
    REQUIRE(specs[0].node_ids[0] == 0);
    REQUIRE(specs[0].node_ids[1] == 1);
    REQUIRE(specs[0].node_ids[2] == 2);
    REQUIRE(specs[0].value == Approx(1.5e-8));
}

TEST_CASE("parseBoundarySegmentSpecs parses insulating boundary", "[boundary]")
{
    nlohmann::json cfg;
    cfg["boundaries"] = nlohmann::json::array({
        {{"name", "left_insulating"}, {"type", "insulating"}, {"node_ids", {0, 3, 6}}}
    });

    const auto specs = parseBoundarySegmentSpecs(cfg);
    REQUIRE(specs.size() == 1);
    REQUIRE(specs[0].name == "left_insulating");
    REQUIRE(specs[0].type == BoundaryType::Insulating);
    REQUIRE(specs[0].node_ids.size() == 3);
    REQUIRE(specs[0].value == Approx(0.0)); // Default value
}

TEST_CASE("parseBoundarySegmentSpecs parses symmetry boundary", "[boundary]")
{
    nlohmann::json cfg;
    cfg["boundaries"] = nlohmann::json::array({
        {{"name", "symmetry_axis"}, {"type", "symmetry"}, {"node_ids", {0, 1}}}
    });

    const auto specs = parseBoundarySegmentSpecs(cfg);
    REQUIRE(specs.size() == 1);
    REQUIRE(specs[0].name == "symmetry_axis");
    REQUIRE(specs[0].type == BoundaryType::Symmetry);
    REQUIRE(specs[0].value == Approx(0.0));
}

TEST_CASE("parseBoundarySegmentSpecs rejects Dirichlet type", "[boundary]")
{
    nlohmann::json cfg;
    cfg["boundaries"] = nlohmann::json::array({
        {{"name", "dir"}, {"type", "dirichlet"}, {"node_ids", {0, 1}}}
    });

    REQUIRE_THROWS_AS(parseBoundarySegmentSpecs(cfg), std::runtime_error);
}

TEST_CASE("parseBoundarySegmentSpecs requires at least 2 node IDs", "[boundary]")
{
    nlohmann::json cfg;
    cfg["boundaries"] = nlohmann::json::array({
        {{"name", "bad"}, {"type", "neumann"}, {"node_ids", {0}}}
    });

    REQUIRE_THROWS_AS(parseBoundarySegmentSpecs(cfg), std::invalid_argument);
}

TEST_CASE("parseBoundarySegmentSpecs requires node_ids field", "[boundary]")
{
    nlohmann::json cfg;
    cfg["boundaries"] = nlohmann::json::array({
        {{"name", "bad"}, {"type", "neumann"}}
    });

    REQUIRE_THROWS_AS(parseBoundarySegmentSpecs(cfg), std::invalid_argument);
}

TEST_CASE("parseBoundarySegmentSpecs rejects non-finite Neumann value", "[boundary]")
{
    nlohmann::json cfg;
    cfg["boundaries"] = nlohmann::json::array({
        {{"name", "bad"}, {"type", "neumann"}, {"node_ids", {0, 1}},
         {"normal_displacement_C_per_m2", std::numeric_limits<double>::infinity()}}
    });

    REQUIRE_THROWS_AS(parseBoundarySegmentSpecs(cfg), std::invalid_argument);
}

TEST_CASE("parseBoundarySegmentSpecs normalises boundary type strings", "[boundary]")
{
    nlohmann::json cfg;
    cfg["boundaries"] = nlohmann::json::array({
        {{"name", "b1"}, {"type", "Neumann"}, {"node_ids", {0, 1}}},
        {{"name", "b2"}, {"type", "INSULATING"}, {"node_ids", {1, 2}}},
        {{"name", "b3"}, {"type", "Symmetry"}, {"node_ids", {2, 3}}}
    });

    const auto specs = parseBoundarySegmentSpecs(cfg);
    REQUIRE(specs.size() == 3);
    REQUIRE(specs[0].type == BoundaryType::Neumann);
    REQUIRE(specs[1].type == BoundaryType::Insulating);
    REQUIRE(specs[2].type == BoundaryType::Symmetry);
}
