#include "vela/boundary/BoundaryCondition.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>

namespace vela {

namespace {

/// Normalise a boundary/contact ``type`` string so that the parser accepts
/// any casing, hyphenated, or underscored variant: ``Metal-Gate`` and
/// ``metal_gate`` both collapse to ``metalgate``.
std::string normalizeTypeKey(std::string_view text)
{
    std::string out;
    out.reserve(text.size());
    for (char ch : text) {
        if (ch == '-' || ch == '_' || ch == ' ' || ch == '\t')
            continue;
        out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
    }
    return out;
}

} // namespace

// ---------------------------------------------------------------------------
// String helpers
// ---------------------------------------------------------------------------

ContactType contactTypeFromString(std::string_view text)
{
    const std::string key = normalizeTypeKey(text);
    if (key == "ohmic")     return ContactType::Ohmic;
    if (key == "dirichlet") return ContactType::Dirichlet;
    if (key == "metalgate") return ContactType::MetalGate;
    if (key == "gate")      return ContactType::MetalGate;
    if (key == "schottky")  return ContactType::Schottky;
    if (key == "floating")  return ContactType::Floating;
    throw std::invalid_argument(
        "contactTypeFromString: unknown contact type '" + std::string(text) + "'. "
        "Expected one of: ohmic, dirichlet, metal_gate, schottky, floating.");
}

BoundaryType boundaryTypeFromString(std::string_view text)
{
    const std::string key = normalizeTypeKey(text);
    if (key == "dirichlet")   return BoundaryType::Dirichlet;
    if (key == "neumann")     return BoundaryType::Neumann;
    if (key == "insulating")  return BoundaryType::Insulating;
    if (key == "symmetry")    return BoundaryType::Symmetry;
    throw std::invalid_argument(
        "boundaryTypeFromString: unknown boundary type '" + std::string(text) + "'. "
        "Expected one of: dirichlet, neumann, insulating, symmetry.");
}

std::string toString(ContactType type)
{
    switch (type) {
        case ContactType::Ohmic:     return "ohmic";
        case ContactType::Dirichlet: return "dirichlet";
        case ContactType::MetalGate: return "metal_gate";
        case ContactType::Schottky:  return "schottky";
        case ContactType::Floating:  return "floating";
    }
    return "unknown";
}

std::string toString(BoundaryType type)
{
    switch (type) {
        case BoundaryType::Dirichlet:  return "dirichlet";
        case BoundaryType::Neumann:    return "neumann";
        case BoundaryType::Insulating: return "insulating";
        case BoundaryType::Symmetry:   return "symmetry";
    }
    return "unknown";
}

// ---------------------------------------------------------------------------
// Parsers
// ---------------------------------------------------------------------------

std::vector<ContactBoundarySpec>
parseContactBoundarySpecs(const nlohmann::json& cfg)
{
    std::vector<ContactBoundarySpec> specs;
    if (!cfg.contains("contacts"))
        return specs;

    const auto& contacts = cfg.at("contacts");
    if (!contacts.is_array()) {
        throw std::invalid_argument(
            "parseContactBoundarySpecs: 'contacts' must be a JSON array.");
    }

    specs.reserve(contacts.size());
    for (const auto& ct : contacts) {
        ContactBoundarySpec spec;
        spec.name = ct.at("name").get<std::string>();
        spec.bias = ct.at("bias").get<Real>();

        if (ct.contains("type")) {
            spec.rawType = ct.at("type").get<std::string>();
            spec.type = contactTypeFromString(spec.rawType);
        } else {
            // Legacy decks without an explicit type are treated as Ohmic so
            // that DD/Gummel/Newton paths see the historical contact-bias
            // semantics unchanged.
            spec.type = ContactType::Ohmic;
        }

        const bool hasFlatband = ct.contains("flatband_voltage");
        const bool hasWorkFunction = ct.contains("work_function_eV");
        if (hasFlatband && hasWorkFunction) {
            throw std::runtime_error(
                "parseContactBoundarySpecs: contact '" + spec.name +
                "' cannot set both flatband_voltage and work_function_eV.");
        }
        if (hasFlatband)
            spec.flatbandVoltage = ct.at("flatband_voltage").get<Real>();
        if (hasWorkFunction)
            spec.workFunction_eV = ct.at("work_function_eV").get<Real>();

        if (ct.contains("barrier_eV"))
            spec.barrier_eV = ct.at("barrier_eV").get<Real>();
        if (ct.contains("surface_recombination_velocity")) {
            spec.surfaceRecombinationVelocity =
                ct.at("surface_recombination_velocity").get<Real>();
        }

        specs.push_back(std::move(spec));
    }
    return specs;
}

std::vector<BoundarySegmentSpec>
parseBoundarySegmentSpecs(const nlohmann::json& cfg)
{
    std::vector<BoundarySegmentSpec> specs;
    if (!cfg.contains("boundaries"))
        return specs;

    const auto& boundaries = cfg.at("boundaries");
    if (!boundaries.is_array()) {
        throw std::invalid_argument(
            "parseBoundarySegmentSpecs: 'boundaries' must be a JSON array.");
    }

    specs.reserve(boundaries.size());
    for (const auto& bd : boundaries) {
        BoundarySegmentSpec spec;
        spec.name = bd.at("name").get<std::string>();
        spec.rawType = bd.at("type").get<std::string>();
        spec.type = boundaryTypeFromString(spec.rawType);

        if (!bd.contains("node_ids")) {
            throw std::invalid_argument(
                "parseBoundarySegmentSpecs: boundary '" + spec.name +
                "' must have a 'node_ids' array.");
        }

        spec.node_ids = bd.at("node_ids").get<std::vector<Index>>();
        if (spec.node_ids.size() < 2) {
            throw std::invalid_argument(
                "parseBoundarySegmentSpecs: boundary '" + spec.name +
                "' must have at least 2 node IDs.");
        }

        // For Neumann boundaries, read the normal displacement value
        if (spec.type == BoundaryType::Neumann) {
            spec.value = bd.value("normal_displacement_C_per_m2", 0.0);
            if (!std::isfinite(spec.value)) {
                throw std::invalid_argument(
                    "parseBoundarySegmentSpecs: boundary '" + spec.name +
                    "' has non-finite normal_displacement_C_per_m2.");
            }
        }

        // Insulating and symmetry are zero Neumann, no value needed
        // Dirichlet via boundaries is not yet implemented
        if (spec.type == BoundaryType::Dirichlet) {
            throw std::runtime_error(
                "parseBoundarySegmentSpecs: boundary '" + spec.name +
                "' has type 'dirichlet' which is not yet implemented via the boundaries array. "
                "Use contacts for Dirichlet boundary conditions.");
        }

        specs.push_back(std::move(spec));
    }
    return specs;
}

Real effectivePoissonDirichletPotential(const ContactBoundarySpec& spec)
{
    Real value = spec.bias;
    if (spec.flatbandVoltage && spec.workFunction_eV) {
        throw std::runtime_error(
            "effectivePoissonDirichletPotential: contact '" + spec.name +
            "' cannot set both flatband_voltage and work_function_eV.");
    }
    if (spec.flatbandVoltage)
        value -= *spec.flatbandVoltage;
    if (spec.workFunction_eV)
        value -= *spec.workFunction_eV; // 1 eV/q == 1 V convention.
    return value;
}

} // namespace vela
