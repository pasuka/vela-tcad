#include "vela/boundary/BoundaryCondition.h"

#include "vela/core/PhysicalConstants.h"

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
        if (ct.contains("electron_barrier_eV"))
            spec.electronBarrier_eV = ct.at("electron_barrier_eV").get<Real>();
        if (ct.contains("hole_barrier_eV"))
            spec.holeBarrier_eV = ct.at("hole_barrier_eV").get<Real>();
        if (ct.contains("surface_recombination_velocity")) {
            spec.surfaceRecombinationVelocity =
                ct.at("surface_recombination_velocity").get<Real>();
        } else if (ct.contains("surface_recombination_velocity_m_per_s")) {
            spec.surfaceRecombinationVelocity =
                ct.at("surface_recombination_velocity_m_per_s").get<Real>();
        }
        if (ct.contains("emission_model")) {
            spec.emissionModel = ct.at("emission_model").get<std::string>();
        }

        // Validate finite numerics where present so configuration mistakes
        // surface here instead of at solve time.
        const auto requireFinite = [&](const std::optional<Real>& v, const char* field) {
            if (v && !std::isfinite(*v))
                throw std::runtime_error(
                    "parseContactBoundarySpecs: contact '" + spec.name +
                    "' has non-finite '" + field + "'.");
        };
        requireFinite(spec.barrier_eV, "barrier_eV");
        requireFinite(spec.electronBarrier_eV, "electron_barrier_eV");
        requireFinite(spec.holeBarrier_eV, "hole_barrier_eV");
        requireFinite(spec.surfaceRecombinationVelocity,
                      "surface_recombination_velocity_m_per_s");

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

// ---------------------------------------------------------------------------
// Schottky / metal-semiconductor contact prototype
// ---------------------------------------------------------------------------

namespace {

bool isFinitePositive(Real x)
{
    return std::isfinite(x) && x > 0.0;
}

} // namespace

Real schottkyElectronBarrier_eV(const ContactBoundarySpec& spec,
                                Real electronAffinity_eV)
{
    if (spec.electronBarrier_eV)
        return *spec.electronBarrier_eV;
    if (spec.barrier_eV)
        return *spec.barrier_eV;
    if (spec.workFunction_eV && std::isfinite(electronAffinity_eV)) {
        const Real phi = *spec.workFunction_eV - electronAffinity_eV;
        if (std::isfinite(phi))
            return phi;
    }
    throw std::invalid_argument(
        "schottkyElectronBarrier_eV: contact '" + spec.name +
        "' must set 'barrier_eV', 'electron_barrier_eV', or "
        "'work_function_eV' (with a known electron affinity) when type='schottky'.");
}

Real schottkyHoleBarrier_eV(const ContactBoundarySpec& spec,
                            Real electronBarrier_eV,
                            Real bandgap_eV)
{
    if (spec.holeBarrier_eV)
        return *spec.holeBarrier_eV;
    if (std::isfinite(bandgap_eV) && bandgap_eV > 0.0) {
        const Real phiP = bandgap_eV - electronBarrier_eV;
        // Clamp to a small positive value so we never exponentiate +inf when
        // a deck specifies an unphysical phi_Bn > Eg.
        return std::max(phiP, 1.0e-3);
    }
    // Fall back to a symmetric prototype assumption.
    return std::max(electronBarrier_eV, 1.0e-3);
}

ContactState computeSchottkyContactState(const ContactBoundarySpec& spec,
                                         Real ni,
                                         Real Nc,
                                         Real Nv,
                                         Real bandgap_eV,
                                         Real electronAffinity_eV,
                                         Real netDoping_m3,
                                         Real temperature_K)
{
    (void)netDoping_m3; // Reserved for a future image-force / doping correction.

    if (!isFinitePositive(temperature_K)) {
        throw std::invalid_argument(
            "computeSchottkyContactState: contact '" + spec.name +
            "' requires a positive finite temperature_K.");
    }
    if (!std::isfinite(ni) || ni < 0.0) {
        throw std::invalid_argument(
            "computeSchottkyContactState: contact '" + spec.name +
            "' has non-finite or negative intrinsic density 'ni'.");
    }

    // Allow the model selector to be empty (default) or "dirichlet_barrier".
    if (!spec.emissionModel.empty() &&
        spec.emissionModel != "dirichlet_barrier") {
        throw std::runtime_error(
            "computeSchottkyContactState: contact '" + spec.name +
            "' uses unsupported emission_model '" + spec.emissionModel +
            "'. Only 'dirichlet_barrier' is implemented in this prototype.");
    }

    const Real Vt = constants::kb * temperature_K / constants::q;
    const Real phiBn = schottkyElectronBarrier_eV(spec, electronAffinity_eV);
    const Real phiBp = schottkyHoleBarrier_eV(spec, phiBn, bandgap_eV);

    ContactState state;
    state.phin = spec.bias;
    state.phip = spec.bias;

    // Electrostatic potential at the contact: pin the metal Fermi level so
    // E_Fm = E_C(surface) - q*phi_Bn.  In the standard Vela convention the
    // intrinsic level is at psi=0 with E_Fi reference.  We approximate the
    // surface offset relative to the metal Fermi level as
    //   psi_contact = bias - (phi_Bn - chi - 0.5*Eg)
    // when both bandgap and affinity are known so the offset is measured
    // relative to mid-gap, otherwise fall back to a flat
    //   psi_contact = bias - phi_Bn (1 eV/q == 1 V)
    // form.  This is intentionally a smoke-level prototype.
    Real psiOffset = phiBn;
    if (std::isfinite(bandgap_eV) && bandgap_eV > 0.0 &&
        std::isfinite(electronAffinity_eV)) {
        // Offset of the conduction band edge below the metal Fermi level,
        // measured from the intrinsic mid-gap reference used elsewhere in
        // Vela (psi=0 corresponds to the intrinsic Fermi level in the bulk).
        psiOffset = phiBn - 0.5 * bandgap_eV;
    }
    if (spec.workFunction_eV && std::isfinite(electronAffinity_eV) &&
        !spec.barrier_eV && !spec.electronBarrier_eV &&
        std::isfinite(bandgap_eV) && bandgap_eV > 0.0) {
        // When the deck pins the metal work function rather than the barrier
        // height directly, use the standard relation Vbi = phi_M - chi - Eg/2
        // so the Poisson Dirichlet value matches the legacy
        // effectivePoissonDirichletPotential semantics for metal gates.
        psiOffset = *spec.workFunction_eV - electronAffinity_eV - 0.5 * bandgap_eV;
    }
    state.psi = spec.bias - psiOffset;

    // Carrier densities.  Prefer the Boltzmann form n = Nc*exp(-phi_Bn/Vt)
    // so the result is independent of the bulk doping.  Fall back to a
    // ni-based estimate when band edge density-of-states is unavailable.
    if (isFinitePositive(Nc)) {
        state.n = Nc * std::exp(-phiBn / Vt);
    } else {
        const Real halfEg = (std::isfinite(bandgap_eV) && bandgap_eV > 0.0)
            ? 0.5 * bandgap_eV : 0.0;
        state.n = ni * std::exp(-(phiBn - halfEg) / Vt);
    }
    if (isFinitePositive(Nv)) {
        state.p = Nv * std::exp(-phiBp / Vt);
    } else {
        const Real halfEg = (std::isfinite(bandgap_eV) && bandgap_eV > 0.0)
            ? 0.5 * bandgap_eV : 0.0;
        state.p = ni * std::exp(-(phiBp - halfEg) / Vt);
    }

    // Numeric guards: keep carrier densities finite and strictly positive so
    // downstream solvers do not degenerate.  Use a small but non-zero floor.
    constexpr Real kMinDensity = 1.0e-30;
    if (!std::isfinite(state.n) || state.n < kMinDensity)
        state.n = kMinDensity;
    if (!std::isfinite(state.p) || state.p < kMinDensity)
        state.p = kMinDensity;

    return state;
}

} // namespace vela


