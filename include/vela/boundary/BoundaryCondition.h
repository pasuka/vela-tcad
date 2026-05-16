#pragma once

#include "vela/core/Types.h"

#include <nlohmann/json_fwd.hpp>
#include <optional>
#include <string>
#include <vector>

namespace vela {

/// Logical type assigned to a contact entry in the device-level deck.
///
/// The enum is intentionally larger than the set of physical models that are
/// currently implemented in the solvers.  Schottky and Floating are part of
/// the schema/parser milestone (M1.1) so future tasks can plug their physics
/// in without touching deck parsing again.
enum class ContactType {
    Ohmic,
    Dirichlet,
    MetalGate,
    Schottky,
    Floating,
};

/// Logical type assigned to a non-contact boundary segment.
enum class BoundaryType {
    Dirichlet,
    Neumann,
    Insulating,
    Symmetry,
};

/// Parsed representation of a single ``contacts[]`` entry.
///
/// The struct keeps both the originally requested ``type`` and a handful of
/// optional model parameters (``flatbandVoltage``, ``workFunction_eV`` and
/// reserved Schottky knobs) so that downstream physics layers can decide how
/// to use them.  Only fields that were present in the deck are populated; the
/// rest remain ``std::nullopt``.
struct ContactBoundarySpec {
    std::string  name;
    ContactType  type = ContactType::Ohmic;
    Real         bias = 0.0;

    std::optional<Real> flatbandVoltage;       ///< Optional flat-band voltage [V].
    std::optional<Real> workFunction_eV;       ///< Metal work-function offset [eV].
    std::optional<Real> barrier_eV;            ///< Schottky barrier height [eV].
    std::optional<Real> surfaceRecombinationVelocity;

    /// Verbatim copy of the spelling found in the deck for diagnostics.
    /// Empty when the legacy untyped form was used.
    std::string rawType;
};

/// Parsed representation of a single boundary-segment entry.
///
/// Supports explicit Neumann, insulating, and symmetry boundary conditions.
/// The boundary segment is defined by a polyline of node IDs.
///
/// For Neumann boundaries:
///   - value represents normal_displacement_C_per_m2 (D.n at the boundary)
///   - Positive value means outward flux (field pointing out of domain)
///   - RHS contribution: value * edge_length / 2 to each endpoint
///
/// For insulating/symmetry boundaries:
///   - Equivalent to zero Neumann (D.n = 0)
///   - No matrix/RHS modification needed
struct BoundarySegmentSpec {
    std::string         name;
    BoundaryType        type = BoundaryType::Dirichlet;
    std::vector<Index>  node_ids;  ///< Polyline defining the boundary segment
    Real                value = 0.0;  ///< Boundary value (interpretation depends on type)
    std::string         rawType;
};

// ---------------------------------------------------------------------------
// String <-> enum helpers
// ---------------------------------------------------------------------------

/// Parse a contact type string with case-insensitive matching that also
/// accepts ``-`` and ``_`` as separators (e.g. ``metal-gate`` and
/// ``metal_gate`` map to ``MetalGate``).  Throws ``std::invalid_argument``
/// for unknown values.
ContactType contactTypeFromString(std::string_view text);

/// Parse a boundary type string with the same normalization rules as
/// ``contactTypeFromString``.  Throws ``std::invalid_argument`` for unknown
/// values.
BoundaryType boundaryTypeFromString(std::string_view text);

/// Canonical (lower-case, single-word) string for a contact type, useful for
/// log/diagnostic output.
std::string toString(ContactType type);
std::string toString(BoundaryType type);

// ---------------------------------------------------------------------------
// Parsers
// ---------------------------------------------------------------------------

/// Parse the top-level ``contacts`` array from a deck JSON object.
///
/// Behaviour:
///   * Each entry must have a string ``name`` and a numeric ``bias``.
///   * The optional ``type`` field is normalised via
///     ``contactTypeFromString``.  Legacy decks that omit ``type`` are
///     treated as ``ContactType::Ohmic`` so the existing Gummel/Newton
///     contact-bias path continues to work unchanged.
///   * ``flatband_voltage`` and ``work_function_eV`` are mutually exclusive
///     and either one shifts the effective Poisson Dirichlet potential by
///     subtracting from the applied bias.  The historical convention
///     ``1 eV/q = 1 V`` is preserved.
///   * Schottky-specific fields (``barrier_eV`` and
///     ``surface_recombination_velocity``) are accepted and stored so the
///     parser stays stable across the M1.x milestones.
std::vector<ContactBoundarySpec>
parseContactBoundarySpecs(const nlohmann::json& cfg);

/// Parse the top-level ``boundaries`` array from a deck JSON object.
///
/// Behaviour:
///   * Each entry must have a string ``name`` and a string ``type``.
///   * The ``type`` field is normalised via ``boundaryTypeFromString``.
///   * The ``node_ids`` field must be an array of at least 2 node indices
///     defining a polyline boundary segment.
///   * For Neumann boundaries, the optional ``normal_displacement_C_per_m2``
///     field specifies the normal displacement (D.n) at the boundary.
///   * Insulating and symmetry boundaries are equivalent to zero Neumann.
///   * Dirichlet boundaries are not yet implemented via this path.
std::vector<BoundarySegmentSpec>
parseBoundarySegmentSpecs(const nlohmann::json& cfg);

/// Compute the effective electrostatic Dirichlet potential for a contact when
/// it is mapped to a Poisson Dirichlet node.  Matches the legacy formulation:
///
///   psi_contact = bias - flatband_voltage - work_function_eV
///
/// Only one of ``flatband_voltage`` / ``work_function_eV`` may be set; the
/// parser already enforces this, but the helper is also safe when called on a
/// hand-constructed spec.
Real effectivePoissonDirichletPotential(const ContactBoundarySpec& spec);

} // namespace vela
