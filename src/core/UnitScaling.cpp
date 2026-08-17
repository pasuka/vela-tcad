#include "vela/core/UnitScaling.h"

#include <nlohmann/json.hpp>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace vela {

PhysicalUnitSystem::PhysicalUnitSystem(Real lengthMPerInternal,
                                       Real concentrationM3PerInternal,
                                       Real sheetDensityM2PerInternal,
                                       Real mobilityM2PerVSPerInternal,
                                       Real velocityMPerSPerInternal,
                                       Real electricFieldVPerMPerInternal,
                                       Real inverseLengthMInvPerInternal,
                                       Real surfaceFieldCoefficientMPerVPerInternal,
                                       Real currentDensityAM2PerInternal)
    : lengthMPerInternal_(lengthMPerInternal),
      concentrationM3PerInternal_(concentrationM3PerInternal),
      sheetDensityM2PerInternal_(sheetDensityM2PerInternal),
      mobilityM2PerVSPerInternal_(mobilityM2PerVSPerInternal),
      velocityMPerSPerInternal_(velocityMPerSPerInternal),
      electricFieldVPerMPerInternal_(electricFieldVPerMPerInternal),
      inverseLengthMInvPerInternal_(inverseLengthMInvPerInternal),
      surfaceFieldCoefficientMPerVPerInternal_(surfaceFieldCoefficientMPerVPerInternal),
      currentDensityAM2PerInternal_(currentDensityAM2PerInternal)
{}

PhysicalUnitSystem PhysicalUnitSystem::legacySI()
{
    return PhysicalUnitSystem(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
}

PhysicalUnitSystem PhysicalUnitSystem::tcadInternal()
{
    return PhysicalUnitSystem(1.0e-6, 1.0e6, 1.0e4, 1.0e-4, 1.0e-2, 1.0e2, 1.0e2, 1.0e-2, 1.0e4);
}

Real PhysicalUnitSystem::chargeAreaFactor() const
{
    return concentrationM3PerInternal_ * areaM2PerInternal();
}

Real PhysicalUnitSystem::chargeLineFactor() const
{
    return sheetDensityM2PerInternal_ * lengthMPerInternal_;
}

Real PhysicalUnitSystem::chargeVolumeFactor() const
{
    return concentrationM3PerInternal_ * volumeM3PerInternal();
}

Real PhysicalUnitSystem::chargeSheetFactor() const
{
    return sheetDensityM2PerInternal_ * areaM2PerInternal();
}

Real PhysicalUnitSystem::fieldFromCoordinateDeltaFactor() const
{
    return (1.0 / lengthMPerInternal_) / electricFieldVPerMPerInternal_;
}

Real PhysicalUnitSystem::continuitySourceIntegralFactor() const
{
    // Match a volumetric particle source integrated over a 2-D control area
    // to the native SG edge-flux line integral.  The latter is converted to
    // physical current per device depth by
    //   currentDensityAM2PerInternal * lengthMPerInternal,
    // with q applied by the caller.  The previous area/mobility expression
    // omitted the concentration and field/current-density unit scales; in the
    // TCAD cm/um convention it understated every continuity source by 1e4.
    return concentrationM3PerInternal_ * areaM2PerInternal()
         / (currentDensityAM2PerInternal_ * lengthMPerInternal_);
}

Real PhysicalUnitSystem::currentPerInternalDepthFactor() const
{
    return 1.0e-6;
}

Real PhysicalUnitSystem::internalCurrentPerDeviceDepthToAPerUm(Real value) const
{
    return value * currentPerInternalDepthFactor();
}

UnitScalingConfig::UnitScalingConfig(UnitScalingMode modeValue)
    : mode(modeValue),
      physicalUnitSystem(modeValue == UnitScalingMode::UnitScaling
                             ? PhysicalUnitSystem::tcadInternal()
                             : PhysicalUnitSystem::legacySI())
{}

Real UnitScalingConfig::lengthToSI(Real value) const
{
    return unitSystem().internalLengthToMeters(value);
}

Real UnitScalingConfig::concentrationToSI(Real value) const
{
    return unitSystem().internalConcentrationToM3(value);
}

Real UnitScalingConfig::sheetDensityToSI(Real value) const
{
    return unitSystem().internalSheetDensityToM2(value);
}

Real UnitScalingConfig::mobilityToSI(Real value) const
{
    return unitSystem().internalMobilityToM2PerVS(value);
}

Real UnitScalingConfig::velocityToSI(Real value) const
{
    return unitSystem().internalVelocityToMPerS(value);
}

Real UnitScalingConfig::electricFieldToSI(Real value) const
{
    return unitSystem().internalElectricFieldToVPerM(value);
}

Real UnitScalingConfig::inverseLengthToSI(Real value) const
{
    return unitSystem().internalInverseLengthToMInv(value);
}

Real UnitScalingConfig::surfaceFieldCoefficientToSI(Real value) const
{
    return unitSystem().internalSurfaceFieldCoefficientToMPerV(value);
}

int parseDeckFormatVersion(const nlohmann::json& cfg)
{
    if (!cfg.is_object() || !cfg.contains("format_version"))
        return 0;

    const auto& value = cfg.at("format_version");
    if (!value.is_number_integer()) {
        throw std::invalid_argument(
            "format_version must be the integer 2.");
    }

    // Compare before narrowing: get<int>() would wrap a wide JSON integer such
    // as 4294967298 onto 2 and accept it as a supported version.
    const bool isSupported = value.is_number_unsigned()
        ? value.get<std::uint64_t>() == 2u
        : value.get<std::int64_t>() == 2;
    if (isSupported)
        return 2;

    throw std::invalid_argument(
        "Unsupported format_version " + value.dump() +
        ". The only supported deck format version is 2, in which every value is "
        "expressed in the TCAD internal units (um, cm^-3, cm^2/(V s), V/cm). "
        "Run the deck migration tool to convert an older deck.");
}

namespace {

// Class A of the format_version 2 contract
// (docs/superpowers/specs/2026-08-17-unit-system-v2-contract.md): the value is
// already interpreted in TCAD internal units, only the key name still claims an
// SI unit. In a version-2 deck the v1 name is rejected and the v2 name is the
// single accepted spelling.
struct KeyRename {
    std::string_view v1;
    std::string_view v2;
};

constexpr KeyRename kExactRenames[] = {
    {"ni", "ni_cm3"},
    {"mun", "mun_cm2_V_s"},
    {"mup", "mup_cm2_V_s"},
    {"donors", "donors_cm3"},
    {"acceptors", "acceptors_cm3"},
    {"depth_m", "depth_um"},
    {"contact_radius", "contact_radius_um"},
    {"surface_recombination_velocity", "surface_recombination_velocity_cm_per_s"},
    {"minimum_field_V_m", "minimum_field_V_per_cm"},
    {"switch_field_V_m", "switch_field_V_per_cm"},
    {"electron_B_V_m", "electron_B_V_per_cm"},
    {"hole_B_V_m", "hole_B_V_per_cm"},
    {"electron_b_low_V_m", "electron_b_low_V_per_cm"},
    {"hole_b_low_V_m", "hole_b_low_V_per_cm"},
    {"electron_b_high_V_m", "electron_b_high_V_per_cm"},
    {"hole_b_high_V_m", "hole_b_high_V_per_cm"},
};

// Suffix families: the mobility, impact-ionization, recombination, doping and
// diagnostics keys all rename by unit suffix only. Longer suffixes are listed
// first so that `_m2_V_s` is not matched as `_m2`.
constexpr KeyRename kSuffixRenames[] = {
    {"_m6_per_s", "_cm6_per_s"},
    {"_m2_V_s", "_cm2_V_s"},
    {"_m_per_V", "_cm_per_V"},
    {"_m_per_s", "_cm_per_s"},
    {"_m_inv", "_cm_inv"},
    {"_V_per_m", "_V_per_cm"},
    {"_m_s", "_cm_s"},
    {"_m3", "_cm3"},
    {"_m2", "_cm2"},
};

// `solver.band_to_band` is Class C: it stays SI and keeps its SI key names.
// `solver.normalization` is spelled in v2 units already.
constexpr std::string_view kUntouchedSubtrees[] = {"band_to_band", "normalization"};

bool endsWith(const std::string& text, std::string_view suffix)
{
    return text.size() > suffix.size() &&
        text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// Maps a key between the two spellings of the rename table. `toV2` selects the
// direction: v1 name to v2 name, or v2 name back to the name the field parsers
// consume.
std::optional<std::string> renamedKey(const std::string& key, bool toV2)
{
    for (const KeyRename& rename : kExactRenames) {
        const std::string_view from = toV2 ? rename.v1 : rename.v2;
        if (key == from)
            return std::string(toV2 ? rename.v2 : rename.v1);
    }
    for (const KeyRename& rename : kSuffixRenames) {
        const std::string_view from = toV2 ? rename.v1 : rename.v2;
        const std::string_view to = toV2 ? rename.v2 : rename.v1;
        if (endsWith(key, from))
            return key.substr(0, key.size() - from.size()) + std::string(to);
    }
    return std::nullopt;
}

bool isUntouchedSubtree(const std::string& key)
{
    for (std::string_view name : kUntouchedSubtrees) {
        if (key == name)
            return true;
    }
    return false;
}

nlohmann::json canonicalizeNode(const nlohmann::json& node)
{
    if (node.is_array()) {
        nlohmann::json out = nlohmann::json::array();
        for (const auto& element : node)
            out.push_back(canonicalizeNode(element));
        return out;
    }
    if (!node.is_object())
        return node;

    nlohmann::json out = nlohmann::json::object();
    for (auto it = node.begin(); it != node.end(); ++it) {
        const std::string& key = it.key();
        if (isUntouchedSubtree(key)) {
            out[key] = it.value();
            continue;
        }

        if (const std::optional<std::string> v2 = renamedKey(key, true)) {
            throw std::invalid_argument(
                "Deck key '" + key + "' was renamed to '" + *v2 +
                "' in format_version 2 and is not accepted; the old name is not "
                "an alias. Run the deck migration tool to convert an older deck.");
        }

        const std::optional<std::string> parsed = renamedKey(key, false);
        const std::string& outKey = parsed ? *parsed : key;
        if (out.contains(outKey)) {
            throw std::invalid_argument(
                "Deck key '" + key + "' collides with '" + outKey +
                "' after the format_version 2 rename; remove one of them.");
        }
        out[outKey] = canonicalizeNode(it.value());
    }
    return out;
}

} // namespace

nlohmann::json canonicalizeDeckKeys(const nlohmann::json& doc, int formatVersion)
{
    if (formatVersion != 2)
        return doc;
    return canonicalizeNode(doc);
}

nlohmann::json canonicalizeDeck(const nlohmann::json& cfg)
{
    return canonicalizeDeckKeys(cfg, parseDeckFormatVersion(cfg));
}

UnitScalingConfig parseUnitScalingConfig(const nlohmann::json& cfg)
{
    const bool hasNormalization = cfg.is_object() && cfg.contains("solver") &&
        cfg.at("solver").is_object() && cfg.at("solver").contains("normalization");

    if (parseDeckFormatVersion(cfg) == 2) {
        if (cfg.contains("scaling")) {
            throw std::invalid_argument(
                "scaling is not accepted with format_version 2; the deck already "
                "uses the TCAD internal units. Move "
                "characteristic_length_um, reference_concentration_cm3 and "
                "reference_mobility_cm2_V_s to solver.normalization.");
        }
        UnitScalingConfig scaling{UnitScalingMode::UnitScaling};
        scaling.deckFormatVersion = 2;
        return scaling;
    }

    if (hasNormalization) {
        throw std::invalid_argument(
            "solver.normalization requires format_version 2.");
    }

    if (!cfg.contains("scaling"))
        return {};

    const auto& scaling = cfg.at("scaling");
    if (!scaling.is_object()) {
        throw std::invalid_argument(
            "scaling must be an object with mode set to 'unit_scaling'.");
    }
    if (scaling.contains("system")) {
        throw std::invalid_argument(
            "scaling.system is not supported; use scaling.mode = 'unit_scaling'.");
    }
    if (!scaling.contains("mode")) {
        throw std::invalid_argument(
            "scaling.mode is required when scaling is present; supported value is 'unit_scaling'.");
    }

    const std::string mode = scaling.at("mode").get<std::string>();
    if (mode == "unit_scaling")
        return UnitScalingConfig{UnitScalingMode::UnitScaling};

    throw std::invalid_argument(
        "Unsupported scaling.mode '" + mode + "'. Supported value is 'unit_scaling'. "
        "Omit scaling to keep legacy SI input behavior.");
}

} // namespace vela
