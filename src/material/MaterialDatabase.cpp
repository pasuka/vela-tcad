#include "vela/material/MaterialDatabase.h"

#include <fstream>
#include <cmath>
#include <nlohmann/json.hpp>
#include <set>
#include <stdexcept>
#include <vector>

namespace vela {

namespace {

void setOptionalReal(const nlohmann::json& j,
                     const char* key,
                     std::optional<Real>& field)
{
    if (j.contains(key))
        field = j.at(key).get<Real>();
}

void setOptionalScaledReal(const nlohmann::json& j,
                           const char* key,
                           std::optional<Real>& field,
                           Real (UnitScalingConfig::*scale)(Real) const,
                           UnitScalingConfig scaling)
{
    if (j.contains(key))
        field = (scaling.*scale)(j.at(key).get<Real>());
}

void rejectUnknownKeys(const nlohmann::json& object,
                       const std::set<std::string>& allowed,
                       const std::string& context)
{
    for (auto it = object.begin(); it != object.end(); ++it) {
        if (allowed.count(it.key()) == 0) {
            throw std::invalid_argument(
                context + ": unexpected key '" + it.key() + "'.");
        }
    }
}

Real finiteInRange(const nlohmann::json& j,
                   const char* key,
                   Real minimum,
                   Real maximum,
                   const std::string& context)
{
    const Real value = j.at(key).get<Real>();
    if (!std::isfinite(value) || value < minimum || value > maximum) {
        throw std::invalid_argument(
            context + "." + key + " must be finite and in [" +
            std::to_string(minimum) + ", " + std::to_string(maximum) + "].");
    }
    return value;
}

bool isTemplatesLdmosMaterialsV1(const nlohmann::json& root)
{
    return root.is_object() &&
        root.value("schema", std::string{}) ==
            "vela.templates_ldmos.materials.v1";
}

void validateTemplatesLdmosMaterialRoot(const nlohmann::json& root)
{
    rejectUnknownKeys(
        root,
        {"schema", "benchmark", "revision", "unit_system", "materials"},
        "MaterialDatabase: vela.templates_ldmos.materials.v1");
    if (!root.contains("benchmark") || !root.at("benchmark").is_string() ||
        root.at("benchmark").get<std::string>().empty() ||
        !root.contains("revision") || !root.at("revision").is_number_integer() ||
        root.at("revision").get<int>() < 1) {
        throw std::invalid_argument(
            "MaterialDatabase: versioned materials require benchmark and positive revision.");
    }
    const auto& units = root.at("unit_system");
    if (!units.is_object()) {
        throw std::invalid_argument(
            "MaterialDatabase: versioned materials unit_system must be an object.");
    }
    rejectUnknownKeys(
        units,
        {"concentration", "mobility", "energy", "temperature",
         "thermal_conductivity", "specific_heat", "mass_density"},
        "MaterialDatabase: versioned materials unit_system");
    const std::vector<std::pair<const char*, const char*>> expected = {
        {"concentration", "cm^-3"},
        {"mobility", "cm^2/(V*s)"},
        {"energy", "eV"},
        {"temperature", "K"},
        {"thermal_conductivity", "W/(m*K)"},
        {"specific_heat", "J/(kg*K)"},
        {"mass_density", "kg/m^3"},
    };
    for (const auto& [key, value] : expected) {
        if (units.value(key, std::string{}) != value) {
            throw std::invalid_argument(
                std::string("MaterialDatabase: versioned materials unit_system.") +
                key + " must be '" + value + "'.");
        }
    }
}

Material templatesLdmosMaterialFromJson(const nlohmann::json& j,
                                        const Material* base,
                                        UnitScalingConfig scaling)
{
    const std::string context = "MaterialDatabase: versioned material";
    rejectUnknownKeys(
        j,
        {"name", "role", "eps_r", "intrinsic_carrier_density_cm3",
         "electron_mobility_cm2_per_V_s", "hole_mobility_cm2_per_V_s",
         "bandgap_eV", "electron_affinity_eV",
         "conduction_band_density_of_states_cm3",
         "valence_band_density_of_states_cm3", "temperature_K",
         "electron_quantum_gamma", "electron_quantum_dos_mass_ratio",
         "electron_quantum_coefficient_mass_ratio",
         "thermal_conductivity_W_per_m_K", "specific_heat_J_per_kg_K",
         "mass_density_kg_per_m3", "provenance"},
        context);
    if (!j.contains("name") || !j.at("name").is_string() ||
        j.at("name").get<std::string>().empty()) {
        throw std::invalid_argument(context + ".name must be a non-empty string.");
    }
    if (!j.contains("role") || !j.at("role").is_string() ||
        (j.at("role") != "transport_semiconductor" &&
         j.at("role") != "dielectric" &&
         j.at("role") != "electrostatic_only_semiconductor_poly")) {
        throw std::invalid_argument(context + ".role is missing or unsupported.");
    }
    if (!j.contains("provenance") || !j.at("provenance").is_string() ||
        j.at("provenance").get<std::string>().empty()) {
        throw std::invalid_argument(context + ".provenance must be a non-empty string.");
    }
    if (!j.contains("eps_r"))
        throw std::invalid_argument(context + ".eps_r is required.");
    Material mat = base != nullptr ? *base : Material{};
    mat.name = j.at("name").get<std::string>();
    const auto& units = scaling.unitSystem();
    mat.eps_r = finiteInRange(j, "eps_r", 1.0, 1000.0, context);
    if (j.contains("intrinsic_carrier_density_cm3")) {
        const Real value = finiteInRange(
            j, "intrinsic_carrier_density_cm3", 0.0, 1.0e30, context);
        mat.ni = units.m3ToInternalConcentration(value * 1.0e6);
    }
    if (j.contains("electron_mobility_cm2_per_V_s")) {
        const Real value = finiteInRange(
            j, "electron_mobility_cm2_per_V_s", 0.0, 1.0e9, context);
        mat.mun = units.m2PerVSToInternalMobility(value * 1.0e-4);
    }
    if (j.contains("hole_mobility_cm2_per_V_s")) {
        const Real value = finiteInRange(
            j, "hole_mobility_cm2_per_V_s", 0.0, 1.0e9, context);
        mat.mup = units.m2PerVSToInternalMobility(value * 1.0e-4);
    }
    if (j.contains("bandgap_eV"))
        mat.bandgap_eV = finiteInRange(j, "bandgap_eV", 0.0, 100.0, context);
    if (j.contains("electron_affinity_eV"))
        mat.electron_affinity_eV = finiteInRange(
            j, "electron_affinity_eV", 0.0, 100.0, context);
    if (j.contains("conduction_band_density_of_states_cm3")) {
        const Real value = finiteInRange(
            j, "conduction_band_density_of_states_cm3", 0.0, 1.0e30, context);
        mat.Nc_m3 = units.m3ToInternalConcentration(value * 1.0e6);
    }
    if (j.contains("valence_band_density_of_states_cm3")) {
        const Real value = finiteInRange(
            j, "valence_band_density_of_states_cm3", 0.0, 1.0e30, context);
        mat.Nv_m3 = units.m3ToInternalConcentration(value * 1.0e6);
    }
    if (j.contains("temperature_K"))
        mat.temperature_K = finiteInRange(j, "temperature_K", 1.0, 1.0e5, context);
    if (j.contains("electron_quantum_gamma"))
        mat.electron_quantum_gamma = finiteInRange(
            j, "electron_quantum_gamma", 0.0, 1.0e6, context);
    if (j.contains("electron_quantum_dos_mass_ratio"))
        mat.electron_quantum_dos_mass_ratio = finiteInRange(
            j, "electron_quantum_dos_mass_ratio", 1.0e-12, 1.0e6, context);
    if (j.contains("electron_quantum_coefficient_mass_ratio"))
        mat.electron_quantum_coefficient_mass_ratio = finiteInRange(
            j, "electron_quantum_coefficient_mass_ratio", 1.0e-12, 1.0e6, context);
    if (j.contains("thermal_conductivity_W_per_m_K"))
        mat.thermal_conductivity_W_per_m_K = finiteInRange(
            j, "thermal_conductivity_W_per_m_K", 0.0, 1.0e7, context);
    if (j.contains("specific_heat_J_per_kg_K"))
        mat.specific_heat_J_per_kg_K = finiteInRange(
            j, "specific_heat_J_per_kg_K", 0.0, 1.0e9, context);
    if (j.contains("mass_density_kg_per_m3"))
        mat.mass_density_kg_per_m3 = finiteInRange(
            j, "mass_density_kg_per_m3", 0.0, 1.0e8, context);
    return mat;
}

Material materialFromJson(const nlohmann::json& j,
                          const Material* base,
                          UnitScalingConfig scaling)
{
    Material mat = base != nullptr ? *base : Material{};
    mat.name = j.at("name").get<std::string>();
    if (j.contains("eps_r")) mat.eps_r = j.at("eps_r").get<Real>();
    if (j.contains("ni")) mat.ni = scaling.concentrationToInternal(j.at("ni").get<Real>());
    if (j.contains("mun")) mat.mun = scaling.mobilityToInternal(j.at("mun").get<Real>());
    if (j.contains("mup")) mat.mup = scaling.mobilityToInternal(j.at("mup").get<Real>());
    setOptionalReal(j, "bandgap_eV", mat.bandgap_eV);
    setOptionalReal(j, "electron_affinity_eV", mat.electron_affinity_eV);
    setOptionalScaledReal(
        j, "Nc_m3", mat.Nc_m3, &UnitScalingConfig::concentrationToInternal, scaling);
    setOptionalScaledReal(
        j, "Nv_m3", mat.Nv_m3, &UnitScalingConfig::concentrationToInternal, scaling);
    setOptionalReal(j, "temperature_K", mat.temperature_K);
    setOptionalReal(j, "electron_quantum_gamma", mat.electron_quantum_gamma);
    setOptionalReal(
        j, "electron_quantum_dos_mass_ratio",
        mat.electron_quantum_dos_mass_ratio);
    setOptionalReal(
        j, "electron_quantum_coefficient_mass_ratio",
        mat.electron_quantum_coefficient_mass_ratio);
    return mat;
}

std::vector<nlohmann::json> materialEntries(const nlohmann::json& root)
{
    if (root.is_array())
        return root.get<std::vector<nlohmann::json>>();
    if (root.is_object() && root.contains("materials")) {
        const auto& materials = root.at("materials");
        if (!materials.is_array())
            throw std::runtime_error(
                "MaterialDatabase: 'materials' must be an array of material objects.");
        return materials.get<std::vector<nlohmann::json>>();
    }
    if (root.is_object()) {
        std::vector<nlohmann::json> entries;
        for (auto it = root.begin(); it != root.end(); ++it) {
            if (!it.value().is_object())
                continue;
            nlohmann::json entry = it.value();
            if (!entry.contains("name"))
                entry["name"] = it.key();
            entries.push_back(std::move(entry));
        }
        return entries;
    }
    throw std::runtime_error(
        "MaterialDatabase: materials JSON must be an array, an object with a "
        "'materials' array, or an object map.");
}

} // namespace

MaterialDatabase::MaterialDatabase()
{
    // Silicon. All units are SI unless the field name explicitly says eV.
    Material si;
    si.name  = "Si";
    si.eps_r = 11.7;
    si.ni    = 1.0e16;   // [m^-3]
    si.mun   = 0.135;    // [m^2/V/s]
    si.mup   = 0.048;    // [m^2/V/s]
    si.bandgap_eV = 1.12;
    si.electron_affinity_eV = 4.05;
    si.Nc_m3 = 2.8e25;
    si.Nv_m3 = 1.04e25;
    si.temperature_K = 300.0;
    db_["Si"] = si;

    // Silicon dioxide (insulator - ni, mun, mup remain 0)
    Material sio2;
    sio2.name  = "SiO2";
    sio2.eps_r = 3.9;
    sio2.bandgap_eV = 9.0;
    sio2.electron_affinity_eV = 0.95;
    sio2.temperature_K = 300.0;
    db_["SiO2"] = sio2;
}
MaterialDatabase::MaterialDatabase(UnitScalingConfig scaling)
    : MaterialDatabase()
{
    if (!scaling.isUnitScaling())
        return;

    const PhysicalUnitSystem& units = scaling.unitSystem();
    for (auto& [_, mat] : db_) {
        mat.ni = units.m3ToInternalConcentration(mat.ni);
        mat.mun = units.m2PerVSToInternalMobility(mat.mun);
        mat.mup = units.m2PerVSToInternalMobility(mat.mup);
        if (mat.Nc_m3)
            mat.Nc_m3 = units.m3ToInternalConcentration(*mat.Nc_m3);
        if (mat.Nv_m3)
            mat.Nv_m3 = units.m3ToInternalConcentration(*mat.Nv_m3);
    }
}

MaterialDatabase::MaterialDatabase(const std::string& jsonPath)
    : MaterialDatabase()
{
    loadJson(jsonPath);
}

MaterialDatabase::MaterialDatabase(const std::string& jsonPath, UnitScalingConfig scaling)
    : MaterialDatabase(scaling)
{
    loadJson(jsonPath, scaling);
}

void MaterialDatabase::loadJson(const std::string& jsonPath)
{
    loadJson(jsonPath, UnitScalingConfig{});
}

void MaterialDatabase::loadJson(const std::string& jsonPath, UnitScalingConfig scaling)
{
    std::ifstream ifs(jsonPath);
    if (!ifs.is_open())
        throw std::runtime_error("MaterialDatabase: cannot open materials file: " + jsonPath);

    try {
        nlohmann::json root;
        ifs >> root;
        // A materials file inherits the deck's format version, so a version-2
        // deck also gets the version-2 key contract on its materials file.
        root = canonicalizeDeckKeys(root, scaling.deckFormatVersion);

        const bool versionedTemplatesLdmos = isTemplatesLdmosMaterialsV1(root);
        if (root.is_object() && root.contains("schema") &&
            !versionedTemplatesLdmos) {
            throw std::invalid_argument(
                "MaterialDatabase: unsupported materials schema '" +
                root.at("schema").get<std::string>() + "'.");
        }
        if (versionedTemplatesLdmos)
            validateTemplatesLdmosMaterialRoot(root);

        for (const nlohmann::json& entry : materialEntries(root)) {
            const std::string name = entry.at("name").get<std::string>();
            const Material* base = nullptr;
            auto it = db_.find(name);
            if (it != db_.end())
                base = &it->second;
            addMaterial(versionedTemplatesLdmos
                ? templatesLdmosMaterialFromJson(entry, base, scaling)
                : materialFromJson(entry, base, scaling));
        }
    } catch (const std::exception& e) {
        throw std::runtime_error(
            "MaterialDatabase: failed to load materials file '" + jsonPath + "': " + e.what());
    }
}

void MaterialDatabase::addMaterial(const Material& mat)
{
    if (mat.name.empty())
        throw std::invalid_argument("MaterialDatabase: material name must not be empty.");
    db_[mat.name] = mat;
}

const Material& MaterialDatabase::getMaterial(const std::string& name) const
{
    auto it = db_.find(name);
    if (it == db_.end())
        throw std::out_of_range("MaterialDatabase: unknown material '" + name + "'");
    return it->second;
}

Material MaterialDatabase::getMaterial(const std::string& name, Real temperature_K) const
{
    return getMaterial(name).atTemperature(temperature_K);
}

bool MaterialDatabase::hasMaterial(const std::string& name) const
{
    return db_.count(name) > 0;
}

} // namespace vela
