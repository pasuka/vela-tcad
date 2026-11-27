#include "vela/material/MaterialDatabase.h"
#include "vela/material/SemiconductorArithmetic.h"
#include "vela/core/ArithmeticDecimal.h"
#include "vela/core/ArithmeticStageInput.h"
#include "vela/core/PhysicalConstants.h"

#include <array>
#include <fstream>
#include <cmath>
#include <nlohmann/json.hpp>
#include <set>
#include <stdexcept>
#include <vector>

namespace vela {

SemiconductorArithmetic evaluateSemiconductorArithmetic(
    const nlohmann::json& material, ArithmeticReal temperature)
{
    if (!(temperature > 0) || !std::isfinite(temperature))
        throw std::invalid_argument("Arithmetic material temperature must be positive and finite");
    using Real = vela::ArithmeticReal;
    const auto& sections = material.at("sections");
    const auto param = [](const nlohmann::json& section, const char* name) {
        return parseArithmeticDecimal(section.at(name).at("literals").at(0).get<std::string>());
    };
    const auto paramOrZero = [&param](const nlohmann::json& section, const char* name) {
        return section.contains(name) ? param(section, name) : Real{0};
    };
    const auto& band = sections.at("Bandgap");
    const auto& electron = sections.at("eDOSMass");
    const auto& hole = sections.at("hDOSMass");
    if (param(electron, "Formula") != 1 || param(hole, "Formula") != 1)
        throw std::invalid_argument("Arithmetic DOS input requires Formula 1");

    const Real tpar = param(band, "Tpar");
    const Real beta = param(band, "beta");
    const Real alpha = param(band, "alpha");
    if (beta + tpar <= 0 || beta + temperature <= 0)
        throw std::invalid_argument("Arithmetic Bandgap input has nonpositive Varshni denominator");
    const Real gapZero = param(band, "Eg0") + alpha * tpar * tpar / (beta + tpar);
    const Real gap = gapZero - alpha * temperature * temperature / (beta + temperature);
    Real offset = 0;
    if (material.at("fold_old_slotboom").get<bool>())
        offset = paramOrZero(band, "dEg0(OldSlotboom)");
    else if (material.value("fold_bennett", false))
        offset = paramOrZero(band, "dEg0(Bennett)");
    if (gap + offset <= 0 || gapZero + offset <= 0)
        throw std::invalid_argument("Arithmetic DOS input has nonpositive effective gap");
    const Real mt = param(electron, "a") * (gapZero + offset) / (gap + offset);
    const Real sixMt = 6 * mt;
    const Real me = std::cbrt(sixMt * sixMt * param(electron, "ml")) + param(electron, "mm");

    Real numerator = param(hole, "a");
    Real denominator = 1;
    Real power = temperature;
    const char* numeratorFields[]{"b", "c", "d", "e"};
    const char* denominatorFields[]{"f", "g", "h", "i"};
    for (int i = 0; i < 4; ++i) {
        numerator += paramOrZero(hole, numeratorFields[i]) * power;
        denominator += paramOrZero(hole, denominatorFields[i]) * power;
        power *= temperature;
    }
    if (numerator < 0 || denominator <= 0)
        throw std::invalid_argument("Arithmetic DOS input has invalid hole mass ratio");
    const Real mh = std::pow(numerator / denominator, Real{2} / 3) + paramOrZero(hole, "mm");
    if (me <= 0 || mh <= 0) throw std::invalid_argument("Arithmetic DOS masses must be positive");

    const auto constants = vela::constants::arithmeticConstants<Real>();
    const auto dos = [temperature, constants](Real mass) {
        const Real x = 2 * constants.dosPi * (mass * constants.m0) * constants.kb * temperature /
                       (constants.h * constants.h);
        return 2 * x * std::sqrt(x) / 1000000;
    };
    const Real nc = dos(me);
    const Real nv = dos(mh);
    const Real vt = vela::constants::thermalVoltage(temperature);
    const Real ni = std::sqrt(nc * nv) * std::exp(-(gap + offset) / (2 * vt));
    const Real affinity = param(band, "Chi0");
    const Real ratio = me / mh;
    const Real reference = affinity + gapZero / 2 + vt / 2 * std::log(ratio * std::sqrt(ratio));
    if (!(nc > 0) || !(nv > 0) || !(ni > 0) || !std::isfinite(nc) ||
        !std::isfinite(nv) || !std::isfinite(ni) || !std::isfinite(reference))
        throw std::invalid_argument("Arithmetic semiconductor outputs must be finite with positive concentrations");
    if (material.value("poly_bound_convention", false)) {
        return {gap + offset, affinity + param(band, "Bgn2Chi") * (gapZero - gap),
                nc, nv, ni, reference};
    }
    return {gap, affinity, nc, nv, ni, reference};
}

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
        {"schema", "benchmark", "revision", "unit_system", "materials", "physical_constants_profile"},
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

std::optional<MaterialDatabase::DecimalSemiconductor> applyDecimalSemiconductor(
    const nlohmann::json& j, Material& mat, UnitScalingConfig scaling)
{
    if (!j.contains("decimal_semiconductor"))
        return std::nullopt;
    if (compiledArithmeticPrecision != ArithmeticPrecision::Extended80)
        throw std::invalid_argument("decimal_semiconductor requires an extended80 build");
    if (j.value("role", std::string{}) == "dielectric")
        throw std::invalid_argument("decimal_semiconductor requires a semiconductor material");
    const auto& input = j.at("decimal_semiconductor");
    const auto temperatureText = input.at("temperature_K").get<std::string>();
    const auto temperature = parseArithmeticDecimal(temperatureText);
    if (temperature < 1 || temperature > 1.0e5 || !mat.temperature_K ||
        parseArithmeticDecimal<Real>(temperatureText) != *mat.temperature_K)
        throw std::invalid_argument("decimal_semiconductor temperature differs from the material deck");
    const auto values = evaluateSemiconductorArithmetic(input, temperature);
    const auto inRange = [](ArithmeticReal value, ArithmeticReal maximum) {
        return value >= 0 && value <= maximum;
    };
    if (!inRange(values.bandgap_eV, 100) || !inRange(values.affinity_eV, 100) ||
        !inRange(values.electronsDOS_cm3, VELA_ARITHMETIC_LITERAL(1e30)) ||
        !inRange(values.holesDOS_cm3, VELA_ARITHMETIC_LITERAL(1e30)) ||
        !inRange(values.intrinsicDensity_cm3, VELA_ARITHMETIC_LITERAL(1e30)))
        throw std::invalid_argument("decimal_semiconductor derived values are outside the material ranges");
    const auto concentrationFactor =
        VELA_ARITHMETIC_LITERAL(1e6) / scaling.unitSystem().concentrationM3PerInternal();
    const MaterialDatabase::DecimalSemiconductor exact{
        temperature, values.bandgap_eV, values.affinity_eV,
        values.electronsDOS_cm3 * concentrationFactor,
        values.holesDOS_cm3 * concentrationFactor,
        values.intrinsicDensity_cm3 * concentrationFactor,
        j.contains("decimal_intrinsic_reference_potential_V")
            ? decimalCompatibilityNumber(j.at("decimal_intrinsic_reference_potential_V"),
                mat.name + ".intrinsic_reference_potential_V", mat.intrinsic_reference_potential_V)
            : static_cast<ArithmeticReal>(mat.intrinsic_reference_potential_V)
    };
    mat.bandgap_eV = static_cast<Real>(exact.bandgap_eV);
    mat.electron_affinity_eV = static_cast<Real>(exact.affinity_eV);
    mat.Nc_m3 = static_cast<Real>(exact.electronsDOS_internal);
    mat.Nv_m3 = static_cast<Real>(exact.holesDOS_internal);
    mat.ni = static_cast<Real>(exact.intrinsicDensity_internal);
    return exact;
}

Material templatesLdmosMaterialFromJson(const nlohmann::json& j,
                                        const Material* base,
                                        UnitScalingConfig scaling,
                                        std::optional<MaterialDatabase::DecimalSemiconductor>& decimal)
{
    const std::string context = "MaterialDatabase: versioned material";
    rejectUnknownKeys(
        j,
        {"name", "role", "eps_r", "intrinsic_carrier_density_cm3",
         "electron_mobility_cm2_per_V_s", "hole_mobility_cm2_per_V_s",
         "bandgap_eV", "electron_affinity_eV", "intrinsic_reference_potential_V",
         "conduction_band_density_of_states_cm3",
         "valence_band_density_of_states_cm3", "temperature_K",
         "electron_quantum_gamma", "electron_quantum_affinity_eV",
         "electron_quantum_dos_mass_ratio",
         "electron_quantum_coefficient_mass_ratio",
         "thermal_conductivity_W_per_m_K", "specific_heat_J_per_kg_K",
         "mass_density_kg_per_m3", "provenance", "decimal_semiconductor",
         "decimal_eps_r", "decimal_intrinsic_reference_potential_V", "decimal_mobility_cm2_per_V_s"},
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
    if (const auto exact = j.find("decimal_mobility_cm2_per_V_s"); exact != j.end()) {
        if (!exact->is_array() || exact->size() != 2)
            throw std::invalid_argument(context + ".decimal_mobility_cm2_per_V_s requires two decimal literals.");
        const std::array keys{"electron_mobility_cm2_per_V_s", "hole_mobility_cm2_per_V_s"};
        for (std::size_t i = 0; i < keys.size(); ++i) {
            if (!exact->at(i).is_string() ||
                parseArithmeticDecimal<double>(exact->at(i).get<std::string>()) != j.at(keys[i]).get<Real>())
                throw std::invalid_argument(context + ".decimal_mobility_cm2_per_V_s differs from the numeric fields.");
        }
    }
    if (j.contains("bandgap_eV"))
        mat.bandgap_eV = finiteInRange(j, "bandgap_eV", 0.0, 100.0, context);
    if (j.contains("intrinsic_reference_potential_V"))
        mat.intrinsic_reference_potential_V = finiteInRange(
            j, "intrinsic_reference_potential_V", -10.0, 10.0, context);
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
    decimal = applyDecimalSemiconductor(j, mat, scaling);
    if (j.contains("electron_quantum_gamma"))
        mat.electron_quantum_gamma = finiteInRange(
            j, "electron_quantum_gamma", 0.0, 1.0e6, context);
    if (j.contains("electron_quantum_affinity_eV"))
        mat.electron_quantum_affinity_eV = finiteInRange(
            j, "electron_quantum_affinity_eV", 0.0, 100.0, context);
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
                          UnitScalingConfig scaling,
                          std::optional<MaterialDatabase::DecimalSemiconductor>& decimal)
{
    Material mat = base != nullptr ? *base : Material{};
    mat.name = j.at("name").get<std::string>();
    if (j.contains("eps_r")) mat.eps_r = j.at("eps_r").get<Real>();
    if (j.contains("ni")) mat.ni = scaling.concentrationToInternal(j.at("ni").get<Real>());
    if (j.contains("mun")) mat.mun = scaling.mobilityToInternal(j.at("mun").get<Real>());
    if (j.contains("mup")) mat.mup = scaling.mobilityToInternal(j.at("mup").get<Real>());
    setOptionalReal(j, "bandgap_eV", mat.bandgap_eV);
    setOptionalReal(j, "electron_affinity_eV", mat.electron_affinity_eV);
    if (j.contains("intrinsic_reference_potential_V"))
        mat.intrinsic_reference_potential_V = finiteInRange(
            j, "intrinsic_reference_potential_V", -10.0, 10.0, "MaterialDatabase");
    setOptionalScaledReal(
        j, "Nc_m3", mat.Nc_m3, &UnitScalingConfig::concentrationToInternal, scaling);
    setOptionalScaledReal(
        j, "Nv_m3", mat.Nv_m3, &UnitScalingConfig::concentrationToInternal, scaling);
    setOptionalReal(j, "temperature_K", mat.temperature_K);
    setOptionalReal(j, "electron_quantum_gamma", mat.electron_quantum_gamma);
    setOptionalReal(j, "electron_quantum_affinity_eV", mat.electron_quantum_affinity_eV);
    setOptionalReal(
        j, "electron_quantum_dos_mass_ratio",
        mat.electron_quantum_dos_mass_ratio);
    setOptionalReal(
        j, "electron_quantum_coefficient_mass_ratio",
        mat.electron_quantum_coefficient_mass_ratio);
    decimal = applyDecimalSemiconductor(j, mat, scaling);
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
    // No-BGN 300 K snapshot emitted from share/vela/materials/models.par.
    // All units are SI unless the field name explicitly says eV.
    Material si;
    si.name  = "Si";
    si.eps_r = 11.7;
    si.mun   = 0.1417;   // [m^2/V/s]
    si.mup   = 0.04705;  // [m^2/V/s]
    si.bandgap_eV = 1.1241592307692307;
    si.electron_affinity_eV = 4.05;
#if defined(VELA_CONSTANTS_SENT_U2022_12) && VELA_CONSTANTS_SENT_U2022_12
    // Same Formula-1 no-BGN material law, regenerated under this profile.
    si.Nc_m3 = 2.856653288716314e25;
    si.Nv_m3 = 3.1046289324748374e25;
#else
    si.Nc_m3 = 2.8566791451911786e25;
    si.Nv_m3 = 3.104657033455896e25;
#endif
    si.temperature_K = 300.0;
    si.ni = std::sqrt(*si.Nc_m3 * *si.Nv_m3) *
        std::exp(-*si.bandgap_eV * constants::q /
                 (2.0 * constants::kb * *si.temperature_K));
    db_["Si"] = si;

    // Silicon dioxide (insulator - ni, mun, mup remain 0)
    Material sio2;
    sio2.name  = "SiO2";
    sio2.eps_r = 3.9;
    sio2.bandgap_eV = 9.0;
    sio2.electron_affinity_eV = 0.9;
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
        constants::requireProfile(root.is_object() ? root.value("physical_constants_profile", std::string("codata2018")) : std::string("codata2018"));
        const bool versionedTemplatesLdmos = isTemplatesLdmosMaterialsV1(root);
        // Versioned material documents declare their own TCAD unit contract.
        // Only unversioned overrides inherit the deck's v1/v2 key spelling;
        // rewriting C3 here would turn its native *_cm3 fields into *_m3.
        if (!versionedTemplatesLdmos)
            root = canonicalizeDeckKeys(root, scaling.deckFormatVersion);
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
            std::optional<DecimalSemiconductor> decimal;
            addMaterial(versionedTemplatesLdmos
                ? templatesLdmosMaterialFromJson(entry, base, scaling, decimal)
                : materialFromJson(entry, base, scaling, decimal));
            if (decimal)
                decimalSemiconductors_.insert_or_assign(name, *decimal);
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
    decimalSemiconductors_.erase(mat.name);
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

const MaterialDatabase::DecimalSemiconductor* MaterialDatabase::decimalSemiconductor(
    const std::string& name, ArithmeticReal temperature_K) const
{
    const auto it = decimalSemiconductors_.find(name);
    return it != decimalSemiconductors_.end() && it->second.temperature_K == temperature_K
        ? &it->second : nullptr;
}

bool MaterialDatabase::hasMaterial(const std::string& name) const
{
    return db_.count(name) > 0;
}

} // namespace vela
