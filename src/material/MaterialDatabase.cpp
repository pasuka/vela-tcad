#include "vela/material/MaterialDatabase.h"

#include <fstream>
#include <nlohmann/json.hpp>
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

Material materialFromJson(const nlohmann::json& j, const Material* base = nullptr)
{
    Material mat = base != nullptr ? *base : Material{};
    mat.name = j.at("name").get<std::string>();
    if (j.contains("eps_r")) mat.eps_r = j.at("eps_r").get<Real>();
    if (j.contains("ni")) mat.ni = j.at("ni").get<Real>();
    if (j.contains("mun")) mat.mun = j.at("mun").get<Real>();
    if (j.contains("mup")) mat.mup = j.at("mup").get<Real>();
    setOptionalReal(j, "bandgap_eV", mat.bandgap_eV);
    setOptionalReal(j, "electron_affinity_eV", mat.electron_affinity_eV);
    setOptionalReal(j, "Nc_m3", mat.Nc_m3);
    setOptionalReal(j, "Nv_m3", mat.Nv_m3);
    setOptionalReal(j, "temperature_K", mat.temperature_K);
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

MaterialDatabase::MaterialDatabase(const std::string& jsonPath)
    : MaterialDatabase()
{
    loadJson(jsonPath);
}

void MaterialDatabase::loadJson(const std::string& jsonPath)
{
    std::ifstream ifs(jsonPath);
    if (!ifs.is_open())
        throw std::runtime_error("MaterialDatabase: cannot open materials file: " + jsonPath);

    try {
        nlohmann::json root;
        ifs >> root;

        for (const nlohmann::json& entry : materialEntries(root)) {
            const std::string name = entry.at("name").get<std::string>();
            const Material* base = nullptr;
            auto it = db_.find(name);
            if (it != db_.end())
                base = &it->second;
            addMaterial(materialFromJson(entry, base));
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
