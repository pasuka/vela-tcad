#pragma once

#include "vela/core/UnitScaling.h"
#include "vela/material/Material.h"
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace vela {

/**
 * @brief Look-up table of built-in semiconductor materials.
 *
 * Pre-loaded with Si and SiO2.  Additional materials can be registered
 * at runtime via addMaterial().
 */
class MaterialDatabase {
public:
    MaterialDatabase();
    explicit MaterialDatabase(const std::string& jsonPath);
    MaterialDatabase(const std::string& jsonPath, UnitScalingConfig scaling);

    /// Load material entries from JSON, adding new entries or overriding built-ins.
    void loadJson(const std::string& jsonPath);
    void loadJson(const std::string& jsonPath, UnitScalingConfig scaling);

    /// Add or overwrite a material entry.
    void addMaterial(const Material& mat);

    /// Retrieve a material by name.
    /// @throws std::out_of_range if the material is not found.
    const Material& getMaterial(const std::string& name) const;

    /// Retrieve a temperature-adjusted material copy.
    Material getMaterial(const std::string& name, Real temperature_K) const;

    /// Check whether a material is registered.
    bool hasMaterial(const std::string& name) const;

private:
    std::unordered_map<std::string, Material> db_;
};

} // namespace vela
