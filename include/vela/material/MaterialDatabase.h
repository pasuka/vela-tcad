#pragma once

#include "vela/material/Material.h"
#include <unordered_map>
#include <string>
#include <stdexcept>

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

    /// Add or overwrite a material entry.
    void addMaterial(const Material& mat);

    /// Retrieve a material by name.
    /// @throws std::out_of_range if the material is not found.
    const Material& getMaterial(const std::string& name) const;

    /// Check whether a material is registered.
    bool hasMaterial(const std::string& name) const;

private:
    std::unordered_map<std::string, Material> db_;
};

} // namespace vela
