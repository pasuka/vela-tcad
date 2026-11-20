#pragma once

#include "vela/core/UnitScaling.h"
#include "vela/material/Material.h"
#include "vela/core/ArithmeticPrecision.h"
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
    struct DecimalSemiconductor {
        ArithmeticReal temperature_K;
        ArithmeticReal bandgap_eV, affinity_eV;
        ArithmeticReal electronsDOS_internal, holesDOS_internal, intrinsicDensity_internal;
        ArithmeticReal boundReferencePotential_V;
    };
    MaterialDatabase();
    explicit MaterialDatabase(UnitScalingConfig scaling);
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

    /// Exact operating-temperature values from the active decimal PAR
    /// binding. Returns null when the bound material or temperature differs.
    const DecimalSemiconductor* decimalSemiconductor(const std::string& name,
                                                     ArithmeticReal temperature_K) const;

    /// Check whether a material is registered.
    bool hasMaterial(const std::string& name) const;

private:
    std::unordered_map<std::string, Material> db_;
    std::unordered_map<std::string, DecimalSemiconductor> decimalSemiconductors_;
};

} // namespace vela
