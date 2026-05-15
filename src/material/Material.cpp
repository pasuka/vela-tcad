#include "vela/material/Material.h"
#include "vela/core/PhysicalConstants.h"
#include <cmath>
#include <stdexcept>

namespace vela {

Material Material::atTemperature(Real targetTemperature_K) const
{
    if (targetTemperature_K <= 0.0)
        throw std::invalid_argument("Material::atTemperature: temperature_K must be positive.");

    const Real referenceTemperature = temperature_K.value_or(constants::T0);
    if (referenceTemperature <= 0.0)
        throw std::invalid_argument("Material::atTemperature: reference temperature_K must be positive.");

    Material adjusted = *this;
    adjusted.temperature_K = targetTemperature_K;
    const Real ratio = targetTemperature_K / referenceTemperature;

    if (ni > 0.0 && bandgap_eV.has_value()) {
        // Intrinsic density ratio from Nc,Nv ~ T^(3/2) and exp(-Eg/(2 kT)).
        const Real exponent = -(*bandgap_eV) * constants::q / (2.0 * constants::kb) *
            (1.0 / targetTemperature_K - 1.0 / referenceTemperature);
        adjusted.ni = ni * std::pow(ratio, 1.5) * std::exp(exponent);
    }

    if (mun > 0.0)
        adjusted.mun = mun * std::pow(ratio, -2.2);
    if (mup > 0.0)
        adjusted.mup = mup * std::pow(ratio, -2.2);
    if (Nc_m3.has_value())
        adjusted.Nc_m3 = *Nc_m3 * std::pow(ratio, 1.5);
    if (Nv_m3.has_value())
        adjusted.Nv_m3 = *Nv_m3 * std::pow(ratio, 1.5);

    return adjusted;
}

} // namespace vela
