#pragma once

#include <stdexcept>
#include <string_view>

namespace vela {

/// Field used by GradQuasiFermi models in elements touching an electrode.
enum class GradQfContactField { ElectricField, QuasiFermiGradient };

inline GradQfContactField gradQfContactFieldFromString(std::string_view value)
{
    if (value == "electric_field")
        return GradQfContactField::ElectricField;
    if (value == "quasi_fermi_gradient")
        return GradQfContactField::QuasiFermiGradient;
    throw std::invalid_argument(
        "solver.grad_qf_contact_field must be electric_field or quasi_fermi_gradient");
}

} // namespace vela
