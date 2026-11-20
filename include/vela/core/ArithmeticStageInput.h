#pragma once

#include "vela/core/ArithmeticDecimal.h"
#include <nlohmann/json.hpp>

namespace vela {

inline ArithmeticReal decimalCompatibilityNumber(const nlohmann::json& input,
                                                  const std::string& label, double compatibility)
{
    if constexpr (compiledArithmeticPrecision != ArithmeticPrecision::Extended80)
        throw std::invalid_argument(label + " requires an extended80 build");
    const auto text = input.get<std::string>();
    if (parseArithmeticDecimal<double>(text) != compatibility)
        throw std::invalid_argument(label + " differs from deck compatibility value");
    return parseArithmeticDecimal(text);
}

inline ArithmeticReal decimalStageControl(const nlohmann::json& input,
                                          const std::string& label, double compatibility)
{
    return decimalCompatibilityNumber(input, "decimal_stage_controls." + label, compatibility);
}

// Read the authoritative decimal form and verify that it describes the deck's
// compatibility value. Callers with legacy state types narrow at their boundary.
inline ArithmeticReal stageControl(const nlohmann::json& deck, const char* group,
                                   const std::string& name, double compatibility)
{
    const auto metadata = deck.find("decimal_stage_controls");
    if (metadata == deck.end())
        return compatibility;
    return decimalStageControl(metadata->at(group).at(name),
                               std::string(group) + "." + name, compatibility);
}

} // namespace vela
