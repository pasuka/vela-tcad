#include "vela/simulation/CurveSweep.h"
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace vela {

std::string toString(CurveSweepMode mode)
{
    switch (mode) {
    case CurveSweepMode::IV:
        return "iv";
    case CurveSweepMode::CVQuasistatic:
        return "cv_quasistatic";
    case CurveSweepMode::BVReverse:
        return "bv_reverse";
    }
    return "iv";
}

CurveSweepMode curveSweepModeFromString(const std::string& mode)
{
    std::string normalized = mode;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    std::replace(normalized.begin(), normalized.end(), '-', '_');

    if (normalized.empty() || normalized == "iv")
        return CurveSweepMode::IV;
    if (normalized == "cv" || normalized == "cv_quasistatic")
        return CurveSweepMode::CVQuasistatic;
    if (normalized == "bv" || normalized == "bv_reverse" || normalized == "reverse_breakdown")
        return CurveSweepMode::BVReverse;

    throw std::invalid_argument(
        "CurveSweep: sweep.mode must be 'iv', 'cv_quasistatic', or 'bv_reverse'.");
}

} // namespace vela
