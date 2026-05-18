#pragma once

#include "vela/core/Types.h"
#include <string>

namespace vela {

struct RegionFixedChargeSpec {
    std::string region;          ///< Region name (matches Region::name)
    Real        fixedCharge = 0; ///< Fixed charge density [m^-3], in units of q
};

struct InterfaceSheetChargeSpec {
    std::string region0;          ///< First region name adjacent to the interface
    std::string region1;          ///< Second region name adjacent to the interface
    Real        sheetCharge = 0;  ///< Legacy total sheet charge density [m^-2], in units of q
    Real        fixedCharge = 0;  ///< Fixed interface charge density [m^-2], in units of q
    Real        trapDensity = 0;  ///< Interface trap density [m^-2], in units of q when occupied
    Real        trapOccupancy = 0; ///< Occupied trap fraction [-]

    Real totalSheetCharge() const { return sheetCharge + fixedCharge + trapDensity * trapOccupancy; }
};

} // namespace vela
