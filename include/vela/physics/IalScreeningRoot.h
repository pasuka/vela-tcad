#pragma once
#include "vela/core/Types.h"
#include <string>

namespace vela {
// Safeguarded Halley is the default; Legacy preserves the original floating-point path.
enum class IalScreeningMethod { Legacy, Newton, Halley, Toms748, U202212SingleNewton };
inline constexpr auto defaultIalScreeningMethod = IalScreeningMethod::Halley;
IalScreeningMethod ialScreeningMethod(const std::string& name);
const char* ialScreeningMethodName(IalScreeningMethod method);
// Compatibility selection returns the observed one-step cutoff, not a stationary root.
Real ialScreeningMinimum(Real mass, Real temperature,
                         IalScreeningMethod method=defaultIalScreeningMethod);
} // namespace vela
