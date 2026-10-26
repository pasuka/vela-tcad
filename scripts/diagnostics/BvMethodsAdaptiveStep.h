#pragma once
#include "../../include/vela/simulation/AdaptiveVoltageStepControl.h"
namespace vela::diagnostics {
using detail::InnerStepLimit;
using detail::OuterStepControl;
using detail::defaultPredictedInnerStep_V;
using detail::stepNumber;
using detail::positiveStepNumber;
using detail::stepObject;
using detail::resolveInnerStepLimit;
using detail::resolveOuterStepControl;
using detail::guardedOuterStep;
}
