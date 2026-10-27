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

struct PredictorStepDecision {
    std::optional<double> guarded;
    std::optional<double> applied;
};
// reportOnly keeps the ordinary outer controller's proposal. The shadow guard
// may report a sub-minimum step without rejecting the unmodified proposal.
inline PredictorStepDecision predictorStepDecision(double proposed,
    const OuterStepControl& control,double innerLimit,double deltaInner,
    double deltaOuter,bool reportOnly)
{
    const auto guarded=guardedOuterStep(proposed,control,innerLimit,deltaInner,deltaOuter);
    return {guarded,reportOnly?std::optional<double>(std::min(proposed,control.maximum)):guarded};
}
}
