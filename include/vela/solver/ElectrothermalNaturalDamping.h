#pragma once
#include "vela/core/Types.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace vela::experimental {
struct NaturalNewtonTrial {
    Real theta=0.,nextAlpha=0.;
    bool decreasing=false;
};
// Error-oriented trial test using the CURRENT factorization for both vectors:
// d=-J_k^-1 F_k, correction=-J_k^-1 F(candidate), in identical column units.
// This implements the NLEQ_ERR corrector test and local shrinking estimate,
// not its cross-iteration predictor or its optional extra correction step.
inline NaturalNewtonTrial electrothermalNaturalTrial(const VectorXd& direction,
                                                     const VectorXd& correction,
                                                     Real alpha,bool projected=false) {
    if(direction.size()!=correction.size() || !direction.allFinite() || !correction.allFinite() ||
       !(alpha>0. && alpha<=1.))throw std::invalid_argument("Invalid natural Newton trial");
    const Real norm=direction.stableNorm();
    if(!(norm>0. && std::isfinite(norm)))throw std::invalid_argument("Zero or nonfinite Newton direction");
    const Real theta=correction.stableNorm()/norm;
    const Real defect=(correction-(1.-alpha)*direction).stableNorm();
    Real next=.5*alpha;
    // A projected path is not the straight Newton path used by this estimate.
    if(!projected && defect>0.)next=std::min(next,.5*norm*alpha*alpha/defect);
    return {theta,next,std::isfinite(theta)&&theta<1.};
}
}
