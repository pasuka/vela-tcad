#pragma once
#include "vela/core/Types.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace vela::experimental {
struct NearSteadyRebase {
    std::array<unsigned,2> changed{};
    Real maxRoundingEstimate_V=0.;
};
// Representation only: preserve psi/T and the physical QF sum. Small physical
// QFs should not remain the cancellation of a reference and a large increment.
inline NearSteadyRebase rebaseSmallElectrothermalQf(VectorXd& x,VectorXd& electrons,VectorXd& holes) {
    if(x.size()!=4*electrons.size() || holes.size()!=electrons.size() ||
       !x.allFinite() || !electrons.allFinite() || !holes.allFinite())
        throw std::invalid_argument("Invalid near-steady referenced state");
    NearSteadyRebase result;
    for(Eigen::Index i=0;i<electrons.size();++i)for(int k=0;k<2;++k) {
        Real& reference=k==0?electrons[i]:holes[i];Real& increment=x[4*i+k+1];
        const Real physical=reference+increment;
        if(reference!=0. && std::abs(physical)<=1e-3) {
            const long double error=(static_cast<long double>(reference)+increment)-physical;
            result.maxRoundingEstimate_V=std::max(result.maxRoundingEstimate_V,static_cast<Real>(std::abs(error)));
            reference=0.;increment=physical;++result.changed[k];
        }
    }
    return result;
}
}
