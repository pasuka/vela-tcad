#pragma once
#include "vela/equation/ElectrothermalAssembler.h"
#include <set>
#include <stdexcept>

namespace vela::experimental {
// Partial derivative F_V at fixed state and QF references. Neutral potential
// shifts one volt per contact volt; equilibrium density is gauge invariant.
// Finite hole exchange therefore has no explicit contact-bias derivative.
inline VectorXd electrothermalContactBiasDerivative(Index nodes,
        const ElectrothermalBoundary& bc,const std::set<Index>& moving) {
    if(moving.empty())throw std::invalid_argument("Tangent requires moving neutral contacts");
    VectorXd derivative=VectorXd::Zero(4*nodes);
    for(Index i:moving) {
        if(i>=nodes || !bc.neutralContactBias_V.contains(i) || bc.potential_V.contains(i) ||
           bc.electronQf_V.contains(i) || bc.holeQf_V.contains(i))
            throw std::invalid_argument("Ambiguous tangent contact boundary");
        derivative[4*i]=derivative[4*i+1]=-1.;
        if(!bc.holeRecombination.contains(i))derivative[4*i+2]=-1.;
    }
    return derivative;
}
}
