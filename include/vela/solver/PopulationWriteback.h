#pragma once
#include "vela/solver/PopulationState.h"
#include <limits>

namespace vela {
// Boltzmann, fixed-prefactor finite map. Correct the population coordinates
// for rounding of physical psi, preserving the intended finite QF change.
// This is not a new tangent or an acceptance criterion.
inline PopulationTrial compensateBoltzmannPsiWriteback(
    const VectorXd& before, const VectorXd& raw, PopulationTrial trial,
    Real potentialScale, Real vt, const std::vector<bool>& active) {
    const int n=static_cast<int>(active.size());
    if(n==0 || before.size()!=3*n || raw.size()!=before.size() ||
       trial.state.size()!=before.size() || !before.allFinite() ||
       !raw.allFinite() || !trial.state.allFinite() ||
       !(trial.alpha>0. && trial.alpha<=1.) ||
       !(potentialScale>0.) || !std::isfinite(potentialScale) ||
       !(vt>0.) || !std::isfinite(vt))
        throw std::invalid_argument("invalid Boltzmann psi writeback inputs");
    if(raw.isZero(0.))return trial;
    for(int i=0;i<n;++i)if(active[i]) {
        // Round physical psi in exactly the same way as the state consumer.
        const Real psiBefore=before[i]*potentialScale;
        const Real psiWritten=trial.state[i]*potentialScale;
        const long double error=static_cast<long double>(psiWritten)-psiBefore-
            static_cast<long double>(trial.alpha)*raw[i]*potentialScale;
        for(int block:{1,2}) {
            const int k=block*n+i;
            const auto& protectedNodes=block==1?trial.protectedElectrons:trial.protectedHoles;
            // Local nonpositive projection retains ownership of its finite map.
            if(std::find(protectedNodes.begin(),protectedNodes.end(),i)!=protectedNodes.end())continue;
            if(!(before[k]>0.))throw std::domain_error("nonpositive Boltzmann writeback parent");
            const long double ideal=static_cast<long double>(before[k])+static_cast<long double>(trial.alpha)*raw[k];
            if(!(ideal>0.))continue; // Original protection/constraint owns this domain.
            const long double value=ideal*std::exp((block==1?error:-error)/vt);
            if(!(value>0.) || !std::isfinite(value) || value>std::numeric_limits<Real>::max())
                throw std::domain_error("unrepresentable compensated Boltzmann population");
            const Real written=static_cast<Real>(value);
            if(!(written>0.) || !std::isfinite(written))
                throw std::domain_error("compensated Boltzmann population underflow");
            trial.state[k]=written;
        }
    }
    return trial;
}
} // namespace vela
