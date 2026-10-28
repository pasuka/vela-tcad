#pragma once
#include "vela/equation/CoupledDDAssembler.h"

namespace vela {
// Explicit 300 K Tri3 BVMethods runtime. One autonomous residual owns
// transport, SRH, avalanche and their jets; no native field/cache substitution.
class BvMethodsCoupledRuntime {
    struct Impl;
    std::unique_ptr<Impl> impl_;
public:
    explicit BvMethodsCoupledRuntime(const CoupledDDAssembler& owner);
    ~BvMethodsCoupledRuntime();
    VectorXd transportResidual(const VectorXd&) const;
    SparseMatrixd transportJacobian(const VectorXd&) const;
    VectorXd residual(const VectorXd&, const CoupledDDBoundaryConditions&) const;
    SparseMatrixd jacobian(const VectorXd&, const CoupledDDBoundaryConditions&, bool approximate=false) const;
    std::vector<CoupledDDCarrierTermDiagnostic> terms(
        const VectorXd&, const CoupledDDBoundaryConditions&, bool includeImpact=true) const;
    VectorXd avalancheResidual(const VectorXd&, const CoupledDDBoundaryConditions&) const;
    SparseMatrixd avalancheJacobian(const VectorXd&, const CoupledDDBoundaryConditions&) const;
    std::vector<CoupledDDPoissonTermDiagnostic> populationPoissonTerms(
        const VectorXd&, const CoupledDDBoundaryConditions&) const;
};
}
