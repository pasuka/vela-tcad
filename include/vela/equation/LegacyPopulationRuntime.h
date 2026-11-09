#pragma once
#include "vela/equation/CoupledStateRuntime.h"

namespace vela {
// Potential/density coordinates for the existing classical SLOT operator.
// QF fields are derived consumers; n/p remain authoritative stored unknowns.
class LegacyPopulationRuntime final : public CoupledStateRuntime {
    struct Impl;
    std::unique_ptr<Impl> impl_;
public:
    explicit LegacyPopulationRuntime(const CoupledDDAssembler&);
    ~LegacyPopulationRuntime() override;
    VectorXd transportResidual(const VectorXd&) const override;
    SparseMatrixd transportJacobian(const VectorXd&) const override;
    VectorXd residual(const VectorXd&,const CoupledDDBoundaryConditions&) const override;
    SparseMatrixd jacobian(const VectorXd&,const CoupledDDBoundaryConditions&,bool=false) const override;
    std::vector<CoupledDDCarrierTermDiagnostic> terms(const VectorXd&,const CoupledDDBoundaryConditions&,bool=true) const override;
    CoupledDDPhysicalDiagnostics physicalDiagnostics(const VectorXd&,const CoupledDDBoundaryConditions&) const override;
    VectorXd avalancheResidual(const VectorXd&,const CoupledDDBoundaryConditions&) const override;
    SparseMatrixd avalancheJacobian(const VectorXd&,const CoupledDDBoundaryConditions&) const override;
    std::vector<CoupledDDPoissonTermDiagnostic> populationPoissonTerms(const VectorXd&,const CoupledDDBoundaryConditions&) const override;
};
}
