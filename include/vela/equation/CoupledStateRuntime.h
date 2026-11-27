#pragma once
#include "vela/equation/CoupledDDAssembler.h"

namespace vela {
// A coordinate runtime must retain the owning assembler's physical model.
class CoupledStateRuntime {
public:
    virtual ~CoupledStateRuntime() = default;
    virtual VectorXd transportResidual(const VectorXd&) const = 0;
    virtual SparseMatrixd transportJacobian(const VectorXd&) const = 0;
    virtual VectorXd residual(const VectorXd&,const CoupledDDBoundaryConditions&) const = 0;
    virtual SparseMatrixd jacobian(const VectorXd&,const CoupledDDBoundaryConditions&,bool=false) const = 0;
    virtual std::vector<CoupledDDCarrierTermDiagnostic> terms(const VectorXd&,const CoupledDDBoundaryConditions&,bool=true) const = 0;
    virtual CoupledDDPhysicalDiagnostics physicalDiagnostics(const VectorXd&,const CoupledDDBoundaryConditions&) const = 0;
    virtual VectorXd avalancheResidual(const VectorXd&,const CoupledDDBoundaryConditions&) const = 0;
    virtual SparseMatrixd avalancheJacobian(const VectorXd&,const CoupledDDBoundaryConditions&) const = 0;
    virtual std::vector<CoupledDDPoissonTermDiagnostic> populationPoissonTerms(const VectorXd&,const CoupledDDBoundaryConditions&) const = 0;
};
}
