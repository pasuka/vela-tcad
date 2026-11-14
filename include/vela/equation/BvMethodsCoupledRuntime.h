#pragma once
#include "vela/equation/CoupledStateRuntime.h"

namespace vela {
// Explicit 300 K Tri3 BVMethods runtime. One autonomous residual owns
// transport, SRH, avalanche and their jets; no native field/cache substitution.
class BvMethodsCoupledRuntime final : public CoupledStateRuntime {
    struct Impl;
    std::unique_ptr<Impl> impl_;
public:
    explicit BvMethodsCoupledRuntime(const CoupledDDAssembler& owner);
    ~BvMethodsCoupledRuntime();
    VectorXd transportResidual(const VectorXd& x) const override {
        return transportResidual(x, BvMethodsTransportPurpose::other);
    }
    VectorXd transportResidual(const VectorXd&, BvMethodsTransportPurpose) const;
    SparseMatrixd transportJacobian(const VectorXd& x) const override { return transportJacobian(x,nullptr); }
    SparseMatrixd transportJacobian(const VectorXd&,const VectorXd* columns) const;
    VectorXd residual(const VectorXd&, const CoupledDDBoundaryConditions&) const;
    SparseMatrixd jacobian(const VectorXd& x, const CoupledDDBoundaryConditions& b, bool approximate=false) const override { return jacobian(x,b,approximate,nullptr,true); }
    SparseMatrixd jacobian(const VectorXd&, const CoupledDDBoundaryConditions&, bool approximate, const VectorXd* columns, bool freezeStatisticalFloor=true) const;
    std::vector<CoupledDDCarrierTermDiagnostic> terms(
        const VectorXd&, const CoupledDDBoundaryConditions&, bool includeImpact=true) const;
    CoupledDDPhysicalDiagnostics physicalDiagnostics(
        const VectorXd&,const CoupledDDBoundaryConditions&) const;
    VectorXd avalancheResidual(const VectorXd&, const CoupledDDBoundaryConditions&) const;
    SparseMatrixd avalancheJacobian(const VectorXd&, const CoupledDDBoundaryConditions&) const;
    std::vector<CoupledDDPoissonTermDiagnostic> populationPoissonTerms(
        const VectorXd&, const CoupledDDBoundaryConditions&) const;
};
}
