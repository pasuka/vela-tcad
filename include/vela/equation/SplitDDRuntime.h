#pragma once
#include "vela/core/Types.h"
#include <nlohmann/json_fwd.hpp>
#include <memory>
#include <string>
#include <vector>
#include <array>

namespace vela {
struct CoupledDDBoundaryConditions;
struct CoupledDDCarrierTermDiagnostic;
struct CoupledDDEdgeFluxDiagnostic;
struct DDSolution;
struct ContactCurrentDetailedResult;

// Explicit plain-PhuMob runtime. Keep the multiprecision implementation out of
// the solver headers. Geometry/physics are immutable; only the coordinate tail
// changes between Newton candidates.
class SplitDDRuntime {
    struct Impl;
    std::unique_ptr<Impl> impl_;
public:
    explicit SplitDDRuntime(nlohmann::json model);
    ~SplitDDRuntime();
    void setLow(const VectorXd& low);
    VectorXd low() const;
    VectorXd candidate(const VectorXd& x,const VectorXd& step,Real alpha,
                       const VectorXd& baseLow);
    VectorXd residual(const VectorXd& x,const CoupledDDBoundaryConditions& bcs) const;
    VectorXd symmetricDifference(const VectorXd& x,const VectorXd& step,
        const CoupledDDBoundaryConditions& bcs,VectorXd& forward,VectorXd& backward);
    VectorXd density(const VectorXd& x,bool electron) const;
    std::vector<CoupledDDCarrierTermDiagnostic> terms(
        const VectorXd& x,const CoupledDDBoundaryConditions& bcs) const;
    std::vector<CoupledDDEdgeFluxDiagnostic> edges(const VectorXd& x) const;
    Real current(const VectorXd& x,const std::string& contact) const;
    ContactCurrentDetailedResult contact(const VectorXd& x,const std::string& name) const;
    void save(const VectorXd& x,DDSolution& solution) const;
    VectorXd restore(const DDSolution& solution);
    // Explicit warm continuation only. Strict checkpoint restore stays strict.
    VectorXd restoreForContinuation(const DDSolution& solution);
    void applyBoundary(VectorXd& x,const CoupledDDBoundaryConditions& bcs);
    const std::string& fingerprint() const;
};
DDSolution predictSplitDDState(const DDSolution& previous,const DDSolution& current,
    Real ratio,const std::array<bool,3>& fields);
}
