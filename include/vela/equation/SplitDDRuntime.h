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
struct MobilityModelConfig;

// Explicit plain-PhuMob runtime. Keep the multiprecision implementation out of
// the solver headers. Geometry/physics are immutable; only the coordinate tail
// changes between Newton candidates.
class SplitDDRuntime {
    struct Impl;
    std::unique_ptr<Impl> impl_;
public:
    explicit SplitDDRuntime(nlohmann::json model,
        std::shared_ptr<const MobilityModelConfig> ialConfig = {});
    ~SplitDDRuntime();
    void setLow(const VectorXd& low);
    VectorXd low() const;
    VectorXd candidate(const VectorXd& x,const VectorXd& step,Real alpha,
                       const VectorXd& baseLow);
    /// Full high/low coordinate change, rounded only after the subtraction.
    VectorXd coordinateDelta(const VectorXd& x,const VectorXd& base,
                             const VectorXd& baseLow) const;
    /// Read-only source-free constant-material branch; excludes Si/HFS rows.
    nlohmann::json constantTransportDirection(const VectorXd& x,const VectorXd& direction,
        const CoupledDDBoundaryConditions& bcs) const;
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
// Build a normalized target without adding the tiny voltage in binary64.
void setSplitBoundaryTarget(CoupledDDBoundaryConditions& bcs,int component,
    Index node,Real builtIn,Real voltage,Real scale,Real voltageLow = 0.);
void shiftSplitBoundaryVoltage(CoupledDDBoundaryConditions& bcs,Index node,
    Real voltageDelta,Real scale,Real previousLow = 0.);
DDSolution predictSplitDDState(const DDSolution& previous,const DDSolution& current,
    Real ratio,const std::array<bool,3>& fields);
}
