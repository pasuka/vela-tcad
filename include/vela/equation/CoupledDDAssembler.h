#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include <unordered_map>
#include <vector>

namespace vela {

struct CoupledDDState {
    VectorXd psi;
    VectorXd phin;
    VectorXd phip;
};

struct CoupledDDBoundaryConditions {
    std::unordered_map<Index, Real> psi;
    std::unordered_map<Index, Real> phin;
    std::unordered_map<Index, Real> phip;
};

class CoupledDDAssembler {
public:
    CoupledDDAssembler(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       double Vt,
                       double taun,
                       double taup);

    VectorXd pack(const CoupledDDState& state) const;
    CoupledDDState unpack(const VectorXd& x) const;

    VectorXd residual(const VectorXd& x,
                      const CoupledDDBoundaryConditions& bcs) const;

    SparseMatrixd finiteDifferenceJacobian(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        Real relativeStep = 1.0e-6) const;

    VectorXd electronDensity(const VectorXd& x) const;
    VectorXd holeDensity(const VectorXd& x) const;

    bool hasPositiveFiniteCarriers(const VectorXd& x) const;
    Index numNodes() const { return mesh_.numNodes(); }
    const std::vector<Real>& intrinsicDensity() const { return ni_; }

private:
    int psiOffset() const { return 0; }
    int phinOffset() const { return static_cast<int>(mesh_.numNodes()); }
    int phipOffset() const { return 2 * static_cast<int>(mesh_.numNodes()); }

    const DeviceMesh& mesh_;
    const MaterialDatabase& matdb_;
    const DopingModel& doping_;
    double Vt_;
    double taun_;
    double taup_;
    std::vector<Real> ni_;

    // Mesh-derived quantities cached at construction time.
    std::vector<std::vector<Index>> edgeCells_;
    std::vector<Real> vol_;
    std::vector<Real> couple_;
};

} // namespace vela
