#pragma once

#include "vela/core/Types.h"
#include "vela/equation/ChargeSpec.h"
#include "vela/equation/DDAssembler.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/Material.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/physics/RecombinationModel.h"
#include <memory>
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

struct CoupledDDCarrierTermDiagnostic {
    Index nodeId = 0;
    Real electronFlux = 0.0;
    Real holeFlux = 0.0;
    Real electronRecombination = 0.0;
    Real holeRecombination = 0.0;
    Real electronImpact = 0.0;
    Real holeImpact = 0.0;
    Real impactElectronSource = 0.0;
    Real impactHoleSource = 0.0;
    Real impactCombinedSource = 0.0;
    Real electronGauge = 0.0;
    Real holeGauge = 0.0;
    Real electronBoundary = 0.0;
    Real holeBoundary = 0.0;
    Real electronResidual = 0.0;
    Real holeResidual = 0.0;
};

struct CoupledDDEdgeFluxDiagnostic {
    Index edgeId = 0;
    Index node0 = 0;
    Index node1 = 0;
    Real x0 = 0.0;
    Real y0 = 0.0;
    Real x1 = 0.0;
    Real y1 = 0.0;
    Real length_m = 0.0;
    Real couple_m = 0.0;
    Real netDopingAvg_m3 = 0.0;
    Real ni0_m3 = 0.0;
    Real ni1_m3 = 0.0;
    Real psi0_V = 0.0;
    Real psi1_V = 0.0;
    Real phin0_V = 0.0;
    Real phin1_V = 0.0;
    Real phip0_V = 0.0;
    Real phip1_V = 0.0;
    Real electricField_V_m = 0.0;
    Real electronMobility_m2_V_s = 0.0;
    Real holeMobility_m2_V_s = 0.0;
    // Signed Scharfetter-Gummel continuity edge flux (added to node0's residual
    // and subtracted from node1's), identical to the residual edge loop.
    Real electronFlux = 0.0;
    Real holeFlux = 0.0;
};

class CoupledDDAssembler {
public:
    CoupledDDAssembler(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       double Vt,
                       double taun,
                       double taup);

    CoupledDDAssembler(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       double Vt,
                       double taun,
                       double taup,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges);

    CoupledDDAssembler(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       double Vt,
                       const MobilityModelConfig& mobilityConfig,
                       const RecombinationModelConfig& recombinationConfig,
                       const BandgapNarrowingConfig& bandgapNarrowingConfig = {},
                       const ImpactIonizationModelConfig& impactIonizationConfig = {},
                       std::vector<RegionFixedChargeSpec> fixedCharges = {},
                       std::vector<InterfaceSheetChargeSpec> sheetCharges = {},
                       DDScalingSpec scaling = {});

    VectorXd pack(const CoupledDDState& state) const;
    CoupledDDState unpack(const VectorXd& x) const;

    VectorXd residual(const VectorXd& x,
                      const CoupledDDBoundaryConditions& bcs) const;

    SparseMatrixd assembleJacobian(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs) const;

    SparseMatrixd finiteDifferenceJacobian(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs,
        Real relativeStep = 1.0e-6) const;

    VectorXd electronDensity(const VectorXd& x) const;
    VectorXd holeDensity(const VectorXd& x) const;
    std::vector<CoupledDDCarrierTermDiagnostic> carrierContinuityTermDiagnostics(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs) const;

    std::vector<CoupledDDEdgeFluxDiagnostic> sgEdgeFluxDiagnostics(
        const VectorXd& x,
        const CoupledDDBoundaryConditions& bcs) const;

    bool hasPositiveFiniteCarriers(const VectorXd& x) const;
    Index numNodes() const { return mesh_.numNodes(); }
    const std::vector<Real>& intrinsicDensity() const { return ni_; }
    bool usesScaledState() const { return scaling_.enabled; }
    Real potentialScale() const { return scaling_.V0; }
    Real concentrationScale() const { return scaling_.C0; }

private:
    int psiOffset() const { return 0; }
    int phinOffset() const { return static_cast<int>(mesh_.numNodes()); }
    int phipOffset() const { return 2 * static_cast<int>(mesh_.numNodes()); }

    const DeviceMesh& mesh_;
    const MaterialDatabase& matdb_;
    const DopingModel& doping_;
    double Vt_;
    MobilityModelConfig mobilityConfig_;
    std::unique_ptr<MobilityModel> mobility_;
    RecombinationModel recombination_;
    ImpactIonizationModelConfig impactIonizationConfig_;
    std::unique_ptr<ImpactIonizationModel> impactIonization_;
    bool impactIonizationEnabled_ = false;
    bool bgnEnabled_ = false;
    std::vector<Real> ni_;
    std::vector<Material> cellMaterials_;

    // Mesh-derived quantities cached at construction time.
    std::vector<std::vector<Index>> edgeCells_;
    std::vector<std::vector<Index>> nodeCells_;
    std::vector<Real> vol_;
    std::vector<Real> couple_;
    VectorXd fixedInterfaceChargeRhs_;
    DDScalingSpec scaling_;
};

} // namespace vela
