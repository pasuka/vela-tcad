#pragma once

#include "vela/core/Types.h"
#include "vela/core/UnitScaling.h"
#include "vela/mesh/DeviceMesh.h"

#include <unordered_map>
#include <string>
#include <vector>
#include <array>
#include <optional>

namespace vela {

/// Sentaurus-compatible electron density-gradient correction for Boltzmann
/// statistics.  The solved energy divided by q is stored in volts.
struct DensityGradientQuantumPotentialConfig {
    bool enabled = false;
    std::string couplingMode = "outer"; ///< "outer" or frozen imported potential.
    /// How the converged quantum potential enters electron transport.
    /// direct_band_edge preserves Vela's historical behavior and matches
    /// Sentaurus Math{DirectQuantumCorrection}. sentaurus_exponential keeps
    /// the classical Fermi argument and multiplies the carrier density by
    /// exp(-Lambda/kT), matching the Sentaurus default approximation.
    std::string transportCoupling = "direct_band_edge";
    /// Optional continuation parameter for sentaurus_exponential: zero is
    /// direct-band-edge coupling and one is the full exponential form.
    Real transportCouplingWeight = 1.0;
    std::string formulation = "potential_based"; ///< Sentaurus default or density-based audit.
    Real gamma = 3.6;
    /// Electron DOS mass divided by the free-electron mass. O-2018.06 Silicon
    /// eDOSMass Formula1 evaluates to 1.0618016171622988 at 300 K.
    Real effectiveMassRatio = 1.0618016171622988;
    /// Optional mass used only by gamma*hbar^2/(6*m*q).  When absent, the
    /// DOS mass above retains the legacy coupled-mass behavior.
    std::optional<Real> coefficientMassRatio;
    /// Experimental complete global material domain.  Sentaurus solves the
    /// global eQuantumPotential equation in insulators with xi=eta=0, but the
    /// default remains the qualified semiconductor-side interface treatment
    /// until discontinuous potential-like interface unknowns are represented.
    bool includeInsulators = false;
    Real insulatorGamma = 1.0;
    Real insulatorEffectiveMassRatio = 0.42;
    /// Optional insulator counterpart of coefficientMassRatio.
    std::optional<Real> insulatorCoefficientMassRatio;
    std::string interfaceBoundary = "homogeneous_neumann";
    ///< "homogeneous_neumann" or Sentaurus Eq. 233 "sentaurus_step".
    Real theta = 0.5;
    /// Fraction of the positive BGN magnitude assigned to electron affinity.
    Real conductionBandNarrowingFraction = 0.5;
    int maxIterations = 30;
    Real relativeTolerance = 1.0e-7;
    Real absoluteTolerance_V = 1.0e-10;
    /// Optional outer fixed-point Lambda-change tolerance.  When absent,
    /// absoluteTolerance_V retains the legacy shared inner/outer behavior.
    std::optional<Real> outerAbsoluteTolerance_V;
    Real damping = 0.5;
    Real maxUpdate_V = 0.1;
    int outerMaxIterations = 20;
    std::string outerAcceleration = "none"; ///< "none" or vector Aitken relaxation.
    Real outerRelaxation = 1.0;
    Real outerRelaxationMin = 0.1;
    Real outerRelaxationMax = 1.5;
    /// Optional output prefix for a decomposition of the initial global
    /// Eq. 231 residual. Empty keeps the diagnostic disabled.
    std::string residualDiagnosticPrefix;
    /// Evaluate diagnostic drives from the input checkpoint instead of the
    /// preceding frozen-quantum DD solve. Intended for fixed-reference audits.
    bool residualDiagnosticUseInitialState = false;
    /// Global Eq. 231 nonlinear spatial operator. The exponential-fitted
    /// form exactly preserves div(grad(w))+|grad(w)|^2/2 = 2*laplace(u)/u
    /// for u=exp(w/2); p1_direct retains the expanded audit form. The
    /// cvfem_full option integrates the complete expanded operator over
    /// median-dual sub-control volumes, following the Charon/DEVSIM CVFEM
    /// architecture rather than the Galerkin P1 test functions. The
    /// p1_lambda_direct keeps the expanded P1 audit operator but switches to
    /// a globally continuous Lambda unknown. gss_potentiallike_fitted retains
    /// the legacy continuous potential-like state and changes only the flux.
    /// conservative_sqrt_fitted uses the exact theta=1/2 sqrt-density weak
    /// form with a common fixed row scaling, preserving flux/reaction balance.
    /// gss_density_fitted uses the GSS sqrt(n) fitted flux together with the
    /// continuous-Lambda/material-side trace contract. sentaurus_box uses
    /// the verified Formula-0 fitted edge operator, signed AverageBox
    /// circumcentric element-vertex measures, and algebraically eliminated
    /// region-side fitted/reaction traces at a shared material vertex.
    std::string globalDiscretization = "p1_direct";
    /// Dimensionless half-jump offsets for each material-side fitted trace at
    /// a transport/insulator vertex.  Zeros preserve the generic operator.
    Real sentaurusInterfaceInsulatorHalfJumpOffset = 0.0;
    Real sentaurusInterfaceSiliconHalfJumpOffset = 0.0;
    Real sentaurusInterfacePolysiliconHalfJumpOffset = 0.0;
    /// Algebraically eliminated region-side reaction traces at shared
    /// transport/nontransport vertices.  Each value multiplies Lambda and
    /// its consistent diagonal Jacobian on that material side.  Ones retain
    /// the generic single-trace control-volume reaction.  They are explicit
    /// sentaurus_box oracle controls, not global material defaults.
    Real sentaurusInterfaceSiliconReactionWeight = 1.0;
    Real sentaurusInterfacePolysiliconReactionWeight = 1.0;
    Real sentaurusInterfaceInsulatorAtSiliconReactionWeight = 1.0;
    Real sentaurusInterfaceInsulatorAtPolysiliconReactionWeight = 1.0;
    /// Constant terms in the algebraically eliminated affine material-side
    /// reaction traces.  Zeros retain a purely multiplicative trace.
    Real sentaurusInterfaceSiliconReactionOffset_V = 0.0;
    Real sentaurusInterfacePolysiliconReactionOffset_V = 0.0;
    Real sentaurusInterfaceInsulatorAtSiliconReactionOffset_V = 0.0;
    Real sentaurusInterfaceInsulatorAtPolysiliconReactionOffset_V = 0.0;
    /// Reaction control-volume closure for a pure insulator/insulator
    /// re-entrant vertex: exactly six incident triangles split 2:4 between
    /// two nontransport materials.  One preserves the generic box measure.
    Real sentaurusInsulatorReentrantCornerReactionWeight = 1.0;
    /// Experimental DEVSIM/Garcia-Asenov oxide interface closure.  The
    /// oxide-side integrated Eq. 231 row receives the WKB penetration source
    /// b_n,ox/x_n at every transport/insulator interface segment.
    std::string oxideBoundary = "none"; ///< "none" or "devsim_wkb".
    Real oxideQuantumMassRatio = 0.14;
    Real oxideBarrierMassRatio = 0.4;
    Real oxideBarrierHeight_V = 3.15;
};

struct DensityGradientQuantumPotentialResult {
    VectorXd potential_V;
    /// Continuous potential-like unknown Phi/q used by the legacy
    /// all-material formulations. Empty for Lambda-primary formulations,
    /// including gss_density_fitted.
    VectorXd potentialLike_V;
    int iterations = 0;
    Real residualInfinityNorm = 0.0;
    bool converged = false;
    Real lastUpdateInfinityNorm_V = 0.0;
    Real potentialInfinityNorm_V = 0.0;
    Index maxUpdateNode = 0;
    Real maxUpdateNodeValue_V = 0.0;
};

/// One semiconductor-to-unsolved-nonmetal interface segment for the
/// Sentaurus Eq. 233 analytic step boundary.  barrierHeight_V is
/// Ec(unsolved)-Ec(solved), sampled at the two interface vertices.
struct DensityGradientStepBoundary {
    Index edgeId = 0;
    Real barrierHeightN0_V = 0.0;
    Real barrierHeightN1_V = 0.0;
    Real barrierEffectiveMassRatio = 1.0;
    Real solvedGamma = 1.0;
    Real theta = 0.5;
    Real alphaDeterminantCubeRoot = 1.0;
    /// Outward normal derivative of Ec-related driving potential on the
    /// solved side.  Eq. 233 supplies dn(Lambda); the transformed u equation
    /// requires dn(drive-Lambda).
    Real normalDrivingPotentialGradient_V_per_m = 0.0;
};

/// Tri3 material data for the global potential-like formulation.  The shared
/// unknown is Phi/q = Lambda-(psi+affinity+DOS drive), while Lambda and the
/// exponential auxiliary variable remain region-side quantities. In
/// particular, materialBandDrive_V must use the cell-side material at shared
/// interface vertices; it is not the ownership-selected output trace.
struct DensityGradientCellMaterial {
    Index cellId = 0;
    bool isTransport = true;
    Real coefficientVm2 = 0.0;
    std::array<Real, 3> materialBandDrive_V{};
    std::array<Real, 3> dynamicDrivingPotential_V{};
    std::array<Real, 3> initialLambda_V{};
};

/// Stable Sentaurus Eq. 233 step function and its derivative.
Real densityGradientStepFunction(Real x);
Real densityGradientStepFunctionDerivative(Real x);

/// Return gamma*hbar^2/(6*m*q), in V*m^2.
Real densityGradientCoefficientVm2(
    const DensityGradientQuantumPotentialConfig& config);

/// Return gamma*hbar^2/(6*m*q) for an explicit material parameter pair.
Real densityGradientCoefficientVm2(Real gamma, Real effectiveMassRatio);

/// WKB penetration length hbar/sqrt(2*m*q*barrier), in metres.
Real densityGradientOxidePenetrationDepthM(
    Real barrierEffectiveMassRatio,
    Real barrierHeight_V);

/// Solve Lambda*u + C*laplacian(u)=0 with
/// u=sqrt(n/Nref) and n=classicalDensity*exp(-Lambda/Vt).
/// Edges to inactive (insulating) nodes are omitted, giving the homogeneous
/// Neumann boundary used at semiconductor/insulator interfaces.  Values on
/// nodes listed in dirichletPotential_V are imposed exactly.
DensityGradientQuantumPotentialResult solveElectronDensityGradientPotential(
    const DeviceMesh& mesh,
    const VectorXd& classicalDensity,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletPotential_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units = PhysicalUnitSystem::legacySI(),
    DensityGradientQuantumPotentialConfig config = {},
    const VectorXd& initialPotential_V = {});

/// Global Eq. 231 solve. Legacy global discretizations use a continuous
/// potential-like interface unknown. gss_density_fitted uses continuous
/// Lambda and reconstructs the potential-like and sqrt(n) traces per cell.
/// nodeOutputBandDrive_V selects the region-side output convention where a
/// legacy potential-like result must be converted back to Lambda.
DensityGradientQuantumPotentialResult
solveElectronDensityGradientPotentialLikeGlobal(
    const DeviceMesh& mesh,
    const std::vector<DensityGradientCellMaterial>& cellMaterials,
    const VectorXd& nodeOutputBandDrive_V,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletLambda_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units = PhysicalUnitSystem::legacySI(),
    DensityGradientQuantumPotentialConfig config = {},
    const VectorXd& initialLambda_V = {},
    const VectorXd& initialPotentialLike_V = {});

/// Residual and partial derivatives for the fitted sentaurus_box equation.
/// A material-drive column adds the same voltage to every material trace of
/// that node and its output trace. Dynamic drives are independent on the
/// transport and insulator sides. The primary unknown is Phi/q in volts.
struct DensityGradientBoxLinearization {
    VectorXd residual;
    SparseMatrixd potentialLike;
    VectorXd materialDriveDiagonal;
    SparseMatrixd transportDrive;
    SparseMatrixd insulatorDrive;
};

DensityGradientBoxLinearization linearizeElectronDensityGradientBox(
    const DeviceMesh& mesh,
    const std::vector<DensityGradientCellMaterial>& cellMaterials,
    const VectorXd& nodeOutputBandDrive_V,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletLambda_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    const DensityGradientQuantumPotentialConfig& config,
    const VectorXd& potentialLike_V,
    bool assembleJacobian = true);

/// Material-resolved Eq. 231 variant.  coefficientVm2 stores
/// gamma*hbar^2/(6*m*q) per node.  Eq. 231 keeps gamma/m outside the
/// divergence, so it scales each material control-volume row rather than the
/// continuous normal alpha-weighted interface flux.
DensityGradientQuantumPotentialResult solveElectronDensityGradientPotential(
    const DeviceMesh& mesh,
    const VectorXd& classicalDensity,
    const VectorXd& coefficientVm2,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletPotential_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units = PhysicalUnitSystem::legacySI(),
    DensityGradientQuantumPotentialConfig config = {},
    const VectorXd& initialPotential_V = {});

DensityGradientQuantumPotentialResult solveElectronDensityGradientPotential(
    const DeviceMesh& mesh,
    const VectorXd& classicalDensity,
    const VectorXd& coefficientVm2,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletPotential_V,
    const std::vector<DensityGradientStepBoundary>& stepBoundaries,
    Real thermalVoltage_V,
    PhysicalUnitSystem units = PhysicalUnitSystem::legacySI(),
    DensityGradientQuantumPotentialConfig config = {},
    const VectorXd& initialPotential_V = {});

} // namespace vela
