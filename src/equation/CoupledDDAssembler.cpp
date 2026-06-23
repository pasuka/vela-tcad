#include "vela/equation/CoupledDDAssembler.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/Bernoulli.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/physics/CarrierStatistics.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace vela {
namespace {

Real bernoulliDerivative(Real x)
{
    const Real ax = std::abs(x);
    if (ax < 1.0e-8)
        return -0.5 + x / 6.0 - x * x * x / 180.0;
    if (x > 500.0)
        return (1.0 - x) * std::exp(-x);
    if (x < -500.0)
        return -1.0;

    const Real em1 = std::expm1(x);
    const Real ex = em1 + 1.0;
    return (em1 - x * ex) / (em1 * em1);
}

Real limitedExp(Real value)
{
    return std::exp(std::clamp(value, -500.0, 500.0));
}

bool transportMobilityDependsOnPotentials(const MobilityModelConfig& config)
{
    if (!config.jacobianFieldDerivatives)
        return false;
    return config.model == "caughey_thomas_field" ||
           config.model == "masetti_field" ||
           config.model == "caughey_thomas_field_surface" ||
           isSurfaceMobilityModel(config);
}

} // namespace

CoupledDDAssembler::CoupledDDAssembler(const DeviceMesh& mesh,
                                       const MaterialDatabase& matdb,
                                       const DopingModel& doping,
                                       double Vt,
                                       double taun,
                                       double taup)
    : CoupledDDAssembler(mesh, matdb, doping, Vt, taun, taup, {}, {})
{}

CoupledDDAssembler::CoupledDDAssembler(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    double Vt,
    double taun,
    double taup,
    std::vector<RegionFixedChargeSpec> fixedCharges,
    std::vector<InterfaceSheetChargeSpec> sheetCharges)
    : CoupledDDAssembler(mesh,
                         matdb,
                         doping,
                         Vt,
                         MobilityModelConfig{},
                         recombinationModelConfig({"srh"}, taun, taup),
                         BandgapNarrowingConfig{},
                         ImpactIonizationModelConfig{},
                         std::move(fixedCharges),
                         std::move(sheetCharges))
{}

CoupledDDAssembler::CoupledDDAssembler(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    double Vt,
    const MobilityModelConfig& mobilityConfig,
    const RecombinationModelConfig& recombinationConfig,
    const BandgapNarrowingConfig& bandgapNarrowingConfig,
    const ImpactIonizationModelConfig& impactIonizationConfig,
    std::vector<RegionFixedChargeSpec> fixedCharges,
    std::vector<InterfaceSheetChargeSpec> sheetCharges,
    DDScalingSpec scaling)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , Vt_(Vt)
    , mobilityConfig_(mobilityConfig)
    , mobility_(makeMobilityModel(mobilityConfig))
    , recombination_(recombinationConfig)
    , impactIonizationConfig_(impactIonizationConfig)
    , impactIonization_(makeImpactIonizationModel(impactIonizationConfig))
    , impactIonizationEnabled_(impactIonizationConfig.model != "none")
    , bgnEnabled_(bandgapNarrowingConfig.model != "none")
    , ni_(detail::buildValidatedEffectiveNodeNi(
          "CoupledDDAssembler",
          mesh,
          matdb,
          doping,
          bandgapNarrowingConfig,
          Vt))
    , cellMaterials_(detail::buildCellMaterials(
          mesh,
          matdb,
          Vt * constants::q / constants::kb))
    , edgeCells_(detail::buildEdgeCellMap(mesh))
    , nodeCells_(detail::buildNodeCellMap(mesh))
    , vol_(detail::computeNodeVolumes(mesh))
    , couple_(detail::computeEdgeCouplings(mesh))
    , fixedInterfaceChargeRhs_(detail::computeFixedAndInterfaceChargeRhs(
          mesh, edgeCells_, fixedCharges, sheetCharges, "CoupledDDAssembler"))
    , scaling_(scaling)
{
    if (scaling_.enabled) {
        const auto isPositiveFinite = [](Real value) {
            return value > 0.0 && std::isfinite(value);
        };
        if (!isPositiveFinite(scaling_.V0) ||
            !isPositiveFinite(scaling_.C0) ||
            !isPositiveFinite(scaling_.mu0) ||
            !isPositiveFinite(scaling_.D0) ||
            !isPositiveFinite(scaling_.L0) ||
            !isPositiveFinite(scaling_.permittivityReference_F_per_m)) {
            throw std::invalid_argument(
                "CoupledDDAssembler: scaling references must be positive and finite when scaling is enabled.");
        }
    }
}

VectorXd CoupledDDAssembler::pack(const CoupledDDState& state) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (state.psi.size() != N || state.phin.size() != N || state.phip.size() != N)
        throw std::invalid_argument("CoupledDDAssembler::pack: state vector size mismatch.");

    VectorXd x(3 * N);
    x.segment(psiOffset(), N) = state.psi;
    x.segment(phinOffset(), N) = state.phin;
    x.segment(phipOffset(), N) = state.phip;
    return x;
}

CoupledDDState CoupledDDAssembler::unpack(const VectorXd& x) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (x.size() != 3 * N)
        throw std::invalid_argument("CoupledDDAssembler::unpack: vector size mismatch.");

    CoupledDDState state;
    state.psi = x.segment(psiOffset(), N);
    state.phin = x.segment(phinOffset(), N);
    state.phip = x.segment(phipOffset(), N);
    return state;
}

VectorXd CoupledDDAssembler::electronDensity(const VectorXd& x) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    VectorXd n(N);
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    for (int i = 0; i < N; ++i)
        n(i) = vela::electronDensity(ni_[static_cast<Index>(i)],
                                     x(psiOffset() + i) * potentialScale,
                                     x(phinOffset() + i) * potentialScale,
                                     Vt_);
    return n;
}

VectorXd CoupledDDAssembler::holeDensity(const VectorXd& x) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    VectorXd p(N);
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    for (int i = 0; i < N; ++i)
        p(i) = vela::holeDensity(ni_[static_cast<Index>(i)],
                                 x(psiOffset() + i) * potentialScale,
                                 x(phipOffset() + i) * potentialScale,
                                 Vt_);
    return p;
}

bool CoupledDDAssembler::hasPositiveFiniteCarriers(const VectorXd& x) const
{
    const VectorXd n = electronDensity(x);
    const VectorXd p = holeDensity(x);
    for (int i = 0; i < n.size(); ++i) {
        if (!std::isfinite(n(i)) || !std::isfinite(p(i)) || n(i) <= 0.0 || p(i) <= 0.0)
            return false;
    }
    return true;
}

VectorXd CoupledDDAssembler::residual(const VectorXd& x,
                                      const CoupledDDBoundaryConditions& bcs) const
{
    const Index Nidx = mesh_.numNodes();
    const int N = static_cast<int>(Nidx);
    if (x.size() != 3 * N)
        throw std::invalid_argument("CoupledDDAssembler::residual: vector size mismatch.");

    // n and p are needed for Poisson source and configured recombination.
    const VectorXd n = electronDensity(x);
    const VectorXd p = holeDensity(x);
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    const VectorXd psi = x.segment(psiOffset(), N) * potentialScale;

    const std::vector<Real> nodeElectricFields = impactIonizationEnabled_
        ? detail::computeNodeElectricFields(psi, mesh_)
        : std::vector<Real>{};
    const bool qfImpact = impactIonizationConfig_.drivingForce == "quasi_fermi_gradient";
    const VectorXd phinPhysical = x.segment(phinOffset(), N) * potentialScale;
    const VectorXd phipPhysical = x.segment(phipOffset(), N) * potentialScale;
    const std::vector<Real> nodeElectronDrivingFields = (impactIonizationEnabled_ && qfImpact)
        ? detail::computeNodeScalarGradientMagnitudes(
            phinPhysical, mesh_)
        : nodeElectricFields;
    const std::vector<Real> nodeHoleDrivingFields = (impactIonizationEnabled_ && qfImpact)
        ? detail::computeNodeScalarGradientMagnitudes(
            phipPhysical, mesh_)
        : nodeElectricFields;
    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";
    const bool sgCurrentAvalanche = impactIonizationEnabled_ &&
        detail::usesEdgeCurrentAvalancheSource(impactIonizationConfig_);
    const std::vector<Real> sgAvalancheSourceIntegrals = sgCurrentAvalanche
        ? detail::sgEdgeCurrentAvalancheSourceIntegrals(
            impactIonizationConfig_,
            *impactIonization_,
            mobilityConfig_,
            *mobility_,
            edgeCells_,
            mesh_,
            doping_,
            cellMaterials_,
            psi,
            phinPhysical,
            phipPhysical,
            n,
            p,
            ni_,
            Vt_)
        : std::vector<Real>{};

    VectorXd r = VectorXd::Zero(3 * N);
    std::vector<bool> hasElectronContribution(static_cast<std::size_t>(N), false);
    std::vector<bool> hasHoleContribution(static_cast<std::size_t>(N), false);

    const auto& edgeCells = edgeCells_;
    const auto& vol = vol_;
    const auto& couple = couple_;

    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30) continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real psi_i = x(psiOffset() + i) * potentialScale;
        const Real psi_j = x(psiOffset() + j) * potentialScale;
        const Real phin_i = x(phinOffset() + i) * potentialScale;
        const Real phin_j = x(phinOffset() + j) * potentialScale;
        const Real phip_i = x(phipOffset() + i) * potentialScale;
        const Real phip_j = x(phipOffset() + j) * potentialScale;
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h);
        const Real electronMobilityField =
            qfMobility ? std::abs((phin_j - phin_i) / h) : electricField;
        const Real holeMobilityField =
            qfMobility ? std::abs((phip_j - phip_i) / h) : electricField;

        const Real eps = detail::edgeEpsilon(edgeCells, mesh_, matdb_, e);
        const Real G = eps * couple[e] / h;
        const Real psiFlux = G * (psi_i - psi_j);
        r(psiOffset() + i) += psiFlux;
        r(psiOffset() + j) -= psiFlux;

        const Real mun = detail::edgeMobility(
            edgeCells, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Electron,
            electronMobilityField,
            &mobilityConfig_,
            &psi);
        if (mun > 0.0) {
            hasElectronContribution[static_cast<std::size_t>(i)] = true;
            hasElectronContribution[static_cast<std::size_t>(j)] = true;

            const Real coef = mun * Vt_ * couple[e] / h;
            const Index idxI = static_cast<Index>(i);
            const Index idxJ = static_cast<Index>(j);
            // With bandgap narrowing the per-node intrinsic density varies, so
            // the ni-gradient drift term must be retained (VariableNi form).
            // Without BGN, reproduce the baseline discretization: the balanced
            // quasi-Fermi form on equal-ni edges (overflow-safe, it never
            // materializes the clamped carrier densities) and the density-based
            // SG flux on material-interface edges where ni differs.
            Real nFlux;
            if (bgnEnabled_) {
                nFlux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                    ni_[idxI], ni_[idxJ], psi_i, psi_j, phin_i, phin_j, Vt_, coef);
            } else if (ni_[idxI] == ni_[idxJ]) {
                nFlux = sgElectronContinuityFluxFromQuasiFermiStable(
                    ni_[idxI], psi_i, psi_j, phin_i, phin_j, Vt_, coef);
            } else {
                nFlux = sgElectronContinuityFlux(n(i), n(j), dpsi, Vt_, coef);
            }
            r(phinOffset() + i) += nFlux;
            r(phinOffset() + j) -= nFlux;
        }

        const Real mup = detail::edgeMobility(
            edgeCells, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Hole,
            holeMobilityField,
            &mobilityConfig_,
            &psi);
        if (mup > 0.0) {
            hasHoleContribution[static_cast<std::size_t>(i)] = true;
            hasHoleContribution[static_cast<std::size_t>(j)] = true;

            const Real coef = mup * Vt_ * couple[e] / h;
            const Index idxI = static_cast<Index>(i);
            const Index idxJ = static_cast<Index>(j);
            // See the electron flux above for the BGN gating rationale.
            Real pFlux;
            if (bgnEnabled_) {
                pFlux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                    ni_[idxI], ni_[idxJ], psi_i, psi_j, phip_i, phip_j, Vt_, coef);
            } else if (ni_[idxI] == ni_[idxJ]) {
                pFlux = sgHoleContinuityFluxFromQuasiFermiStable(
                    ni_[idxI], psi_i, psi_j, phip_i, phip_j, Vt_, coef);
            } else {
                pFlux = sgHoleContinuityFlux(p(i), p(j), dpsi, Vt_, coef);
            }
            r(phipOffset() + i) += pFlux;
            r(phipOffset() + j) -= pFlux;
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        r(psiOffset() + ii) -= constants::q *
            (p(ii) - n(ii) + doping_.netDoping(i)) * vol[i];

        const Real ni = ni_[i];
        if (ni <= 0.0)
            continue;

        if (recombination_.srhEnabled() || recombination_.augerEnabled()) {
            // Compute n*p - ni^2 via the identity n*p = ni^2 * exp((phip-phin)/Vt),
            // i.e. n*p - ni^2 = ni^2 * expm1((phip-phin)/Vt). This avoids
            // catastrophic cancellation when phip ~= phin (near equilibrium), where
            // the naive form n*p - ni^2 is dominated by floating-point rounding
            // in exp(+u) * exp(-u) != 1.
            const Real dPhi =
                (x(phipOffset() + ii) - x(phinOffset() + ii)) * potentialScale;
            const Real excessProduct = (ni > 0.0)
                ? ni * ni * std::expm1(dPhi / Vt_)
                : n(ii) * p(ii);
            const Real R = recombination_.totalRateFromExcessProduct(
                excessProduct, n(ii), p(ii), ni);
            if (R != 0.0) {
                r(phinOffset() + ii) += R * vol[i];
                r(phipOffset() + ii) += R * vol[i];
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }

        if (impactIonizationEnabled_ && sgCurrentAvalanche) {
            const Real source = sgAvalancheSourceIntegrals[i];
            if (source != 0.0) {
                r(phinOffset() + ii) -= source;
                r(phipOffset() + ii) -= source;
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        } else if (impactIonizationEnabled_) {
            const Real G = detail::impactIonizationGenerationRate(
                impactIonizationConfig_,
                *impactIonization_,
                mobilityConfig_,
                *mobility_,
                nodeCells_,
                mesh_,
                doping_,
                cellMaterials_,
                i,
                nodeElectricFields[i],
                nodeElectronDrivingFields[i],
                nodeHoleDrivingFields[i],
                n(ii),
                p(ii));
            if (G != 0.0) {
                r(phinOffset() + ii) -= G * vol[i];
                r(phipOffset() + ii) -= G * vol[i];
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }
    }

    // Insulating nodes such as SiO2 (mun = mup = ni = 0) can have no
    // transport or recombination contribution in the continuity equations.
    // Pin those otherwise unconstrained quasi-Fermi unknowns to avoid zero
    // residual/Jacobian rows. Explicit boundary conditions below take
    // precedence over this internal gauge constraint.
    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        if (!hasElectronContribution[static_cast<std::size_t>(ii)])
            r(phinOffset() + ii) = x(phinOffset() + ii);
        if (!hasHoleContribution[static_cast<std::size_t>(ii)])
            r(phipOffset() + ii) = x(phipOffset() + ii);
    }

    for (int i = 0; i < N; ++i)
        r(psiOffset() + i) -= fixedInterfaceChargeRhs_(i);

    if (scaling_.enabled) {
        const Real poissonScale =
            scaling_.permittivityReference_F_per_m * scaling_.V0;
        const Real continuityScale = scaling_.C0 * scaling_.D0;
        for (int i = 0; i < N; ++i) {
            r(psiOffset() + i) /= poissonScale;
            r(phinOffset() + i) /= continuityScale;
            r(phipOffset() + i) /= continuityScale;
        }
    }

    if (scaling_.enabled) {
        for (Index i = 0; i < Nidx; ++i) {
            const int ii = static_cast<int>(i);
            if (!hasElectronContribution[static_cast<std::size_t>(ii)])
                r(phinOffset() + ii) = x(phinOffset() + ii);
            if (!hasHoleContribution[static_cast<std::size_t>(ii)])
                r(phipOffset() + ii) = x(phipOffset() + ii);
        }
    }

    // Boundary-condition maps are independent so multi-terminal MOS callers can
    // pin electrostatic potential, electron quasi-Fermi potential, and hole
    // quasi-Fermi potential on the source/drain/body nodes selected by contact
    // name without relying on contact order.
    for (const auto& [node, value] : bcs.psi)
        r(psiOffset() + static_cast<int>(node)) = x(psiOffset() + static_cast<int>(node)) - value;
    for (const auto& [node, value] : bcs.phin)
        r(phinOffset() + static_cast<int>(node)) = x(phinOffset() + static_cast<int>(node)) - value;
    for (const auto& [node, value] : bcs.phip)
        r(phipOffset() + static_cast<int>(node)) = x(phipOffset() + static_cast<int>(node)) - value;

    return r;
}

std::vector<CoupledDDCarrierTermDiagnostic>
CoupledDDAssembler::carrierContinuityTermDiagnostics(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs) const
{
    const Index Nidx = mesh_.numNodes();
    const int N = static_cast<int>(Nidx);
    if (x.size() != 3 * N)
        throw std::invalid_argument(
            "CoupledDDAssembler::carrierContinuityTermDiagnostics: vector size mismatch.");

    std::vector<CoupledDDCarrierTermDiagnostic> terms(Nidx);
    for (Index i = 0; i < Nidx; ++i)
        terms[i].nodeId = i;

    const VectorXd n = electronDensity(x);
    const VectorXd p = holeDensity(x);
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    const VectorXd psi = x.segment(psiOffset(), N) * potentialScale;

    const std::vector<Real> nodeElectricFields = impactIonizationEnabled_
        ? detail::computeNodeElectricFields(psi, mesh_)
        : std::vector<Real>{};
    const bool qfImpact = impactIonizationConfig_.drivingForce == "quasi_fermi_gradient";
    const VectorXd phinPhysical = x.segment(phinOffset(), N) * potentialScale;
    const VectorXd phipPhysical = x.segment(phipOffset(), N) * potentialScale;
    const std::vector<Real> nodeElectronDrivingFields = (impactIonizationEnabled_ && qfImpact)
        ? detail::computeNodeScalarGradientMagnitudes(phinPhysical, mesh_)
        : nodeElectricFields;
    const std::vector<Real> nodeHoleDrivingFields = (impactIonizationEnabled_ && qfImpact)
        ? detail::computeNodeScalarGradientMagnitudes(phipPhysical, mesh_)
        : nodeElectricFields;
    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";
    const bool sgCurrentAvalanche = impactIonizationEnabled_ &&
        detail::usesEdgeCurrentAvalancheSource(impactIonizationConfig_);
    const detail::SgAvalancheSourceComponentIntegrals sgAvalancheSourceComponents =
        sgCurrentAvalanche
        ? detail::sgEdgeCurrentAvalancheSourceComponentIntegrals(
            impactIonizationConfig_,
            *impactIonization_,
            mobilityConfig_,
            *mobility_,
            edgeCells_,
            mesh_,
            doping_,
            cellMaterials_,
            psi,
            phinPhysical,
            phipPhysical,
            n,
            p,
            ni_,
            Vt_)
        : detail::SgAvalancheSourceComponentIntegrals{};

    std::vector<bool> hasElectronContribution(static_cast<std::size_t>(N), false);
    std::vector<bool> hasHoleContribution(static_cast<std::size_t>(N), false);

    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30)
            continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real psi_i = x(psiOffset() + i) * potentialScale;
        const Real psi_j = x(psiOffset() + j) * potentialScale;
        const Real phin_i = x(phinOffset() + i) * potentialScale;
        const Real phin_j = x(phinOffset() + j) * potentialScale;
        const Real phip_i = x(phipOffset() + i) * potentialScale;
        const Real phip_j = x(phipOffset() + j) * potentialScale;
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h);
        const Real electronMobilityField =
            qfMobility ? std::abs((phin_j - phin_i) / h) : electricField;
        const Real holeMobilityField =
            qfMobility ? std::abs((phip_j - phip_i) / h) : electricField;

        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Electron,
            electronMobilityField,
            &mobilityConfig_,
            &psi);
        if (mun > 0.0) {
            hasElectronContribution[static_cast<std::size_t>(i)] = true;
            hasElectronContribution[static_cast<std::size_t>(j)] = true;
            const Real coef = mun * Vt_ * couple_[e] / h;
            const Real nFlux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni_[static_cast<Index>(i)],
                ni_[static_cast<Index>(j)],
                psi_i,
                psi_j,
                phin_i,
                phin_j,
                Vt_,
                coef,
                bgnEnabled_);
            terms[static_cast<Index>(i)].electronFlux += nFlux;
            terms[static_cast<Index>(j)].electronFlux -= nFlux;
        }

        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Hole,
            holeMobilityField,
            &mobilityConfig_,
            &psi);
        if (mup > 0.0) {
            hasHoleContribution[static_cast<std::size_t>(i)] = true;
            hasHoleContribution[static_cast<std::size_t>(j)] = true;
            const Real coef = mup * Vt_ * couple_[e] / h;
            const Real pFlux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                ni_[static_cast<Index>(i)],
                ni_[static_cast<Index>(j)],
                psi_i,
                psi_j,
                phip_i,
                phip_j,
                Vt_,
                coef,
                bgnEnabled_);
            terms[static_cast<Index>(i)].holeFlux += pFlux;
            terms[static_cast<Index>(j)].holeFlux -= pFlux;
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        const Real ni = ni_[i];
        if (ni <= 0.0)
            continue;

        if (recombination_.srhEnabled() || recombination_.augerEnabled()) {
            const Real dPhi =
                (x(phipOffset() + ii) - x(phinOffset() + ii)) * potentialScale;
            const Real excessProduct = (ni > 0.0)
                ? ni * ni * std::expm1(dPhi / Vt_)
                : n(ii) * p(ii);
            const Real R = recombination_.totalRateFromExcessProduct(
                excessProduct, n(ii), p(ii), ni);
            if (R != 0.0) {
                const Real contribution = R * vol_[i];
                terms[i].electronRecombination += contribution;
                terms[i].holeRecombination += contribution;
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }

        if (impactIonizationEnabled_ && sgCurrentAvalanche) {
            const Real source = sgAvalancheSourceComponents.combined[i];
            const Real contribution = -source;
            if (contribution != 0.0) {
                terms[i].electronImpact += contribution;
                terms[i].holeImpact += contribution;
                terms[i].impactElectronSource += sgAvalancheSourceComponents.electron[i];
                terms[i].impactHoleSource += sgAvalancheSourceComponents.hole[i];
                terms[i].impactCombinedSource += source;
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        } else if (impactIonizationEnabled_) {
            const Real G = detail::impactIonizationGenerationRate(
                impactIonizationConfig_,
                *impactIonization_,
                mobilityConfig_,
                *mobility_,
                nodeCells_,
                mesh_,
                doping_,
                cellMaterials_,
                i,
                nodeElectricFields[i],
                nodeElectronDrivingFields[i],
                nodeHoleDrivingFields[i],
                n(ii),
                p(ii));
            if (G != 0.0) {
                const Real contribution = -G * vol_[i];
                terms[i].electronImpact += contribution;
                terms[i].holeImpact += contribution;
                terms[i].impactCombinedSource += G * vol_[i];
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }
    }

    if (scaling_.enabled) {
        const Real continuityScale = scaling_.C0 * scaling_.D0;
        for (auto& term : terms) {
            term.electronFlux /= continuityScale;
            term.holeFlux /= continuityScale;
            term.electronRecombination /= continuityScale;
            term.holeRecombination /= continuityScale;
            term.electronImpact /= continuityScale;
            term.holeImpact /= continuityScale;
            term.impactElectronSource /= continuityScale;
            term.impactHoleSource /= continuityScale;
            term.impactCombinedSource /= continuityScale;
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        if (!hasElectronContribution[static_cast<std::size_t>(ii)])
            terms[i].electronGauge = x(phinOffset() + ii);
        if (!hasHoleContribution[static_cast<std::size_t>(ii)])
            terms[i].holeGauge = x(phipOffset() + ii);
    }

    for (auto& term : terms) {
        term.electronResidual = term.electronFlux
            + term.electronRecombination
            + term.electronImpact
            + term.electronGauge;
        term.holeResidual = term.holeFlux
            + term.holeRecombination
            + term.holeImpact
            + term.holeGauge;
    }

    auto applyElectronBoundary = [&](Index node, Real value) {
        CoupledDDCarrierTermDiagnostic& term = terms[node];
        term.electronFlux = 0.0;
        term.electronRecombination = 0.0;
        term.electronImpact = 0.0;
        term.impactElectronSource = 0.0;
        term.impactHoleSource = 0.0;
        term.impactCombinedSource = 0.0;
        term.electronGauge = 0.0;
        term.electronBoundary = x(phinOffset() + static_cast<int>(node)) - value;
        term.electronResidual = term.electronBoundary;
    };
    auto applyHoleBoundary = [&](Index node, Real value) {
        CoupledDDCarrierTermDiagnostic& term = terms[node];
        term.holeFlux = 0.0;
        term.holeRecombination = 0.0;
        term.holeImpact = 0.0;
        term.impactElectronSource = 0.0;
        term.impactHoleSource = 0.0;
        term.impactCombinedSource = 0.0;
        term.holeGauge = 0.0;
        term.holeBoundary = x(phipOffset() + static_cast<int>(node)) - value;
        term.holeResidual = term.holeBoundary;
    };
    for (const auto& [node, value] : bcs.phin)
        applyElectronBoundary(node, value);
    for (const auto& [node, value] : bcs.phip)
        applyHoleBoundary(node, value);

    return terms;
}

std::vector<CoupledDDEdgeFluxDiagnostic>
CoupledDDAssembler::sgEdgeFluxDiagnostics(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs) const
{
    (void)bcs;
    const Index Nidx = mesh_.numNodes();
    const int N = static_cast<int>(Nidx);
    if (x.size() != 3 * N)
        throw std::invalid_argument(
            "CoupledDDAssembler::sgEdgeFluxDiagnostics: vector size mismatch.");

    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    const VectorXd psi = x.segment(psiOffset(), N) * potentialScale;

    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";
    const Real continuityScale =
        scaling_.enabled ? (scaling_.C0 * scaling_.D0) : 1.0;

    std::vector<CoupledDDEdgeFluxDiagnostic> edges;
    edges.reserve(static_cast<std::size_t>(mesh_.numEdges()));

    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30)
            continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Index idxI = static_cast<Index>(i);
        const Index idxJ = static_cast<Index>(j);
        const Real psi_i = x(psiOffset() + i) * potentialScale;
        const Real psi_j = x(psiOffset() + j) * potentialScale;
        const Real phin_i = x(phinOffset() + i) * potentialScale;
        const Real phin_j = x(phinOffset() + j) * potentialScale;
        const Real phip_i = x(phipOffset() + i) * potentialScale;
        const Real phip_j = x(phipOffset() + j) * potentialScale;
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h);
        const Real electronMobilityField =
            qfMobility ? std::abs((phin_j - phin_i) / h) : electricField;
        const Real holeMobilityField =
            qfMobility ? std::abs((phip_j - phip_i) / h) : electricField;

        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Electron, electronMobilityField, &mobilityConfig_, &psi);
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Hole, holeMobilityField, &mobilityConfig_, &psi);

        Real nFlux = 0.0;
        if (mun > 0.0) {
            const Real coef = mun * Vt_ * couple_[e] / h;
            nFlux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni_[idxI], ni_[idxJ], psi_i, psi_j, phin_i, phin_j, Vt_, coef,
                bgnEnabled_);
        }
        Real pFlux = 0.0;
        if (mup > 0.0) {
            const Real coef = mup * Vt_ * couple_[e] / h;
            pFlux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                ni_[idxI], ni_[idxJ], psi_i, psi_j, phip_i, phip_j, Vt_, coef,
                bgnEnabled_);
        }

        const Node& node0 = mesh_.getNode(edge.n0);
        const Node& node1 = mesh_.getNode(edge.n1);
        CoupledDDEdgeFluxDiagnostic record;
        record.edgeId = e;
        record.node0 = idxI;
        record.node1 = idxJ;
        record.x0 = node0.x;
        record.y0 = node0.y;
        record.x1 = node1.x;
        record.y1 = node1.y;
        record.length_m = h;
        record.couple_m = couple_[e];
        record.netDopingAvg_m3 =
            0.5 * (doping_.netDoping(idxI) + doping_.netDoping(idxJ));
        record.ni0_m3 = ni_[idxI];
        record.ni1_m3 = ni_[idxJ];
        record.psi0_V = psi_i;
        record.psi1_V = psi_j;
        record.phin0_V = phin_i;
        record.phin1_V = phin_j;
        record.phip0_V = phip_i;
        record.phip1_V = phip_j;
        record.electricField_V_m = electricField;
        record.electronMobility_m2_V_s = mun;
        record.holeMobility_m2_V_s = mup;
        record.electronFlux = nFlux / continuityScale;
        record.holeFlux = pFlux / continuityScale;
        edges.push_back(record);
    }

    return edges;
}

SparseMatrixd CoupledDDAssembler::assembleJacobian(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs) const
{
    const Index Nidx = mesh_.numNodes();
    const int N = static_cast<int>(Nidx);
    if (x.size() != 3 * N)
        throw std::invalid_argument("CoupledDDAssembler::assembleJacobian: vector size mismatch.");

    const VectorXd n = electronDensity(x);
    const VectorXd p = holeDensity(x);
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    const VectorXd psi = x.segment(psiOffset(), N) * potentialScale;
    const VectorXd phinState = x.segment(phinOffset(), N) * potentialScale;
    const VectorXd phipState = x.segment(phipOffset(), N) * potentialScale;
    const std::vector<Real> nodeElectricFields = impactIonizationEnabled_
        ? detail::computeNodeElectricFields(psi, mesh_)
        : std::vector<Real>{};
    const bool qfImpact = impactIonizationConfig_.drivingForce == "quasi_fermi_gradient";
    const std::vector<Real> nodeElectronDrivingFields = (impactIonizationEnabled_ && qfImpact)
        ? detail::computeNodeScalarGradientMagnitudes(
            x.segment(phinOffset(), N) * potentialScale, mesh_)
        : nodeElectricFields;
    const std::vector<Real> nodeHoleDrivingFields = (impactIonizationEnabled_ && qfImpact)
        ? detail::computeNodeScalarGradientMagnitudes(
            x.segment(phipOffset(), N) * potentialScale, mesh_)
        : nodeElectricFields;
    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";
    const std::vector<bool> contactNodes = detail::contactNodeMask(mesh_);
    const bool transportMobilityDerivative = transportMobilityDependsOnPotentials(mobilityConfig_);
    const bool sgCurrentAvalanche = impactIonizationEnabled_ &&
        detail::usesEdgeCurrentAvalancheSource(impactIonizationConfig_);

    std::vector<bool> constrainedRows(static_cast<std::size_t>(3 * N), false);
    for (const auto& [node, value] : bcs.psi) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(psiOffset() + static_cast<int>(node))] = true;
    }
    for (const auto& [node, value] : bcs.phin) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(phinOffset() + static_cast<int>(node))] = true;
    }
    for (const auto& [node, value] : bcs.phip) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(phipOffset() + static_cast<int>(node))] = true;
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(N) * 27);
    auto scaleDerivative = [&](int row, int, Real value) {
        if (!scaling_.enabled)
            return value;
        const Real rowScale = (row < N)
            ? (scaling_.permittivityReference_F_per_m * scaling_.V0)
            : (scaling_.C0 * scaling_.D0);
        return value * scaling_.V0 / rowScale;
    };

    auto add = [&](int row, int col, Real value) {
        if (value != 0.0 && !constrainedRows[static_cast<std::size_t>(row)])
            triplets.emplace_back(row, col, scaleDerivative(row, col, value));
    };

    auto addGauge = [&](int row, int col, Real value) {
        if (value != 0.0 && !constrainedRows[static_cast<std::size_t>(row)])
            triplets.emplace_back(row, col, value);
    };

    std::vector<bool> hasElectronContribution(static_cast<std::size_t>(N), false);
    std::vector<bool> hasHoleContribution(static_cast<std::size_t>(N), false);

    struct EdgeAvalancheNodeSources {
        Real node0 = 0.0;
        Real node1 = 0.0;
    };

    // Scharfetter-Gummel edge-current avalanche source for one edge, evaluated
    // directly from the endpoint potentials. This mirrors
    // detail::sgEdgeCurrentAvalancheSourceRecords and returns the directionally
    // partitioned nodal source used by the continuity residuals.
    auto edgeAvalancheNodeSources =
        [&](Index e, int i, int j, Real h,
            Real psi_i, Real psi_j,
            auto&& phinAt,
            auto&& phipAt) -> EdgeAvalancheNodeSources {
        EdgeAvalancheNodeSources sources;
        const Real couple_e = couple_[e];
        if (h <= 1.0e-30 || couple_e <= 0.0)
            return sources;
        const Index idxI = static_cast<Index>(i);
        const Index idxJ = static_cast<Index>(j);
        const Real niI = ni_[idxI];
        const Real niJ = ni_[idxJ];
        const Real phin_i = phinAt(idxI);
        const Real phin_j = phinAt(idxJ);
        const Real phip_i = phipAt(idxI);
        const Real phip_j = phipAt(idxJ);
        const Real n_i = niI * limitedExp((psi_i - phin_i) / Vt_);
        const Real n_j = niJ * limitedExp((psi_j - phin_j) / Vt_);
        const Real p_i = niI * limitedExp((phip_i - psi_i) / Vt_);
        const Real p_j = niJ * limitedExp((phip_j - psi_j) / Vt_);
        const Real electronQf_i = detail::electronQfForAvalancheGradient(
            psi_i, phin_i, n_i, niI, Vt_, impactIonizationConfig_);
        const Real electronQf_j = detail::electronQfForAvalancheGradient(
            psi_j, phin_j, n_j, niJ, Vt_, impactIonizationConfig_);
        const Real holeQf_i = detail::holeQfForAvalancheGradient(
            psi_i, phip_i, p_i, niI, Vt_, impactIonizationConfig_);
        const Real holeQf_j = detail::holeQfForAvalancheGradient(
            psi_j, phip_j, p_j, niJ, Vt_, impactIonizationConfig_);
        const Real electricField = std::abs((psi_j - psi_i) / h);
        const Real electronQfField = std::abs((electronQf_j - electronQf_i) / h);
        const Real holeQfField = std::abs((holeQf_j - holeQf_i) / h);
        const Real electronCoefficientField = detail::edgeHighFieldDrivingField(
            qfImpact, electronQfField, electricField, edgeCells_, mesh_, e, contactNodes);
        const Real holeCoefficientField = detail::edgeHighFieldDrivingField(
            qfImpact, holeQfField, electricField, edgeCells_, mesh_, e, contactNodes);
        const Real electronMobilityField = qfMobility ? electronQfField : electricField;
        const Real holeMobilityField = qfMobility ? holeQfField : electricField;
        const Real nAvg = 0.5 * (n_i + n_j);
        const Real pAvg = 0.5 * (p_i + p_j);
        const Real electronImpactField = detail::electronAvalancheDrivingField(
            impactIonizationConfig_, electronCoefficientField, electricField, nAvg);
        const Real holeImpactField = detail::holeAvalancheDrivingField(
            impactIonizationConfig_, holeCoefficientField, electricField, pAvg);
        const Real edgeArea =
            0.5 * h * couple_e * impactIonizationConfig_.sourceGeometryScale;

        Real electronSource = 0.0;
        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Electron, electronMobilityField, &mobilityConfig_, &psi);
        if (mun > 0.0) {
            const Real fluxN = std::abs(sgElectronContinuityFluxFromQuasiFermiVariableNi(
                niI, niJ, psi_i, psi_j, phin_i, phin_j, Vt_, mun * Vt_ / h, bgnEnabled_));
            const Real alphaN = impactIonization_->electronCoefficient(electronImpactField);
            electronSource = alphaN * fluxN * edgeArea;
        }

        Real holeSource = 0.0;
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Hole, holeMobilityField, &mobilityConfig_, &psi);
        if (mup > 0.0) {
            const Real fluxP = std::abs(sgHoleContinuityFluxFromQuasiFermiVariableNi(
                niI, niJ, psi_i, psi_j, phip_i, phip_j, Vt_, mup * Vt_ / h, bgnEnabled_));
            const Real alphaP = impactIonization_->holeCoefficient(holeImpactField);
            holeSource = alphaP * fluxP * edgeArea;
        }

        const detail::EdgeAvalancheDirectionalWeights weights =
            detail::edgeAvalancheDirectionalWeights(
                edgeCells_,
                mesh_,
                e,
                [&](Index node) {
                    const int nodeIndex = static_cast<int>(node);
                    const Real psiNode = node == idxI ? psi_i : (node == idxJ ? psi_j : psi(nodeIndex));
                    const Real phinNode = phinAt(node);
                    const Real nNode = ni_[node] * limitedExp((psiNode - phinNode) / Vt_);
                    return detail::electronQfForAvalancheGradient(
                        psiNode, phinNode, nNode, ni_[node], Vt_, impactIonizationConfig_);
                },
                [&](Index node) {
                    const int nodeIndex = static_cast<int>(node);
                    const Real psiNode = node == idxI ? psi_i : (node == idxJ ? psi_j : psi(nodeIndex));
                    const Real phipNode = phipAt(node);
                    const Real pNode = ni_[node] * limitedExp((phipNode - psiNode) / Vt_);
                    return detail::holeQfForAvalancheGradient(
                        psiNode, phipNode, pNode, ni_[node], Vt_, impactIonizationConfig_);
                });
        sources.node0 = weights.electronNode0 * electronSource + weights.holeNode0 * holeSource;
        sources.node1 = weights.electronNode1 * electronSource + weights.holeNode1 * holeSource;
        return sources;
    };

    auto edgeElectronTransportFlux =
        [&](Index e, int i, int j, Real h,
            Real psi_i, Real psi_j, Real phin_i, Real phin_j) -> Real {
        const Real couple_e = couple_[e];
        if (h <= 1.0e-30 || couple_e <= 0.0)
            return 0.0;
        const Index idxI = static_cast<Index>(i);
        const Index idxJ = static_cast<Index>(j);
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h);
        const Real electronMobilityField = qfMobility
            ? std::abs((phin_j - phin_i) / h)
            : electricField;
        VectorXd psiForSurface;
        const VectorXd* psiForMobility = &psi;
        if (isSurfaceMobilityModel(mobilityConfig_)) {
            psiForSurface = psi;
            psiForSurface(i) = psi_i;
            psiForSurface(j) = psi_j;
            psiForMobility = &psiForSurface;
        }
        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Electron, electronMobilityField, &mobilityConfig_, psiForMobility);
        if (mun <= 0.0)
            return 0.0;
        const Real coef = mun * Vt_ * couple_e / h;
        if (bgnEnabled_) {
            return sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni_[idxI], ni_[idxJ], psi_i, psi_j, phin_i, phin_j, Vt_, coef);
        }
        if (ni_[idxI] == ni_[idxJ]) {
            return sgElectronContinuityFluxFromQuasiFermiStable(
                ni_[idxI], psi_i, psi_j, phin_i, phin_j, Vt_, coef);
        }
        const Real n_i = ni_[idxI] * limitedExp((psi_i - phin_i) / Vt_);
        const Real n_j = ni_[idxJ] * limitedExp((psi_j - phin_j) / Vt_);
        return sgElectronContinuityFlux(n_i, n_j, dpsi, Vt_, coef);
    };

    auto edgeHoleTransportFlux =
        [&](Index e, int i, int j, Real h,
            Real psi_i, Real psi_j, Real phip_i, Real phip_j) -> Real {
        const Real couple_e = couple_[e];
        if (h <= 1.0e-30 || couple_e <= 0.0)
            return 0.0;
        const Index idxI = static_cast<Index>(i);
        const Index idxJ = static_cast<Index>(j);
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h);
        const Real holeMobilityField = qfMobility
            ? std::abs((phip_j - phip_i) / h)
            : electricField;
        VectorXd psiForSurface;
        const VectorXd* psiForMobility = &psi;
        if (isSurfaceMobilityModel(mobilityConfig_)) {
            psiForSurface = psi;
            psiForSurface(i) = psi_i;
            psiForSurface(j) = psi_j;
            psiForMobility = &psiForSurface;
        }
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Hole, holeMobilityField, &mobilityConfig_, psiForMobility);
        if (mup <= 0.0)
            return 0.0;
        const Real coef = mup * Vt_ * couple_e / h;
        if (bgnEnabled_) {
            return sgHoleContinuityFluxFromQuasiFermiVariableNi(
                ni_[idxI], ni_[idxJ], psi_i, psi_j, phip_i, phip_j, Vt_, coef);
        }
        if (ni_[idxI] == ni_[idxJ]) {
            return sgHoleContinuityFluxFromQuasiFermiStable(
                ni_[idxI], psi_i, psi_j, phip_i, phip_j, Vt_, coef);
        }
        const Real p_i = ni_[idxI] * limitedExp((phip_i - psi_i) / Vt_);
        const Real p_j = ni_[idxJ] * limitedExp((phip_j - psi_j) / Vt_);
        return sgHoleContinuityFlux(p_i, p_j, dpsi, Vt_, coef);
    };

    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30) continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real psi_i = x(psiOffset() + i) * potentialScale;
        const Real psi_j = x(psiOffset() + j) * potentialScale;
        const Real phin_i = x(phinOffset() + i) * potentialScale;
        const Real phin_j = x(phinOffset() + j) * potentialScale;
        const Real phip_i = x(phipOffset() + i) * potentialScale;
        const Real phip_j = x(phipOffset() + j) * potentialScale;
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h);
        const Real electronMobilityField =
            qfMobility ? std::abs((phin_j - phin_i) / h) : electricField;
        const Real holeMobilityField =
            qfMobility ? std::abs((phip_j - phip_i) / h) : electricField;
        const Real u = dpsi / Vt_;
        const Real Bu = bernoulli(u);
        const Real dBu = bernoulliDerivative(u);
        const Real Bminus = bernoulli(-u);
        const Real dBminusDu = -bernoulliDerivative(-u);

        const Real eps = detail::edgeEpsilon(edgeCells_, mesh_, matdb_, e);
        const Real G = eps * couple_[e] / h;
        add(psiOffset() + i, psiOffset() + i,  G);
        add(psiOffset() + i, psiOffset() + j, -G);
        add(psiOffset() + j, psiOffset() + i, -G);
        add(psiOffset() + j, psiOffset() + j,  G);

        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Electron,
            electronMobilityField,
            &mobilityConfig_,
            &psi);
        if (mun > 0.0) {
            hasElectronContribution[static_cast<std::size_t>(i)] = true;
            hasElectronContribution[static_cast<std::size_t>(j)] = true;

            if (transportMobilityDerivative) {
                const Real vals[4] = {psi_i, psi_j, phin_i, phin_j};
                const int cols[4] = {
                    psiOffset() + i, psiOffset() + j,
                    phinOffset() + i, phinOffset() + j,
                };
                for (int k = 0; k < 4; ++k) {
                    const Real step = 1.0e-6 * std::max(1.0, std::abs(vals[k]));
                    Real vp[4] = {vals[0], vals[1], vals[2], vals[3]};
                    Real vm[4] = {vals[0], vals[1], vals[2], vals[3]};
                    vp[k] += step;
                    vm[k] -= step;
                    const Real fp = edgeElectronTransportFlux(
                        e, i, j, h, vp[0], vp[1], vp[2], vp[3]);
                    const Real fm = edgeElectronTransportFlux(
                        e, i, j, h, vm[0], vm[1], vm[2], vm[3]);
                    const Real dF = (fp - fm) / (2.0 * step);
                    add(phinOffset() + i, cols[k], dF);
                    add(phinOffset() + j, cols[k], -dF);
                }
            } else {
                const Real coef = mun * Vt_ * couple_[e] / h;
                Real dF_dpsi_i = 0.0;
                Real dF_dpsi_j = 0.0;
                Real dF_dphin_i = 0.0;
                Real dF_dphin_j = 0.0;
                const Index idxI = static_cast<Index>(i);
                const Index idxJ = static_cast<Index>(j);
                const Real eta = u + std::log(ni_[idxJ] / ni_[idxI]);
                const Real Bplus = bernoulli(eta);
                const Real Bminus = bernoulli(-eta);
                const Real dBplus = bernoulliDerivative(eta);
                const Real dBminusArg = bernoulliDerivative(-eta);
                const Real niI = ni_[idxI];
                const Real niJ = ni_[idxJ];
                const Real nI = niI * limitedExp((psi_i - phin_i) / Vt_);
                const Real nJ = niJ * limitedExp((psi_j - phin_j) / Vt_);
                dF_dpsi_i = coef / Vt_ * ((dBminusArg + Bminus) * nI + dBplus * nJ);
                dF_dpsi_j = coef / Vt_ * (-dBminusArg * nI - (dBplus + Bplus) * nJ);
                dF_dphin_i = coef * Bminus * (-nI / Vt_);
                dF_dphin_j = coef * Bplus * ( nJ / Vt_);

                add(phinOffset() + i, psiOffset() + i, dF_dpsi_i);
                add(phinOffset() + i, psiOffset() + j, dF_dpsi_j);
                add(phinOffset() + i, phinOffset() + i, dF_dphin_i);
                add(phinOffset() + i, phinOffset() + j, dF_dphin_j);
                add(phinOffset() + j, psiOffset() + i, -dF_dpsi_i);
                add(phinOffset() + j, psiOffset() + j, -dF_dpsi_j);
                add(phinOffset() + j, phinOffset() + i, -dF_dphin_i);
                add(phinOffset() + j, phinOffset() + j, -dF_dphin_j);
            }
        }

        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Hole,
            holeMobilityField,
            &mobilityConfig_,
            &psi);
        if (mup > 0.0) {
            hasHoleContribution[static_cast<std::size_t>(i)] = true;
            hasHoleContribution[static_cast<std::size_t>(j)] = true;

            if (transportMobilityDerivative) {
                const Real vals[4] = {psi_i, psi_j, phip_i, phip_j};
                const int cols[4] = {
                    psiOffset() + i, psiOffset() + j,
                    phipOffset() + i, phipOffset() + j,
                };
                for (int k = 0; k < 4; ++k) {
                    const Real step = 1.0e-6 * std::max(1.0, std::abs(vals[k]));
                    Real vp[4] = {vals[0], vals[1], vals[2], vals[3]};
                    Real vm[4] = {vals[0], vals[1], vals[2], vals[3]};
                    vp[k] += step;
                    vm[k] -= step;
                    const Real fp = edgeHoleTransportFlux(
                        e, i, j, h, vp[0], vp[1], vp[2], vp[3]);
                    const Real fm = edgeHoleTransportFlux(
                        e, i, j, h, vm[0], vm[1], vm[2], vm[3]);
                    const Real dF = (fp - fm) / (2.0 * step);
                    add(phipOffset() + i, cols[k], dF);
                    add(phipOffset() + j, cols[k], -dF);
                }
            } else {
                const Real coef = mup * Vt_ * couple_[e] / h;
                Real dF_dpsi_i = 0.0;
                Real dF_dpsi_j = 0.0;
                Real dF_dphip_i = 0.0;
                Real dF_dphip_j = 0.0;
                const Index idxI = static_cast<Index>(i);
                const Index idxJ = static_cast<Index>(j);
                const Real eta = u + std::log(ni_[idxI] / ni_[idxJ]);
                const Real Bplus = bernoulli(eta);
                const Real Bminus = bernoulli(-eta);
                const Real dBplus = bernoulliDerivative(eta);
                const Real dBminusArg = bernoulliDerivative(-eta);
                const Real niI = ni_[idxI];
                const Real niJ = ni_[idxJ];
                const Real pI = niI * limitedExp((phip_i - psi_i) / Vt_);
                const Real pJ = niJ * limitedExp((phip_j - psi_j) / Vt_);
                dF_dpsi_i = coef / Vt_ * (-(dBplus + Bplus) * pI - dBminusArg * pJ);
                dF_dpsi_j = coef / Vt_ * (dBplus * pI + (dBminusArg + Bminus) * pJ);
                dF_dphip_i = coef * Bplus * ( pI / Vt_);
                dF_dphip_j = coef * Bminus * (-pJ / Vt_);

                add(phipOffset() + i, psiOffset() + i, dF_dpsi_i);
                add(phipOffset() + i, psiOffset() + j, dF_dpsi_j);
                add(phipOffset() + i, phipOffset() + i, dF_dphip_i);
                add(phipOffset() + i, phipOffset() + j, dF_dphip_j);
                add(phipOffset() + j, psiOffset() + i, -dF_dpsi_i);
                add(phipOffset() + j, psiOffset() + j, -dF_dpsi_j);
                add(phipOffset() + j, phipOffset() + i, -dF_dphip_i);
                add(phipOffset() + j, phipOffset() + j, -dF_dphip_j);
            }
        }

        if (sgCurrentAvalanche) {
            // Finite-difference the directionally partitioned nodal avalanche
            // source with respect to the six endpoint potentials. This captures
            // the flux (carrier density), driving-field (alpha), field-dependent
            // mobility, and qF-gradient partition derivatives together.
            auto phinAt = [&](Index node) { return phinState(static_cast<int>(node)); };
            auto phipAt = [&](Index node) { return phipState(static_cast<int>(node)); };
            const EdgeAvalancheNodeSources base = edgeAvalancheNodeSources(
                e, i, j, h, psi_i, psi_j, phinAt, phipAt);
            std::vector<Index> qfStencilNodes = {static_cast<Index>(i), static_cast<Index>(j)};
            if (e < edgeCells_.size()) {
                for (const Index cellId : edgeCells_[e]) {
                    const Cell& cell = mesh_.getCell(cellId);
                    for (const Index node : cell.node_ids) {
                        if (std::find(qfStencilNodes.begin(), qfStencilNodes.end(), node) ==
                            qfStencilNodes.end()) {
                            qfStencilNodes.push_back(node);
                        }
                    }
                }
            }

            std::vector<int> cols;
            std::vector<Real> dS0;
            std::vector<Real> dS1;
            bool anyNonzero = false;
            auto appendDerivative = [&](int col, const EdgeAvalancheNodeSources& sp,
                                        const EdgeAvalancheNodeSources& sm, Real step) {
                cols.push_back(col);
                dS0.push_back((sp.node0 - sm.node0) / (2.0 * step));
                dS1.push_back((sp.node1 - sm.node1) / (2.0 * step));
                if (dS0.back() != 0.0 || dS1.back() != 0.0)
                    anyNonzero = true;
            };

            const Real psiVals[2] = {psi_i, psi_j};
            const int psiCols[2] = {psiOffset() + i, psiOffset() + j};
            for (int k = 0; k < 2; ++k) {
                const Real step = 1.0e-7 * std::max(1.0, std::abs(psiVals[k]));
                Real psiP[2] = {psi_i, psi_j};
                Real psiM[2] = {psi_i, psi_j};
                psiP[k] += step;
                psiM[k] -= step;
                const EdgeAvalancheNodeSources sp = edgeAvalancheNodeSources(
                    e, i, j, h, psiP[0], psiP[1], phinAt, phipAt);
                const EdgeAvalancheNodeSources sm = edgeAvalancheNodeSources(
                    e, i, j, h, psiM[0], psiM[1], phinAt, phipAt);
                appendDerivative(psiCols[k], sp, sm, step);
            }

            for (const Index node : qfStencilNodes) {
                const int nodeIndex = static_cast<int>(node);
                const Real phinValue = phinState(nodeIndex);
                const Real step = 1.0e-7 * std::max(1.0, std::abs(phinValue));
                auto phinPlusAt = [&](Index queryNode) {
                    return queryNode == node ? phinValue + step
                                             : phinState(static_cast<int>(queryNode));
                };
                auto phinMinusAt = [&](Index queryNode) {
                    return queryNode == node ? phinValue - step
                                             : phinState(static_cast<int>(queryNode));
                };
                const EdgeAvalancheNodeSources sp = edgeAvalancheNodeSources(
                    e, i, j, h, psi_i, psi_j, phinPlusAt, phipAt);
                const EdgeAvalancheNodeSources sm = edgeAvalancheNodeSources(
                    e, i, j, h, psi_i, psi_j, phinMinusAt, phipAt);
                appendDerivative(phinOffset() + nodeIndex, sp, sm, step);
            }

            for (const Index node : qfStencilNodes) {
                const int nodeIndex = static_cast<int>(node);
                const Real phipValue = phipState(nodeIndex);
                const Real step = 1.0e-7 * std::max(1.0, std::abs(phipValue));
                auto phipPlusAt = [&](Index queryNode) {
                    return queryNode == node ? phipValue + step
                                             : phipState(static_cast<int>(queryNode));
                };
                auto phipMinusAt = [&](Index queryNode) {
                    return queryNode == node ? phipValue - step
                                             : phipState(static_cast<int>(queryNode));
                };
                const EdgeAvalancheNodeSources sp = edgeAvalancheNodeSources(
                    e, i, j, h, psi_i, psi_j, phinAt, phipPlusAt);
                const EdgeAvalancheNodeSources sm = edgeAvalancheNodeSources(
                    e, i, j, h, psi_i, psi_j, phinAt, phipMinusAt);
                appendDerivative(phipOffset() + nodeIndex, sp, sm, step);
            }

            if (base.node0 != 0.0 || base.node1 != 0.0 || anyNonzero) {
                const int node0Rows[2] = {phinOffset() + i, phipOffset() + i};
                const int node1Rows[2] = {phinOffset() + j, phipOffset() + j};
                for (int row : node0Rows)
                    for (std::size_t k = 0; k < cols.size(); ++k)
                        add(row, cols[k], -dS0[k]);
                for (int row : node1Rows)
                    for (std::size_t k = 0; k < cols.size(); ++k)
                        add(row, cols[k], -dS1[k]);
                hasElectronContribution[static_cast<std::size_t>(i)] = true;
                hasElectronContribution[static_cast<std::size_t>(j)] = true;
                hasHoleContribution[static_cast<std::size_t>(i)] = true;
                hasHoleContribution[static_cast<std::size_t>(j)] = true;
            }
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        const Real dni_dpsi = n(ii) / Vt_;
        const Real dni_dphin = -n(ii) / Vt_;
        const Real dpi_dpsi = -p(ii) / Vt_;
        const Real dpi_dphip = p(ii) / Vt_;

        add(psiOffset() + ii, psiOffset() + ii,
            -constants::q * (dpi_dpsi - dni_dpsi) * vol_[i]);
        add(psiOffset() + ii, phinOffset() + ii,
            -constants::q * (-dni_dphin) * vol_[i]);
        add(psiOffset() + ii, phipOffset() + ii,
            -constants::q * dpi_dphip * vol_[i]);

        const Real ni = ni_[i];
        if (ni <= 0.0)
            continue;

        if (recombination_.srhEnabled() || recombination_.augerEnabled()) {
            const Real dPhi =
                (x(phipOffset() + ii) - x(phinOffset() + ii)) * potentialScale;
            const Real np = (ni > 0.0) ? ni * ni * std::exp(dPhi / Vt_) : n(ii) * p(ii);
            const Real excessProduct = (ni > 0.0)
                ? ni * ni * std::expm1(dPhi / Vt_)
                : n(ii) * p(ii);
            const Real R = recombination_.totalRateFromExcessProduct(
                excessProduct, n(ii), p(ii), ni);
            if (R != 0.0) {
                const auto deriv = recombination_.totalRateDerivativesFromExcessProduct(
                    excessProduct, n(ii), p(ii), ni);

                const Real dExcess_dphin = -np / Vt_;
                const Real dExcess_dphip =  np / Vt_;
                const Real dR_dpsi = deriv.dRateDn * dni_dpsi + deriv.dRateDp * dpi_dpsi;
                const Real dR_dphin = deriv.dRateDn * dni_dphin
                                     + deriv.dRateDExcess * dExcess_dphin;
                const Real dR_dphip = deriv.dRateDp * dpi_dphip
                                     + deriv.dRateDExcess * dExcess_dphip;

                add(phinOffset() + ii, psiOffset() + ii, dR_dpsi * vol_[i]);
                add(phinOffset() + ii, phinOffset() + ii, dR_dphin * vol_[i]);
                add(phinOffset() + ii, phipOffset() + ii, dR_dphip * vol_[i]);
                add(phipOffset() + ii, psiOffset() + ii, dR_dpsi * vol_[i]);
                add(phipOffset() + ii, phinOffset() + ii, dR_dphin * vol_[i]);
                add(phipOffset() + ii, phipOffset() + ii, dR_dphip * vol_[i]);

                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }

        if (impactIonizationEnabled_ && sgCurrentAvalanche) {
            // The SG edge-current avalanche probe has nonlocal edge derivatives.
            // Omit source derivatives here rather than applying the legacy
            // node-local approximation to the wrong discretization.
        } else if (impactIonizationEnabled_) {
            const Real electronImpactField = detail::electronAvalancheDrivingField(
                impactIonizationConfig_, nodeElectronDrivingFields[i], nodeElectricFields[i], n(ii));
            const Real holeImpactField = detail::holeAvalancheDrivingField(
                impactIonizationConfig_, nodeHoleDrivingFields[i], nodeElectricFields[i], p(ii));
            const Real alphaN = impactIonization_->electronCoefficient(electronImpactField);
            const Real alphaP = impactIonization_->holeCoefficient(holeImpactField);
            const Real G = detail::impactIonizationGenerationRate(
                impactIonizationConfig_,
                *impactIonization_,
                mobilityConfig_,
                *mobility_,
                nodeCells_,
                mesh_,
                doping_,
                cellMaterials_,
                i,
                nodeElectricFields[i],
                nodeElectronDrivingFields[i],
                nodeHoleDrivingFields[i],
                n(ii),
                p(ii));
            if (G != 0.0) {
                // Local carrier-density derivatives are included in the analytic Jacobian.
                // Driving-field and mobility derivatives are intentionally omitted because
                // both use edge/node max operations; finite-difference Jacobian remains
                // available for exact derivatives of configured avalanche runs.
                Real electronFactor = 0.0;
                Real holeFactor = 0.0;
                if (impactIonizationConfig_.generation == "current_density") {
                    const Real mun = detail::nodeMobility(
                        nodeCells_,
                        mesh_,
                        doping_,
                        *mobility_,
                        cellMaterials_,
                        i,
                        CarrierType::Electron,
                        electronImpactField);
                    const Real mup = detail::nodeMobility(
                        nodeCells_,
                        mesh_,
                        doping_,
                        *mobility_,
                        cellMaterials_,
                        i,
                        CarrierType::Hole,
                        holeImpactField);
                    electronFactor = alphaN * mun * std::abs(electronImpactField);
                    holeFactor = alphaP * mup * std::abs(holeImpactField);
                } else {
                    const Real denominator = alphaN * n(ii) + alphaP * p(ii);
                    const Real velocity = (denominator != 0.0) ? (G / denominator) : 0.0;
                    electronFactor = velocity * alphaN;
                    holeFactor = velocity * alphaP;
                }
                const Real dG_dpsi = electronFactor * dni_dpsi + holeFactor * dpi_dpsi;
                const Real dG_dphin = electronFactor * dni_dphin;
                const Real dG_dphip = holeFactor * dpi_dphip;

                add(phinOffset() + ii, psiOffset() + ii, -dG_dpsi * vol_[i]);
                add(phinOffset() + ii, phinOffset() + ii, -dG_dphin * vol_[i]);
                add(phinOffset() + ii, phipOffset() + ii, -dG_dphip * vol_[i]);
                add(phipOffset() + ii, psiOffset() + ii, -dG_dpsi * vol_[i]);
                add(phipOffset() + ii, phinOffset() + ii, -dG_dphin * vol_[i]);
                add(phipOffset() + ii, phipOffset() + ii, -dG_dphip * vol_[i]);

                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        if (!hasElectronContribution[static_cast<std::size_t>(ii)])
            addGauge(phinOffset() + ii, phinOffset() + ii, 1.0);
        if (!hasHoleContribution[static_cast<std::size_t>(ii)])
            addGauge(phipOffset() + ii, phipOffset() + ii, 1.0);
    }

    for (const auto& [node, value] : bcs.psi) {
        (void)value;
        triplets.emplace_back(psiOffset() + static_cast<int>(node),
                              psiOffset() + static_cast<int>(node), 1.0);
    }
    for (const auto& [node, value] : bcs.phin) {
        (void)value;
        triplets.emplace_back(phinOffset() + static_cast<int>(node),
                              phinOffset() + static_cast<int>(node), 1.0);
    }
    for (const auto& [node, value] : bcs.phip) {
        (void)value;
        triplets.emplace_back(phipOffset() + static_cast<int>(node),
                              phipOffset() + static_cast<int>(node), 1.0);
    }

    SparseMatrixd J(3 * N, 3 * N);
    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}

SparseMatrixd CoupledDDAssembler::finiteDifferenceJacobian(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    Real relativeStep) const
{
    if (relativeStep <= 0.0)
        throw std::invalid_argument(
            "CoupledDDAssembler::finiteDifferenceJacobian: relativeStep must be positive.");

    const int M = x.size();
    const VectorXd r0 = residual(x, bcs);
    std::vector<Eigen::Triplet<double>> triplets;
    // Heuristic: for a 2-D mesh with M = 3*N unknowns and sparse connectivity
    // (average ~6-7 neighbours per node), the Jacobian has O(N * 7 * 9) ~= 63*N
    // non-zeros.  Reserving M * 7 avoids the M^2 allocation that would OOM for
    // any realistically-sized mesh while still keeping reallocations rare.
    triplets.reserve(static_cast<std::size_t>(M) * 7);

    // Minimum absolute perturbation to prevent h == 0 when x(col) == 0.
    constexpr Real minAbsStep = 1.0e-15;

    for (int col = 0; col < M; ++col) {
        const Real h = std::max(relativeStep * std::max(1.0, std::abs(x(col))), minAbsStep);
        VectorXd xp = x;
        xp(col) += h;
        const VectorXd rp = residual(xp, bcs);
        const VectorXd dr = (rp - r0) / h;
        for (int row = 0; row < M; ++row) {
            if (dr(row) != 0.0)
                triplets.emplace_back(row, col, dr(row));
        }
    }

    SparseMatrixd J(M, M);
    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}

} // namespace vela
