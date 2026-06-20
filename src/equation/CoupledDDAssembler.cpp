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
        detail::usesDensityGradientAvalancheCurrent(impactIonizationConfig_);
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
            const Real nFlux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni_[idxI],
                ni_[idxJ],
                psi_i,
                psi_j,
                phin_i,
                phin_j,
                Vt_,
                coef);
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
            const Real pFlux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                ni_[idxI],
                ni_[idxJ],
                psi_i,
                psi_j,
                phip_i,
                phip_j,
                Vt_,
                coef);
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
        detail::usesDensityGradientAvalancheCurrent(impactIonizationConfig_);
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
                coef);
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
                coef);
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
    const bool sgCurrentAvalanche = impactIonizationEnabled_ &&
        detail::usesDensityGradientAvalancheCurrent(impactIonizationConfig_);

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

            if (sgCurrentAvalanche) {
                const Real electronCoefficientField = qfImpact
                    ? std::abs((phin_j - phin_i) / h)
                    : electricField;
                const Real nAvg = 0.5 * (n(i) + n(j));
                const Real electronImpactField = detail::electronAvalancheDrivingField(
                    impactIonizationConfig_, electronCoefficientField, electricField, nAvg);
                const Real alphaN = impactIonization_->electronCoefficient(electronImpactField);
                const Real flux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                    ni_[idxI],
                    ni_[idxJ],
                    psi_i,
                    psi_j,
                    phin_i,
                    phin_j,
                    Vt_,
                    coef);
                if (alphaN > 0.0 && flux != 0.0) {
                    const Real edgeAreaToCouple = 0.5 * h;
                    const Real sourceFactor =
                        -0.5 * alphaN * edgeAreaToCouple * (flux > 0.0 ? 1.0 : -1.0);
                    const int rows[] = {
                        phinOffset() + i,
                        phinOffset() + j,
                        phipOffset() + i,
                        phipOffset() + j,
                    };
                    for (int row : rows) {
                        add(row, psiOffset() + i, sourceFactor * dF_dpsi_i);
                        add(row, psiOffset() + j, sourceFactor * dF_dpsi_j);
                        add(row, phinOffset() + i, sourceFactor * dF_dphin_i);
                        add(row, phinOffset() + j, sourceFactor * dF_dphin_j);
                    }
                }
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

            if (sgCurrentAvalanche) {
                const Real holeCoefficientField = qfImpact
                    ? std::abs((phip_j - phip_i) / h)
                    : electricField;
                const Real pAvg = 0.5 * (p(i) + p(j));
                const Real holeImpactField = detail::holeAvalancheDrivingField(
                    impactIonizationConfig_, holeCoefficientField, electricField, pAvg);
                const Real alphaP = impactIonization_->holeCoefficient(holeImpactField);
                const Real flux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                    ni_[idxI],
                    ni_[idxJ],
                    psi_i,
                    psi_j,
                    phip_i,
                    phip_j,
                    Vt_,
                    coef);
                if (alphaP > 0.0 && flux != 0.0) {
                    const Real edgeAreaToCouple = 0.5 * h;
                    const Real sourceFactor =
                        -0.5 * alphaP * edgeAreaToCouple * (flux > 0.0 ? 1.0 : -1.0);
                    const int rows[] = {
                        phinOffset() + i,
                        phinOffset() + j,
                        phipOffset() + i,
                        phipOffset() + j,
                    };
                    for (int row : rows) {
                        add(row, psiOffset() + i, sourceFactor * dF_dpsi_i);
                        add(row, psiOffset() + j, sourceFactor * dF_dpsi_j);
                        add(row, phipOffset() + i, sourceFactor * dF_dphip_i);
                        add(row, phipOffset() + j, sourceFactor * dF_dphip_j);
                    }
                }
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
