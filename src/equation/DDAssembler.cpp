#include "vela/equation/DDAssembler.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/ScharfetterGummel.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace vela {

namespace {

VectorXd electronQuasiFermiFromDensity(const VectorXd& psi,
                                       const VectorXd& n,
                                       const std::vector<Real>& ni,
                                       Real Vt,
                                       const DDScalingSpec& scaling)
{
    VectorXd phin = psi;
    for (int i = 0; i < psi.size(); ++i) {
        const Real n_si = scaling.enabled ? n(i) * scaling.C0 : n(i);
        if (ni[static_cast<Index>(i)] > 0.0 && n_si > 0.0)
            phin(i) = psi(i) - Vt * std::log(n_si / ni[static_cast<Index>(i)]);
    }
    return phin;
}

VectorXd holeQuasiFermiFromDensity(const VectorXd& psi,
                                   const VectorXd& p,
                                   const std::vector<Real>& ni,
                                   Real Vt,
                                   const DDScalingSpec& scaling)
{
    VectorXd phip = psi;
    for (int i = 0; i < psi.size(); ++i) {
        const Real p_si = scaling.enabled ? p(i) * scaling.C0 : p(i);
        if (ni[static_cast<Index>(i)] > 0.0 && p_si > 0.0)
            phip(i) = psi(i) + Vt * std::log(p_si / ni[static_cast<Index>(i)]);
    }
    return phip;
}

}

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

DDAssembler::DDAssembler(const DeviceMesh&       mesh,
                         const MaterialDatabase& matdb,
                         const DopingModel&      doping,
                         double                  Vt,
                         double                  taun,
                         double                  taup)
    : DDAssembler(mesh, matdb, doping, Vt, taun, taup, {}, {})
{}

DDAssembler::DDAssembler(const DeviceMesh&       mesh,
                         const MaterialDatabase& matdb,
                         const DopingModel&      doping,
                         double                  Vt,
                         double                  taun,
                         double                  taup,
                         std::vector<RegionFixedChargeSpec> fixedCharges,
                         std::vector<InterfaceSheetChargeSpec> sheetCharges)
    : DDAssembler(mesh,
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

DDAssembler::DDAssembler(const DeviceMesh&               mesh,
                         const MaterialDatabase&         matdb,
                         const DopingModel&              doping,
                         double                          Vt,
                         const MobilityModelConfig&      mobilityConfig,
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
          "DDAssembler",
          mesh,
          matdb,
          doping,
          bandgapNarrowingConfig,
          Vt))
    , edgeCells_(detail::buildEdgeCellMap(mesh))
    , nodeCells_(detail::buildNodeCellMap(mesh))
    , vol_(detail::computeNodeVolumes(mesh))
    , couple_(detail::computeEdgeCouplings(mesh))
    , fixedInterfaceChargeRhs_(detail::computeFixedAndInterfaceChargeRhs(
          mesh, edgeCells_, fixedCharges, sheetCharges, "DDAssembler"))
    , scaling_(scaling)
    , A_(static_cast<int>(mesh.numNodes()),
         static_cast<int>(mesh.numNodes()))
    , b_(VectorXd::Zero(static_cast<int>(mesh.numNodes())))
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
                "DDAssembler: scaling references must be positive and finite when scaling is enabled.");
        }
    }
}

// ---------------------------------------------------------------------------
// Dirichlet boundary conditions
// ---------------------------------------------------------------------------

void DDAssembler::applyDirichlet(const std::unordered_map<Index, Real>& bcs)
{
    // The caller chooses which contact nodes belong to each scalar equation:
    // MOS source/drain/body contacts pass carrier-density Dirichlet maps, and
    // the electrostatic solve passes the corresponding potential map.
    detail::applyDirichletBC(A_, b_, bcs);
}

// ---------------------------------------------------------------------------
// Poisson with carriers
// ---------------------------------------------------------------------------

void DDAssembler::assemblePoissonWithCarriers(const VectorXd& n,
                                              const VectorXd& p,
                                              const VectorXd& psi)
{
    const Index N = mesh_.numNodes();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = VectorXd::Zero(static_cast<int>(N));

    // Edge flux terms
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (h < 1.0e-30) continue;

        const Real eps = detail::edgeEpsilon(edgeCells_, mesh_, matdb_, e);
        const Real G   = eps * couple_[e] / h;

        auto i = static_cast<int>(edge.n0);
        auto j = static_cast<int>(edge.n1);

        const Real matrixScale = scaling_.enabled
            ? (1.0 / scaling_.permittivityReference_F_per_m)
            : 1.0;
        triplets.emplace_back(i, i,  G * matrixScale);
        triplets.emplace_back(j, j,  G * matrixScale);
        triplets.emplace_back(i, j, -G * matrixScale);
        triplets.emplace_back(j, i, -G * matrixScale);
    }

    A_.setFromTriplets(triplets.begin(), triplets.end());

    // Linearised carrier charge:
    //   q*(n+p)/Vt added to diagonal; adjusted RHS keeps the equation
    //   equivalent to the nonlinear form around the current psi.
    for (Index i = 0; i < N; ++i) {
        const int  ii     = static_cast<int>(i);
        const Real ni_v   = scaling_.enabled ? n(ii) * scaling_.C0 : n(ii);
        const Real pi_v   = scaling_.enabled ? p(ii) * scaling_.C0 : p(ii);
        const Real psi_v  = scaling_.enabled ? psi(ii) * scaling_.V0 : psi(ii);
        const Real vol_i  = vol_[i];

        const Real diagCarrier = constants::q * (ni_v + pi_v) / Vt_ * vol_i;
        const Real matrixScale = scaling_.enabled
            ? (1.0 / scaling_.permittivityReference_F_per_m)
            : 1.0;
        A_.coeffRef(ii, ii) += diagCarrier * matrixScale;

        const Real rhs_si = constants::q *
                 (pi_v - ni_v + doping_.netDoping(i)) * vol_i
                 + diagCarrier * psi_v;
        b_(ii) = scaling_.enabled
            ? rhs_si / (scaling_.permittivityReference_F_per_m * scaling_.V0)
            : rhs_si;
    }

    if (scaling_.enabled)
        b_ += fixedInterfaceChargeRhs_ /
            (scaling_.permittivityReference_F_per_m * scaling_.V0);
    else
        b_ += fixedInterfaceChargeRhs_;
}

// ---------------------------------------------------------------------------
// Electron continuity
// ---------------------------------------------------------------------------

void DDAssembler::assembleElectronContinuity(const VectorXd& psi,
                                             const VectorXd& n_old,
                                             const VectorXd& p_old)
{
    const Index N = mesh_.numNodes();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = VectorXd::Zero(static_cast<int>(N));

    const Real temperature_K = Vt_ * constants::q / constants::kb;
    const std::vector<Material> cellMaterials =
        detail::buildCellMaterials(mesh_, matdb_, temperature_K);
    const VectorXd psiForMobility = scaling_.enabled ? (psi * scaling_.V0) : psi;
    const VectorXd phinForMobility =
        electronQuasiFermiFromDensity(psiForMobility, n_old, ni_, Vt_, scaling_);
    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";

    // SG matrix entries from all edges
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (h < 1.0e-30) continue;

        const Real psi0 = scaling_.enabled
            ? psi(static_cast<int>(edge.n0)) * scaling_.V0
            : psi(static_cast<int>(edge.n0));
        const Real psi1 = scaling_.enabled
            ? psi(static_cast<int>(edge.n1)) * scaling_.V0
            : psi(static_cast<int>(edge.n1));
        const Real electricField = std::abs((psi1 - psi0) / h);
        const Real electronMobilityField = qfMobility
            ? std::abs((phinForMobility(static_cast<int>(edge.n1)) -
                        phinForMobility(static_cast<int>(edge.n0))) / h)
            : electricField;
        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials, e, CarrierType::Electron,
            electronMobilityField,
            &mobilityConfig_,
            &psiForMobility);
        if (mun <= 0.0) continue; // skip insulator edges

        const Real coef = scaling_.enabled
            ? (mun / scaling_.mu0) * couple_[e] / h
            : mun * Vt_ * couple_[e] / h;

        auto i = static_cast<int>(edge.n0);
        auto j = static_cast<int>(edge.n1);

        const Real dpsi = scaling_.enabled
            ? (psi(j) - psi(i))
            : (psi(j) - psi(i));
        const SGEdgeWeights weights = scaling_.enabled
            ? sgEdgeWeights(dpsi, 1.0)
            : sgEdgeWeights(dpsi, Vt_);

        // Electron continuity flux from i to j:
        //   F_nij = coef * (B(-u) * n_i - B(+u) * n_j)
        triplets.emplace_back(i, i,  coef * weights.b_minus);
        triplets.emplace_back(j, j,  coef * weights.b_plus);
        triplets.emplace_back(i, j, -coef * weights.b_plus);
        triplets.emplace_back(j, i, -coef * weights.b_minus);
    }

    A_.setFromTriplets(triplets.begin(), triplets.end());

    const VectorXd psi_si = psiForMobility;
    const std::vector<Real> nodeElectricFields = impactIonizationEnabled_
        ? detail::computeNodeElectricFields(psi_si, mesh_)
        : std::vector<Real>{};
    const bool qfImpact = impactIonizationConfig_.drivingForce == "quasi_fermi_gradient";
    const VectorXd phipForImpact =
        holeQuasiFermiFromDensity(psi_si, p_old, ni_, Vt_, scaling_);
    const std::vector<Real> nodeElectronDrivingFields = (impactIonizationEnabled_ && qfImpact)
        ? detail::computeNodeScalarGradientMagnitudes(phinForMobility, mesh_)
        : nodeElectricFields;
    const std::vector<Real> nodeHoleDrivingFields = (impactIonizationEnabled_ && qfImpact)
        ? detail::computeNodeScalarGradientMagnitudes(phipForImpact, mesh_)
        : nodeElectricFields;
    const bool sgCurrentAvalanche = impactIonizationEnabled_ &&
        detail::usesDensityGradientAvalancheCurrent(impactIonizationConfig_);
    VectorXd n_physical = n_old;
    VectorXd p_physical = p_old;
    if (scaling_.enabled) {
        n_physical *= scaling_.C0;
        p_physical *= scaling_.C0;
    }
    const std::vector<Real> sgAvalancheSourceIntegrals = sgCurrentAvalanche
        ? detail::sgEdgeCurrentAvalancheSourceIntegrals(
            impactIonizationConfig_,
            *impactIonization_,
            mobilityConfig_,
            *mobility_,
            edgeCells_,
            mesh_,
            doping_,
            cellMaterials,
            psi_si,
            phinForMobility,
            phipForImpact,
            n_physical,
            p_physical,
            ni_,
            Vt_)
        : std::vector<Real>{};

    // Recombination source term linearised w.r.t. n.
    // Positive source derivatives move to the LHS diagonal; constants move to RHS.
    for (Index i = 0; i < N; ++i) {
        const int  ii    = static_cast<int>(i);
        const Real ni_i  = ni_[i];
        const Real n_v   = scaling_.enabled ? n_old(ii) * scaling_.C0 : n_old(ii);
        const Real p_v   = scaling_.enabled ? p_old(ii) * scaling_.C0 : p_old(ii);
        const Real vol_i = vol_[i];

        if (ni_i <= 0.0)
            continue;

        const RecombinationLinearization linearization =
            recombination_.electronLinearization(n_v, p_v, ni_i);
        if (scaling_.enabled) {
            A_.coeffRef(ii, ii) += linearization.diagonal * vol_i / scaling_.D0;
            b_(ii) += linearization.rhs * vol_i / (scaling_.C0 * scaling_.D0);
        } else {
            A_.coeffRef(ii, ii) += linearization.diagonal * vol_i;
            b_(ii) += linearization.rhs * vol_i;
        }
        if (impactIonizationEnabled_ && sgCurrentAvalanche) {
            const Real gen = sgAvalancheSourceIntegrals[i];
            b_(ii) += scaling_.enabled ? (gen / (scaling_.C0 * scaling_.D0)) : gen;
        } else if (impactIonizationEnabled_) {
            const Real gen = detail::impactIonizationGenerationRate(
                impactIonizationConfig_,
                *impactIonization_,
                mobilityConfig_,
                *mobility_,
                nodeCells_,
                mesh_,
                doping_,
                cellMaterials,
                i,
                nodeElectricFields[i],
                nodeElectronDrivingFields[i],
                nodeHoleDrivingFields[i],
                n_v,
                p_v) * vol_i;
            b_(ii) += scaling_.enabled ? (gen / (scaling_.C0 * scaling_.D0)) : gen;
        }
    }

    // Guard: if any diagonal is still zero (insulator node with all edges
    // skipped due to zero mobility), pin that node to n = 0 to avoid
    // a singular system.
    for (Index i = 0; i < N; ++i) {
        const int ii = static_cast<int>(i);
        if (A_.coeff(ii, ii) == 0.0) {
            A_.coeffRef(ii, ii) = 1.0;
            b_(ii) = 0.0;
        }
    }
}

// ---------------------------------------------------------------------------
// Hole continuity
// ---------------------------------------------------------------------------

void DDAssembler::assembleHoleContinuity(const VectorXd& psi,
                                         const VectorXd& n_old,
                                         const VectorXd& p_old)
{
    const Index N = mesh_.numNodes();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = VectorXd::Zero(static_cast<int>(N));

    const Real temperature_K = Vt_ * constants::q / constants::kb;
    const std::vector<Material> cellMaterials =
        detail::buildCellMaterials(mesh_, matdb_, temperature_K);
    const VectorXd psiForMobility = scaling_.enabled ? (psi * scaling_.V0) : psi;
    const VectorXd phipForMobility =
        holeQuasiFermiFromDensity(psiForMobility, p_old, ni_, Vt_, scaling_);
    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";

    // SG matrix entries for holes
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (h < 1.0e-30) continue;

        const Real psi0 = scaling_.enabled
            ? psi(static_cast<int>(edge.n0)) * scaling_.V0
            : psi(static_cast<int>(edge.n0));
        const Real psi1 = scaling_.enabled
            ? psi(static_cast<int>(edge.n1)) * scaling_.V0
            : psi(static_cast<int>(edge.n1));
        const Real electricField = std::abs((psi1 - psi0) / h);
        const Real holeMobilityField = qfMobility
            ? std::abs((phipForMobility(static_cast<int>(edge.n1)) -
                        phipForMobility(static_cast<int>(edge.n0))) / h)
            : electricField;
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials, e, CarrierType::Hole,
            holeMobilityField,
            &mobilityConfig_,
            &psiForMobility);
        if (mup <= 0.0) continue; // skip insulator edges

        const Real coef = scaling_.enabled
            ? (mup / scaling_.mu0) * couple_[e] / h
            : mup * Vt_ * couple_[e] / h;

        auto i = static_cast<int>(edge.n0);
        auto j = static_cast<int>(edge.n1);

        const Real dpsi = scaling_.enabled
            ? (psi(j) - psi(i))
            : (psi(j) - psi(i));
        const SGEdgeWeights weights = scaling_.enabled
            ? sgEdgeWeights(dpsi, 1.0)
            : sgEdgeWeights(dpsi, Vt_);

        // Hole continuity flux from i to j:
        //   F_pij = coef * (B(+u) * p_i - B(-u) * p_j)
        triplets.emplace_back(i, i,  coef * weights.b_plus);
        triplets.emplace_back(j, j,  coef * weights.b_minus);
        triplets.emplace_back(i, j, -coef * weights.b_minus);
        triplets.emplace_back(j, i, -coef * weights.b_plus);
    }

    A_.setFromTriplets(triplets.begin(), triplets.end());

    const VectorXd psi_si = psiForMobility;
    const std::vector<Real> nodeElectricFields = impactIonizationEnabled_
        ? detail::computeNodeElectricFields(psi_si, mesh_)
        : std::vector<Real>{};
    const bool qfImpact = impactIonizationConfig_.drivingForce == "quasi_fermi_gradient";
    const VectorXd phinForImpact =
        electronQuasiFermiFromDensity(psi_si, n_old, ni_, Vt_, scaling_);
    const std::vector<Real> nodeElectronDrivingFields = (impactIonizationEnabled_ && qfImpact)
        ? detail::computeNodeScalarGradientMagnitudes(phinForImpact, mesh_)
        : nodeElectricFields;
    const std::vector<Real> nodeHoleDrivingFields = (impactIonizationEnabled_ && qfImpact)
        ? detail::computeNodeScalarGradientMagnitudes(phipForMobility, mesh_)
        : nodeElectricFields;
    const bool sgCurrentAvalanche = impactIonizationEnabled_ &&
        detail::usesDensityGradientAvalancheCurrent(impactIonizationConfig_);
    VectorXd n_physical = n_old;
    VectorXd p_physical = p_old;
    if (scaling_.enabled) {
        n_physical *= scaling_.C0;
        p_physical *= scaling_.C0;
    }
    const std::vector<Real> sgAvalancheSourceIntegrals = sgCurrentAvalanche
        ? detail::sgEdgeCurrentAvalancheSourceIntegrals(
            impactIonizationConfig_,
            *impactIonization_,
            mobilityConfig_,
            *mobility_,
            edgeCells_,
            mesh_,
            doping_,
            cellMaterials,
            psi_si,
            phinForImpact,
            phipForMobility,
            n_physical,
            p_physical,
            ni_,
            Vt_)
        : std::vector<Real>{};

    // Recombination source term linearised w.r.t. p.
    // Positive source derivatives move to the LHS diagonal; constants move to RHS.
    for (Index i = 0; i < N; ++i) {
        const int  ii    = static_cast<int>(i);
        const Real ni_i  = ni_[i];
        const Real n_v   = scaling_.enabled ? n_old(ii) * scaling_.C0 : n_old(ii);
        const Real p_v   = scaling_.enabled ? p_old(ii) * scaling_.C0 : p_old(ii);
        const Real vol_i = vol_[i];

        if (ni_i <= 0.0)
            continue;

        const RecombinationLinearization linearization =
            recombination_.holeLinearization(n_v, p_v, ni_i);
        if (scaling_.enabled) {
            A_.coeffRef(ii, ii) += linearization.diagonal * vol_i / scaling_.D0;
            b_(ii) += linearization.rhs * vol_i / (scaling_.C0 * scaling_.D0);
        } else {
            A_.coeffRef(ii, ii) += linearization.diagonal * vol_i;
            b_(ii) += linearization.rhs * vol_i;
        }
        if (impactIonizationEnabled_ && sgCurrentAvalanche) {
            const Real gen = sgAvalancheSourceIntegrals[i];
            b_(ii) += scaling_.enabled ? (gen / (scaling_.C0 * scaling_.D0)) : gen;
        } else if (impactIonizationEnabled_) {
            const Real gen = detail::impactIonizationGenerationRate(
                impactIonizationConfig_,
                *impactIonization_,
                mobilityConfig_,
                *mobility_,
                nodeCells_,
                mesh_,
                doping_,
                cellMaterials,
                i,
                nodeElectricFields[i],
                nodeElectronDrivingFields[i],
                nodeHoleDrivingFields[i],
                n_v,
                p_v) * vol_i;
            b_(ii) += scaling_.enabled ? (gen / (scaling_.C0 * scaling_.D0)) : gen;
        }
    }

    // Guard: pin insulator nodes (zero-diagonal) to p = 0
    for (Index i = 0; i < N; ++i) {
        const int ii = static_cast<int>(i);
        if (A_.coeff(ii, ii) == 0.0) {
            A_.coeffRef(ii, ii) = 1.0;
            b_(ii) = 0.0;
        }
    }
}

} // namespace vela
