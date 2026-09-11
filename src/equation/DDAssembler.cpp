#include "vela/equation/DDAssembler.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/ScharfetterGummel.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <utility>
#include <vector>

namespace vela {

namespace {

VectorXd electronQuasiFermiFromDensity(const VectorXd& psi,
                                       const VectorXd& n,
                                       const std::vector<Real>& ni,
                                       const std::vector<Real>& Nc,
                                       Real Vt,
                                       const DDScalingSpec& scaling,
                                       const CarrierStatisticsConfig& statistics)
{
    VectorXd phin = psi;
    for (int i = 0; i < psi.size(); ++i) {
        const Real n_si = scaling.enabled ? n(i) * scaling.C0 : n(i);
        const Index node = static_cast<Index>(i);
        if (ni[node] > 0.0 && n_si > 0.0)
            phin(i) = electronQuasiFermiPotential(
                ni[node], Nc[node], psi(i), n_si, Vt, statistics);
    }
    return phin;
}

VectorXd holeQuasiFermiFromDensity(const VectorXd& psi,
                                   const VectorXd& p,
                                   const std::vector<Real>& ni,
                                   const std::vector<Real>& Nv,
                                   Real Vt,
                                   const DDScalingSpec& scaling,
                                   const CarrierStatisticsConfig& statistics)
{
    VectorXd phip = psi;
    for (int i = 0; i < psi.size(); ++i) {
        const Real p_si = scaling.enabled ? p(i) * scaling.C0 : p(i);
        const Index node = static_cast<Index>(i);
        if (ni[node] > 0.0 && p_si > 0.0)
            phip(i) = holeQuasiFermiPotential(
                ni[node], Nv[node], psi(i), p_si, Vt, statistics);
    }
    return phip;
}

Real reducedFermiEnergyFromDensity(Real density, Real densityOfStates)
{
    if (!(density > 0.0) || !(densityOfStates > 0.0))
        return -500.0;
    const Real ratio = density / densityOfStates;
    if (ratio > 0.0 && std::isfinite(ratio))
        return inverseFermiDiracHalf(ratio);
    // Preserve the nondegenerate eta ~= log(n/Nc) limit when the ratio
    // underflows even though both operands are positive and finite.
    return std::log(density) - std::log(densityOfStates);
}

Real fermiDensityDerivativeEtaFromDensity(Real density, Real densityOfStates)
{
    if (!(density > 0.0) || !(densityOfStates > 0.0))
        return 0.0;
    const Real eta = reducedFermiEnergyFromDensity(density, densityOfStates);
    const Real derivative = densityOfStates * fermiDiracHalfDerivative(eta);
    // Deep in the nondegenerate tail the analytic derivative can underflow;
    // dn/deta -> n is the correct limiting value.
    return derivative > 0.0 && std::isfinite(derivative) ? derivative : density;
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
                         DDScalingSpec scaling,
                         CarrierStatisticsConfig carrierStatistics)
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
    , impactIonizationCoupled_(
          impactIonizationConfig.model != "none" &&
          impactIonizationConfig.couplingMode == "self_consistent")
    , carrierStatistics_(std::move(carrierStatistics))
    , ni_(detail::buildValidatedEffectiveNodeNi(
          "DDAssembler",
          mesh,
          matdb,
          doping,
          bandgapNarrowingConfig,
          Vt))
    , Nc_(detail::buildNodeDensityOfStates(
          mesh, matdb, Vt * constants::q / constants::kb, true))
    , Nv_(detail::buildNodeDensityOfStates(
          mesh, matdb, Vt * constants::q / constants::kb, false))
    , cellMaterials_(detail::buildCellMaterials(
          mesh, matdb, Vt * constants::q / constants::kb))
    , edgeCells_(detail::buildEdgeCellMap(mesh))
    , nodeCells_(detail::buildNodeCellMap(mesh))
    , vol_(detail::computeNodeVolumes(mesh))
    , poissonChargeVol_(detail::computePoissonChargeVolumes(
          mesh, cellMaterials_, scaling.poissonChargeVolumePolicy))
    , couple_(detail::computeTransportEdgeCouplings(mesh))
    , fixedInterfaceChargeRhs_(detail::computeFixedAndInterfaceChargeRhs(
          mesh,
          edgeCells_,
          fixedCharges,
          sheetCharges,
          "DDAssembler",
          scaling.enabled ? scaling.chargeAreaFactor : 1.0,
          scaling.enabled ? scaling.chargeLineFactor : 1.0))
    , scaling_(scaling)
    , A_(static_cast<int>(mesh.numNodes()),
         static_cast<int>(mesh.numNodes()))
    , b_(VectorXd::Zero(static_cast<int>(mesh.numNodes())))
{
    (void)usesFermiDirac(carrierStatistics_);
    if (usesFermiDirac(carrierStatistics_)) {
        for (Index node = 0; node < mesh_.numNodes(); ++node) {
            if (ni_[node] > 0.0 && (!(Nc_[node] > 0.0) || !(Nv_[node] > 0.0))) {
                throw std::invalid_argument(
                    "DDAssembler: Fermi-Dirac statistics require Nc_m3 and Nv_m3 "
                    "for every transport node.");
            }
        }
    }
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

        if (!std::isfinite(G)) {
            std::ostringstream message;
            message << "DDAssembler::assemblePoissonWithCarriers: non-finite edge "
                    << "coefficient at edge " << e << " (nodes " << edge.n0
                    << ',' << edge.n1 << ", eps_F_per_m=" << eps
                    << ", couple=" << couple_[e] << ", length=" << h << ')';
            throw std::runtime_error(message.str());
        }

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
        const Real vol_i  = poissonChargeVol_[i];
        const Real chargeAreaFactor = scaling_.enabled ? scaling_.chargeAreaFactor : 1.0;

        Real electronDerivativeEta = ni_v;
        Real holeDerivativeEta = pi_v;
        if (usesFermiDirac(carrierStatistics_)) {
            electronDerivativeEta = fermiDensityDerivativeEtaFromDensity(ni_v, Nc_[i]);
            holeDerivativeEta = fermiDensityDerivativeEtaFromDensity(pi_v, Nv_[i]);
        }
        const Real diagCarrier = constants::q *
            (electronDerivativeEta + holeDerivativeEta) / Vt_ * vol_i * chargeAreaFactor;
        const Real matrixScale = scaling_.enabled
            ? (1.0 / scaling_.permittivityReference_F_per_m)
            : 1.0;
        if (!std::isfinite(diagCarrier) || !std::isfinite(matrixScale)) {
            std::ostringstream message;
            message << "DDAssembler::assemblePoissonWithCarriers: non-finite carrier "
                    << "diagonal at node " << i << " (n=" << ni_v
                    << ", p=" << pi_v << ", dn_deta=" << electronDerivativeEta
                    << ", dp_deta=" << holeDerivativeEta << ", Vt=" << Vt_
                    << ", volume=" << vol_i << ", charge_area_factor="
                    << chargeAreaFactor << ", matrix_scale=" << matrixScale
                    << ", concentration_scale=" << scaling_.C0 << ')';
            throw std::runtime_error(message.str());
        }
        A_.coeffRef(ii, ii) += diagCarrier * matrixScale;

        const Real rhs_si = constants::q *
                 (pi_v - ni_v + doping_.netDoping(i)) * vol_i * chargeAreaFactor
                 + diagCarrier * psi_v;
        if (!std::isfinite(rhs_si)) {
            std::ostringstream message;
            message << "DDAssembler::assemblePoissonWithCarriers: non-finite Poisson "
                    << "RHS at node " << i << " (n=" << ni_v << ", p=" << pi_v
                    << ", net_doping=" << doping_.netDoping(i) << ", volume="
                    << vol_i << ", charge_area_factor=" << chargeAreaFactor
                    << ", carrier_diagonal=" << diagCarrier << ", psi=" << psi_v
                    << ')';
            throw std::runtime_error(message.str());
        }
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
        electronQuasiFermiFromDensity(
            psiForMobility, n_old, ni_, Nc_, Vt_, scaling_, carrierStatistics_);
    if (mobilityConfig_.model == "ialmob") {
        const VectorXd phipForMobility = holeQuasiFermiFromDensity(
            psiForMobility,p_old,ni_,Nv_,Vt_,scaling_,carrierStatistics_);
        const Real densityScale=scaling_.enabled?scaling_.C0:1.;
        updateIalTransportState(mobilityConfig_,mesh_,doping_,psiForMobility,
            n_old*densityScale,p_old*densityScale,phinForMobility,phipForMobility);
    }
    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";
    const Real fieldFactor = scaling_.enabled ? scaling_.fieldFromCoordinateDeltaFactor : 1.0;
    const bool vectorQfMobility = qfMobility &&
        mobilityConfig_.highFieldGradientDiscretization == "transport_cell_vector";
    const std::vector<Real> electronVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials, phinForMobility, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> contactElectricMobilityFields =
        mobilityConfig_.contactElectricFieldFallback
        ? detail::contactCellElectricFieldEdgeMagnitudesForMobility(
              mobilityConfig_, mesh_, edgeCells_, cellMaterials,
              psiForMobility, fieldFactor)
        : std::vector<Real>{};
    const Real sourceIntegralFactor = scaling_.enabled
        ? scaling_.unitSystem.continuitySourceIntegralFactor()
        : 1.0;

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
        const Real electricField = std::abs((psi1 - psi0) / h) * fieldFactor;
        const Real electronMobilityField = detail::mobilityHighFieldDrivingField(
            mobilityConfig_, e,
            vectorQfMobility ? electronVectorMobilityFields[e]
                : std::abs((phinForMobility(static_cast<int>(edge.n1)) -
                            phinForMobility(static_cast<int>(edge.n0))) / h) * fieldFactor,
            electricField, contactElectricMobilityFields);
        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials, e, CarrierType::Electron,
            electronMobilityField,
            &mobilityConfig_,
            &psiForMobility);
        if (mun <= 0.0) continue; // skip insulator edges

        Real coef = scaling_.enabled
            ? (mun / scaling_.mu0) * fieldFactor * couple_[e] / h
            : mun * Vt_ * fieldFactor * couple_[e] / h;

        auto i = static_cast<int>(edge.n0);
        auto j = static_cast<int>(edge.n1);

        SGEdgeWeights weights;
        if (usesFermiDirac(carrierStatistics_)) {
            const Real n0 = scaling_.enabled ? n_old(i) * scaling_.C0 : n_old(i);
            const Real n1 = scaling_.enabled ? n_old(j) * scaling_.C0 : n_old(j);
            const Real eta0 = reducedFermiEnergyFromDensity(n0, Nc_[edge.n0]);
            const Real eta1 = reducedFermiEnergyFromDensity(n1, Nc_[edge.n1]);
            const Real factor = sgGeneralizedEinsteinFactor(n0, n1, eta0, eta1);
            const Real driftPotential = (psi1 - psi0) + Vt_ * std::log(
                (ni_[edge.n1] / Nc_[edge.n1]) /
                (ni_[edge.n0] / Nc_[edge.n0]));
            weights = sgEdgeWeights(driftPotential, Vt_ * factor);
            coef *= factor;
        } else {
            const Real dpsi = psi(j) - psi(i);
            weights = scaling_.enabled
                ? sgEdgeWeights(dpsi, 1.0)
                : sgEdgeWeights(dpsi, Vt_);
        }

        // Electron continuity flux from i to j:
        //   F_nij = coef * (B(-u) * n_i - B(+u) * n_j)
        triplets.emplace_back(i, i,  coef * weights.b_minus);
        triplets.emplace_back(j, j,  coef * weights.b_plus);
        triplets.emplace_back(i, j, -coef * weights.b_plus);
        triplets.emplace_back(j, i, -coef * weights.b_minus);
    }

    A_.setFromTriplets(triplets.begin(), triplets.end());

    const VectorXd psi_si = psiForMobility;
    const std::vector<Real> nodeElectricFields =
        impactIonizationCoupled_
        ? detail::computeNodeElectricFields(psi_si, mesh_, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> bandToBandSourceIntegrals =
        recombination_.bandToBandEnabled()
        ? detail::bandToBandGenerationNodeSourceIntegrals(
              recombination_.bandToBand(), scaling_.unitSystem, mesh_,
              cellMaterials, psi_si, fieldFactor)
        : std::vector<Real>{};
    VectorXd n_physical = n_old;
    VectorXd p_physical = p_old;
    if (scaling_.enabled) {
        n_physical *= scaling_.C0;
        p_physical *= scaling_.C0;
    }
    const bool qfImpact =
        detail::usesQuasiFermiAvalancheDrivingForce(impactIonizationConfig_);
    const VectorXd phipForImpact =
        holeQuasiFermiFromDensity(
            psi_si, p_old, ni_, Nv_, Vt_, scaling_, carrierStatistics_);
    const std::vector<Real> nodeElectronDrivingFields = (impactIonizationCoupled_ && qfImpact)
        ? detail::computeElectronAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig_, mesh_, nodeCells_, psi_si, phinForMobility,
            n_physical, ni_, Vt_, fieldFactor)
        : nodeElectricFields;
    const std::vector<Real> nodeHoleDrivingFields = (impactIonizationCoupled_ && qfImpact)
        ? detail::computeHoleAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig_, mesh_, nodeCells_, psi_si, phipForImpact,
            p_physical, ni_, Vt_, fieldFactor)
        : nodeElectricFields;
    const bool sgCurrentAvalanche = impactIonizationCoupled_ &&
        detail::usesEdgeCurrentAvalancheSource(impactIonizationConfig_);
    const std::vector<Real> sgAvalancheSourceIntegrals = sgCurrentAvalanche
        ? detail::currentDensityAvalancheSourceIntegrals(
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
            Vt_,
            fieldFactor,
            carrierStatistics_,
            Nc_,
            Nv_)
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

        Real recombinationNi = ni_i;
        if (usesFermiDirac(carrierStatistics_)) {
            const Real product = equilibriumCarrierProduct(
                n_v, p_v, ni_i, Nc_[i], Nv_[i], Vt_, carrierStatistics_);
            recombinationNi = std::sqrt(std::max<Real>(product, 0.0));
        }
        const RecombinationLinearization linearization =
            recombination_.electronLinearization(
                n_v, p_v, recombinationNi,
                recombination_.srhDopingConcentration(
                    doping_.donors(i), doping_.acceptors(i)));
        if (scaling_.enabled) {
            A_.coeffRef(ii, ii) +=
                linearization.diagonal * vol_i * sourceIntegralFactor / scaling_.D0;
            b_(ii) += linearization.rhs * vol_i * sourceIntegralFactor /
                (scaling_.C0 * scaling_.D0);
        } else {
            A_.coeffRef(ii, ii) += linearization.diagonal * vol_i;
            b_(ii) += linearization.rhs * vol_i;
        }
        if (recombination_.bandToBandEnabled()) {
            const Real source =
                bandToBandSourceIntegrals[i] * sourceIntegralFactor;
            b_(ii) += scaling_.enabled
                ? source / (scaling_.C0 * scaling_.D0)
                : source;
        }
        if (impactIonizationCoupled_ && sgCurrentAvalanche) {
            const Real gen = sgAvalancheSourceIntegrals[i] * sourceIntegralFactor;
            b_(ii) += scaling_.enabled ? (gen / (scaling_.C0 * scaling_.D0)) : gen;
        } else if (impactIonizationCoupled_) {
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
            b_(ii) += scaling_.enabled ? (gen * sourceIntegralFactor / (scaling_.C0 * scaling_.D0)) : gen;
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
        holeQuasiFermiFromDensity(
            psiForMobility, p_old, ni_, Nv_, Vt_, scaling_, carrierStatistics_);
    if (mobilityConfig_.model == "ialmob") {
        const VectorXd phinForMobility = electronQuasiFermiFromDensity(
            psiForMobility,n_old,ni_,Nc_,Vt_,scaling_,carrierStatistics_);
        const Real densityScale=scaling_.enabled?scaling_.C0:1.;
        updateIalTransportState(mobilityConfig_,mesh_,doping_,psiForMobility,
            n_old*densityScale,p_old*densityScale,phinForMobility,phipForMobility);
    }
    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";
    const Real fieldFactor = scaling_.enabled ? scaling_.fieldFromCoordinateDeltaFactor : 1.0;
    const bool vectorQfMobility = qfMobility &&
        mobilityConfig_.highFieldGradientDiscretization == "transport_cell_vector";
    const std::vector<Real> holeVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials, phipForMobility, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> contactElectricMobilityFields =
        mobilityConfig_.contactElectricFieldFallback
        ? detail::contactCellElectricFieldEdgeMagnitudesForMobility(
              mobilityConfig_, mesh_, edgeCells_, cellMaterials,
              psiForMobility, fieldFactor)
        : std::vector<Real>{};
    const Real sourceIntegralFactor = scaling_.enabled
        ? scaling_.unitSystem.continuitySourceIntegralFactor()
        : 1.0;

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
        const Real electricField = std::abs((psi1 - psi0) / h) * fieldFactor;
        const Real holeMobilityField = detail::mobilityHighFieldDrivingField(
            mobilityConfig_, e,
            vectorQfMobility ? holeVectorMobilityFields[e]
                : std::abs((phipForMobility(static_cast<int>(edge.n1)) -
                            phipForMobility(static_cast<int>(edge.n0))) / h) * fieldFactor,
            electricField, contactElectricMobilityFields);
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials, e, CarrierType::Hole,
            holeMobilityField,
            &mobilityConfig_,
            &psiForMobility);
        if (mup <= 0.0) continue; // skip insulator edges

        Real coef = scaling_.enabled
            ? (mup / scaling_.mu0) * fieldFactor * couple_[e] / h
            : mup * Vt_ * fieldFactor * couple_[e] / h;

        auto i = static_cast<int>(edge.n0);
        auto j = static_cast<int>(edge.n1);

        SGEdgeWeights weights;
        if (usesFermiDirac(carrierStatistics_)) {
            const Real p0 = scaling_.enabled ? p_old(i) * scaling_.C0 : p_old(i);
            const Real p1 = scaling_.enabled ? p_old(j) * scaling_.C0 : p_old(j);
            const Real eta0 = reducedFermiEnergyFromDensity(p0, Nv_[edge.n0]);
            const Real eta1 = reducedFermiEnergyFromDensity(p1, Nv_[edge.n1]);
            const Real factor = sgGeneralizedEinsteinFactor(p0, p1, eta0, eta1);
            const Real driftPotential = (psi1 - psi0) + Vt_ * std::log(
                (ni_[edge.n0] / Nv_[edge.n0]) /
                (ni_[edge.n1] / Nv_[edge.n1]));
            weights = sgEdgeWeights(driftPotential, Vt_ * factor);
            coef *= factor;
        } else {
            const Real dpsi = psi(j) - psi(i);
            weights = scaling_.enabled
                ? sgEdgeWeights(dpsi, 1.0)
                : sgEdgeWeights(dpsi, Vt_);
        }

        // Hole continuity flux from i to j:
        //   F_pij = coef * (B(+u) * p_i - B(-u) * p_j)
        triplets.emplace_back(i, i,  coef * weights.b_plus);
        triplets.emplace_back(j, j,  coef * weights.b_minus);
        triplets.emplace_back(i, j, -coef * weights.b_minus);
        triplets.emplace_back(j, i, -coef * weights.b_plus);
    }

    A_.setFromTriplets(triplets.begin(), triplets.end());

    const VectorXd psi_si = psiForMobility;
    const std::vector<Real> nodeElectricFields =
        impactIonizationCoupled_
        ? detail::computeNodeElectricFields(psi_si, mesh_, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> bandToBandSourceIntegrals =
        recombination_.bandToBandEnabled()
        ? detail::bandToBandGenerationNodeSourceIntegrals(
              recombination_.bandToBand(), scaling_.unitSystem, mesh_,
              cellMaterials, psi_si, fieldFactor)
        : std::vector<Real>{};
    VectorXd n_physical = n_old;
    VectorXd p_physical = p_old;
    if (scaling_.enabled) {
        n_physical *= scaling_.C0;
        p_physical *= scaling_.C0;
    }
    const bool qfImpact =
        detail::usesQuasiFermiAvalancheDrivingForce(impactIonizationConfig_);
    const VectorXd phinForImpact =
        electronQuasiFermiFromDensity(
            psi_si, n_old, ni_, Nc_, Vt_, scaling_, carrierStatistics_);
    const std::vector<Real> nodeElectronDrivingFields = (impactIonizationCoupled_ && qfImpact)
        ? detail::computeElectronAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig_, mesh_, nodeCells_, psi_si, phinForImpact,
            n_physical, ni_, Vt_, fieldFactor)
        : nodeElectricFields;
    const std::vector<Real> nodeHoleDrivingFields = (impactIonizationCoupled_ && qfImpact)
        ? detail::computeHoleAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig_, mesh_, nodeCells_, psi_si, phipForMobility,
            p_physical, ni_, Vt_, fieldFactor)
        : nodeElectricFields;
    const bool sgCurrentAvalanche = impactIonizationCoupled_ &&
        detail::usesEdgeCurrentAvalancheSource(impactIonizationConfig_);
    const std::vector<Real> sgAvalancheSourceIntegrals = sgCurrentAvalanche
        ? detail::currentDensityAvalancheSourceIntegrals(
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
            Vt_,
            fieldFactor,
            carrierStatistics_,
            Nc_,
            Nv_)
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

        Real recombinationNi = ni_i;
        if (usesFermiDirac(carrierStatistics_)) {
            const Real product = equilibriumCarrierProduct(
                n_v, p_v, ni_i, Nc_[i], Nv_[i], Vt_, carrierStatistics_);
            recombinationNi = std::sqrt(std::max<Real>(product, 0.0));
        }
        const RecombinationLinearization linearization =
            recombination_.holeLinearization(
                n_v, p_v, recombinationNi,
                recombination_.srhDopingConcentration(
                    doping_.donors(i), doping_.acceptors(i)));
        if (scaling_.enabled) {
            A_.coeffRef(ii, ii) +=
                linearization.diagonal * vol_i * sourceIntegralFactor / scaling_.D0;
            b_(ii) += linearization.rhs * vol_i * sourceIntegralFactor /
                (scaling_.C0 * scaling_.D0);
        } else {
            A_.coeffRef(ii, ii) += linearization.diagonal * vol_i;
            b_(ii) += linearization.rhs * vol_i;
        }
        if (recombination_.bandToBandEnabled()) {
            const Real source =
                bandToBandSourceIntegrals[i] * sourceIntegralFactor;
            b_(ii) += scaling_.enabled
                ? source / (scaling_.C0 * scaling_.D0)
                : source;
        }
        if (impactIonizationCoupled_ && sgCurrentAvalanche) {
            const Real gen = sgAvalancheSourceIntegrals[i] * sourceIntegralFactor;
            b_(ii) += scaling_.enabled ? (gen / (scaling_.C0 * scaling_.D0)) : gen;
        } else if (impactIonizationCoupled_) {
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
            b_(ii) += scaling_.enabled ? (gen * sourceIntegralFactor / (scaling_.C0 * scaling_.D0)) : gen;
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
