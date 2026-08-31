#include "vela/post/ContactCurrent.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/Bernoulli.h"
#include "vela/discretization/ElementQfGradient.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/equation/AssemblerUtils.h"
#include <unordered_set>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace vela {
namespace {

Real validatedThermalVoltage(Real temperature_K)
{
    if (temperature_K <= 0.0)
        throw std::invalid_argument("ContactCurrent: temperature_K must be positive.");
    return constants::kb * temperature_K / constants::q;
}

struct NeumaierSum {
    Real sum = 0.0;
    Real correction = 0.0;

    void add(Real value)
    {
        const Real next = sum + value;
        if (std::abs(sum) >= std::abs(value))
            correction += (sum - next) + value;
        else
            correction += (value - next) + sum;
        sum = next;
    }

    Real value() const { return sum + correction; }
};

} // namespace

ContactCurrent::ContactCurrent(const DeviceMesh& mesh,
                               const MaterialDatabase& matdb,
                               const DopingModel& doping,
                               MobilityModelConfig mobilityConfig,
                               Real temperature_K,
                               DDScalingSpec scaling,
                               BandgapNarrowingConfig bandgapNarrowingConfig,
                               CarrierStatisticsConfig carrierStatistics,
                               DensityGradientQuantumPotentialConfig
                                   electronQuantumPotential)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , edgeCells_(detail::buildEdgeCellMap(mesh))
    , couple_(detail::computeTransportEdgeCouplings(mesh))
    , mobilityConfig_(mobilityConfig)
    , mobility_(makeMobilityModel(mobilityConfig))
    , thermalVoltage_(validatedThermalVoltage(temperature_K))
    , scaling_(scaling)
    , ni_(detail::buildValidatedEffectiveNodeNi(
          "ContactCurrent",
          mesh,
          matdb,
          doping,
          bandgapNarrowingConfig,
          validatedThermalVoltage(temperature_K)))
    , Nc_(detail::buildNodeDensityOfStates(mesh, matdb, temperature_K, true))
    , Nv_(detail::buildNodeDensityOfStates(mesh, matdb, temperature_K, false))
    , carrierStatistics_(std::move(carrierStatistics))
    , electronQuantumPotentialConfig_(std::move(electronQuantumPotential))
{}


ContactCurrentResult ContactCurrent::compute(const DDSolution& solution,
                                             const std::string& contactName) const
{
    return computeDetailed(solution, contactName).totals;
}

ContactCurrentResult ContactCurrent::compute(
    const DDSolution& solution,
    const std::string& contactName,
    const ContactCurrentEdgeOverrides& overrides) const
{
    return computeDetailed(solution, contactName, overrides).totals;
}

ContactCurrentDetailedResult ContactCurrent::computeDetailed(
    const DDSolution& solution,
    const std::string& contactName) const
{
    static const ContactCurrentEdgeOverrides noOverrides;
    return computeDetailed(solution, contactName, noOverrides);
}

ContactCurrentResult ContactCurrent::computeFromResidual(
    const CoupledDDAssembler& assembler,
    const VectorXd& x,
    const std::string& contactName) const
{
    const Contact* contact = nullptr;
    for (const Contact& candidate : mesh_.contacts()) {
        if (candidate.name == contactName) {
            contact = &candidate;
            break;
        }
    }
    if (contact == nullptr)
        throw std::invalid_argument("ContactCurrent: unknown contact '" + contactName + "'.");

    const int N = static_cast<int>(assembler.numNodes());
    const int phinOffset = N;
    const int phipOffset = 2 * N;
    const VectorXd r = assembler.residual(x, CoupledDDBoundaryConditions{});
    const Real continuityScale = assembler.continuityResidualScale();
    const Real currentLineFactor = scaling_.enabled
        ? scaling_.currentDensityLineIntegralFactor
        : 1.0;

    ContactCurrentResult result;
    for (Index node : contact->node_ids) {
        const int row = static_cast<int>(node);
        result.electronCurrent -= constants::q * r(phinOffset + row) * continuityScale * currentLineFactor;
        result.holeCurrent -= constants::q * r(phipOffset + row) * continuityScale * currentLineFactor;
    }
    result.totalCurrent = result.electronCurrent - result.holeCurrent;
    return result;
}
ContactCurrentDetailedResult ContactCurrent::computeDetailed(
    const DDSolution& solution,
    const std::string& contactName,
    const ContactCurrentEdgeOverrides& overrides) const
{
    const Contact* contact = nullptr;
    for (const Contact& candidate : mesh_.contacts()) {
        if (candidate.name == contactName) {
            contact = &candidate;
            break;
        }
    }
    if (contact == nullptr)
        throw std::invalid_argument("ContactCurrent: unknown contact '" + contactName + "'.");

    std::unordered_set<Index> contactNodes(contact->node_ids.begin(), contact->node_ids.end());
    const Real temperature_K = thermalVoltage_ * constants::q / constants::kb;
    const std::vector<Material> cellMaterials =
        detail::buildCellMaterials(mesh_, matdb_, temperature_K);
    const Real fieldFactor = scaling_.enabled
        ? scaling_.fieldFromCoordinateDeltaFactor : 1.0;
    const bool hasReferencedElectronQf =
        solution.phinIncrement.size() == solution.phin.size();
    const bool hasReferencedHoleQf =
        solution.phipIncrement.size() == solution.phip.size();
    VectorXd electronQf = hasReferencedElectronQf
        ? solution.phinIncrement : solution.phin;
    VectorXd holeQf = hasReferencedHoleQf
        ? solution.phipIncrement : solution.phip;
    if (hasReferencedElectronQf) {
        for (int i = 0; i < electronQf.size(); ++i) {
            electronQf(i) = static_cast<Real>(
                static_cast<long double>(
                    solution.electronQuasiFermiReferenceAt(i)) +
                static_cast<long double>(solution.phinIncrement(i)));
        }
    }
    if (hasReferencedHoleQf) {
        for (int i = 0; i < holeQf.size(); ++i) {
            holeQf(i) = static_cast<Real>(
                static_cast<long double>(
                    solution.holeQuasiFermiReferenceAt(i)) +
                static_cast<long double>(solution.phipIncrement(i)));
        }
    }
    const bool vectorQfMobility =
        mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient" &&
        mobilityConfig_.highFieldGradientDiscretization == "transport_cell_vector";
    const std::vector<Real> electronVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials, electronQf, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> holeVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials, holeQf, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> contactElectricMobilityFields =
        mobilityConfig_.contactElectricFieldFallback
        ? detail::contactCellElectricFieldEdgeMagnitudesForMobility(
              mobilityConfig_, mesh_, edgeCells_, cellMaterials,
              solution.psi, fieldFactor)
        : std::vector<Real>{};

    ContactCurrentDetailedResult detailed;
    NeumaierSum electronCurrentCompensated;
    NeumaierSum holeCurrentCompensated;
    long double electronCurrentLongDouble = 0.0L;
    long double holeCurrentLongDouble = 0.0L;
    if (mobilityConfig_.carrierCurrentDiscretization ==
        "element_qf_gradient") {
        if (contact->edge_node_ids.empty()) {
            throw std::invalid_argument(
                "ContactCurrent: element_qf_gradient requires exact "
                "contact edge_node_ids.");
        }
        const Real currentLineFactor = scaling_.enabled
            ? scaling_.currentDensityLineIntegralFactor : 1.0;
        for (const auto& boundaryNodes : contact->edge_node_ids) {
            Index edgeId = mesh_.numEdges();
            for (Index candidate = 0; candidate < mesh_.numEdges(); ++candidate) {
                const Edge& edge = mesh_.getEdge(candidate);
                if ((edge.n0 == boundaryNodes[0] && edge.n1 == boundaryNodes[1]) ||
                    (edge.n0 == boundaryNodes[1] && edge.n1 == boundaryNodes[0])) {
                    edgeId = candidate;
                    break;
                }
            }
            if (edgeId == mesh_.numEdges())
                throw std::invalid_argument(
                    "ContactCurrent: exact contact edge is absent from mesh topology.");

            Index transportCell = mesh_.numCells();
            for (const Index cellId : edgeCells_[edgeId]) {
                if (detail::isTransportMaterial(cellMaterials[cellId])) {
                    transportCell = cellId;
                    break;
                }
            }
            if (transportCell == mesh_.numCells())
                continue;
            const Cell& cell = mesh_.getCell(transportCell);
            if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
                throw std::invalid_argument(
                    "ContactCurrent: element_qf_gradient supports Tri3 transport cells only.");

            bool valid = false;
            Real area = 0.0;
            const Point2 psiGradient = detail::cellScalarGradient(
                mesh_, cell,
                [&](Index node) { return solution.psi(static_cast<int>(node)); },
                valid, area);
            const Point2 electronGradient = detail::cellScalarGradient(
                mesh_, cell,
                [&](Index node) { return electronQf(static_cast<int>(node)); },
                valid, area);
            const Point2 holeGradient = detail::cellScalarGradient(
                mesh_, cell,
                [&](Index node) { return holeQf(static_cast<int>(node)); },
                valid, area);
            if (!valid)
                continue;

            Real averageElectronDensity = 0.0;
            Real averageHoleDensity = 0.0;
            Real mobilityDoping = 0.0;
            for (const Index node : cell.node_ids) {
                averageElectronDensity += solution.n(static_cast<int>(node));
                averageHoleDensity += solution.p(static_cast<int>(node));
                mobilityDoping +=
                    mobilityConfig_.dopingConcentrationBasis == "total_impurity"
                    ? doping_.totalImpurity(node) : doping_.netDoping(node);
            }
            averageElectronDensity /= 3.0;
            averageHoleDensity /= 3.0;
            mobilityDoping /= 3.0;
            if (mobilityConfig_.dopingConcentrationBasis ==
                "cell_reconstructed_total_impurity") {
                mobilityDoping = detail::cellAverageTotalImpurity(
                    mesh_, doping_, transportCell);
            }
            const Real electricField = psiGradient.norm() * fieldFactor;
            const bool qfDrive =
                mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";
            const Real mun = mobility_->electronMobility(
                cellMaterials[transportCell], mobilityDoping,
                averageElectronDensity, averageHoleDensity,
                qfDrive ? electronGradient.norm() * fieldFactor : electricField);
            const Real mup = mobility_->holeMobility(
                cellMaterials[transportCell], mobilityDoping,
                averageElectronDensity, averageHoleDensity,
                qfDrive ? holeGradient.norm() * fieldFactor : electricField);
            const auto current = elementQfGradientCellCurrent(
                mesh_, cell,
                [&](Index node) { return electronQf(static_cast<int>(node)); },
                [&](Index node) { return holeQf(static_cast<int>(node)); },
                [&](Index node) { return solution.n(static_cast<int>(node)); },
                [&](Index node) { return solution.p(static_cast<int>(node)); },
                mun, mup, fieldFactor);
            if (!current.valid)
                continue;

            const Node& node0 = mesh_.getNode(boundaryNodes[0]);
            const Node& node1 = mesh_.getNode(boundaryNodes[1]);
            const Point2 tangent{node1.x - node0.x, node1.y - node0.y};
            const Real length = tangent.norm();
            if (length <= 1.0e-30)
                continue;
            Point2 outward{tangent.y() / length, -tangent.x() / length};
            Point2 centroid = Point2::Zero();
            for (const Index node : cell.node_ids) {
                const Node& point = mesh_.getNode(node);
                centroid += Point2{point.x, point.y};
            }
            centroid /= 3.0;
            const Point2 midpoint{
                0.5 * (node0.x + node1.x),
                0.5 * (node0.y + node1.y)};
            if (outward.dot(centroid - midpoint) > 0.0)
                outward = -outward;

            const Real electronCurrent = constants::q *
                current.electronParticleCurrent.dot(outward) * length *
                currentLineFactor;
            const Real holeCurrent = constants::q *
                current.holeParticleCurrent.dot(outward) * length *
                currentLineFactor;
            detailed.totals.electronCurrent += electronCurrent;
            detailed.totals.holeCurrent += holeCurrent;
            electronCurrentCompensated.add(electronCurrent);
            holeCurrentCompensated.add(holeCurrent);
            electronCurrentLongDouble += electronCurrent;
            holeCurrentLongDouble += holeCurrent;

            ContactCurrentEdgeDiagnostic edgeDiag;
            edgeDiag.edgeId = edgeId;
            edgeDiag.node0 = boundaryNodes[0];
            edgeDiag.node1 = boundaryNodes[1];
            edgeDiag.edgeLength_m =
                scaling_.unitSystem.internalLengthToMeters(length);
            edgeDiag.edgeCouple_m = edgeDiag.edgeLength_m;
            edgeDiag.electronUsedQuasiFermi = true;
            edgeDiag.holeUsedQuasiFermi = true;
            edgeDiag.mun = mun;
            edgeDiag.mup = mup;
            edgeDiag.electronCurrent = electronCurrent;
            edgeDiag.holeCurrent = holeCurrent;
            edgeDiag.totalCurrent = electronCurrent - holeCurrent;
            detailed.edges.push_back(std::move(edgeDiag));
        }
        detailed.totals.totalCurrent =
            detailed.totals.electronCurrent - detailed.totals.holeCurrent;
        detailed.precision.electronCurrentCompensated =
            electronCurrentCompensated.value();
        detailed.precision.holeCurrentCompensated =
            holeCurrentCompensated.value();
        detailed.precision.totalCurrentCompensated =
            detailed.precision.electronCurrentCompensated -
            detailed.precision.holeCurrentCompensated;
        detailed.precision.electronCurrentLongDoubleReference =
            static_cast<Real>(electronCurrentLongDouble);
        detailed.precision.holeCurrentLongDoubleReference =
            static_cast<Real>(holeCurrentLongDouble);
        detailed.precision.totalCurrentLongDoubleReference =
            static_cast<Real>(electronCurrentLongDouble - holeCurrentLongDouble);
        return detailed;
    }
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const bool n0OnContact = contactNodes.count(edge.n0) > 0;
        const bool n1OnContact = contactNodes.count(edge.n1) > 0;
        if (n0OnContact == n1OnContact)
            continue;
        if (edge.length < 1.0e-30 || couple_[e] <= 0.0)
            continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);

        // The solver returns physical DDSolution fields in the active internal unit system.
        const Real psi_i = solution.psi(i);
        const Real psi_j = solution.psi(j);
        const Real n_i = solution.n(i);
        const Real n_j = solution.n(j);
        const Real p_i = solution.p(i);
        const Real p_j = solution.p(j);
        const Real dpsi = psi_j - psi_i;
        const Real edgeLength = edge.length;

        const Real currentLineFactor = scaling_.enabled
            ? scaling_.currentDensityLineIntegralFactor
            : 1.0;
        const Real electricField = std::abs(dpsi / edgeLength) * fieldFactor;
        const Real phin_i = hasReferencedElectronQf
            ? solution.phinIncrement(i) : solution.phin(i);
        const Real phin_j = hasReferencedElectronQf
            ? solution.phinIncrement(j) : solution.phin(j);
        const Real phip_i = hasReferencedHoleQf
            ? solution.phipIncrement(i) : solution.phip(i);
        const Real phip_j = hasReferencedHoleQf
            ? solution.phipIncrement(j) : solution.phip(j);
        const Real electronReference_i = hasReferencedElectronQf
            ? solution.electronQuasiFermiReferenceAt(i) : 0.0;
        const Real electronReference_j = hasReferencedElectronQf
            ? solution.electronQuasiFermiReferenceAt(j) : 0.0;
        const Real holeReference_i = hasReferencedHoleQf
            ? solution.holeQuasiFermiReferenceAt(i) : 0.0;
        const Real holeReference_j = hasReferencedHoleQf
            ? solution.holeQuasiFermiReferenceAt(j) : 0.0;
        const long double electronQfDropLong = hasReferencedElectronQf
            ? (static_cast<long double>(electronReference_j) -
               static_cast<long double>(electronReference_i)) +
              (static_cast<long double>(phin_j) -
               static_cast<long double>(phin_i))
            : static_cast<long double>(solution.phin(j)) -
              static_cast<long double>(solution.phin(i));
        const long double holeQfDropLong = hasReferencedHoleQf
            ? (static_cast<long double>(holeReference_j) -
               static_cast<long double>(holeReference_i)) +
              (static_cast<long double>(phip_j) -
               static_cast<long double>(phip_i))
            : static_cast<long double>(solution.phip(j)) -
              static_cast<long double>(solution.phip(i));
        const Real electronQfDrop = static_cast<Real>(electronQfDropLong);
        const Real holeQfDrop = static_cast<Real>(holeQfDropLong);
        const bool sentaurusExponential =
            electronQuantumPotentialConfig_.enabled &&
            usesFermiDirac(carrierStatistics_) &&
            electronQuantumPotentialConfig_.transportCoupling ==
                "sentaurus_exponential";
        const Real exponentialWeight = sentaurusExponential
            ? electronQuantumPotentialConfig_.transportCouplingWeight : 0.0;
        const Real quantum_i =
            solution.electronQuantumPotential.size() == solution.psi.size()
                ? solution.electronQuantumPotential(i) : 0.0;
        const Real quantum_j =
            solution.electronQuantumPotential.size() == solution.psi.size()
                ? solution.electronQuantumPotential(j) : 0.0;
        const Real electronPsi_i = psi_i - electronReference_i
            - (1.0 - exponentialWeight) * quantum_i;
        const Real electronPsi_j = psi_j - electronReference_j
            - (1.0 - exponentialWeight) * quantum_j;
        const Real holePsi_i = psi_i - holeReference_i;
        const Real holePsi_j = psi_j - holeReference_j;
        Real phip_i_forHole = phip_i;
        Real phip_j_forHole = phip_j;
        bool holeQfDropOverrideApplied = false;
        const auto holeDropIt = overrides.holeQuasiFermiDropByEdge.find(e);
        if (holeDropIt != overrides.holeQuasiFermiDropByEdge.end()
            && std::isfinite(holeDropIt->second)) {
            phip_j_forHole = phip_i_forHole + holeDropIt->second;
            holeQfDropOverrideApplied = true;
        }
        const long double holeQfDropForFluxLong = holeQfDropOverrideApplied
            ? static_cast<long double>(holeDropIt->second)
            : holeQfDropLong;
        const Real holeQfDropForFlux =
            static_cast<Real>(holeQfDropForFluxLong);
        const Real electronMobilityField = detail::mobilityHighFieldDrivingField(
            mobilityConfig_, e,
            vectorQfMobility ? electronVectorMobilityFields[e]
                             : std::abs(electronQfDrop / edgeLength) * fieldFactor,
            electricField, contactElectricMobilityFields);
        const Real holeMobilityField = detail::mobilityHighFieldDrivingField(
            mobilityConfig_, e,
            vectorQfMobility ? holeVectorMobilityFields[e]
                             : std::abs(holeQfDrop / edgeLength) * fieldFactor,
            electricField, contactElectricMobilityFields);

        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials, e, CarrierType::Electron,
            electronMobilityField,
            &mobilityConfig_,
            &solution.psi);
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials, e, CarrierType::Hole,
            holeMobilityField,
            &mobilityConfig_,
            &solution.psi);

        // SG fluxes in physical units.  Mirror CoupledDDAssembler residual:
        // use the cancellation-free quasi-Fermi balanced form, including the
        // variable-ni generalization needed for BGN/effective-ni edges.  The
        // density-based form B(-u)*n0 - B(+u)*n1 does not cancel flat
        // quasi-Fermi levels when ni varies across the edge.
        const Index idxI = edge.n0;
        const Index idxJ = edge.n1;
        const Real ni_i = ni_[idxI];
        const Real ni_j = ni_[idxJ];
        const bool fermiDirac = usesFermiDirac(carrierStatistics_);
        Real electronGeneralizedFactor = 1.0;
        Real electronGeneralizedArgument = dpsi / thermalVoltage_;
        Real holeGeneralizedFactor = 1.0;
        Real holeGeneralizedArgument = dpsi / thermalVoltage_;
        Real electronContinuityFlux01 = 0.0;
        long double electronContinuityFluxLongDouble01 = 0.0L;
        Real electronFlux01 = 0.0;
        if (mun > 0.0) {
            const Real coef = mun * thermalVoltage_ * fieldFactor / edgeLength;
            if (fermiDirac) {
                const Real etaI = (electronPsi_i - phin_i) / thermalVoltage_
                    + std::log(ni_i / Nc_[idxI]);
                const Real etaJ = (electronPsi_j - phin_j) / thermalVoltage_
                    + std::log(ni_j / Nc_[idxJ]);
                const Real driftPotential =
                    (psi_j - psi_i) -
                    (1.0 - exponentialWeight) * (quantum_j - quantum_i)
                    + thermalVoltage_ * std::log(
                        (ni_j / Nc_[idxJ]) / (ni_i / Nc_[idxI]));
                const Real fermiDensityI = Nc_[idxI] * fermiDiracHalf(etaI);
                const Real fermiDensityJ = Nc_[idxJ] * fermiDiracHalf(etaJ);
                electronGeneralizedFactor = sgGeneralizedEinsteinFactor(
                    fermiDensityI, fermiDensityJ, etaI, etaJ);
                electronGeneralizedArgument = driftPotential /
                    (thermalVoltage_ * electronGeneralizedFactor);
                electronContinuityFlux01 =
                    sgElectronFermiDiracQuantumContinuityFlux(
                    n_i, n_j, fermiDensityI, fermiDensityJ,
                    etaI, etaJ, driftPotential,
                    0.0, electronQfDrop, thermalVoltage_, coef);
                electronContinuityFluxLongDouble01 =
                    sgElectronFermiDiracQuantumContinuityFluxLongDouble(
                        n_j, fermiDensityI, fermiDensityJ,
                        etaI, etaJ, driftPotential,
                        electronQfDropLong, thermalVoltage_, coef);
            } else {
                const Real electronPsiFromNode0Qf = electronPsi_i - phin_i;
                const Real electronPsiFromNode0QfAtJ =
                    electronPsi_j - phin_j + electronQfDrop;
                electronContinuityFlux01 = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                    ni_i, ni_j,
                    electronPsiFromNode0Qf, electronPsiFromNode0QfAtJ,
                    0.0, electronQfDrop, thermalVoltage_, coef);
            }
            if (!fermiDirac)
                electronContinuityFluxLongDouble01 = electronContinuityFlux01;
            // sgElectronFlux = -sgElectronContinuityFlux by definition.
            electronFlux01 = -electronContinuityFlux01;
        }
        Real holeContinuityFlux01 = 0.0;
        long double holeContinuityFluxLongDouble01 = 0.0L;
        Real holeFlux01 = 0.0;
        if (mup > 0.0) {
            const Real coef = mup * thermalVoltage_ * fieldFactor / edgeLength;
            if (fermiDirac) {
                const Real etaI = (phip_i_forHole - holePsi_i) / thermalVoltage_
                    + std::log(ni_i / Nv_[idxI]);
                const Real etaJ = (phip_j_forHole - holePsi_j) / thermalVoltage_
                    + std::log(ni_j / Nv_[idxJ]);
                const Real driftPotential = psi_j - psi_i
                    + thermalVoltage_ * std::log(
                        (ni_i / Nv_[idxI]) / (ni_j / Nv_[idxJ]));
                holeGeneralizedFactor = sgGeneralizedEinsteinFactor(
                    p_i, p_j, etaI, etaJ);
                holeGeneralizedArgument = driftPotential /
                    (thermalVoltage_ * holeGeneralizedFactor);
                holeContinuityFlux01 = sgHoleFermiDiracContinuityFlux(
                    p_i, p_j, etaI, etaJ, driftPotential,
                    0.0, holeQfDropForFlux,
                    thermalVoltage_, coef);
                holeContinuityFluxLongDouble01 =
                    sgHoleFermiDiracContinuityFluxLongDouble(
                        p_i, p_j, etaI, etaJ, driftPotential,
                        holeQfDropForFluxLong, thermalVoltage_, coef);
            } else {
                const Real holePsiFromNode0Qf = holePsi_i - phip_i_forHole;
                const Real holePsiFromNode0QfAtJ =
                    holePsi_j - phip_j_forHole + holeQfDropForFlux;
                holeContinuityFlux01 = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                    ni_i, ni_j,
                    holePsiFromNode0Qf, holePsiFromNode0QfAtJ,
                    0.0, holeQfDropForFlux,
                    thermalVoltage_, coef);
            }
            if (!fermiDirac)
                holeContinuityFluxLongDouble01 = holeContinuityFlux01;
            holeFlux01 = -holeContinuityFlux01;
        }

        // Algebraic SG split: J = J_drift + J_diffusion.
        const SGEdgeWeights weights = sgEdgeWeights(dpsi, thermalVoltage_);
        const Real bAvg = 0.5 * (weights.b_plus + weights.b_minus);
        Real electronDriftFlux01 = 0.0;
        Real electronDiffusionFlux01 = 0.0;
        Real holeDriftFlux01 = 0.0;
        Real holeDiffusionFlux01 = 0.0;
        if (fermiDirac) {
            if (mun > 0.0) {
                const Real bMinus = bernoulli(-electronGeneralizedArgument);
                const Real bPlus = bernoulli(electronGeneralizedArgument);
                const Real coefficient = mun * thermalVoltage_ * fieldFactor /
                    edgeLength * electronGeneralizedFactor;
                electronDriftFlux01 = -coefficient * 0.5 * (bMinus - bPlus)
                    * (n_i + n_j);
                electronDiffusionFlux01 = -coefficient * 0.5 * (bMinus + bPlus)
                    * (n_i - n_j);
            }
            if (mup > 0.0) {
                const Real bPlus = bernoulli(holeGeneralizedArgument);
                const Real bMinus = bernoulli(-holeGeneralizedArgument);
                const Real coefficient = mup * thermalVoltage_ * fieldFactor /
                    edgeLength * holeGeneralizedFactor;
                holeDriftFlux01 = -coefficient * 0.5 * (bPlus - bMinus)
                    * (p_i + p_j);
                holeDiffusionFlux01 = -coefficient * 0.5 * (bPlus + bMinus)
                    * (p_i - p_j);
            }
        } else {
            electronDriftFlux01 = (mun > 0.0)
                ? mun * (dpsi / edgeLength) * fieldFactor * (0.5 * (n_i + n_j))
                : 0.0;
            electronDiffusionFlux01 = (mun > 0.0)
                ? mun * (thermalVoltage_ / edgeLength) * fieldFactor * bAvg * (n_i - n_j)
                : 0.0;
            holeDriftFlux01 = (mup > 0.0)
                ? mup * (dpsi / edgeLength) * fieldFactor * (0.5 * (p_i + p_j))
                : 0.0;
            holeDiffusionFlux01 = (mup > 0.0)
                ? mup * (thermalVoltage_ / edgeLength) * fieldFactor * bAvg * (p_j - p_i)
                : 0.0;
        }

        const Real outwardSign = n0OnContact ? 1.0 : -1.0;
        // Current density in the active internal unit system times edge length gives current per internal device depth.
        const Real electronCurrent = constants::q * outwardSign * electronFlux01 * couple_[e] * currentLineFactor;
        const Real electronDriftCurrent = constants::q * outwardSign * electronDriftFlux01 * couple_[e] * currentLineFactor;
        const Real electronDiffusionCurrent = constants::q * outwardSign * electronDiffusionFlux01 * couple_[e] * currentLineFactor;
        const Real holeCurrent = constants::q * outwardSign * holeFlux01 * couple_[e] * currentLineFactor;
        const Real holeDriftCurrent = constants::q * outwardSign * holeDriftFlux01 * couple_[e] * currentLineFactor;
        const Real holeDiffusionCurrent = constants::q * outwardSign * holeDiffusionFlux01 * couple_[e] * currentLineFactor;
        const long double electronCurrentLongDoubleEdge =
            -static_cast<long double>(constants::q) *
            static_cast<long double>(outwardSign) *
            electronContinuityFluxLongDouble01 *
            static_cast<long double>(couple_[e]) *
            static_cast<long double>(currentLineFactor);
        const long double holeCurrentLongDoubleEdge =
            -static_cast<long double>(constants::q) *
            static_cast<long double>(outwardSign) *
            holeContinuityFluxLongDouble01 *
            static_cast<long double>(couple_[e]) *
            static_cast<long double>(currentLineFactor);

        detailed.totals.electronCurrent += electronCurrent;
        detailed.totals.electronDriftCurrent += electronDriftCurrent;
        detailed.totals.electronDiffusionCurrent += electronDiffusionCurrent;
        detailed.totals.holeCurrent += holeCurrent;
        detailed.totals.holeDriftCurrent += holeDriftCurrent;
        detailed.totals.holeDiffusionCurrent += holeDiffusionCurrent;
        electronCurrentCompensated.add(electronCurrent);
        holeCurrentCompensated.add(holeCurrent);
        electronCurrentLongDouble += electronCurrentLongDoubleEdge;
        holeCurrentLongDouble += holeCurrentLongDoubleEdge;

        ContactCurrentEdgeDiagnostic edgeDiag;
        edgeDiag.edgeId = e;
        edgeDiag.node0 = edge.n0;
        edgeDiag.node1 = edge.n1;
        edgeDiag.edgeLength_m = scaling_.unitSystem.internalLengthToMeters(edge.length);
        edgeDiag.edgeCouple_m = scaling_.unitSystem.internalLengthToMeters(couple_[e]);
        edgeDiag.outwardSign = outwardSign;
        edgeDiag.bernoulliU = dpsi / thermalVoltage_;
        edgeDiag.bernoulliBplus = weights.b_plus;
        edgeDiag.bernoulliBminus = weights.b_minus;
        edgeDiag.electronUsedQuasiFermi = true;
        edgeDiag.holeUsedQuasiFermi = true;
        edgeDiag.psi0 = psi_i;
        edgeDiag.psi1 = psi_j;
        edgeDiag.phin0 = phin_i + electronReference_i;
        edgeDiag.phin1 = phin_j + electronReference_j;
        edgeDiag.phip0 = phip_i_forHole + holeReference_i;
        edgeDiag.phip1 = phip_j_forHole + holeReference_j;
        edgeDiag.holeQfDropOverrideApplied = holeQfDropOverrideApplied;
        edgeDiag.n0 = n_i;
        edgeDiag.n1 = n_j;
        edgeDiag.p0 = p_i;
        edgeDiag.p1 = p_j;
        edgeDiag.ni0 = ni_i;
        edgeDiag.ni1 = ni_j;
        edgeDiag.mun = mun;
        edgeDiag.mup = mup;
        edgeDiag.electronContinuityFlux = electronContinuityFlux01;
        edgeDiag.holeContinuityFlux = holeContinuityFlux01;
        edgeDiag.electronCurrent = electronCurrent;
        edgeDiag.electronCurrentLongDoubleReference =
            static_cast<Real>(electronCurrentLongDoubleEdge);
        edgeDiag.electronDriftCurrent = electronDriftCurrent;
        edgeDiag.electronDiffusionCurrent = electronDiffusionCurrent;
        edgeDiag.holeCurrent = holeCurrent;
        edgeDiag.holeCurrentLongDoubleReference =
            static_cast<Real>(holeCurrentLongDoubleEdge);
        edgeDiag.holeDriftCurrent = holeDriftCurrent;
        edgeDiag.holeDiffusionCurrent = holeDiffusionCurrent;
        edgeDiag.totalCurrent = electronCurrent - holeCurrent;
        detailed.edges.push_back(std::move(edgeDiag));
    }

    // Sign convention: electronCurrent and holeCurrent accumulate
    //   q * (carrier-particle inflow into the contact from the device) * couple.
    // With the electron carrier charge being -q, the contribution of electrons
    // to the conventional current supplied into the contact from the external
    // circuit is -(particle inflow), so the total terminal current is
    //   I_total = I_electron - I_hole.
    // Using `I_electron + I_hole` (as previously) double-adds the volume
    // recombination integral into the terminal current and breaks the
    // Kirchhoff balance |I_anode| = |I_cathode| for a two-terminal device.
    detailed.totals.totalCurrent = detailed.totals.electronCurrent - detailed.totals.holeCurrent;
    detailed.precision.electronCurrentCompensated =
        electronCurrentCompensated.value();
    detailed.precision.holeCurrentCompensated = holeCurrentCompensated.value();
    detailed.precision.totalCurrentCompensated =
        detailed.precision.electronCurrentCompensated -
        detailed.precision.holeCurrentCompensated;
    detailed.precision.electronCurrentLongDoubleReference =
        static_cast<Real>(electronCurrentLongDouble);
    detailed.precision.holeCurrentLongDoubleReference =
        static_cast<Real>(holeCurrentLongDouble);
    detailed.precision.totalCurrentLongDoubleReference =
        static_cast<Real>(electronCurrentLongDouble - holeCurrentLongDouble);
    return detailed;
}


ContactCurrentResult ContactCurrent::compute(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const DDSolution& solution,
    const std::string& contactName,
    const MobilityModelConfig& mobilityConfig,
    Real temperature_K,
    DDScalingSpec scaling,
    const BandgapNarrowingConfig& bandgapNarrowingConfig,
    const CarrierStatisticsConfig& carrierStatistics,
    const DensityGradientQuantumPotentialConfig& electronQuantumPotential)
{
    return ContactCurrent(mesh, matdb, doping, mobilityConfig, temperature_K, scaling,
                          bandgapNarrowingConfig, carrierStatistics,
                          electronQuantumPotential)
        .compute(solution, contactName);
}

} // namespace vela
