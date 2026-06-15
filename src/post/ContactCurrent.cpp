#include "vela/post/ContactCurrent.h"
#include "vela/core/PhysicalConstants.h"
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

} // namespace

ContactCurrent::ContactCurrent(const DeviceMesh& mesh,
                               const MaterialDatabase& matdb,
                               const DopingModel& doping,
                               MobilityModelConfig mobilityConfig,
                               Real temperature_K,
                               DDScalingSpec scaling,
                               BandgapNarrowingConfig bandgapNarrowingConfig)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , edgeCells_(detail::buildEdgeCellMap(mesh))
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
{}


ContactCurrentResult ContactCurrent::compute(const DDSolution& solution,
                                             const std::string& contactName) const
{
    return computeDetailed(solution, contactName).totals;
}

ContactCurrentDetailedResult ContactCurrent::computeDetailed(
    const DDSolution& solution,
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

    std::unordered_set<Index> contactNodes(contact->node_ids.begin(), contact->node_ids.end());
    const Real temperature_K = thermalVoltage_ * constants::q / constants::kb;
    const std::vector<Material> cellMaterials =
        detail::buildCellMaterials(mesh_, matdb_, temperature_K);

    ContactCurrentDetailedResult detailed;
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const bool n0OnContact = contactNodes.count(edge.n0) > 0;
        const bool n1OnContact = contactNodes.count(edge.n1) > 0;
        if (n0OnContact == n1OnContact)
            continue;
        if (edge.length < 1.0e-30 || edge.couple <= 0.0)
            continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);

        // The solver returns physical DDSolution fields, so current post-processing
        // always operates in SI units regardless of the optional scaling config.
        const Real psi_i = solution.psi(i);
        const Real psi_j = solution.psi(j);
        const Real n_i = solution.n(i);
        const Real n_j = solution.n(j);
        const Real p_i = solution.p(i);
        const Real p_j = solution.p(j);
        const Real dpsi = psi_j - psi_i;
        const Real edgeLength = edge.length;

        const Real electricField = std::abs(dpsi / edgeLength);

        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials, e, CarrierType::Electron,
            electricField,
            &mobilityConfig_,
            &solution.psi);
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials, e, CarrierType::Hole,
            electricField,
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
        const Real phin_i = solution.phin(i);
        const Real phin_j = solution.phin(j);
        const Real phip_i = solution.phip(i);
        const Real phip_j = solution.phip(j);

        Real electronContinuityFlux01 = 0.0;
        Real electronFlux01 = 0.0;
        if (mun > 0.0) {
            const Real coef = mun * thermalVoltage_ / edgeLength;
            electronContinuityFlux01 = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni_i, ni_j, psi_i, psi_j, phin_i, phin_j, thermalVoltage_, coef);
            // sgElectronFlux = -sgElectronContinuityFlux by definition.
            electronFlux01 = -electronContinuityFlux01;
        }
        Real holeContinuityFlux01 = 0.0;
        Real holeFlux01 = 0.0;
        if (mup > 0.0) {
            const Real coef = mup * thermalVoltage_ / edgeLength;
            holeContinuityFlux01 = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                ni_i, ni_j, psi_i, psi_j, phip_i, phip_j, thermalVoltage_, coef);
            holeFlux01 = -holeContinuityFlux01;
        }

        // Algebraic SG split: J = J_drift + J_diffusion.
        const SGEdgeWeights weights = sgEdgeWeights(dpsi, thermalVoltage_);
        const Real bAvg = 0.5 * (weights.b_plus + weights.b_minus);
        const Real electronDriftFlux01 = (mun > 0.0)
            ? mun * (dpsi / edgeLength) * (0.5 * (n_i + n_j))
            : 0.0;
        const Real electronDiffusionFlux01 = (mun > 0.0)
            ? mun * (thermalVoltage_ / edgeLength) * bAvg * (n_i - n_j)
            : 0.0;
        const Real holeDriftFlux01 = (mup > 0.0)
            ? mup * (dpsi / edgeLength) * (0.5 * (p_i + p_j))
            : 0.0;
        const Real holeDiffusionFlux01 = (mup > 0.0)
            ? mup * (thermalVoltage_ / edgeLength) * bAvg * (p_j - p_i)
            : 0.0;

        const Real outwardSign = n0OnContact ? 1.0 : -1.0;
        // Current density [A/m^2] * edge.couple [m] = [A/m]
        const Real electronCurrent = constants::q * outwardSign * electronFlux01 * edge.couple;
        const Real electronDriftCurrent = constants::q * outwardSign * electronDriftFlux01 * edge.couple;
        const Real electronDiffusionCurrent = constants::q * outwardSign * electronDiffusionFlux01 * edge.couple;
        const Real holeCurrent = constants::q * outwardSign * holeFlux01 * edge.couple;
        const Real holeDriftCurrent = constants::q * outwardSign * holeDriftFlux01 * edge.couple;
        const Real holeDiffusionCurrent = constants::q * outwardSign * holeDiffusionFlux01 * edge.couple;

        detailed.totals.electronCurrent += electronCurrent;
        detailed.totals.electronDriftCurrent += electronDriftCurrent;
        detailed.totals.electronDiffusionCurrent += electronDiffusionCurrent;
        detailed.totals.holeCurrent += holeCurrent;
        detailed.totals.holeDriftCurrent += holeDriftCurrent;
        detailed.totals.holeDiffusionCurrent += holeDiffusionCurrent;

        ContactCurrentEdgeDiagnostic edgeDiag;
        edgeDiag.edgeId = e;
        edgeDiag.node0 = edge.n0;
        edgeDiag.node1 = edge.n1;
        edgeDiag.edgeLength_m = edge.length;
        edgeDiag.edgeCouple_m = edge.couple;
        edgeDiag.outwardSign = outwardSign;
        edgeDiag.bernoulliU = dpsi / thermalVoltage_;
        edgeDiag.bernoulliBplus = weights.b_plus;
        edgeDiag.bernoulliBminus = weights.b_minus;
        edgeDiag.electronUsedQuasiFermi = true;
        edgeDiag.holeUsedQuasiFermi = true;
        edgeDiag.psi0 = psi_i;
        edgeDiag.psi1 = psi_j;
        edgeDiag.phin0 = phin_i;
        edgeDiag.phin1 = phin_j;
        edgeDiag.phip0 = phip_i;
        edgeDiag.phip1 = phip_j;
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
        edgeDiag.electronDriftCurrent = electronDriftCurrent;
        edgeDiag.electronDiffusionCurrent = electronDiffusionCurrent;
        edgeDiag.holeCurrent = holeCurrent;
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
    const BandgapNarrowingConfig& bandgapNarrowingConfig)
{
    return ContactCurrent(mesh, matdb, doping, mobilityConfig, temperature_K, scaling, bandgapNarrowingConfig)
        .compute(solution, contactName);
}

} // namespace vela
