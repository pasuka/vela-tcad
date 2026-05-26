#include "vela/post/ContactCurrent.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/equation/AssemblerUtils.h"
#include <unordered_set>
#include <cmath>
#include <stdexcept>
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
                               DDScalingSpec scaling)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , edgeCells_(detail::buildEdgeCellMap(mesh))
    , mobilityConfig_(mobilityConfig)
    , mobility_(makeMobilityModel(mobilityConfig))
    , thermalVoltage_(validatedThermalVoltage(temperature_K))
    , scaling_(scaling)
{}


ContactCurrentResult ContactCurrent::compute(const DDSolution& solution,
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

    ContactCurrentResult result;
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

        // SG fluxes in physical units
        const Real electronFlux01 = (mun > 0.0)
            ? sgElectronFlux(n_i, n_j, dpsi, thermalVoltage_, mun, edgeLength)
            : 0.0;
        const Real holeFlux01 = (mup > 0.0)
            ? sgHoleFlux(p_i, p_j, dpsi, thermalVoltage_, mup, edgeLength)
            : 0.0;

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
        result.electronCurrent += constants::q * outwardSign * electronFlux01 * edge.couple;
        result.electronDriftCurrent += constants::q * outwardSign * electronDriftFlux01 * edge.couple;
        result.electronDiffusionCurrent += constants::q * outwardSign * electronDiffusionFlux01 * edge.couple;
        result.holeCurrent += constants::q * outwardSign * holeFlux01 * edge.couple;
        result.holeDriftCurrent += constants::q * outwardSign * holeDriftFlux01 * edge.couple;
        result.holeDiffusionCurrent += constants::q * outwardSign * holeDiffusionFlux01 * edge.couple;
    }

    result.totalCurrent = result.electronCurrent + result.holeCurrent;
    return result;
}


ContactCurrentResult ContactCurrent::compute(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const DDSolution& solution,
    const std::string& contactName,
    const MobilityModelConfig& mobilityConfig,
    Real temperature_K,
    DDScalingSpec scaling)
{
    return ContactCurrent(mesh, matdb, doping, mobilityConfig, temperature_K, scaling).compute(solution, contactName);
}

} // namespace vela
