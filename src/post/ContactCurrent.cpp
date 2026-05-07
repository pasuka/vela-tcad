#include "vela/post/ContactCurrent.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/equation/AssemblerUtils.h"
#include <unordered_set>
#include <stdexcept>

namespace vela {

ContactCurrent::ContactCurrent(const DeviceMesh& mesh,
                               const MaterialDatabase& matdb,
                               const DopingModel& doping,
                               MobilityModelConfig mobilityConfig)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , edgeCells_(detail::buildEdgeCellMap(mesh))
    , mobility_(makeMobilityModel(mobilityConfig))
{
}

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

    ContactCurrentResult result;
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const bool n0OnContact = contactNodes.count(edge.n0) > 0;
        const bool n1OnContact = contactNodes.count(edge.n1) > 0;
        if (n0OnContact == n1OnContact)
            continue;
        if (edge.length < 1.0e-30 || edge.couple <= 0.0)
            continue;

        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, matdb_, doping_, *mobility_, e, CarrierType::Electron);
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, matdb_, doping_, *mobility_, e, CarrierType::Hole);

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real dpsi = solution.psi(j) - solution.psi(i);

        const Real electronFlux01 = (mun > 0.0)
            ? sgElectronFlux(solution.n(i), solution.n(j), dpsi,
                             constants::Vt_300, mun, edge.length)
            : 0.0;
        const Real holeFlux01 = (mup > 0.0)
            ? sgHoleFlux(solution.p(i), solution.p(j), dpsi,
                         constants::Vt_300, mup, edge.length)
            : 0.0;

        const Real outwardSign = n0OnContact ? 1.0 : -1.0;
        result.electronCurrent += constants::q * outwardSign * electronFlux01 * edge.couple;
        result.holeCurrent += constants::q * outwardSign * holeFlux01 * edge.couple;
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
    const MobilityModelConfig& mobilityConfig)
{
    return ContactCurrent(mesh, matdb, doping, mobilityConfig).compute(solution, contactName);
}

} // namespace vela
