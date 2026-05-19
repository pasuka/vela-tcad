#include "vela/post/StoredCharge.h"
#include "vela/mesh/BoxGeometryBuilder.h"
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace vela {
namespace {
std::unordered_map<Index, Real> selectedNodeVolumes(
    const DeviceMesh& mesh,
    const std::vector<std::string>& regionNames)
{
    std::unordered_map<Index, Real> volumes;
    if (regionNames.empty()) {
        for (Index i = 0; i < mesh.numNodes(); ++i)
            volumes[i] = mesh.getNode(i).volume;
        return volumes;
    }

    std::unordered_set<std::string> wanted(regionNames.begin(), regionNames.end());
    std::unordered_set<std::string> found;
    for (const Region& region : mesh.regions()) {
        if (wanted.count(region.name) == 0)
            continue;

        found.insert(region.name);
        for (Index cellId : region.cell_ids) {
            const Cell& cell = mesh.getCell(cellId);
            if (cell.type != CellType::Tri3 || cell.node_ids.size() < 3)
                continue;

            const Real area = BoxGeometryBuilder::triangleArea(
                mesh.getNode(cell.node_ids[0]),
                mesh.getNode(cell.node_ids[1]),
                mesh.getNode(cell.node_ids[2]));
            if (area <= 0.0)
                continue;

            const Real nodeShare = area / 3.0;
            for (Index nodeId : cell.node_ids)
                volumes[nodeId] += nodeShare;
        }
    }

    for (const std::string& name : wanted) {
        if (found.count(name) == 0)
            throw std::invalid_argument("StoredCharge: unknown region '" + name + "'.");
    }
    if (volumes.empty())
        throw std::invalid_argument("StoredCharge: selected regions contain no supported cell volume.");

    return volumes;
}
} // namespace

StoredCharge::StoredCharge(const DeviceMesh& mesh)
    : mesh_(mesh)
{}

StoredChargeResult StoredCharge::compute(const DDSolution& solution,
                                         const StoredChargeConfig& config) const
{
    if (solution.n.size() != static_cast<int>(mesh_.numNodes()) ||
        solution.p.size() != static_cast<int>(mesh_.numNodes())) {
        throw std::invalid_argument("StoredCharge: solution size does not match mesh.");
    }
    if (!config.perMeter && config.depth_m <= 0.0)
        throw std::invalid_argument("StoredCharge: depth_m must be positive for total charge.");

    const std::unordered_map<Index, Real> selectedVolumes = selectedNodeVolumes(mesh_, config.regions);

    Real chargePerMeter = 0.0;
    for (const auto& [nodeId, volume] : selectedVolumes) {
        const Real mobileSum = solution.n(static_cast<int>(nodeId)) + solution.p(static_cast<int>(nodeId));
        chargePerMeter += constants::q * mobileSum * volume;
    }

    StoredChargeResult result;
    result.perMeter = config.perMeter;
    result.charge = config.perMeter ? chargePerMeter : chargePerMeter * config.depth_m;
    return result;
}

} // namespace vela
