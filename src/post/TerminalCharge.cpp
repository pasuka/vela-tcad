#include "vela/post/TerminalCharge.h"
#include "vela/mesh/BoxGeometryBuilder.h"
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace vela {
namespace {

const Contact* findContact(const DeviceMesh& mesh, const std::string& name)
{
    for (const Contact& contact : mesh.contacts()) {
        if (contact.name == name)
            return &contact;
    }
    return nullptr;
}

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
            throw std::invalid_argument("TerminalCharge: unknown region '" + name + "'.");
    }
    if (volumes.empty())
        throw std::invalid_argument("TerminalCharge: selected regions contain no supported cell volume.");

    return volumes;
}

bool withinContactRadius(const DeviceMesh& mesh,
                         const Contact& contact,
                         Index nodeId,
                         Real radius)
{
    if (radius <= 0.0) {
        for (Index contactNodeId : contact.node_ids) {
            if (contactNodeId == nodeId)
                return true;
        }
        return false;
    }

    const Node& node = mesh.getNode(nodeId);
    for (Index contactNodeId : contact.node_ids) {
        const Node& contactNode = mesh.getNode(contactNodeId);
        const Real dx = node.x - contactNode.x;
        const Real dy = node.y - contactNode.y;
        if (std::hypot(dx, dy) <= radius)
            return true;
    }
    return false;
}

} // namespace

TerminalCharge::TerminalCharge(const DeviceMesh& mesh, const DopingModel& doping)
    : mesh_(mesh)
    , doping_(doping)
{}

TerminalChargeResult TerminalCharge::compute(const DDSolution& solution,
                                             const TerminalChargeConfig& config) const
{
    if (solution.n.size() != static_cast<int>(mesh_.numNodes()) ||
        solution.p.size() != static_cast<int>(mesh_.numNodes())) {
        throw std::invalid_argument("TerminalCharge: solution size does not match mesh.");
    }
    if (!config.perMeter && config.depth_m <= 0.0)
        throw std::invalid_argument("TerminalCharge: depth_m must be positive for total charge.");

    const Contact* contact = nullptr;
    if (!config.contact.empty()) {
        contact = findContact(mesh_, config.contact);
        if (contact == nullptr)
            throw std::invalid_argument("TerminalCharge: unknown contact '" + config.contact + "'.");
    }

    std::unordered_map<Index, Real> selectedVolumes = selectedNodeVolumes(mesh_, config.regions);

    if (contact != nullptr) {
        for (auto it = selectedVolumes.begin(); it != selectedVolumes.end();) {
            if (!withinContactRadius(mesh_, *contact, it->first, config.contactRadius))
                it = selectedVolumes.erase(it);
            else
                ++it;
        }
    }

    Real chargePerMeter = 0.0;
    for (const auto& [nodeId, volume] : selectedVolumes) {
        Real rho = 0.0;
        if (config.includeMobileCharge)
            rho += solution.p(static_cast<int>(nodeId)) - solution.n(static_cast<int>(nodeId));
        if (config.includeIonizedDopants)
            rho += doping_.netDoping(nodeId);
        chargePerMeter += constants::q * rho * volume;
    }

    TerminalChargeResult result;
    result.perMeter = config.perMeter;
    result.charge = config.perMeter ? chargePerMeter : chargePerMeter * config.depth_m;
    return result;
}

TerminalChargeResult TerminalCharge::compute(const DeviceMesh& mesh,
                                             const DopingModel& doping,
                                             const DDSolution& solution,
                                             const TerminalChargeConfig& config)
{
    return TerminalCharge(mesh, doping).compute(solution, config);
}

} // namespace vela
