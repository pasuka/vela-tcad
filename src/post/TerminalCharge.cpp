#include "vela/post/TerminalCharge.h"
#include <cmath>
#include <stdexcept>
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

std::unordered_set<Index> regionNodeSet(const DeviceMesh& mesh,
                                        const std::vector<std::string>& regionNames)
{
    std::unordered_set<std::string> wanted(regionNames.begin(), regionNames.end());
    std::unordered_set<Index> nodes;
    for (const Region& region : mesh.regions()) {
        if (!wanted.empty() && wanted.count(region.name) == 0)
            continue;
        for (Index cellId : region.cell_ids) {
            const Cell& cell = mesh.getCell(cellId);
            for (Index nodeId : cell.node_ids)
                nodes.insert(nodeId);
        }
    }
    return nodes;
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

    std::unordered_set<Index> selected = regionNodeSet(mesh_, config.regions);
    if (selected.empty() && config.regions.empty()) {
        for (Index i = 0; i < mesh_.numNodes(); ++i)
            selected.insert(i);
    }

    if (contact != nullptr) {
        for (auto it = selected.begin(); it != selected.end();) {
            if (!withinContactRadius(mesh_, *contact, *it, config.contactRadius))
                it = selected.erase(it);
            else
                ++it;
        }
    }

    Real chargePerMeter = 0.0;
    for (Index nodeId : selected) {
        const Node& node = mesh_.getNode(nodeId);
        Real rho = 0.0;
        if (config.includeMobileCharge)
            rho += solution.p(static_cast<int>(nodeId)) - solution.n(static_cast<int>(nodeId));
        if (config.includeIonizedDopants)
            rho += doping_.netDoping(nodeId);
        chargePerMeter += constants::q * rho * node.volume;
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
