#include "vela/physics/DopingModel.h"
#include <stdexcept>
#include <unordered_map>

namespace vela {

DopingModel::DopingModel(Index numNodes)
    : donors_(numNodes, 0.0)
    , acceptors_(numNodes, 0.0)
{}

void DopingModel::setNodeDoping(Index nodeId, Real donors, Real acceptors)
{
    if (nodeId >= donors_.size())
        throw std::out_of_range("DopingModel: node id out of range: "
                                + std::to_string(nodeId));
    donors_[nodeId]    = donors;
    acceptors_[nodeId] = acceptors;
}

Real DopingModel::donors(Index nodeId) const
{
    if (nodeId >= donors_.size())
        throw std::out_of_range("DopingModel: node id out of range: "
                                + std::to_string(nodeId));
    return donors_[nodeId];
}

Real DopingModel::acceptors(Index nodeId) const
{
    if (nodeId >= acceptors_.size())
        throw std::out_of_range("DopingModel: node id out of range: "
                                + std::to_string(nodeId));
    return acceptors_[nodeId];
}

Real DopingModel::netDoping(Index nodeId) const
{
    return donors(nodeId) - acceptors(nodeId);
}

DopingModel DopingModel::fromMeshAndRegions(const DeviceMesh&     mesh,
                                             const nlohmann::json& dopingArray)
{
    // Build region-name → {donors, acceptors} map from JSON
    struct RegionDoping { Real donors = 0.0; Real acceptors = 0.0; };
    std::unordered_map<std::string, RegionDoping> regionSpec;

    for (const auto& entry : dopingArray) {
        RegionDoping rd;
        rd.donors    = entry.at("donors").get<Real>();
        rd.acceptors = entry.at("acceptors").get<Real>();
        regionSpec[entry.at("region").get<std::string>()] = rd;
    }

    // Build region-id → doping lookup
    std::unordered_map<Index, RegionDoping> regionIdSpec;
    for (Index r = 0; r < mesh.numRegions(); ++r) {
        const auto& reg = mesh.getRegion(r);
        auto it = regionSpec.find(reg.name);
        if (it != regionSpec.end()) {
            regionIdSpec[reg.id] = it->second;
        }
    }

    // Accumulate per-node doping (average if node shared between regions)
    const Index N = mesh.numNodes();
    std::vector<Real> sumDon(N, 0.0);
    std::vector<Real> sumAcc(N, 0.0);
    std::vector<int>  count(N, 0);

    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& cell = mesh.getCell(c);
        auto it = regionIdSpec.find(cell.region_id);
        if (it == regionIdSpec.end()) continue;

        const RegionDoping& rd = it->second;
        for (Index nid : cell.node_ids) {
            sumDon[nid] += rd.donors;
            sumAcc[nid] += rd.acceptors;
            ++count[nid];
        }
    }

    DopingModel model(N);
    for (Index i = 0; i < N; ++i) {
        if (count[i] > 0) {
            model.setNodeDoping(i,
                sumDon[i] / count[i],
                sumAcc[i] / count[i]);
        }
    }
    return model;
}

} // namespace vela
