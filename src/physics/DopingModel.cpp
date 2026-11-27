#include "vela/physics/DopingModel.h"
#include "vela/core/ArithmeticDecimal.h"
#include <cmath>
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
    if (hasDecimalNodeDoping())
        decimalDoping_[nodeId] = {donors, acceptors,
            static_cast<ArithmeticReal>(donors) - static_cast<ArithmeticReal>(acceptors)};
}

void DopingModel::setNodeDopingFromDecimal(Index nodeId, std::string_view donors, std::string_view acceptors,
                                         ArithmeticReal concentrationFactor)
{
    if (compiledArithmeticPrecision != ArithmeticPrecision::Extended80)
        throw std::invalid_argument("node_doping_decimal_input requires an extended80 runner");
    const auto nd = parseArithmeticDecimal(donors);
    const auto na = parseArithmeticDecimal(acceptors);
    const auto scaledNd = nd * concentrationFactor;
    const auto scaledNa = na * concentrationFactor;
    const Real storedNd = static_cast<Real>(scaledNd);
    const Real storedNa = static_cast<Real>(scaledNa);
    const auto net = (nd - na) * concentrationFactor;
    if (!std::isfinite(storedNd) || !std::isfinite(storedNa) || !std::isfinite(net))
        throw std::invalid_argument("Nodal doping decimal operand overflows internal units");
    if (!hasDecimalNodeDoping()) {
        decimalDoping_.resize(donors_.size());
        for (Index i = 0; i < numNodes(); ++i)
            decimalDoping_[i] = {donors_[i], acceptors_[i],
                static_cast<ArithmeticReal>(donors_[i]) - static_cast<ArithmeticReal>(acceptors_[i])};
    }
    setNodeDoping(nodeId, storedNd, storedNa);
    decimalDoping_[nodeId] = {scaledNd, scaledNa, net};
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

ArithmeticReal DopingModel::netDopingForAssembly(Index nodeId) const
{
    return hasDecimalNodeDoping() ? decimalDoping_.at(nodeId).net : netDoping(nodeId);
}

std::array<ArithmeticReal,2> DopingModel::dopingOperandsForAssembly(Index nodeId) const
{
    if (!hasDecimalNodeDoping())
        return {donors(nodeId), acceptors(nodeId)};
    const auto& value = decimalDoping_.at(nodeId);
    return {value.donors, value.acceptors};
}

Real DopingModel::totalImpurity(Index nodeId) const
{
    return donors(nodeId) + acceptors(nodeId);
}

DopingModel DopingModel::fromMeshAndRegions(
    const DeviceMesh&                    mesh,
    const std::vector<RegionDopingSpec>& specs)
{
    // Build region-name -> {donors, acceptors} map from specs
    std::unordered_map<std::string, const RegionDopingSpec*> regionSpec;
    for (const auto& s : specs)
        regionSpec[s.region] = &s;

    // Build region-id -> spec lookup
    std::unordered_map<Index, const RegionDopingSpec*> regionIdSpec;
    for (Index r = 0; r < mesh.numRegions(); ++r) {
        const auto& reg = mesh.getRegion(r);
        auto it = regionSpec.find(reg.name);
        if (it != regionSpec.end())
            regionIdSpec[reg.id] = it->second;
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

        const RegionDopingSpec& rd = *it->second;
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
