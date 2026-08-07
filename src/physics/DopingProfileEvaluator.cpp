#include "vela/physics/DopingProfileEvaluator.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace vela {
namespace {

std::unordered_map<Index, std::string> buildRegionNameById(const DeviceMesh& mesh)
{
    std::unordered_map<Index, std::string> names;
    names.reserve(mesh.numRegions());
    for (const auto& region : mesh.regions()) {
        names.emplace(region.id, region.name);
    }
    return names;
}

std::vector<std::vector<Index>> buildNodeRegionMembership(const DeviceMesh& mesh)
{
    std::vector<std::vector<Index>> membership(mesh.numNodes());
    for (const auto& cell : mesh.cells()) {
        for (const auto nodeId : cell.node_ids) {
            membership[nodeId].push_back(cell.region_id);
        }
    }
    for (auto& nodeRegions : membership) {
        std::sort(nodeRegions.begin(), nodeRegions.end());
        nodeRegions.erase(std::unique(nodeRegions.begin(), nodeRegions.end()), nodeRegions.end());
    }
    return membership;
}

bool nodeInTargetRegion(const std::vector<Index>& nodeRegionIds,
                        const std::unordered_map<Index, std::string>& regionNameById,
                        const std::string& targetRegion)
{
    for (const auto regionId : nodeRegionIds) {
        const auto regionIt = regionNameById.find(regionId);
        if (regionIt != regionNameById.end() && regionIt->second == targetRegion) {
            return true;
        }
    }
    return false;
}

} // namespace

DopingModel DopingProfileEvaluator::evaluate(const DeviceMesh& mesh, const DeviceIr2D& ir)
{
    DopingModel model(mesh.numNodes());
    const auto regionNameById = buildRegionNameById(mesh);
    const auto nodeMembership = buildNodeRegionMembership(mesh);

    std::vector<const DopingProfileIr*> profiles;
    profiles.reserve(ir.dopingProfiles.size());
    for (const auto& profile : ir.dopingProfiles) {
        profiles.push_back(&profile);
    }
    std::stable_sort(profiles.begin(), profiles.end(), [](const auto* lhs, const auto* rhs) {
        return lhs->priority < rhs->priority;
    });

    for (Index nodeId = 0; nodeId < mesh.numNodes(); ++nodeId) {
        Real donors = 0.0;
        Real acceptors = 0.0;
        const auto& node = mesh.getNode(nodeId);

        for (const auto* profile : profiles) {
            if (!nodeInTargetRegion(nodeMembership[nodeId], regionNameById, profile->targetRegion)) {
                continue;
            }
            if (profile->kind == DopingProfileKind::Constant) {
                donors += profile->donors_cm3;
                acceptors += profile->acceptors_cm3;
                continue;
            }
            if (profile->kind == DopingProfileKind::Gaussian) {
                if (!std::isfinite(profile->gaussianPeak_cm3) ||
                    !std::isfinite(profile->gaussianValueAtDepth_cm3) ||
                    !std::isfinite(profile->gaussianSigmaXUm) ||
                    profile->gaussianPeak_cm3 <= 0.0 ||
                    profile->gaussianValueAtDepth_cm3 < 0.0 ||
                    profile->gaussianValueAtDepth_cm3 >= profile->gaussianPeak_cm3 ||
                    profile->gaussianSigmaXUm <= 0.0) {
                    throw std::invalid_argument(
                        "DopingProfileEvaluator: invalid gaussian parameters for profile '" +
                        profile->name + "'.");
                }
                const Real dx = (node.x - profile->gaussianPeakPosUm.x_um) / profile->gaussianSigmaXUm;
                const Real concentration = profile->gaussianPeak_cm3 * std::exp(-0.5 * dx * dx);
                if (profile->gaussianActsOnDonors) {
                    donors += concentration;
                } else {
                    acceptors += concentration;
                }
            }
        }
        model.setNodeDoping(nodeId, donors, acceptors);
    }

    return model;
}

} // namespace vela
