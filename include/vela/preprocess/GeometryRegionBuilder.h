#pragma once

#include "vela/preprocess/DeviceIr2D.h"
#include <string>
#include <vector>

namespace vela {

struct RegionTopologyTag {
    Index id = 0;
    std::string name;
    std::string material;
};

struct BoundaryTopologyTag {
    std::string ref;
    Index regionId = 0;
    Point2D startUm{};
    Point2D endUm{};
    bool isExternal = true;
    std::string adjacentRegion;
    std::string contactName;
};

struct GeometryRegionTopology {
    std::vector<RegionTopologyTag> regions;
    std::vector<BoundaryTopologyTag> boundaries;
};

class GeometryRegionBuilder {
public:
    GeometryRegionTopology build(const DeviceIr2D& ir) const;
};

} // namespace vela

