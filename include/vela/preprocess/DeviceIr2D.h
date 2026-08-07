#pragma once

#include "vela/core/Types.h"
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace vela {

struct Point2D {
    Real x_um = 0.0;
    Real y_um = 0.0;
};

enum class GeometryPrimitiveKind {
    Rectangle,
    Polygon,
};

struct GeometryPrimitiveIr {
    GeometryPrimitiveKind kind = GeometryPrimitiveKind::Rectangle;
    std::string name;
    std::string region;
    std::string material;
    std::vector<Point2D> points;
};

struct RegionIr {
    std::string name;
    std::string material;
    std::vector<std::size_t> primitiveIndices;
};

struct ContactIr {
    std::string name;
    std::string ownerRegion;
    std::vector<std::string> boundaryRefs;
};

enum class DopingProfileKind {
    Constant,
    Gaussian,
};

struct DopingProfileIr {
    std::string name;
    DopingProfileKind kind = DopingProfileKind::Constant;
    std::string targetRegion;
    int priority = 0;
    Real donors_cm3 = 0.0;
    Real acceptors_cm3 = 0.0;
    Real gaussianPeak_cm3 = 0.0;
    Real gaussianValueAtDepth_cm3 = 0.0;
    Point2D gaussianPeakPosUm{};
    Real gaussianSigmaXUm = 0.0;
    bool gaussianActsOnDonors = true;
};

struct MeshControlIr {
    std::optional<Real> globalTargetSizeUm;
    std::unordered_map<std::string, Real> regionTargetSizeUm;
    bool refineByDopingGradient = false;
    Real dopingGradientThresholdCm3PerUm = 0.0;
    Real minSizeUm = 0.0;
};

struct DeviceIr2D {
    std::vector<GeometryPrimitiveIr> geometry;
    std::vector<RegionIr> regions;
    std::vector<ContactIr> contacts;
    std::vector<DopingProfileIr> dopingProfiles;
    MeshControlIr meshControl;

    RegionIr* findRegion(std::string_view regionName)
    {
        for (auto& region : regions) {
            if (region.name == regionName) {
                return &region;
            }
        }
        return nullptr;
    }

    const RegionIr* findRegion(std::string_view regionName) const
    {
        for (const auto& region : regions) {
            if (region.name == regionName) {
                return &region;
            }
        }
        return nullptr;
    }
};

} // namespace vela
