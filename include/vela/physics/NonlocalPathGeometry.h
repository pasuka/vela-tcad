#pragma once
#include "vela/physics/NonlocalPathBandToBand.h"
#include <array>
#include <map>

namespace vela::experimental {
enum class TunnelingBoundaryKind { Reflecting, Contact, Open };
struct TunnelingTriangle {
    std::array<Index, 3> nodes;
    Index material = 0, region = 0;
    bool semiconductor = true;
};
struct TunnelingBoundaryEdge {
    std::array<Index, 2> nodes;
    TunnelingBoundaryKind kind = TunnelingBoundaryKind::Open;
};
struct TunnelingGeometryMesh {
    std::uint64_t revision = 0;
    Index siliconMaterial = 0;
    std::vector<Point2> positions_m;
    std::vector<TunnelingTriangle> cells;
    // Every exterior edge requires an explicit condition. Interior Si/oxide
    // interfaces are discovered from the cells, never supplied as fit tables.
    std::vector<TunnelingBoundaryEdge> boundaries;
};
struct TunnelingGeometryOptions {
    bool withinSelectedRegions = true;
    Index maximumSegments = 4096;
    Index maximumReflections = 32;
};
struct TunnelingPathSegment {
    TunnelingEndpoint start, end;
    Real begin_m = 0, end_m = 0;
};
struct TunnelingGeometryResult {
    TunnelingPath path;
    std::vector<TunnelingPathSegment> segments;
    Index reflections = 0;
    Real initialValenceGradient_J_per_m = 0;
    Real sourceChord_m = 0; // Nonzero only for the private cell-seeded policy.
};

// Immutable geometry only; no state-dependent path cache and no rate kernel.
class NonlocalPathGeometry {
public:
    explicit NonlocalPathGeometry(TunnelingGeometryMesh mesh);
    TunnelingGeometryResult traceFromPoint(const TunnelingBandState&,
        const NonlocalPathBandToBandConfig&, Index cell, const Point2& position_m,
        const TunnelingGeometryOptions& = {}) const;
    // One deterministic seed per active Si vertex, including contacts.
    std::vector<TunnelingGeometryResult> traceFromNodes(const TunnelingBandState&,
        const NonlocalPathBandToBandConfig&, const TunnelingGeometryOptions& = {}) const;
    // Retain complete face-to-face segments for energy-window integration;
    // this differs from traceFromPoint's single-energy turning endpoint.
    std::vector<TunnelingGeometryResult> traceFromCells(const TunnelingBandState&,
        const NonlocalPathBandToBandConfig&, const TunnelingGeometryOptions& = {}) const;
private:
    using EdgeKey = std::pair<Index, Index>;
    struct CellGeometry {
        std::array<Point2, 3> gradients; // Barycentric gradients [m^-1].
        Real twiceArea_m2;
    };
    TunnelingGeometryMesh mesh_;
    std::vector<CellGeometry> geometry_;
    std::map<EdgeKey, std::vector<Index>> edges_;
    std::map<EdgeKey, TunnelingBoundaryKind> boundaries_;
    std::vector<std::vector<Index>> incident_;
    std::array<Real, 3> weights(Index, const Point2&) const;
    TunnelingEndpoint support(Index, const Point2&) const;
    Point2 gradient(Index, const TunnelingBandState&, bool valence) const;
    bool enters(Index, const Point2&, const Point2&) const;
    void validate(const TunnelingBandState&, const NonlocalPathBandToBandConfig&,
                  const TunnelingGeometryOptions&) const;
    TunnelingGeometryResult trace(const TunnelingBandState&, const NonlocalPathBandToBandConfig&,
        Index, const Point2&, Point2, std::uint64_t, const TunnelingGeometryOptions&) const;
};
const char* tunnelingTerminationName(TunnelingTermination);
nlohmann::json tunnelingGeometryAudit(const TunnelingGeometryResult&);
} // namespace vela::experimental
