#pragma once

#include "vela/mesh/DelaunayBox.h"
#include <set>

namespace vela {

// Unweighted Tri3 element-intersection AverageBox construction. For each
// endpoint, intersect its neighbor bisector half-planes with the edge's
// perpendicular bisector. Assign the resulting face to each cell's side of
// the edge, then average the two endpoint face lengths. The side is an
// unbounded half-plane, NOT the finite triangle: obtuse boxes extend across
// original cell boundaries. No triangulation or physical model is changed.
inline DelaunayBoxResult buildAverageCellBox(const DeviceMesh& mesh)
{
    using Wide = long double;
    struct Part { Index cell; std::size_t local; int side; };
    std::map<std::pair<Index, Index>, std::vector<Part>> edges;
    std::vector<std::set<Index>> neighbors(mesh.numNodes());
    DelaunayBoxResult result;
    result.coefficients.resize(mesh.numCells());
    for (Index cid = 0; cid < mesh.numCells(); ++cid) {
        const auto& c = mesh.getCell(cid);
        if (c.id != cid || c.type != CellType::Tri3 || c.node_ids.size() != 3)
            throw std::invalid_argument("average_intersection requires indexed Tri3 cells.");
        for (Index id : c.node_ids) {
            if (id >= mesh.numNodes() || !std::isfinite(mesh.getNode(id).x) ||
                !std::isfinite(mesh.getNode(id).y))
                throw std::invalid_argument("Invalid AverageBox vertex.");
        }
        for (std::size_t k = 0; k < 3; ++k) {
            const Index i = std::min(c.node_ids[k], c.node_ids[(k+1)%3]);
            const Index j = std::max(c.node_ids[k], c.node_ids[(k+1)%3]);
            const auto& a = mesh.getNode(i);
            const auto& b = mesh.getNode(j);
            const auto& o = mesh.getNode(c.node_ids[(k+2)%3]);
            const Wide ex = Wide(b.x)-a.x, ey = Wide(b.y)-a.y;
            const Wide cross = ex*(Wide(o.y)-a.y)-ey*(Wide(o.x)-a.x);
            if (!std::isfinite(cross) || std::abs(cross) <= 1e-30L)
                throw std::invalid_argument("AverageBox rejects degenerate cells.");
            neighbors[i].insert(j);
            neighbors[j].insert(i);
            edges[{i,j}].push_back({cid,k,cross>0 ? 1 : -1});
        }
    }
    for (const auto& [key, parts] : edges) {
        if (parts.size()>2 || (parts.size()==2 && parts[0].side==parts[1].side))
            throw std::invalid_argument("AverageBox rejects nonmanifold or overlapping edge fans.");
        const auto [i,j] = key;
        const auto& a = mesh.getNode(i);
        const auto& b = mesh.getNode(j);
        const Wide ex = Wide(b.x)-a.x, ey = Wide(b.y)-a.y;
        // x(t) = midpoint + t*(-ey,ex); t is already face length / edge length.
        std::array<std::pair<Wide,Wide>,2> intervals;
        for (int endpoint = 0; endpoint < 2; ++endpoint) {
            const Index v = endpoint==0 ? i : j;
            const auto& origin = mesh.getNode(v);
            const Wide sign = endpoint==0 ? 1 : -1;
            Wide lo = -std::numeric_limits<Wide>::infinity();
            Wide hi = std::numeric_limits<Wide>::infinity();
            for (Index n : neighbors[v]) {
                if (n==i || n==j) continue;
                const auto& p = mesh.getNode(n);
                const Wide nx = Wide(p.x)-origin.x, ny = Wide(p.y)-origin.y;
                const Wide slope = -ey*nx+ex*ny;
                // Translate before multiplying: avoid subtracting squared
                // absolute coordinates near short edges at large offsets.
                const Wide bound = (nx*(nx-sign*ex)+ny*(ny-sign*ey))/2;
                if (slope>0) hi=std::min(hi,bound/slope);
                else if (slope<0) lo=std::max(lo,bound/slope);
                else if (bound<0) { lo=1; hi=0; break; }
            }
            intervals[endpoint]={lo,hi};
        }
        for (const auto& part : parts) {
            Wide coefficient=0;
            for (const auto& [lo,hi] : intervals) {
                const Wide low=part.side>0 ? std::max(lo,Wide(0)) : lo;
                const Wide high=part.side>0 ? hi : std::min(hi,Wide(0));
                coefficient+=std::max(Wide(0),high-low)/2;
            }
            if (!std::isfinite(coefficient) || coefficient>std::numeric_limits<Real>::max())
                throw std::invalid_argument("Unbounded or nonfinite AverageBox face.");
            result.coefficients[part.cell][part.local]=static_cast<Real>(coefficient);
        }
    }
    return result;
}
} // namespace vela
