#pragma once
#include "vela/mesh/DeviceMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

namespace vela {
struct DelaunayBoxResult {
    std::vector<std::array<Real, 3>> coefficients;
    Index transferredEdges = 0;
};

// Unweighted 2-D Delaunay edge transfer. Negative local contributions may
// move only to the other Tri3 cell in the same region. Global non-Delaunay
// edges and obtuse region boundaries need a different intersection algorithm.
inline DelaunayBoxResult buildDelaunayCellBox(const DeviceMesh& mesh)
{
    struct Part { Index cell; std::size_t local; long double value; };
    std::map<std::pair<Index, Index>, std::vector<Part>> edges;
    DelaunayBoxResult result; result.coefficients.resize(mesh.numCells());
    for (Index cid = 0; cid < mesh.numCells(); ++cid) {
        const auto& c = mesh.getCell(cid);
        if (c.id != cid || c.type != CellType::Tri3 || c.node_ids.size() != 3)
            throw std::invalid_argument("Delaunay cell box requires indexed Tri3 cells.");
        for (std::size_t k = 0; k < 3; ++k) {
            const Index i=c.node_ids[k], j=c.node_ids[(k+1)%3];
            const auto& a=mesh.getNode(i);const auto& b=mesh.getNode(j);
            const auto& o=mesh.getNode(c.node_ids[(k+2)%3]);
            const long double ux=(long double)a.x-o.x, uy=(long double)a.y-o.y;
            const long double vx=(long double)b.x-o.x, vy=(long double)b.y-o.y;
            const long double cross=std::abs(ux*vy-uy*vx);
            if (!(cross > 1e-30L)) throw std::invalid_argument("Delaunay cell box rejects degenerate cells.");
            const long double value=(ux*vx+uy*vy)/(2*cross);
            if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite Delaunay cell coefficient.");
            edges[{std::min(i,j),std::max(i,j)}].push_back({cid,k,value});
        }
    }
    for (auto& [edge, parts] : edges) {
        if (parts.size()>2) throw std::invalid_argument("Delaunay cell box rejects nonmanifold edges.");
        long double sum=0,scale=1;
        for (const auto& p : parts) {sum+=p.value;scale=std::max(scale,std::abs(p.value));}
        const long double tol=128*std::numeric_limits<double>::epsilon()*scale;
        if (sum < -tol) throw std::invalid_argument("Delaunay cell box rejects globally non-Delaunay edges.");
        auto negative=std::min_element(parts.begin(),parts.end(),[](const Part& a,const Part& b){return a.value<b.value;});
        if (negative->value < 0) {
            const bool transferable=parts.size()==2 && mesh.getCell(parts[0].cell).region_id==mesh.getCell(parts[1].cell).region_id;
            if (!transferable && negative->value < -tol)
                throw std::invalid_argument("Delaunay cell box rejects obtuse region-boundary edges.");
            if (transferable) {
                for (auto& p : parts) p.value=(&p==&*negative)?0:std::max(sum,0.L);
                ++result.transferredEdges;
            } else negative->value=0;
        }
        for (const auto& p : parts) result.coefficients[p.cell][p.local]=static_cast<Real>(p.value);
    }
    return result;
}
} // namespace vela
