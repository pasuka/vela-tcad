#include "vela/physics/NonlocalPathGeometry.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace vela::experimental {
namespace {
constexpr Real tolerance = 128 * std::numeric_limits<Real>::epsilon();
void demand(bool ok, const char* message) { if (!ok) throw std::invalid_argument(message); }
Real norm(const Point2& x) { return std::hypot(x.x(), x.y()); }
Real cross(const Point2& a, const Point2& b) { return a.x() * b.y() - a.y() * b.x(); }
auto edge(Index a, Index b) { return std::pair<Index, Index>{std::min(a, b), std::max(a, b)}; }
bool selected(const NonlocalPathBandToBandConfig& c, Index region) {
    return std::find(c.regions.begin(), c.regions.end(), region) != c.regions.end();
}
nlohmann::json endpointJson(const TunnelingEndpoint& p) {
    return {{"cell", p.cell}, {"material", p.material}, {"region", p.region}, {"weights", p.weights}};
}
}
NonlocalPathGeometry::NonlocalPathGeometry(TunnelingGeometryMesh mesh) : mesh_(std::move(mesh)) {
    demand(!mesh_.positions_m.empty() && !mesh_.cells.empty(), "empty tunneling geometry");
    incident_.resize(mesh_.positions_m.size());
    std::set<std::array<Index, 3>> uniqueCells;
    for (const auto& p : mesh_.positions_m) demand(p.allFinite(), "nonfinite tunneling coordinate");
    for (Index k = 0; k < mesh_.cells.size(); ++k) {
        const auto& c = mesh_.cells[k];
        const auto& v = c.nodes;
        demand(v[0] < incident_.size() && v[1] < incident_.size() && v[2] < incident_.size() &&
            v[0] != v[1] && v[1] != v[2] && v[0] != v[2], "invalid tunneling triangle nodes");
        auto canonical = v;
        std::sort(canonical.begin(), canonical.end());
        demand(uniqueCells.insert(canonical).second, "duplicate tunneling triangle");
        const Point2 a = mesh_.positions_m[v[1]] - mesh_.positions_m[v[0]];
        const Point2 b = mesh_.positions_m[v[2]] - mesh_.positions_m[v[0]];
        const Real h = std::max(norm(a), norm(b));
        demand(std::isfinite(h) && h > 0, "degenerate tunneling triangle");
        const Point2 u = a / h, w = b / h;
        const Real det = cross(u, w);
        demand(std::isfinite(det) && std::abs(det) > tolerance, "ill-conditioned tunneling triangle");
        const Point2 g1(w.y() / (det * h), -w.x() / (det * h));
        const Point2 g2(-u.y() / (det * h), u.x() / (det * h));
        const Real area = h * h * std::abs(det);
        demand(std::isfinite(area) && area > 0 && g1.allFinite() && g2.allFinite(), "unrepresentable triangle geometry");
        geometry_.push_back({{-g1 - g2, g1, g2}, area});
        for (Index i = 0; i < 3; ++i) {
            incident_[v[i]].push_back(k);
            auto& adjacent = edges_[edge(v[(i + 1) % 3], v[(i + 2) % 3])];
            adjacent.push_back(k);
            demand(adjacent.size() <= 2, "nonmanifold tunneling edge");
        }
    }
    for (const auto& b : mesh_.boundaries) {
        const EdgeKey key = edge(b.nodes[0], b.nodes[1]);
        demand(edges_.contains(key) && edges_.at(key).size() == 1, "boundary must identify an exterior edge");
        demand(b.kind == TunnelingBoundaryKind::Reflecting || b.kind == TunnelingBoundaryKind::Contact ||
            b.kind == TunnelingBoundaryKind::Open, "unknown tunneling boundary kind");
        demand(boundaries_.emplace(key, b.kind).second, "duplicate tunneling boundary");
    }
    for (const auto& [key, cells] : edges_) {
        demand(cells.size() != 1 || boundaries_.contains(key), "unclassified exterior tunneling boundary");
        if (cells.size() == 2) {
            const Point2 a = mesh_.positions_m[key.first];
            const Point2 direction = mesh_.positions_m[key.second] - a;
            std::array<Real, 2> sides{};
            for (Index i = 0; i < 2; ++i) for (Index node : mesh_.cells[cells[i]].nodes)
                if (node != key.first && node != key.second) sides[i] = cross(direction / norm(direction), mesh_.positions_m[node] - a);
            demand(sides[0] != 0 && sides[1] != 0 && std::signbit(sides[0]) != std::signbit(sides[1]),
                "adjacent tunneling triangles overlap on the same side");
        }
    }
}
std::array<Real, 3> NonlocalPathGeometry::weights(Index cell, const Point2& p) const {
    const auto& c = mesh_.cells.at(cell);
    const Point2 delta = p - mesh_.positions_m[c.nodes[0]];
    const Real b = geometry_[cell].gradients[1].dot(delta);
    const Real d = geometry_[cell].gradients[2].dot(delta);
    return {1 - b - d, b, d};
}
TunnelingEndpoint NonlocalPathGeometry::support(Index cell, const Point2& p) const {
    const auto& c = mesh_.cells.at(cell);
    auto w = weights(cell, p);
    if (std::any_of(w.begin(),w.end(),[](Real v){return !std::isfinite(v) || v < -tolerance || v > 1+tolerance;}))
        throw std::invalid_argument("path left its support triangle: "+nlohmann::json{{"cell",cell},{"position_m",{p.x(),p.y()}},{"weights",w}}.dump());
    Real total = 0;
    for (auto& v : w) {
        demand(std::isfinite(v) && v >= -tolerance && v <= 1 + tolerance, "path left its support triangle");
        v = std::clamp(v, Real(0), Real(1)); total += v;
    }
    TunnelingEndpoint result{cell, c.material, c.region, {}};
    for (Index i = 0; i < 3; ++i) if (w[i] > 0) result.weights.emplace_back(c.nodes[i], w[i] / total);
    return result;
}
Point2 NonlocalPathGeometry::gradient(Index cell, const TunnelingBandState& s, bool valence) const {
    const auto& c = mesh_.cells.at(cell);
    const auto value = [&](Index i) { return valence ? s.nodes[i].valence_J : s.nodes[i].conduction_J; };
    // Difference first: a constant band-energy gauge cannot generate a field.
    return geometry_[cell].gradients[1] * (value(c.nodes[1]) - value(c.nodes[0])) +
           geometry_[cell].gradients[2] * (value(c.nodes[2]) - value(c.nodes[0]));
}
bool NonlocalPathGeometry::enters(Index cell, const Point2& p, const Point2& direction) const {
    const auto w = weights(cell, p);
    for (Index i = 0; i < 3; ++i) {
        if (w[i] < -tolerance || w[i] > 1 + tolerance) return false;
        if (w[i] <= tolerance && geometry_[cell].gradients[i].dot(direction) <
            -tolerance * norm(geometry_[cell].gradients[i])) return false;
    }
    return true;
}
void NonlocalPathGeometry::validate(const TunnelingBandState& s,
    const NonlocalPathBandToBandConfig& c, const TunnelingGeometryOptions& options) const {
    s.validate(); c.validate();
    demand(s.identity.mesh == mesh_.revision, "tunneling geometry mesh revision mismatch");
    demand(s.nodes.size() == mesh_.positions_m.size(), "tunneling geometry/state node count mismatch");
    demand(options.maximumSegments > 0, "positive tunneling traversal budget required");
    for (Index region : c.regions) demand(std::any_of(mesh_.cells.begin(), mesh_.cells.end(), [&](const auto& cell) {
        return cell.semiconductor && cell.material == mesh_.siliconMaterial && cell.region == region;
    }), "selected tunneling region is absent or is not Si");
    for (const auto& cell : mesh_.cells) if (cell.semiconductor && cell.material == mesh_.siliconMaterial) {
        for (Index node : cell.nodes) {
            const auto& n = s.nodes[node];
            demand(n.semiconductor && n.material == cell.material && (n.region == cell.region ||
                std::find(n.additionalRegions.begin(), n.additionalRegions.end(), cell.region) != n.additionalRegions.end()),
                "Si cell and actual band-state membership differ");
        }
    }
}
TunnelingGeometryResult NonlocalPathGeometry::traceFromPoint(const TunnelingBandState& s,
    const NonlocalPathBandToBandConfig& c, Index cell, const Point2& point,
    const TunnelingGeometryOptions& options) const {
    validate(s, c, options);
    demand(cell < mesh_.cells.size() && point.allFinite(), "invalid tunneling seed");
    const auto& triangle = mesh_.cells[cell];
    demand(triangle.semiconductor && triangle.material == mesh_.siliconMaterial && selected(c, triangle.region),
        "tunneling seed must be in selected Si");
    support(cell, point).validate(s);
    return trace(s, c, cell, point, -gradient(cell, s, true), 0, options);
}
std::vector<TunnelingGeometryResult> NonlocalPathGeometry::traceFromNodes(const TunnelingBandState& s,
    const NonlocalPathBandToBandConfig& c, const TunnelingGeometryOptions& options) const {
    validate(s, c, options);
    std::vector<TunnelingGeometryResult> results;
    for (Index node = 0; node < incident_.size(); ++node) {
        Index first = mesh_.cells.size();
        Point2 g = Point2::Zero();
        Real area = 0;
        for (Index k : incident_[node]) {
            const auto& cell = mesh_.cells[k];
            if (!cell.semiconductor || cell.material != mesh_.siliconMaterial) continue;
            if (selected(c, cell.region)) first = std::min(first, k);
            if (!options.withinSelectedRegions || selected(c, cell.region)) {
                g += gradient(k, s, true) * geometry_[k].twiceArea_m2;
                area += geometry_[k].twiceArea_m2;
            }
        }
        if (first == mesh_.cells.size()) continue;
        // Deterministic area-weighted P1 recovery, not claimed as native recovery.
        results.push_back(trace(s, c, first, mesh_.positions_m[node], -g / area, node, options));
    }
    return results;
}
TunnelingGeometryResult NonlocalPathGeometry::trace(const TunnelingBandState& s,
    const NonlocalPathBandToBandConfig& c, Index cell, const Point2& start,
    Point2 direction, std::uint64_t id, const TunnelingGeometryOptions& options) const {
    TunnelingGeometryResult result;
    auto& path = result.path;
    path.id = id; path.identity = s.identity;
    path.points.push_back({start, support(cell, start), false});
    long double startValence = 0;
    for (const auto& [node, w] : path.points[0].support.weights) startValence += static_cast<long double>(w) * s.nodes[node].valence_J;
    const long double energy = startValence - c.valenceOffset_J;
    path.energyLower_J = path.energyUpper_J = static_cast<Real>(energy);
    auto finish = [&](TunnelingTermination reason) {
        path.termination = reason;
        path.validate(s, c.maximumPathLength_m);
        return result;
    };
    const Real magnitude = norm(direction);
    result.initialValenceGradient_J_per_m = magnitude;
    if (!std::isfinite(magnitude)) return finish(TunnelingTermination::NumericalFailure);
    if (magnitude == 0) return finish(TunnelingTermination::ZeroDirection);
    direction /= magnitude;
    Point2 point = start;
    Real travelled = 0;
    for (Index iteration = 0; iteration < options.maximumSegments; ++iteration) {
        const auto& triangle = mesh_.cells[cell];
        const auto w = weights(cell, point);
        long double barrier = c.conductionOffset_J - energy;
        for (Index i = 0; i < 3; ++i) barrier += static_cast<long double>(w[i]) * s.nodes[triangle.nodes[i]].conduction_J;
        if (barrier <= 0) return finish(travelled == 0 ? TunnelingTermination::NoBarrier : TunnelingTermination::Complete);
        const Real slope = gradient(cell, s, false).dot(direction);
        Real distance = std::numeric_limits<Real>::infinity();
        std::array<Real, 3> exits;
        for (Index i = 0; i < 3; ++i) {
            const Real dw = geometry_[cell].gradients[i].dot(direction);
            exits[i] = dw < -tolerance * norm(geometry_[cell].gradients[i]) ? (w[i] <= tolerance ? Real(0) : w[i] / -dw) :
                std::numeric_limits<Real>::infinity();
            distance = std::min(distance, exits[i]);
        }
        if (!std::isfinite(distance)) return finish(TunnelingTermination::NumericalFailure);
        const Real remaining = std::max(Real(0), c.maximumPathLength_m - travelled);
        const Real root = slope < 0 ? static_cast<Real>(-barrier / slope) : std::numeric_limits<Real>::infinity();
        const bool complete = root >= 0 && root <= distance * (1 + tolerance) && root <= remaining * (1 + tolerance);
        const bool lengthLimit = !complete && remaining <= distance;
        const Real step = complete ? std::min({root, distance, remaining}) : std::min(distance, remaining);
        if (step > 0) {
            Point2 next = point + step * direction;
            // A crossing classified as a vertex must use that vertex's stored
            // coordinates. Repeated ray additions otherwise accumulate normal
            // drift along a mesh line. This is the existing 128-epsilon local
            // geometric classification, not a spatial probe or wider support
            // tolerance; roots and length truncations retain their own points.
            if (!complete && !lengthLimit) {
                const Real h = std::max({norm(mesh_.positions_m[triangle.nodes[1]]-mesh_.positions_m[triangle.nodes[0]]),
                    norm(mesh_.positions_m[triangle.nodes[2]]-mesh_.positions_m[triangle.nodes[0]]),
                    norm(mesh_.positions_m[triangle.nodes[2]]-mesh_.positions_m[triangle.nodes[1]])});
                for (Index vertex:triangle.nodes) if (norm(next-mesh_.positions_m[vertex])<=tolerance*h) {
                    next=mesh_.positions_m[vertex];break;
                }
            }
            if (!next.allFinite() || norm(next - point) == 0) return finish(TunnelingTermination::NumericalFailure);
            const Real arc = (next - point).norm();
            if (travelled + arc == travelled) return finish(TunnelingTermination::NumericalFailure);
            auto a = support(cell, point), b = support(cell, next);
            result.segments.push_back({a, b, travelled, travelled + arc});
            path.integrationIntervals_m.emplace_back(travelled, travelled + arc);
            path.points.push_back({next, b, false});
            point = next; travelled += arc;
        }
        if (complete) return finish(TunnelingTermination::Complete);
        if (lengthLimit) return finish(TunnelingTermination::LengthLimit);
        // At edges and vertices choose the unique outgoing geometric sector.
        // Along a shared edge, canonical smallest cell ID prevents duplicates.
        std::set<Index> candidates;
        const auto at = weights(cell, point);
        for (Index i = 0; i < 3; ++i) if (at[i] > tolerance)
            candidates.insert(incident_[triangle.nodes[i]].begin(), incident_[triangle.nodes[i]].end());
        Index nextCell = mesh_.cells.size();
        // A ray on an interface belongs to selected Si when that sector exists;
        // an oxide or excluded cell with a smaller ID must not steal the ray.
        for (Index candidate : candidates) if (candidate != cell && enters(candidate, point, direction)) {
            const auto& cc = mesh_.cells[candidate];
            if (cc.semiconductor && cc.material == mesh_.siliconMaterial &&
                (!options.withinSelectedRegions || selected(c, cc.region))) { nextCell = candidate; break; }
        }
        for (Index candidate : candidates) if (candidate != cell && enters(candidate, point, direction)) {
            if (nextCell == mesh_.cells.size()) nextCell = candidate;
            break;
        }
        if (nextCell != mesh_.cells.size() && mesh_.cells[nextCell].semiconductor) {
            const auto& next = mesh_.cells[nextCell];
            if (next.material != mesh_.siliconMaterial) return finish(TunnelingTermination::MaterialMismatch);
            if (options.withinSelectedRegions && !selected(c, next.region)) return finish(TunnelingTermination::RegionExcluded);
            cell = nextCell;
            continue;
        }
        // Include the whole vertex fan: an exterior corner may be formed by
        // edges belonging to different triangles. Collinear subdivisions of a
        // planar boundary have one normal and are not ambiguous corners.
        candidates.insert(cell);
        std::map<EdgeKey, std::pair<Point2, TunnelingBoundaryKind>> blocked;
        for (Index k : candidates) {
            const auto& cc = mesh_.cells[k];
            if (!cc.semiconductor || cc.material != mesh_.siliconMaterial) continue;
            const auto ww = weights(k, point);
            if (std::any_of(ww.begin(), ww.end(), [](Real v) { return v < -tolerance || v > 1 + tolerance; })) continue;
            for (Index i = 0; i < 3; ++i) {
                const Point2 normal = geometry_[k].gradients[i] / norm(geometry_[k].gradients[i]);
                if (std::abs(ww[i]) > tolerance || normal.dot(direction) >= -tolerance) continue;
                const EdgeKey key = edge(cc.nodes[(i + 1) % 3], cc.nodes[(i + 2) % 3]);
                const auto& adjacent = edges_.at(key);
                if (adjacent.size() == 1) blocked.emplace(key, std::make_pair(normal, boundaries_.at(key)));
                else {
                    const Index other = adjacent[0] == k ? adjacent[1] : adjacent[0];
                    if (!mesh_.cells[other].semiconductor)
                        blocked.emplace(key, std::make_pair(normal, TunnelingBoundaryKind::Reflecting));
                }
            }
        }
        if (blocked.empty()) return finish(TunnelingTermination::NumericalFailure);
        const auto [normal, kind] = blocked.begin()->second;
        for (const auto& [key, boundary] : blocked) {
            (void)key;
            if (norm(boundary.first - normal) > tolerance || boundary.second != kind)
                return finish(TunnelingTermination::AmbiguousCorner);
        }
        if (kind == TunnelingBoundaryKind::Contact) return finish(TunnelingTermination::ContactBoundary);
        if (kind == TunnelingBoundaryKind::Open) return finish(TunnelingTermination::DomainBoundary);
        if (result.reflections >= options.maximumReflections) return finish(TunnelingTermination::ReflectionLimit);
        direction -= 2 * direction.dot(normal) * normal;
        direction /= norm(direction);
        path.points.back().reflection = true;
        ++result.reflections;
    }
    return finish(TunnelingTermination::TraversalLimit);
}
std::vector<TunnelingGeometryResult> NonlocalPathGeometry::traceFromCells(
    const TunnelingBandState& s, const NonlocalPathBandToBandConfig& c,
    const TunnelingGeometryOptions& options) const {
    validate(s,c,options);
    demand(c.valenceOffset_J==0 && c.conductionOffset_J==0,"cell diagnostic requires zero band offsets");
    Real globalEc=std::numeric_limits<Real>::infinity();
    for (const auto& n:s.nodes) if(n.semiconductor && n.material==mesh_.siliconMaterial)
        globalEc=std::min(globalEc,n.conduction_J);
    std::vector<TunnelingGeometryResult> results;
    for(Index seed=0;seed<mesh_.cells.size();++seed) {
        const auto& sc=mesh_.cells[seed];
        if(!sc.semiconductor || sc.material!=mesh_.siliconMaterial || !selected(c,sc.region)) continue;
        TunnelingGeometryResult r;auto& p=r.path;p.id=seed;p.identity=s.identity;
        Point2 center=Point2::Zero();for(Index n:sc.nodes) center+=mesh_.positions_m[n]/3;
        Point2 direction=-gradient(seed,s,true);const Real magnitude=norm(direction);
        r.initialValenceGradient_J_per_m=magnitude;
        auto finish=[&](TunnelingTermination reason) {
            p.termination=reason;p.validate(s,c.maximumPathLength_m);results.push_back(std::move(r));
        };
        if(magnitude==0) {p.points.push_back({center,support(seed,center),false});finish(TunnelingTermination::ZeroDirection);continue;}
        demand(std::isfinite(magnitude),"nonfinite cell gradient");direction/=magnitude;
        Real lo=-std::numeric_limits<Real>::infinity(),hi=-lo;Index lf=0,face=0;
        const auto cw=weights(seed,center);
        for(Index i=0;i<3;++i) {
            const Real slope=geometry_[seed].gradients[i].dot(direction);
            if(slope>0 && -cw[i]/slope>lo) {lo=-cw[i]/slope;lf=i;}
            if(slope<0 && -cw[i]/slope<hi) {hi=-cw[i]/slope;face=i;}
        }
        // Project onto the original face; repeated ray additions must not drift
        // out of barycentric support. No spatial probing or tolerance widening.
        auto onFace=[&](Index cell,Index f,const Point2& pos) {
            const auto& ns=mesh_.cells[cell].nodes;Index a=ns[(f+1)%3],b=ns[(f+2)%3];
            const Point2 d=mesh_.positions_m[b]-mesh_.positions_m[a];
            const Real w=(pos-mesh_.positions_m[a]).dot(d)/d.squaredNorm();
            demand(std::isfinite(w) && w>=-tolerance && w<=1+tolerance,"cell face support out of bounds");
            const Real v=std::clamp(w,Real(0),Real(1));
            return TunnelingPathPoint{(1-v)*mesh_.positions_m[a]+v*mesh_.positions_m[b],
                {cell,mesh_.cells[cell].material,mesh_.cells[cell].region,{{a,1-v},{b,v}}},false};
        };
        demand(lo<0 && hi>0 && std::isfinite(hi-lo),"invalid centroid chord");
        r.sourceChord_m=hi-lo;
        p.points.push_back(onFace(seed,lf,center+lo*direction));
        if(r.sourceChord_m>c.maximumPathLength_m) {
            const auto exit=onFace(seed,face,center+hi*direction);
            long double upper=0,entryEc=0,exitEc=0;
            for(const auto& [n,w]:p.points.front().support.weights) {
                upper+=static_cast<long double>(w)*s.nodes[n].valence_J;
                entryEc+=static_cast<long double>(w)*s.nodes[n].conduction_J;
            }
            for(const auto& [n,w]:exit.support.weights) exitEc+=static_cast<long double>(w)*s.nodes[n].conduction_J;
            // Ec is P1 on the entire source chord. If both endpoints exceed
            // its highest available energy, no conduction turn can lie in
            // the shorter retained domain. Otherwise the chord is unresolved.
            if(std::min(entryEc,exitEc)>=upper) {
                p.energyLower_J=p.energyUpper_J=static_cast<Real>(upper);
                finish(TunnelingTermination::NoBarrier);
            } else finish(TunnelingTermination::LengthLimit);
            continue;
        }
        p.points.push_back(onFace(seed,face,center+hi*direction));
        Real travelled=(p.points[1].position_m-p.points[0].position_m).norm();
        p.integrationIntervals_m.emplace_back(0,travelled);
        r.segments.push_back({p.points[0].support,p.points[1].support,0,travelled});
        auto band=[&](const TunnelingEndpoint& e,bool ev) {
            long double v=0;for(const auto& [n,w]:e.weights) v+=static_cast<long double>(w)*(ev?s.nodes[n].valence_J:s.nodes[n].conduction_J);
            return static_cast<Real>(v);
        };
        const Real lower=std::max(globalEc,band(p.points.back().support,true));
        const Real upper=band(p.points.front().support,true);
        Index cell=seed;TunnelingTermination termination=TunnelingTermination::TraversalLimit;
        for(Index iteration=0;iteration<options.maximumSegments;++iteration) {
            if(lower>=upper) {termination=TunnelingTermination::NoBarrier;break;}
            if(band(p.points.back().support,false)<=lower) {termination=TunnelingTermination::Complete;break;}
            const auto& ns=mesh_.cells[cell].nodes;
            const auto key=edge(ns[(face+1)%3],ns[(face+2)%3]);const auto& adjacent=edges_.at(key);
            if(adjacent.size()==2) {
                const Index other=adjacent[0]==cell?adjacent[1]:adjacent[0];const auto& oc=mesh_.cells[other];
                if(oc.semiconductor) {
                    if(oc.material!=mesh_.siliconMaterial) {termination=TunnelingTermination::MaterialMismatch;break;}
                    if(options.withinSelectedRegions && !selected(c,oc.region)) {termination=TunnelingTermination::RegionExcluded;break;}
                    cell=other;
                } else {
                    if(r.reflections>=options.maximumReflections) {termination=TunnelingTermination::ReflectionLimit;break;}
                    const Point2 normal=geometry_[cell].gradients[face].normalized();
                    direction-=2*direction.dot(normal)*normal;++r.reflections;p.points.back().reflection=true;
                }
            } else {
                const auto kind=boundaries_.at(key);
                if(kind==TunnelingBoundaryKind::Contact) {termination=TunnelingTermination::ContactBoundary;break;}
                if(kind==TunnelingBoundaryKind::Open) {termination=TunnelingTermination::DomainBoundary;break;}
                if(r.reflections>=options.maximumReflections) {termination=TunnelingTermination::ReflectionLimit;break;}
                const Point2 normal=geometry_[cell].gradients[face].normalized();
                direction-=2*direction.dot(normal)*normal;++r.reflections;p.points.back().reflection=true;
            }
            const Point2 point=p.points.back().position_m;const auto w=weights(cell,point);
            Real step=std::numeric_limits<Real>::infinity();Index ties=0;
            for(Index i=0;i<3;++i) {
                const Real slope=geometry_[cell].gradients[i].dot(direction);
                if(slope>=-tolerance*norm(geometry_[cell].gradients[i])) continue;
                const Real candidate=std::max(Real(0),w[i])/-slope;
                if(std::isfinite(step) && std::abs(candidate-step)<=tolerance*std::max(candidate,step)) ++ties;
                else if(candidate<step) {step=candidate;face=i;ties=1;}
            }
            if(ties>1) {termination=TunnelingTermination::AmbiguousCorner;break;}
            if(!std::isfinite(step) || step<=0 || travelled+step==travelled) {termination=TunnelingTermination::NumericalFailure;break;}
            if(travelled+step>c.maximumPathLength_m) {termination=TunnelingTermination::LengthLimit;break;}
            auto next=onFace(cell,face,point+step*direction);
            const Real arc=(next.position_m-point).norm();
            r.segments.push_back({support(cell,point),next.support,travelled,travelled+arc});
            p.integrationIntervals_m.emplace_back(travelled,travelled+arc);
            p.points.push_back(std::move(next));travelled+=arc;
        }
        p.energyUpper_J=upper;p.energyLower_J=std::max(lower,band(p.points[0].support,false));
        // A reflected final point need not be the lowest conduction band.
        Real minimumEc=std::numeric_limits<Real>::infinity();
        for(const auto& point:p.points) minimumEc=std::min(minimumEc,band(point.support,false));
        p.energyLower_J=std::min(upper,std::max(lower,minimumEc));
        finish(termination);
    }
    return results;
}
const char* tunnelingTerminationName(TunnelingTermination t) {
    switch (t) {
    case TunnelingTermination::Complete: return "complete";
    case TunnelingTermination::LengthLimit: return "length_limit";
    case TunnelingTermination::InsulatorBoundary: return "insulator_boundary";
    case TunnelingTermination::ZeroDirection: return "zero_direction";
    case TunnelingTermination::ReflectionLimit: return "reflection_limit";
    case TunnelingTermination::NumericalFailure: return "numerical_failure";
    case TunnelingTermination::DomainBoundary: return "domain_boundary";
    case TunnelingTermination::ContactBoundary: return "contact_boundary";
    case TunnelingTermination::RegionExcluded: return "region_excluded";
    case TunnelingTermination::MaterialMismatch: return "material_mismatch";
    case TunnelingTermination::TraversalLimit: return "traversal_limit";
    case TunnelingTermination::AmbiguousCorner: return "ambiguous_corner";
    case TunnelingTermination::NoBarrier: return "no_barrier";
    }
    throw std::invalid_argument("unknown tunneling termination");
}
nlohmann::json tunnelingGeometryAudit(const TunnelingGeometryResult& r) {
    nlohmann::json points = nlohmann::json::array(), segments = nlohmann::json::array();
    for (const auto& p : r.path.points) points.push_back({{"position_m", {p.position_m.x(), p.position_m.y()}},
        {"support", endpointJson(p.support)}, {"reflection", p.reflection}});
    for (const auto& s : r.segments) segments.push_back({{"start", endpointJson(s.start)}, {"end", endpointJson(s.end)},
        {"begin_m", s.begin_m}, {"end_m", s.end_m}});
    const auto& i = r.path.identity;
    return {{"schema", "vela.experimental.nonlocalpath.geometry.v1"}, {"path_id", r.path.id},
        {"termination", tunnelingTerminationName(r.path.termination)}, {"reflections", r.reflections},
        {"initial_valence_gradient_J_per_m", r.initialValenceGradient_J_per_m},
        {"geometrically_complete", r.path.termination == TunnelingTermination::Complete},
        {"requires_boundary_resolution", r.path.termination == TunnelingTermination::AmbiguousCorner},
        {"resource_failure", r.path.termination == TunnelingTermination::NumericalFailure ||
            r.path.termination == TunnelingTermination::TraversalLimit || r.path.termination == TunnelingTermination::ReflectionLimit},
        {"state_identity", {i.mesh, i.material, i.parameters, i.contacts, i.population, i.temperature, i.trial}},
        {"energy_J", r.path.energyLower_J}, {"points", points}, {"segments", segments}};
}
} // namespace vela::experimental
