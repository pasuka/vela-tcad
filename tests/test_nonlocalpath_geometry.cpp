#include "vela/physics/NonlocalPathGeometry.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

using namespace vela;
using namespace vela::experimental;
using Catch::Approx;
namespace {
struct Fixture { TunnelingGeometryMesh mesh; TunnelingBandState state; NonlocalPathBandToBandConfig config; };
Fixture rectangle(Index nx = 8, Index ny = 4, Real width = 4, Real height = 1,
                  Real gap = .7, Real gy = 0, Real quadratic = 0, Real angle = 0,
                  Real lengthScale = 1, Real energyScale = 1, Real gauge = 0) {
    Fixture f;
    f.config.regions = {1};
    f.config.ab = IndirectABParameters{4e20, 1.9e9, 1.6e-19};
    f.config.phononEnergy_J = .037 * 1.6e-19;
    f.config.maximumPathLength_m = 20 * lengthScale;
    Eigen::Matrix2d rotation;
    rotation << std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle);
    for (Index j = 0; j <= ny; ++j) for (Index i = 0; i <= nx; ++i) {
        const Real x = width * i / nx, y = height * j / ny;
        f.mesh.positions_m.push_back(rotation * Point2(x, y) * lengthScale);
        TunnelingBandNode n;
        n.region = 1;
        n.valence_J = energyScale * (-x - gy * y - quadratic * x * x) + gauge;
        n.conduction_J = n.valence_J + energyScale * gap;
        n.conductionDOS_per_m3 = n.valenceDOS_per_m3 = 1e25;
        f.state.nodes.push_back(n);
    }
    for (Index j = 0; j < ny; ++j) for (Index i = 0; i < nx; ++i) {
        const Index a = j * (nx + 1) + i, b = a + 1, d = a + nx + 1, e = d + 1;
        f.mesh.cells.push_back({{a, b, e}, 0, 1, true});
        f.mesh.cells.push_back({{a, e, d}, 0, 1, true});
    }
    std::map<std::pair<Index, Index>, Index> edges;
    for (const auto& c : f.mesh.cells) for (Index k = 0; k < 3; ++k) ++edges[std::minmax(c.nodes[k], c.nodes[(k + 1) % 3])];
    for (const auto& [key, count] : edges) if (count == 1)
        f.mesh.boundaries.push_back({{key.first, key.second}, TunnelingBoundaryKind::Reflecting});
    return f;
}
Real length(const TunnelingGeometryResult& r) { return r.segments.empty() ? 0 : r.segments.back().end_m; }
Point2 end(const TunnelingGeometryResult& r) { return r.path.points.back().position_m; }
void consistent(const TunnelingGeometryResult& r, const Fixture& f) {
    REQUIRE_NOTHROW(r.path.validate(f.state, f.config.maximumPathLength_m));
    Real position = 0;
    for (const auto& segment : r.segments) {
        REQUIRE(segment.begin_m == position);
        REQUIRE(segment.end_m > segment.begin_m);
        REQUIRE(segment.start.cell == segment.end.cell);
        REQUIRE_NOTHROW(segment.start.validate(f.state));
        REQUIRE_NOTHROW(segment.end.validate(f.state));
        position = segment.end_m;
    }
}
void assignRegions(Fixture& f, Real boundaryX) {
    for (auto& c : f.mesh.cells) {
        Real x = 0; for (Index n : c.nodes) x += f.mesh.positions_m[n].x() / 3;
        c.region = x < boundaryX ? 1 : 2;
    }
    std::vector<std::set<Index>> memberships(f.state.nodes.size());
    for (const auto& c : f.mesh.cells) for (Index n : c.nodes) memberships[n].insert(c.region);
    for (Index n = 0; n < memberships.size(); ++n) {
        auto& node = f.state.nodes[n];
        node.region = *memberships[n].begin();
        for (Index r : memberships[n]) if (r != node.region) node.additionalRegions.push_back(r);
    }
}
}
TEST_CASE("NLP N2 linear bands give analytic endpoints and conservative crossing intervals", "[nonlocalpath][geometry]") {
    auto f = rectangle();
    const auto result = NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.1, .03));
    REQUIRE(result.path.termination == TunnelingTermination::Complete);
    REQUIRE(end(result).x() == Approx(.8).margin(1e-13));
    REQUIRE(end(result).y() == Approx(.03).margin(1e-13));
    REQUIRE(length(result) == Approx(.7).epsilon(1e-13));
    REQUIRE(result.segments.size() > 1);
    consistent(result, f);
}
TEST_CASE("NLP N2 geometry is invariant under energy gauge rotation orientation and SI rescaling", "[nonlocalpath][geometry]") {
    for (Real angle : {0., .37, 1.91, -1.2}) for (bool reverse : {false, true}) {
        auto f = rectangle(8, 4, 4, 1, .7, .2, 0, angle, 1e-6, 1.602176634e-19, 3.1e-19);
        if (reverse) for (auto& c : f.mesh.cells) std::swap(c.nodes[0], c.nodes[2]);
        Eigen::Matrix2d rotation;
        rotation << std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle);
        const Point2 start = rotation * Point2(.1, .03) * 1e-6;
        const auto result = NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, start);
        REQUIRE(result.path.termination == TunnelingTermination::Complete);
        const Point2 expected = start + rotation * Point2(1, .2) * (.7e-6 / 1.04);
        REQUIRE((end(result) - expected).norm() < 2e-19);
        REQUIRE(length(result) / 1e-6 == Approx(.7 / std::sqrt(1.04)).epsilon(2e-13));
        consistent(result, f);
    }
}
TEST_CASE("NLP N2 offsets shift endpoint energy without rotating the starting direction", "[nonlocalpath][geometry]") {
    auto f = rectangle();
    f.config.conductionOffset_J = .1;
    f.config.valenceOffset_J = .2;
    const auto r = NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.1, .03));
    REQUIRE(r.path.termination == TunnelingTermination::Complete);
    REQUIRE(end(r).x() == Approx(1.1).epsilon(1e-13));
    REQUIRE(end(r).y() == Approx(.03));
    REQUIRE(r.path.energyLower_J == Approx(-.3));
}
TEST_CASE("NLP N2 planar Neumann reflection satisfies analytic specular path", "[nonlocalpath][geometry]") {
    auto f = rectangle(8, 4, 4, 1, 2.3, .5);
    const auto r = NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 1, Point2(.2, .2));
    REQUIRE(r.path.termination == TunnelingTermination::Complete);
    REQUIRE(r.reflections == 1);
    REQUIRE(end(r).x() == Approx(2.2).epsilon(2e-13));
    REQUIRE(end(r).y() == Approx(.8).epsilon(2e-13));
    REQUIRE(length(r) == Approx(std::sqrt(5.)).epsilon(2e-13));
    consistent(r, f);
    TunnelingGeometryOptions options; options.maximumReflections = 0;
    REQUIRE(NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 1, Point2(.2, .2), options).path.termination ==
        TunnelingTermination::ReflectionLimit);
}
TEST_CASE("NLP N2 oxide interface reflects and never acquires carrier endpoint support", "[nonlocalpath][geometry]") {
    auto f = rectangle(8, 8, 4, 2, 2.3, .5);
    for (auto& c : f.mesh.cells) {
        Real y = 0; for (Index n : c.nodes) y += f.mesh.positions_m[n].y() / 3;
        if (y > 1) { c.semiconductor = false; c.material = 1; }
    }
    const auto r = NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 1, Point2(.2, .2));
    REQUIRE(r.path.termination == TunnelingTermination::Complete);
    REQUIRE(r.reflections == 1);
    REQUIRE(end(r).x() == Approx(2.2).epsilon(2e-13));
    REQUIRE(end(r).y() == Approx(.8).epsilon(2e-13));
    for (const auto& segment : r.segments) REQUIRE(f.mesh.cells[segment.start.cell].semiconductor);
    consistent(r, f);
}
TEST_CASE("NLP N2 selected region traversal and material restrictions are explicit", "[nonlocalpath][geometry]") {
    auto f = rectangle(8, 4, 4, 1, 2.3);
    assignRegions(f, 1);
    const NonlocalPathGeometry geometry(f.mesh);
    auto r = geometry.traceFromPoint(f.state, f.config, 0, Point2(.1, .03));
    REQUIRE(r.path.termination == TunnelingTermination::RegionExcluded);
    REQUIRE(end(r).x() == Approx(1));
    TunnelingGeometryOptions options; options.withinSelectedRegions = false;
    r = geometry.traceFromPoint(f.state, f.config, 0, Point2(.1, .03), options);
    REQUIRE(r.path.termination == TunnelingTermination::Complete);
    REQUIRE(end(r).x() == Approx(2.4));
    consistent(r, f);
    for (auto& c : f.mesh.cells) if (c.region == 2) c.material = 1;
    REQUIRE(NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.1, .03), options).path.termination ==
        TunnelingTermination::MaterialMismatch);
}
TEST_CASE("NLP N2 mirror symmetry and resolvable near-grazing directions preserve endpoints", "[nonlocalpath][geometry]") {
    auto f = rectangle(8, 4, 4, 1, 2.3, .5);
    for (auto& p : f.mesh.positions_m) p.y() = -p.y();
    const auto reflected = NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 1, Point2(.2, -.2));
    REQUIRE(reflected.path.termination == TunnelingTermination::Complete);
    REQUIRE(reflected.reflections == 1);
    REQUIRE(end(reflected).x() == Approx(2.2).epsilon(2e-13));
    REQUIRE(end(reflected).y() == Approx(-.8).epsilon(2e-13));
    for (Real gy : {1e-8, -1e-8}) {
        f = rectangle(8, 4, 4, 1, .7, gy);
        const auto r = NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.1, .03));
        REQUIRE(r.path.termination == TunnelingTermination::Complete);
        REQUIRE(end(r).x() == Approx(.8).epsilon(2e-13));
        REQUIRE(end(r).y() == Approx(.03 + .7 * gy).margin(1e-14));
    }
}
TEST_CASE("NLP N2 endpoint at region boundary takes precedence and absent barrier is explicit", "[nonlocalpath][geometry]") {
    auto f = rectangle(8, 4, 4, 1, .9);
    assignRegions(f, 1);
    const auto r = NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.1, .03));
    REQUIRE(r.path.termination == TunnelingTermination::Complete);
    REQUIRE(end(r).x() == Approx(1).epsilon(2e-13));
    f.config.conductionOffset_J = -1;
    REQUIRE(NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.1, .03)).path.termination ==
        TunnelingTermination::NoBarrier);
}
TEST_CASE("NLP N2 grazing selected Si sector wins over lower numbered excluded cells", "[nonlocalpath][geometry]") {
    auto f = rectangle(4, 4, 2, 2, .7);
    assignRegions(f, 1);
    f.config.regions = {2};
    for (Index n = 0; n < f.state.nodes.size(); ++n) {
        f.state.nodes[n].valence_J = -f.mesh.positions_m[n].y();
        f.state.nodes[n].conduction_J = f.state.nodes[n].valence_J + .7;
    }
    const auto paths = NonlocalPathGeometry(f.mesh).traceFromNodes(f.state, f.config);
    const auto seed = std::find_if(paths.begin(), paths.end(), [](const auto& r) { return r.path.id == 2; });
    REQUIRE(seed != paths.end()); // x=1, y=0, on selected/excluded interface.
    REQUIRE(seed->path.termination == TunnelingTermination::Complete);
    REQUIRE(end(*seed).x() == Approx(1));
    REQUIRE(end(*seed).y() == Approx(.7));
    for (const auto& segment : seed->segments) REQUIRE(segment.start.region == 2);
    for (auto& c : f.mesh.cells) if (c.region == 1) { c.semiconductor = false; c.material = 1; }
    const auto oxidePaths = NonlocalPathGeometry(f.mesh).traceFromNodes(f.state, f.config);
    const auto oxideSeed = std::find_if(oxidePaths.begin(), oxidePaths.end(), [](const auto& r) { return r.path.id == 2; });
    REQUIRE(oxideSeed != oxidePaths.end());
    REQUIRE(oxideSeed->path.termination == TunnelingTermination::Complete);
    REQUIRE(oxideSeed->reflections == 0); // Tangential travel does not reflect.
}
TEST_CASE("NLP N2 length and traversal resource limits retain failed path evidence", "[nonlocalpath][geometry]") {
    auto f = rectangle();
    f.config.maximumPathLength_m = .3;
    const auto r = NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.1, .03));
    REQUIRE(r.path.termination == TunnelingTermination::LengthLimit);
    REQUIRE(length(r) == Approx(.3));
    REQUIRE(end(r).x() == Approx(.4));
    const auto audit = tunnelingGeometryAudit(r);
    REQUIRE(audit.at("termination") == "length_limit");
    REQUIRE_FALSE(audit.at("geometrically_complete").get<bool>());
    REQUIRE_FALSE(audit.at("resource_failure").get<bool>());
    REQUIRE(audit.at("segments").size() == r.segments.size());
    f.config.maximumPathLength_m = 20;
    TunnelingGeometryOptions options; options.maximumSegments = 1;
    REQUIRE(NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.1, .03), options).path.termination ==
        TunnelingTermination::TraversalLimit);
    f.config.maximumPathLength_m = .7;
    REQUIRE(NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.1, .03)).path.termination ==
        TunnelingTermination::Complete);
}
TEST_CASE("NLP N2 zero direction contacts open boundaries and corners have deterministic reasons", "[nonlocalpath][geometry]") {
    auto f = rectangle(1, 1, 1, 1, 10);
    for (auto& n : f.state.nodes) { n.valence_J = 0; n.conduction_J = 1; }
    auto r = NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.2, .1));
    REQUIRE(r.path.termination == TunnelingTermination::ZeroDirection);
    REQUIRE(r.path.points.size() == 1);
    consistent(r, f);
    f = rectangle(1, 1, 1, 1, 10);
    for (auto& b : f.mesh.boundaries) b.kind = TunnelingBoundaryKind::Contact;
    REQUIRE(NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.2, .1)).path.termination ==
        TunnelingTermination::ContactBoundary);
    for (auto& b : f.mesh.boundaries) b.kind = TunnelingBoundaryKind::Open;
    REQUIRE(NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.2, .1)).path.termination ==
        TunnelingTermination::DomainBoundary);
    f = rectangle(1, 1, 1, 1, 10, 1);
    REQUIRE(NonlocalPathGeometry(f.mesh).traceFromPoint(f.state, f.config, 0, Point2(.2, .2)).path.termination ==
        TunnelingTermination::AmbiguousCorner);
}
TEST_CASE("NLP N2 automatic vertex seeds are unique through grazing edges and mesh vertices", "[nonlocalpath][geometry]") {
    auto f = rectangle(12, 6);
    const auto paths = NonlocalPathGeometry(f.mesh).traceFromNodes(f.state, f.config);
    REQUIRE(paths.size() == f.mesh.positions_m.size());
    std::set<std::uint64_t> ids;
    Index complete = 0;
    for (const auto& r : paths) {
        REQUIRE(ids.insert(r.path.id).second);
        REQUIRE(r.path.termination != TunnelingTermination::TraversalLimit);
        REQUIRE(r.path.termination != TunnelingTermination::NumericalFailure);
        if (r.path.termination == TunnelingTermination::Complete) {
            ++complete;
            REQUIRE(length(r) == Approx(.7).epsilon(2e-13));
        }
        consistent(r, f);
    }
    REQUIRE(complete > paths.size() / 2);
}
TEST_CASE("NLP N2 curved band endpoint error converges under nested mesh refinement", "[nonlocalpath][geometry]") {
    const Real exact = (std::sqrt(1 + .8 * .7) - 1) / .4;
    Real previous = 1;
    for (Index nx : {4, 8, 16, 32, 64}) {
        auto f = rectangle(nx, 2, 2, 1, .7, 0, .2);
        const auto paths = NonlocalPathGeometry(f.mesh).traceFromNodes(f.state, f.config);
        const auto& r = paths[nx + 1]; // x=0,y=0.5; gradient is horizontal.
        REQUIRE(r.path.termination == TunnelingTermination::Complete);
        const Real error = std::abs(end(r).x() - exact);
        REQUIRE(error < previous);
        REQUIRE(error <= .2 * std::pow(2. / nx, 2) / 4);
        previous = error;
    }
}
TEST_CASE("NLP N2 rejects invalid meshes stale identity and inconsistent memberships", "[nonlocalpath][geometry]") {
    auto f = rectangle();
    auto bad = f.mesh; bad.boundaries.pop_back(); REQUIRE_THROWS(NonlocalPathGeometry(bad));
    bad = f.mesh; bad.cells[0].nodes[1] = bad.cells[0].nodes[0]; REQUIRE_THROWS(NonlocalPathGeometry(bad));
    bad = f.mesh; bad.cells.push_back(bad.cells[0]); REQUIRE_THROWS(NonlocalPathGeometry(bad));
    bad = f.mesh; bad.boundaries.push_back(bad.boundaries[0]); REQUIRE_THROWS(NonlocalPathGeometry(bad));
    bad = f.mesh;
    const auto v = bad.cells[0].nodes;
    bad.positions_m[v[2]] = .5 * (bad.positions_m[v[0]] + bad.positions_m[v[1]]);
    REQUIRE_THROWS(NonlocalPathGeometry(bad));
    TunnelingGeometryMesh overlap;
    overlap.positions_m = {Point2(0,0), Point2(1,0), Point2(.2,1), Point2(.8,1)};
    overlap.cells = {{{0,1,2}}, {{0,1,3}}};
    overlap.boundaries = {{{0,2}}, {{1,2}}, {{0,3}}, {{1,3}}};
    REQUIRE_THROWS(NonlocalPathGeometry(overlap));
    const NonlocalPathGeometry geometry(f.mesh);
    auto state = f.state; ++state.identity.mesh;
    REQUIRE_THROWS(geometry.traceFromNodes(state, f.config));
    state = f.state; state.nodes[0].region = 7;
    REQUIRE_THROWS(geometry.traceFromNodes(state, f.config));
    REQUIRE_THROWS(geometry.traceFromPoint(f.state, f.config, 0, Point2(-1, -1)));
    auto config = f.config; config.regions = {99};
    REQUIRE_THROWS(geometry.traceFromNodes(f.state, config));
}
