#include "vela/physics/NonlocalPathSourceAssembly.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace vela::experimental {
namespace {
void check(bool ok, const char* text) { if (!ok) throw std::invalid_argument(text); }
bool active(const NonlocalPathBandToBandConfig& c, Index r) {
    return std::find(c.regions.begin(), c.regions.end(), r) != c.regions.end();
}
Real sampleFermi(const TunnelingEndpoint& e, const TunnelingBandState& s, bool electron) {
    long double result = 0;
    for (const auto& [n,w] : e.weights)
        result += static_cast<long double>(w) * (electron ? s.nodes[n].electronFermi_J : s.nodes[n].holeFermi_J);
    return static_cast<Real>(result);
}
TunnelingContactLedger contactLedger(const NonlocalPairSources& p, const TunnelingBandState& s) {
    long double ef=0,ec=0,hf=0,hc=0;
    for (Index i=0;i<s.nodes.size();++i) {
        (s.nodes[i].contact ? ec : ef) += p.electronIntegrals_per_m_s[i];
        (s.nodes[i].contact ? hc : hf) += p.holeIntegrals_per_m_s[i];
    }
    return {static_cast<Real>(ef),static_cast<Real>(ec),static_cast<Real>(hf),static_cast<Real>(hc)};
}
}
TunnelingSourceBoxes tunnelingSourceBoxes(const TunnelingGeometryMesh& m,
    const NonlocalPathBandToBandConfig& c, TunnelingBoxPolicy policy) {
    c.validate();
    // Reuse N2's topology/coordinate validation, including nonmanifold edges.
    const NonlocalPathGeometry validated(m);
    check(policy==TunnelingBoxPolicy::Barycentric || policy==TunnelingBoxPolicy::Delaunay, "unknown source box policy");
    struct Part { Index cell, edge; long double coefficient; };
    std::map<std::pair<Index,Index>,std::vector<Part>> edges;
    std::vector<std::array<long double,3>> coefficients(m.cells.size());
    for (Index k=0;k<m.cells.size();++k) {
        const auto& ns=m.cells[k].nodes;
        for (Index j=0;j<3;++j) {
            const Point2 u=m.positions_m[ns[j]]-m.positions_m[ns[(j+2)%3]];
            const Point2 v=m.positions_m[ns[(j+1)%3]]-m.positions_m[ns[(j+2)%3]];
            const long double cross=std::abs(static_cast<long double>(u.x())*v.y()-static_cast<long double>(u.y())*v.x());
            const long double dot=static_cast<long double>(u.x())*v.x()+static_cast<long double>(u.y())*v.y();
            edges[std::minmax(ns[j],ns[(j+1)%3])].push_back({k,j,dot/(2*cross)});
        }
    }
    if (policy==TunnelingBoxPolicy::Delaunay) for (const auto& [edge,parts] : edges) {
        long double total=0,scale=1;
        for (const auto& p:parts) { total+=p.coefficient;scale=std::max(scale,std::abs(p.coefficient)); }
        const long double tolerance=128*std::numeric_limits<Real>::epsilon()*scale;
        check(total>=-tolerance,"non-Delaunay source boxes");
        const auto low=std::min_element(parts.begin(),parts.end(),[](const Part& a,const Part& b){return a.coefficient<b.coefficient;});
        const bool same=parts.size()==2 && m.cells[parts[0].cell].region==m.cells[parts[1].cell].region &&
            m.cells[parts[0].cell].material==m.cells[parts[1].cell].material;
        if (low->coefficient<0) {
            check(same || low->coefficient>=-tolerance,"obtuse material/exterior boundary source box");
            for (const auto& p:parts) coefficients[p.cell][p.edge]=same ?
                (p.cell==low->cell ? 0 : std::max(total,0.L)) : std::max(p.coefficient,0.L);
        } else for (const auto& p:parts) coefficients[p.cell][p.edge]=p.coefficient;
    }
    TunnelingSourceBoxes result;
    result.seedAreas_m2.assign(m.positions_m.size(),0);result.fieldAreas_m2=result.seedAreas_m2;
    result.cellSeedAreas_m2.assign(m.cells.size(),0);
    for (Index k=0;k<m.cells.size();++k) {
        const auto& cell=m.cells[k];if (!cell.semiconductor || cell.material!=m.siliconMaterial) continue;
        const auto& ns=cell.nodes;
        const Point2 u=m.positions_m[ns[1]]-m.positions_m[ns[0]],v=m.positions_m[ns[2]]-m.positions_m[ns[0]];
        const long double area=std::abs(static_cast<long double>(u.x())*v.y()-static_cast<long double>(u.y())*v.x())/2;
        std::array<long double,3> products;
        for (Index j=0;j<3;++j) products[j]=coefficients[k][j]*(m.positions_m[ns[(j+1)%3]]-m.positions_m[ns[j]]).squaredNorm();
        for (Index j=0;j<3;++j) {
            const Real a=static_cast<Real>(policy==TunnelingBoxPolicy::Barycentric ? area/3 : (products[j]+products[(j+2)%3])/4);
            check(std::isfinite(a) && a>=0,"invalid source box area");
            result.fieldAreas_m2[ns[j]]+=a;
            if (active(c,cell.region)) result.seedAreas_m2[ns[j]]+=a;
            if (active(c,cell.region)) result.cellSeedAreas_m2[k]+=a;
        }
    }
    return result;
}
TunnelingSourceAssembly assembleNonlocalPathSources(const TunnelingGeometryMesh& m,
    const TunnelingBandState& s, const NonlocalPathBandToBandConfig& c,
    const TunnelingWindowInput& windows, TunnelingBoxPolicy policy, const TunnelingGeometryOptions& options) {
    s.validate();c.validate();
    check(windows.identity==s.identity && windows.nodes.size()==s.nodes.size(),"stale or incomplete source window input");
    for (const auto& w:windows.nodes) check(std::isfinite(w.valenceAvailableEnergy_J) && w.valenceAvailableEnergy_J>=0 &&
        std::isfinite(w.conductionAvailableEnergy_J) && w.conductionAvailableEnergy_J>=0,"invalid source window energy");
    TunnelingSourceAssembly result;
    result.constants=c.constants;
    result.sourceMesh=m;
    result.boxes=tunnelingSourceBoxes(m,c,policy);
    auto paths=NonlocalPathGeometry(m).traceFromNodes(s,c,options);
    std::vector<NonlocalPairEvent> events;
    result.resolved=true;
    for (auto& geometry:paths) {
        const Index seed=geometry.path.id;
        TunnelingSourceRecord record;record.seedArea_m2=result.boxes.seedAreas_m2.at(seed);
        const auto termination=geometry.path.termination;
        if (termination==TunnelingTermination::Complete) {
            try {
                const auto p=tunnelingProfileFromGeometry(geometry,s,c);
                const auto& w=windows.nodes[seed];
                const TunnelingReservoirs r{sampleFermi(geometry.path.points.back().support,s,true),
                    sampleFermi(geometry.path.points.front().support,s,false),w.valenceAvailableEnergy_J,w.conductionAvailableEnergy_J};
                record.kernel=evaluatePhononTunneling(p,r,c);
                // Integrate in log coordinates so narrowing a point rate cannot
                // silently discard a representable integrated event.
                record.signedEvent_per_m_s=record.kernel->sign==0 || record.seedArea_m2==0 ? 0 :
                    static_cast<Real>(record.kernel->sign*std::exp(static_cast<long double>(record.kernel->logAbsoluteGeneration)+std::log(static_cast<long double>(record.seedArea_m2))));
                check(std::isfinite(record.signedEvent_per_m_s),"nonfinite integrated path event");
                events.push_back({geometry.path.id,record.signedEvent_per_m_s,
                    geometry.path.points.front().support,geometry.path.points.back().support});
                record.resolved=true;record.disposition="evaluated";
            } catch (const std::invalid_argument& e) {record.disposition=e.what();}
        } else if (termination==TunnelingTermination::LengthLimit || termination==TunnelingTermination::RegionExcluded ||
                   termination==TunnelingTermination::ZeroDirection) {
            record.resolved=true;record.disposition="explicit_geometric_zero";
        } else record.disposition="unresolved_path";
        result.resolved &= record.resolved;
        record.geometry=std::move(geometry);result.records.push_back(std::move(record));
    }
    result.sources=distributePairEvents(s,events);
    result.contacts=contactLedger(result.sources,s);
    return result;
}
void TunnelingSourceAssembly::requireResolved(const TunnelingBandState& s) const {
    sources.validate(s);
    check(sourceMesh.revision==s.identity.mesh && sourceMesh.positions_m.size()==s.nodes.size(),"stale source mesh");
    check(resolved,"partial nonlocal source audit cannot enter actual residual");
    check(boxes.seedAreas_m2.size()==s.nodes.size() && boxes.fieldAreas_m2.size()==s.nodes.size(),"source box size mismatch");
    std::vector<NonlocalPairEvent> events;
    std::set<Index> covered;
    for (const auto& r:records) {
        check(r.resolved && r.geometry.path.identity==s.identity,"unresolved or stale source path");
        if (cellSeeded) {
            check(covered.insert(r.geometry.path.id).second,"duplicate source cell record");
            const auto& cell=sourceMesh.cells.at(r.geometry.path.id);
            check(cell.semiconductor && cell.material==sourceMesh.siliconMaterial &&
                std::find(selectedRegions.begin(),selectedRegions.end(),cell.region)!=selectedRegions.end(),"unexpected source cell record");
            events.insert(events.end(),r.collectionEvents.begin(),r.collectionEvents.end());
        } else if (r.kernel) events.push_back({r.geometry.path.id,r.signedEvent_per_m_s,
            r.geometry.path.points.front().support,r.geometry.path.points.back().support});
    }
    if (cellSeeded) for (Index k=0;k<sourceMesh.cells.size();++k) {
        const auto& cell=sourceMesh.cells[k];
        if (cell.semiconductor && cell.material==sourceMesh.siliconMaterial &&
            std::find(selectedRegions.begin(),selectedRegions.end(),cell.region)!=selectedRegions.end())
            check(covered.contains(k),"missing source cell record");
    }
    const auto expected=distributePairEvents(s,events);
    check(expected.electronIntegrals_per_m_s==sources.electronIntegrals_per_m_s &&
        expected.holeIntegrals_per_m_s==sources.holeIntegrals_per_m_s,"source record/ledger mismatch");
}
nlohmann::json tunnelingSourceAudit(const TunnelingSourceAssembly& a, const TunnelingBandState& s) {
    nlohmann::json records=nlohmann::json::array(),counts=nlohmann::json::object();
    for (const auto& r:a.records) {
        nlohmann::json v{{"geometry",tunnelingGeometryAudit(r.geometry)},{"seed_area_m2",r.seedArea_m2},
            {"resolved",r.resolved},{"disposition",r.disposition},{"event_per_m_s",r.signedEvent_per_m_s},
            {"energy_samples",r.energySamples},{"minimum_policy_stops",r.minimumPolicyStops},{"depth_policy_stops",r.depthPolicyStops}};
        const std::string key=r.kernel || r.disposition=="cell_evaluated"?"evaluated":r.disposition=="length_excluded_energy_window"?"length_excluded_energy_window":r.disposition=="empty_energy_window"?"empty_energy_window":r.disposition=="explicit_geometric_zero"?
            tunnelingTerminationName(r.geometry.path.termination):"unresolved";
        counts[key]=counts.value(key,Index(0))+1;
        if (r.kernel) v["kernel"]={{"action",r.kernel->action},{"log_abs_generation",r.kernel->logAbsoluteGeneration},
            {"generation_per_m3_s",r.kernel->generation_per_m3_s},{"valence_window_factor",r.kernel->valenceWindowFactor},
            {"conduction_window_factor",r.kernel->conductionWindowFactor},{"underflow",r.kernel->underflow},{"sign",r.kernel->sign}};
        records.push_back(std::move(v));
    }
    const auto& l=a.contacts;
    long double absoluteEvents=0;for(const auto& event:a.sources.ledger) absoluteEvents+=std::abs(static_cast<long double>(event.signedRate_per_m_s));
    return {{"schema","vela.experimental.nonlocalpath.source_audit.v1"},{"resolved",a.resolved},
        {"native_window_mapping_qualified",false},{"derivatives_evaluated",false},{"counts",counts},
        {"source_policy",a.cellSeeded?"experimental_cell_fraction_first_crossing":"vertex_canonical_n3"},
        {"independent_quadrature_accuracy_qualified",false},
        {"electron_integrals_per_m_s",a.sources.electronIntegrals_per_m_s},{"hole_integrals_per_m_s",a.sources.holeIntegrals_per_m_s},
        {"electron_generation_per_m3_s",a.sources.volumeRates(true,a.boxes.fieldAreas_m2)},
        {"hole_generation_per_m3_s",a.sources.volumeRates(false,a.boxes.fieldAreas_m2)},
        {"field_areas_m2",a.boxes.fieldAreas_m2},{"pair_imbalance",a.sources.balanceError()},
        {"absolute_event_integral_per_m_s",static_cast<Real>(absoluteEvents)},
        {"contacts",{{"free_electron_per_m_s",l.freeElectron_per_m_s},{"contact_electron_per_m_s",l.contactElectron_per_m_s},
                     {"free_hole_per_m_s",l.freeHole_per_m_s},{"contact_hole_per_m_s",l.contactHole_per_m_s}}},{"records",records}};
}
} // namespace vela::experimental
