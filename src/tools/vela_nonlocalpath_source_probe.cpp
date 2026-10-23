// Private N4 source audit; intentionally absent from installed CLI dispatch.
#include "vela/physics/NonlocalPathSourceAssembly.h"
#include "vela/physics/NonlocalPathCellSource.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <set>
#include "NonlocalPathResidualProbe.h"

using namespace vela;
using namespace vela::experimental;
namespace {
Index index(const nlohmann::json& value) {
    if (!value.is_number_unsigned()) throw std::invalid_argument("unsigned geometry index required");
    return value.get<Index>();
}
template<std::size_t N> std::array<Index, N> indices(const nlohmann::json& values) {
    if (!values.is_array() || values.size() != N) throw std::invalid_argument("invalid geometry index shape");
    std::array<Index, N> result;
    for (std::size_t i = 0; i < N; ++i) result[i] = index(values.at(i));
    return result;
}
}
int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::invalid_argument("usage: source_probe input.json audit.json");
        std::ifstream input(argv[1]);
        if (!input) throw std::invalid_argument("cannot open geometry input");
        const auto j = nlohmann::json::parse(input);
        if (j.at("schema") != "vela.experimental.nonlocalpath.source_input.v1")
            throw std::invalid_argument("unsupported geometry input schema");
        TunnelingGeometryMesh mesh;
        TunnelingBandState state;
        mesh.revision = index(j.at("mesh_revision"));
        state.identity.mesh = mesh.revision;
        state.identity.trial = index(j.at("trial_revision"));
        for (const auto& p : j.at("positions_m")) {
            if (!p.is_array() || p.size() != 2) throw std::invalid_argument("invalid coordinate shape");
            mesh.positions_m.emplace_back(p.at(0).get<Real>(), p.at(1).get<Real>());
        }
        for (const auto& c : j.at("cells")) mesh.cells.push_back({indices<3>(c.at("nodes")),
            index(c.at("material")), index(c.at("region")), c.at("semiconductor").get<bool>()});
        for (const auto& b : j.at("boundaries")) {
            const auto kind = b.at("kind").get<std::string>();
            if (kind != "reflecting" && kind != "contact" && kind != "open")
                throw std::invalid_argument("unknown boundary kind");
            mesh.boundaries.push_back({indices<2>(b.at("nodes")),
                kind == "reflecting" ? TunnelingBoundaryKind::Reflecting :
                kind == "contact" ? TunnelingBoundaryKind::Contact : TunnelingBoundaryKind::Open});
        }
        const auto& bands = j.at("bands");
        if (!bands.is_array() || bands.size() != mesh.positions_m.size())
            throw std::invalid_argument("band-state size mismatch");
        std::vector<std::set<Index>> memberships(bands.size());
        for (const auto& c : mesh.cells) if (c.semiconductor && c.material == mesh.siliconMaterial)
            for (Index n : c.nodes) memberships.at(n).insert(c.region);
        for (Index i = 0; i < bands.size(); ++i) {
            TunnelingBandNode n;
            n.semiconductor = !memberships[i].empty();
            if (n.semiconductor) {
                n.region = *memberships[i].begin();
                for (Index r : memberships[i]) if (r != n.region) n.additionalRegions.push_back(r);
            }
            const auto& b = bands[i];
            n.conduction_J = b.at("conduction_J").get<Real>();
            n.valence_J = b.at("valence_J").get<Real>();
            n.electrons_per_m3 = b.at("electrons_per_m3").get<Real>();
            n.holes_per_m3 = b.at("holes_per_m3").get<Real>();
            n.conductionDOS_per_m3 = b.at("conduction_DOS_per_m3").get<Real>();
            n.valenceDOS_per_m3 = b.at("valence_DOS_per_m3").get<Real>();
            n.temperature_K = b.at("temperature_K").get<Real>();
            n.electronFermi_J = b.at("electron_Fermi_J").get<Real>();
            n.holeFermi_J = b.at("hole_Fermi_J").get<Real>();
            state.nodes.push_back(n);
        }
        for (const auto& b : mesh.boundaries) if (b.kind == TunnelingBoundaryKind::Contact)
            for (Index n : b.nodes) state.nodes.at(n).contact = true;
        if(j.contains("contact_nodes")) for(const auto& n:j.at("contact_nodes")) state.nodes.at(index(n)).contact=true;
        const auto config = nonlocalPathConfigFromJson(j.at("config"));
        TunnelingGeometryOptions options;
        if (j.contains("options")) {
            const auto& o = j.at("options");
            if (!o.is_object()) throw std::invalid_argument("options must be an object");
            for (auto it = o.begin(); it != o.end(); ++it)
                if (it.key() != "within_selected_regions" && it.key() != "maximum_segments" && it.key() != "maximum_reflections")
                    throw std::invalid_argument("unknown geometry option");
            for (const auto* key : {"maximum_segments", "maximum_reflections"})
                if (o.contains(key) && !o.at(key).is_number_unsigned()) throw std::invalid_argument("unsigned traversal limits required");
            options.withinSelectedRegions = o.value("within_selected_regions", true);
            options.maximumSegments = o.value("maximum_segments", Index(4096));
            options.maximumReflections = o.value("maximum_reflections", Index(32));
        }
        TunnelingWindowInput windows; windows.identity=state.identity;
        for (const auto& b:bands) windows.nodes.push_back({b.at("valence_window_J").get<Real>(),b.at("conduction_window_J").get<Real>()});
        const auto policy=j.at("box_policy").get<std::string>();
        if (policy!="barycentric" && policy!="delaunay") throw std::invalid_argument("unknown source box policy");
        const auto sourcePolicy=j.value("source_policy",std::string("vertex_canonical_n3"));
        TunnelingSourceAssembly result;
        if(sourcePolicy=="experimental_cell_fraction_first_crossing") {
            const auto& p=j.at("cell_integration_policy");
            for(auto it=p.begin();it!=p.end();++it) if(it.key()!="variation_reference_per_m2_s")
                throw std::invalid_argument("unknown private cell integration option");
            TunnelingCellIntegrationPolicy integration;
            integration.variationReference_per_m2_s=p.at("variation_reference_per_m2_s").get<Real>();
            if(policy!="delaunay") throw std::invalid_argument("cell diagnostic requires Delaunay boxes");
            result=assembleNonlocalPathCellSources(mesh,state,config,integration,options);
        } else if(sourcePolicy=="vertex_canonical_n3") result=assembleNonlocalPathSources(mesh,state,config,windows,
            policy=="barycentric" ? TunnelingBoxPolicy::Barycentric : TunnelingBoxPolicy::Delaunay,options);
        else throw std::invalid_argument("unknown private source policy");
        if(result.resolved) result.requireResolved(state);
        auto audit=tunnelingSourceAudit(result,state);
        if(j.contains("actual_residual")) {
            result.requireResolved(state);
            audit["actual_residual"]=nonlocalResidualProbe(j,mesh,state,result);
        }
        audit["window_input_provenance"]=sourcePolicy=="experimental_cell_fraction_first_crossing"
            ? nlohmann::json("autonomous source faces and global silicon band extrema") : j.at("window_input_provenance");
        std::ofstream output(argv[2]);
        if (!output) throw std::invalid_argument("cannot open source audit output");
        output<<audit.dump();output.close();
        if (!output) throw std::runtime_error("source audit write failed");
        std::cout<<nlohmann::json{{"resolved",result.resolved},{"counts",audit.at("counts")},
            {"pair_imbalance",result.sources.balanceError()}}.dump()<<'\n';
        return result.resolved ? 0 : 4;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 2;
    }
}
