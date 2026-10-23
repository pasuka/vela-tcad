// Private geometry audit; intentionally absent from installed CLI dispatch.
#include "vela/physics/NonlocalPathGeometry.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <set>

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
        if (argc != 3) throw std::invalid_argument("usage: geometry_probe input.json audit.json");
        std::ifstream input(argv[1]);
        if (!input) throw std::invalid_argument("cannot open geometry input");
        const auto j = nlohmann::json::parse(input);
        if (j.at("schema") != "vela.experimental.nonlocalpath.geometry_input.v1")
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
            // QF/derivatives are not used by geometry; no source/RHS claim.
            state.nodes.push_back(n);
        }
        for (const auto& b : mesh.boundaries) if (b.kind == TunnelingBoundaryKind::Contact)
            for (Index n : b.nodes) state.nodes.at(n).contact = true;
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
        const NonlocalPathGeometry geometry(std::move(mesh));
        const auto results = geometry.traceFromNodes(state, config, options);
        nlohmann::json paths = nlohmann::json::array(), counts = nlohmann::json::object();
        bool failure = false;
        for (const auto& r : results) {
            paths.push_back(tunnelingGeometryAudit(r));
            const auto name = tunnelingTerminationName(r.path.termination);
            counts[name] = counts.value(name, Index(0)) + 1;
            failure |= r.path.termination == TunnelingTermination::NumericalFailure ||
                       r.path.termination == TunnelingTermination::TraversalLimit ||
                       r.path.termination == TunnelingTermination::ReflectionLimit;
        }
        std::ofstream output(argv[2]);
        if (!output) throw std::invalid_argument("cannot open geometry audit output");
        output << nlohmann::json{{"schema", "vela.experimental.nonlocalpath.batch_audit.v1"},
            {"scope", "geometry_only_no_source_or_RHS_binding"}, {"counts", counts}, {"paths", paths}}.dump();
        output.close();
        if (!output) throw std::runtime_error("geometry audit write failed");
        std::cout << nlohmann::json{{"paths", results.size()}, {"counts", counts}}.dump() << '\n';
        return failure ? 3 : 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 2;
    }
}
