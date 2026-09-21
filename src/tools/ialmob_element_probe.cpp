// Fixed-state element mobility replay; does not run a coupled solve.
#include "vela/equation/IalElementMobility.h"
#include "IalProbeParameters.h"
#include <fstream>
#include <iostream>
#include <map>
using json=nlohmann::json;
int main(int argc,char** argv) {
    try {
        if (argc!=2) throw std::invalid_argument("Usage: ialmob_element_probe INPUT.json");
        std::ifstream stream(argv[1]);const auto input=json::parse(stream);
        for (const auto& [k,v]:input.items())
            if (k!="temperature_K"&&k!="electron_parameters_cm"&&k!="hole_parameters_cm"&&k!="elements_SI")
                throw std::invalid_argument("Unsupported element probe input: "+k);
        if (input.at("temperature_K").get<double>()!=300.) throw std::invalid_argument("IALMob requires 300 K");
        std::map<int,vela::IalMobility> electron,hole;
        for (const auto& [key,params]:input.at("electron_parameters_cm").items())
            electron.emplace(std::stoi(key),vela::IalMobility(vela::probe_detail::parameters(params,true),true));
        for (const auto& [key,params]:input.at("hole_parameters_cm").items())
            hole.emplace(std::stoi(key),vela::IalMobility(vela::probe_detail::parameters(params,false),false));
        json rows=json::array();
        for (const auto& cell:input.at("elements_SI")) {
            vela::IalElementGeometry g;
            g.coordinates_m=cell.at("coordinates_m");
            g.interfaceDistance_m=cell.at("interface_distance_m");
            g.vertexMeasure_m2=cell.at("vertex_measure_m2");
            g.partialBoundaryLayer=cell.at("partial_boundary_layer");
            g.boundaryTangent=cell.at("boundary_tangent");
            g.touchesEffectiveElectrode=cell.at("touches_effective_electrode");
            std::array<vela::IalElementVertexState,3> states;
            std::array<const vela::IalMobility*,3> em,hm;
            const auto& vertices=cell.at("vertices");
            if (vertices.size()!=3) throw std::invalid_argument("IALMob requires three vertices");
            for (int k=0;k<3;++k) {
                const auto& v=vertices.at(k);const int family=v.at("orientation_family");
                em[k]=&electron.at(family);hm[k]=&hole.at(family);
                states[k]={v.at("potential_V"),v.at("electron_qf_V"),v.at("hole_qf_V"),
                    v.at("donors_m3"),v.at("acceptors_m3"),v.at("electrons_m3"),v.at("holes_m3"),
                    v.at("electron_response_m3_per_V"),v.at("hole_response_m3_per_V")};
            }
            const auto r=vela::evaluateIalElementMobility(g,states,em,hm);
            rows.push_back({{"cell_id",cell.at("cell_id")},{"electron_low_m2_per_Vs",r.electronLowField.value},
                {"hole_low_m2_per_Vs",r.holeLowField.value},{"electron_m2_per_Vs",r.electron.value},
                {"hole_m2_per_Vs",r.hole.value},{"electron_derivative_per_V",r.electron.derivative},
                {"hole_derivative_per_V",r.hole.derivative}});
        }
        std::cout<<json({{"schema","vela.ialmob.element_probe.v1"},{"elements",rows}}).dump(2)<<'\n';
        return 0;
    } catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
