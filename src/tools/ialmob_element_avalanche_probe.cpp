// Diagnostic source/Jacobian replay with explicit geometry; no coupled solve.
#include "vela/equation/IalElementAvalanche.h"
#include "IalProbeParameters.h"
#include <fstream>
#include <iostream>
#include <map>
using json=nlohmann::json;
int main(int argc,char** argv) {
    try {
        if (argc!=2) throw std::invalid_argument("Usage: ialmob_element_avalanche_probe INPUT.json");
        std::ifstream stream(argv[1]);const auto input=json::parse(stream);
        for (const auto& [k,v]:input.items())
            if (k!="temperature_K"&&k!="electron_parameters_cm"&&k!="hole_parameters_cm"&&k!="elements_SI"&&k!="high_field"&&k!="reference_density_m3")
                throw std::invalid_argument("Unsupported element probe input: "+k);
        if (input.at("temperature_K").get<double>()!=300.) throw std::invalid_argument("IALMob requires 300 K");
        vela::IalElementAvalancheOptions options;
        options.highFieldMobility=input.value("high_field",true);
        if (input.value("reference_density_m3",0.)!=0.) throw std::invalid_argument("Avalanche probe requires zero reference density");
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
            options.mobilityVertex=cell.at("avalanche_mobility_vertex").get<std::size_t>();
            options.sourceMeasure_m2=cell.at("avalanche_source_measure_m2");
            const auto r=vela::evaluateIalElementAvalanche(g,states,em,hm,options);
            const auto serialize=[](const vela::IalElementAvalancheCarrier& c) {
                json shares=json::array();
                for (const auto& s:c.source_per_m_s)
                    shares.push_back({{"value",s.value},{"derivative_per_V",s.derivative}});
                return json{{"density_m3",c.density_m3.value},
                    {"mobility_m2_per_Vs",c.mobility_m2_per_Vs.value},
                    {"current_field_V_per_m",c.currentField_V_per_m.value},
                    {"ionization_field_V_per_m",c.ionizationField_V_per_m.value},
                    {"alpha_per_m",c.alpha_per_m.value},
                    {"generation_m3_per_s",c.generation_m3_per_s.value},
                    {"generation_derivative_per_V",c.generation_m3_per_s.derivative},
                    {"source_per_m_s",shares}};
            };
            rows.push_back({{"cell_id",cell.at("cell_id")},
                {"electron",serialize(r.electron)},{"hole",serialize(r.hole)}});
        }
        std::cout<<json({{"schema","vela.ialmob.element_avalanche_probe.v1"},{"elements",rows}}).dump(2)<<'\n';
        return 0;
    } catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
