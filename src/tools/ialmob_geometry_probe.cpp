// Geometry-only probe: no reference field/label inputs and no Newton solve.
#include "vela/physics/IalInterfaceGeometry.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>
using json=nlohmann::json;
int main(int argc,char** argv) {
    try {
        if (argc!=2) throw std::invalid_argument("Usage: ialmob_geometry_probe INPUT.json");
        std::ifstream stream(argv[1]);const auto input=json::parse(stream);
        for (const auto& [k,v]:input.items())
            if (k!="coordinates_m"&&k!="segments"&&k!="crystal_x"&&k!="crystal_y")
                throw std::invalid_argument("Unsupported geometry input: "+k);
        const auto xy=input.at("coordinates_m").get<std::vector<std::array<vela::Real,2>>>();
        std::vector<vela::IalInterfaceSegment> segments;
        for (const auto& s:input.at("segments"))
            segments.push_back({s.at("node0"),s.at("node1"),s.at("semiconductor_point_m")});
        const auto results=vela::buildIalInterfaceGeometry(xy,segments,
            input.at("crystal_x"),input.at("crystal_y"));
        json rows=json::array();
        for (std::size_t i=0;i<results.size();++i) {
            const auto& r=results[i];rows.push_back({{"node_id",i},{"distance_m",r.distance_m},
                {"orientation_family",r.orientationFamily},{"nearest_interface_vertex",r.nearestInterfaceVertex},
                {"on_interface",r.onInterface}});
        }
        std::cout<<json({{"schema","vela.ialmob.geometry_probe.v1"},{"nodes",rows}}).dump(2)<<'\n';
        return 0;
    } catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
