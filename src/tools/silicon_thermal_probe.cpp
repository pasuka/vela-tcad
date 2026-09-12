#include "vela/physics/SiliconThermalPhysics.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
using json=nlohmann::json;
int main(int argc,char**argv){try{
    if(argc!=2)throw std::invalid_argument("Usage: silicon_thermal_probe INPUT.json");
    std::ifstream f(argv[1]);auto input=json::parse(f);
    if(input.size()!=1||!input.contains("states_SI"))throw std::invalid_argument("Expected states_SI only; audited Siliconc100 model");
    vela::SiliconThermalPhysics model;json rows=json::array();
    for(const auto&s:input.at("states_SI")){
        for(const auto&[k,v]:s.items())if(k!="id"&&k!="potential_V"&&k!="electron_qf_V"&&k!="hole_qf_V"&&k!="temperature_K"&&k!="donors_m3"&&k!="acceptors_m3")throw std::invalid_argument("Unsupported state key: "+k);
        auto r=model.evaluate({s.at("potential_V"),s.at("electron_qf_V"),s.at("hole_qf_V"),s.at("temperature_K"),s.at("donors_m3"),s.at("acceptors_m3")});
        json row={{"id",s.at("id")},{"bgn_eV",r.bandgapNarrowing_eV}};
        const auto put=[&](const char*name,const vela::ThermalQuantity&q){row[name]={{"value",q.value},{"derivative",q.derivative}};};
        put("bandgap_eV",r.bandgap_eV);put("affinity_eV",r.affinity_eV);put("Nc_m3",r.Nc_m3);put("Nv_m3",r.Nv_m3);put("ni_m3",r.ni_m3);put("effective_ni_m3",r.effectiveNi_m3);
        put("conduction_band_eV",r.conductionBand_eV);put("valence_band_eV",r.valenceBand_eV);put("electrons_m3",r.electrons_m3);put("holes_m3",r.holes_m3);put("srh_m3_per_s",r.srhRate_m3_per_s);put("auger_m3_per_s",r.augerRate_m3_per_s);
        rows.push_back(row);
    }
    std::cout<<json({{"scope","local audited silicon temperature physics; no coupled solve"},{"reference_potential_V",model.referencePotential_V()},{"derivative_order",{"psi_V","electron_qf_V","hole_qf_V","temperature_K"}},{"results",rows}}).dump(2)<<'\n';
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
