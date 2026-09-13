#include "vela/physics/SiliconThermalPhysics.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/core/PhysicalConstants.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <cmath>
using json=nlohmann::json;
int main(int argc,char**argv){try{
    if(argc!=2)throw std::invalid_argument("Usage: silicon_thermal_probe INPUT.json");
    std::ifstream f(argv[1]);auto input=json::parse(f);
    if(!input.contains("states_SI") || input.size()!=(input.contains("auger_with_generation")?2:1))
        throw std::invalid_argument("Expected states_SI and optional auger_with_generation; audited Siliconc100 model");
    vela::SiliconThermalParameters parameters;parameters.augerWithGeneration=input.value("auger_with_generation",false);
    vela::SiliconThermalPhysics model(parameters);json rows=json::array();
    for(const auto&s:input.at("states_SI")){
        for(const auto&[k,v]:s.items())if(k!="id"&&k!="potential_V"&&k!="electron_qf_V"&&k!="hole_qf_V"&&k!="temperature_K"&&k!="donors_m3"&&k!="acceptors_m3"&&k!="electron_qf_reference_V"&&k!="hole_qf_reference_V"&&k!="reference_conduction_band_eV"&&k!="reference_valence_band_eV")throw std::invalid_argument("Unsupported state key: "+k);
        auto r=model.evaluate({s.at("potential_V"),s.at("electron_qf_V"),s.at("hole_qf_V"),s.at("temperature_K"),s.at("donors_m3"),s.at("acceptors_m3"),s.value("electron_qf_reference_V",0.),s.value("hole_qf_reference_V",0.)});
        json row={{"id",s.at("id")},{"bgn_eV",r.bandgapNarrowing_eV}};
        const auto put=[&](const char*name,const vela::ThermalQuantity&q){row[name]={{"value",q.value},{"derivative",q.derivative}};};
        put("bandgap_eV",r.bandgap_eV);put("affinity_eV",r.affinity_eV);put("Nc_m3",r.Nc_m3);put("Nv_m3",r.Nv_m3);put("ni_m3",r.ni_m3);put("effective_ni_m3",r.effectiveNi_m3);
        put("conduction_band_eV",r.conductionBand_eV);put("valence_band_eV",r.valenceBand_eV);put("electrons_m3",r.electrons_m3);put("holes_m3",r.holes_m3);put("srh_m3_per_s",r.srhRate_m3_per_s);put("auger_m3_per_s",r.augerRate_m3_per_s);
        put("auger_electron_m6_per_s",r.augerElectron_m6_per_s);put("auger_hole_m6_per_s",r.augerHole_m6_per_s);
        put("electron_lifetime_s",r.electronLifetime_s);put("hole_lifetime_s",r.holeLifetime_s);
        if(s.contains("reference_conduction_band_eV")!=s.contains("reference_valence_band_eV"))
            throw std::invalid_argument("Both reference bands are required");
        if(s.contains("reference_conduction_band_eV")){
            const double ec=s.at("reference_conduction_band_eV"),ev=s.at("reference_valence_band_eV");
            if(!std::isfinite(ec)||!std::isfinite(ev))throw std::invalid_argument("Nonfinite reference band");
            const double vt=vela::constants::kb/vela::constants::q*s.at("temperature_K").get<double>();
            const double en=static_cast<double>(-static_cast<long double>(s.value("electron_qf_reference_V",0.))-s.at("electron_qf_V").get<double>()-ec)/vt;
            const double ep=static_cast<double>(static_cast<long double>(s.value("hole_qf_reference_V",0.))+s.at("hole_qf_V").get<double>()+ev)/vt;
            row["electrons_at_native_band_m3"]=r.Nc_m3.value*vela::fermiDiracHalf(en);
            row["holes_at_native_band_m3"]=r.Nv_m3.value*vela::fermiDiracHalf(ep);
        }
        rows.push_back(row);
    }
    std::cout<<json({{"scope","local audited silicon temperature physics; no coupled solve"},{"auger_with_generation",parameters.augerWithGeneration},{"reference_potential_V",model.referencePotential_V()},{"derivative_order",{"psi_V","electron_qf_V","hole_qf_V","temperature_K"}},{"results",rows}}).dump(2)<<'\n';
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
