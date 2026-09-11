// Isothermal local-state diagnostic. This does not run a coupled solve.
#include "vela/physics/IalMobility.h"
#include "IalProbeParameters.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>

using json=nlohmann::json;

namespace {


void keys(const json& object,const std::initializer_list<std::string>& allowed)
{
    if (!object.is_object()) throw std::invalid_argument("Expected JSON object");
    for (const auto& [key,value]:object.items())
        if (std::find(allowed.begin(),allowed.end(),key)==allowed.end())
            throw std::invalid_argument("Unsupported field: "+key);
}
}

int main(int argc,char** argv)
{
    if (argc!=2) {
        std::cerr<<"Usage: ialmob_probe INPUT.json (results on stdout)\n";
        return 2;
    }
    try {
        std::ifstream input(argv[1]);
        if (!input) throw std::runtime_error("Cannot open probe input");
        const auto cfg=json::parse(input);
        keys(cfg,{"temperature_K","carrier","parameters_cm","states_SI","derivatives"});
        if (cfg.at("temperature_K").get<double>()!=300.)
            throw std::invalid_argument("IALMob probe supports only 300 K");
        const auto carrier=cfg.at("carrier").get<std::string>();
        if (carrier!="electron" && carrier!="hole")
            throw std::invalid_argument("carrier must be electron or hole");
        const auto parameters=vela::probe_detail::parameters(
            cfg.value("parameters_cm",json::object()),carrier=="electron");
        const vela::IalMobility model(parameters,carrier=="electron");
        const auto& states=cfg.at("states_SI");
        if (!states.is_array() || states.empty()) throw std::invalid_argument("states_SI must be a nonempty array");
        json results=json::array();
        for (const auto& state:states) {
            keys(state,{"id","donors_m3","acceptors_m3","electrons_m3","holes_m3",
                        "normalField_V_per_m","interfaceDistance_m","crystal_normal"});
            const vela::IalMobilityState inputState{state.at("donors_m3"),state.at("acceptors_m3"),
                state.at("electrons_m3"),state.at("holes_m3"),state.at("normalField_V_per_m"),
                state.at("interfaceDistance_m")};
            const auto result=model.evaluate(inputState);
            json row={{"id",state.at("id")},{"mobility_m2_per_Vs",result.mobility_m2_per_Vs},
                      {"coulomb3d_m2_per_Vs",result.coulomb3d_m2_per_Vs},
                      {"coulomb2d_m2_per_Vs",result.coulomb2d_m2_per_Vs},
                      {"coulomb_m2_per_Vs",result.coulomb_m2_per_Vs},
                      {"phonon_m2_per_Vs",result.phonon_m2_per_Vs},
                      {"roughness_m2_per_Vs",result.roughness_m2_per_Vs}};
            if (cfg.value("derivatives",false))
                row["mobility_derivatives_SI"]=model.evaluateWithDerivatives(inputState).derivative_SI;
            if (state.contains("crystal_normal"))
                row["orientation_family"]=vela::IalMobility::orientationFamily(
                    state.at("crystal_normal").get<std::array<vela::Real,3>>());
            results.push_back(row);
        }
        // Infinite component mobility means no scattering; JSON uses null.
        std::cout<<json({{"schema","vela.ialmob.local_probe.v1"},{"temperature_K",300.},
            {"carrier",carrier},{"scope","FullPhuMob, PhononCombination=1; no stress, thickness or high-field correction; orientation label does not select parameters"},
            {"infinite_component_representation","null"},{"results",results}}).dump(2)<<'\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<"IALMob probe: "<<error.what()<<'\n';
        return 1;
    }
}
