// Local-temperature state diagnostic. This does not run a coupled solve.
#include "vela/physics/IalMobility.h"
#include "vela/physics/IalHighFieldMobility.h"
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
        keys(cfg,{"temperature_K","carrier","parameters_cm","states_SI","derivatives","high_field"});
        const auto temperature=cfg.at("temperature_K").get<double>();
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
                        "normalField_V_per_m","interfaceDistance_m","crystal_normal","temperature_K","drivingField_V_per_m"});
            const vela::IalMobilityState inputState{state.at("donors_m3"),state.at("acceptors_m3"),
                state.at("electrons_m3"),state.at("holes_m3"),state.at("normalField_V_per_m"),
                state.at("interfaceDistance_m"),state.value("temperature_K",temperature)};
            const auto result=model.evaluate(inputState);
            json row={{"id",state.at("id")},{"mobility_m2_per_Vs",result.mobility_m2_per_Vs},
                      {"coulomb3d_m2_per_Vs",result.coulomb3d_m2_per_Vs},
                      {"coulomb2d_m2_per_Vs",result.coulomb2d_m2_per_Vs},
                      {"coulomb_m2_per_Vs",result.coulomb_m2_per_Vs},
                      {"phonon_m2_per_Vs",result.phonon_m2_per_Vs},
                      {"roughness_m2_per_Vs",result.roughness_m2_per_Vs}};
            row["temperature_K"]=inputState.temperature_K;
            if (cfg.value("derivatives",false)) {
                const auto differential=model.evaluateWithDerivatives(inputState);
                row["mobility_derivatives_SI"]=differential.derivative_SI;
                row["mobility_temperature_derivative_m2_per_Vs_K"]=differential.temperatureDerivative_m2_per_Vs_K;
            }
            if (cfg.contains("high_field")) {
                const auto& hf=cfg.at("high_field");
                keys(hf,{"vsat300_m_per_s","beta300","vsat_exponent","beta_exponent"});
                const vela::IalHighFieldParameters p{hf.at("vsat300_m_per_s"),hf.at("beta300"),
                    hf.at("vsat_exponent"),hf.at("beta_exponent")};
                const auto limited=vela::evaluateIalHighFieldMobility(result.mobility_m2_per_Vs,
                    state.at("drivingField_V_per_m"),inputState.temperature_K,p);
                row["high_field_mobility_m2_per_Vs"]=limited.mobility_m2_per_Vs;
                if (cfg.value("derivatives",false))
                    row["high_field_temperature_derivative_m2_per_Vs_K"]=limited.temperatureDerivative_m2_per_Vs_K+
                        limited.lowFieldDerivative*model.evaluateWithDerivatives(inputState).temperatureDerivative_m2_per_Vs_K;
            } else if (state.contains("drivingField_V_per_m")) {
                throw std::invalid_argument("drivingField_V_per_m requires explicit high_field parameters");
            }
            if (state.contains("crystal_normal"))
                row["orientation_family"]=vela::IalMobility::orientationFamily(
                    state.at("crystal_normal").get<std::array<vela::Real,3>>());
            results.push_back(row);
        }
        // Infinite component mobility means no scattering; JSON uses null.
        std::cout<<json({{"schema","vela.ialmob.local_probe.v1"},{"temperature_K",temperature},
            {"carrier",carrier},{"scope","FullPhuMob, PhononCombination=1; no stress or thickness correction; optional explicit alpha=0 high-field diagnostic; orientation label does not select parameters"},
            {"infinite_component_representation","null"},{"results",results}}).dump(2)<<'\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<"IALMob probe: "<<error.what()<<'\n';
        return 1;
    }
}
