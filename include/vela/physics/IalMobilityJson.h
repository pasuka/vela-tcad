#pragma once
#include "vela/physics/IalMobility.h"
#include <nlohmann/json.hpp>
#include <map>
#include <stdexcept>
namespace vela::ial_json {
inline const std::map<std::string,vela::Real vela::IalMobilityParameters::*> fields{
    {"mumax",&vela::IalMobilityParameters::muMax},
    {"mumin",&vela::IalMobilityParameters::muMin},
    {"theta",&vela::IalMobilityParameters::theta},
    {"alpha",&vela::IalMobilityParameters::alpha},
    {"n_ref",&vela::IalMobilityParameters::nRef},
    {"mass",&vela::IalMobilityParameters::mass},
    {"other_mass",&vela::IalMobilityParameters::otherMass},
    {"nref_D",&vela::IalMobilityParameters::nRefD},
    {"nref_A",&vela::IalMobilityParameters::nRefA},
    {"cref_D",&vela::IalMobilityParameters::cRefD},
    {"cref_A",&vela::IalMobilityParameters::cRefA},
    {"ndop_ref",&vela::IalMobilityParameters::nDopRef},
    {"nsc_ref",&vela::IalMobilityParameters::nScRef},
    {"S",&vela::IalMobilityParameters::S},
    {"p",&vela::IalMobilityParameters::transitionP},
    {"B",&vela::IalMobilityParameters::B},
    {"C",&vela::IalMobilityParameters::C},
    {"lambda",&vela::IalMobilityParameters::lambda},
    {"k",&vela::IalMobilityParameters::k},
    {"delta",&vela::IalMobilityParameters::delta},
    {"eta",&vela::IalMobilityParameters::eta},
    {"lambda_sr",&vela::IalMobilityParameters::lambdaSr},
    {"A",&vela::IalMobilityParameters::A},
    {"alpha_sr",&vela::IalMobilityParameters::alphaSr},
    {"nu",&vela::IalMobilityParameters::nu},
    {"N1",&vela::IalMobilityParameters::N1},
    {"N2",&vela::IalMobilityParameters::N2},
    {"l_crit",&vela::IalMobilityParameters::lCrit},
    {"l_crit_c",&vela::IalMobilityParameters::lCritC},
    {"D1_inv",&vela::IalMobilityParameters::d1Inv},
    {"D2_inv",&vela::IalMobilityParameters::d2Inv},
    {"nu0_inv",&vela::IalMobilityParameters::nu0Inv},
    {"nu1_inv",&vela::IalMobilityParameters::nu1Inv},
    {"nu2_inv",&vela::IalMobilityParameters::nu2Inv},
    {"D1_acc",&vela::IalMobilityParameters::d1Acc},
    {"D2_acc",&vela::IalMobilityParameters::d2Acc},
    {"nu0_acc",&vela::IalMobilityParameters::nu0Acc},
    {"nu1_acc",&vela::IalMobilityParameters::nu1Acc},
    {"nu2_acc",&vela::IalMobilityParameters::nu2Acc}
};
inline IalMobilityParameters parameters(const nlohmann::json& overrides,bool electron) {
    auto p=IalMobility::siliconDefaults(electron);
    if (!overrides.is_object()) throw std::invalid_argument("parameters_cm must be an object");
    for (const auto& [key,value]:overrides.items()) {
        const auto found=fields.find(key);
        if (found==fields.end()) throw std::invalid_argument("Unsupported parameter: "+key);
        p.*(found->second)=value.get<Real>();
    }
    return p;
}
} // namespace vela::ial_json
