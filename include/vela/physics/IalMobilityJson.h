#pragma once
#include "vela/core/ArithmeticStageInput.h"
#include "vela/physics/IalMobility.h"
#include <nlohmann/json.hpp>
#include <map>
#include <stdexcept>
namespace vela::ial_json {
inline IalScreeningMethod screeningMethod(const nlohmann::json& input) {
    return ialScreeningMethod(input.value("screening_method",
        std::string(ialScreeningMethodName(defaultIalScreeningMethod))));
}
// The old electrothermal diagnostic key remains an explicit alias. Reject
// conflicting declarations rather than report one method while using another.
inline IalScreeningMethod electrothermalScreeningMethod(const nlohmann::json& input) {
    const auto mobility=input.value("mobility_SI",nlohmann::json::object());
    const auto options=mobility.is_object()?mobility.value("ialmob",nlohmann::json::object()):nlohmann::json::object();
    const auto method=screeningMethod(options);
    if (!input.contains("diagnostic_ialmob_screening_method")) return method;
    const auto alias=ialScreeningMethod(input.at("diagnostic_ialmob_screening_method").get<std::string>());
    if (options.contains("screening_method") && method!=alias)
        throw std::invalid_argument("Conflicting IALMob screening_method and diagnostic alias");
    return alias;
}
template<class Scalar>
inline const std::map<std::string,Scalar BasicIalMobilityParameters<Scalar>::*> coefficientFields{
    {"mumax",&BasicIalMobilityParameters<Scalar>::muMax},
    {"mumin",&BasicIalMobilityParameters<Scalar>::muMin},
    {"theta",&BasicIalMobilityParameters<Scalar>::theta},
    {"alpha",&BasicIalMobilityParameters<Scalar>::alpha},
    {"n_ref",&BasicIalMobilityParameters<Scalar>::nRef},
    {"mass",&BasicIalMobilityParameters<Scalar>::mass},
    {"other_mass",&BasicIalMobilityParameters<Scalar>::otherMass},
    {"nref_D",&BasicIalMobilityParameters<Scalar>::nRefD},
    {"nref_A",&BasicIalMobilityParameters<Scalar>::nRefA},
    {"cref_D",&BasicIalMobilityParameters<Scalar>::cRefD},
    {"cref_A",&BasicIalMobilityParameters<Scalar>::cRefA},
    {"ndop_ref",&BasicIalMobilityParameters<Scalar>::nDopRef},
    {"nsc_ref",&BasicIalMobilityParameters<Scalar>::nScRef},
    {"S",&BasicIalMobilityParameters<Scalar>::S},
    {"p",&BasicIalMobilityParameters<Scalar>::transitionP},
    {"B",&BasicIalMobilityParameters<Scalar>::B},
    {"C",&BasicIalMobilityParameters<Scalar>::C},
    {"lambda",&BasicIalMobilityParameters<Scalar>::lambda},
    {"k",&BasicIalMobilityParameters<Scalar>::k},
    {"delta",&BasicIalMobilityParameters<Scalar>::delta},
    {"eta",&BasicIalMobilityParameters<Scalar>::eta},
    {"lambda_sr",&BasicIalMobilityParameters<Scalar>::lambdaSr},
    {"A",&BasicIalMobilityParameters<Scalar>::A},
    {"alpha_sr",&BasicIalMobilityParameters<Scalar>::alphaSr},
    {"nu",&BasicIalMobilityParameters<Scalar>::nu},
    {"N1",&BasicIalMobilityParameters<Scalar>::N1},
    {"N2",&BasicIalMobilityParameters<Scalar>::N2},
    {"l_crit",&BasicIalMobilityParameters<Scalar>::lCrit},
    {"l_crit_c",&BasicIalMobilityParameters<Scalar>::lCritC},
    {"D1_inv",&BasicIalMobilityParameters<Scalar>::d1Inv},
    {"D2_inv",&BasicIalMobilityParameters<Scalar>::d2Inv},
    {"nu0_inv",&BasicIalMobilityParameters<Scalar>::nu0Inv},
    {"nu1_inv",&BasicIalMobilityParameters<Scalar>::nu1Inv},
    {"nu2_inv",&BasicIalMobilityParameters<Scalar>::nu2Inv},
    {"alpha1_inv",&BasicIalMobilityParameters<Scalar>::alpha1Inv},
    {"alpha2_inv",&BasicIalMobilityParameters<Scalar>::alpha2Inv},
    {"alpha1_acc",&BasicIalMobilityParameters<Scalar>::alpha1Acc},
    {"alpha2_acc",&BasicIalMobilityParameters<Scalar>::alpha2Acc},
    {"D1_acc",&BasicIalMobilityParameters<Scalar>::d1Acc},
    {"D2_acc",&BasicIalMobilityParameters<Scalar>::d2Acc},
    {"nu0_acc",&BasicIalMobilityParameters<Scalar>::nu0Acc},
    {"nu1_acc",&BasicIalMobilityParameters<Scalar>::nu1Acc},
    {"nu2_acc",&BasicIalMobilityParameters<Scalar>::nu2Acc}
};
inline const auto& fields = coefficientFields<Real>;
inline IalMobilityParameters parameters(const nlohmann::json& overrides,bool electron,
                                        const nlohmann::json* decimals = nullptr) {
    auto p=IalMobility::siliconDefaults(electron);
    if (!overrides.is_object()) throw std::invalid_argument("parameters_cm must be an object");
    if (decimals && (!decimals->is_object() || decimals->size() != overrides.size()))
        throw std::invalid_argument("IALMob decimal coefficients must cover the compatibility fields");
    for (const auto& [key,value]:overrides.items()) {
        const auto found=fields.find(key);
        if (found==fields.end()) throw std::invalid_argument("Unsupported parameter: "+key);
        const auto compatibility = value.get<Real>();
        if (decimals)
            decimalCompatibilityNumber(decimals->at(key), "IALMob decimal_coefficients." + key, compatibility);
        // The low-field model and its screening cache still use Real.
        p.*(found->second)=compatibility;
    }
    return p;
}
// Read the same parameter family used by the compatibility kernel, with the
// decimal values from PAR overriding the corresponding rounded JSON fields.
inline BasicIalMobilityParameters<long double> extendedParameters(
    const nlohmann::json& overrides, bool electron, const nlohmann::json& decimals) {
    if (!overrides.is_object() || !decimals.is_object() || decimals.size() != overrides.size())
        throw std::invalid_argument("IALMob decimal coefficients must cover the compatibility fields");
    auto p = ialSiliconDefaults<long double>(electron);
    const auto& fields = coefficientFields<long double>;
    for (const auto& [key,value] : overrides.items()) {
        const auto found = fields.find(key);
        if (found == fields.end()) throw std::invalid_argument("Unsupported parameter: " + key);
        p.*(found->second) = decimalCompatibilityNumber(
            decimals.at(key), "IALMob decimal_coefficients." + key, value.get<Real>());
    }
    return p;
}
} // namespace vela::ial_json
