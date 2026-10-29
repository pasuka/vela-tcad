#pragma once

#include "vela/physics/experimental/BvMethodsPopulationSg.h"

namespace vela::experimental {

// Diagnostic constitutive contract, with no production configuration dispatch.
// Coefficients are dimensionless element-edge box coefficients; areas average
// vertex mobilities only. Fluxes are signed particle rates per cm of depth/s.
struct Tri3TransportResult {
    AvJet elementMobility_cm2_V_s;
    std::array<AvJet,3> edgeFlux_per_cm_s;
    std::array<AvJet,3> nodeFlux_per_cm_s;
};

inline Tri3TransportResult evaluateBvMethodsTri3Transport(
    const std::array<double,3>& vertexAreas_cm2,
    const std::array<double,3>& edgeCoefficients,
    const std::array<AvJet,3>& density_cm3,
    const std::array<AvJet,3>& quasiFermi_V,
    const std::array<AvJet,3>& vertexMobility_cm2_V_s,
    double thermalVoltage_V,bool electron,
    const std::array<AvJet,3>& einsteinFactor={AvJet(1),AvJet(1),AvJet(1)},
    const std::array<AvJet,3>* directArguments=nullptr,
    const std::array<std::size_t,3>& nodeOrder={0,1,2})
{
    if(!(thermalVoltage_V>0.))throw std::invalid_argument("Transport requires positive thermal voltage");
    Tri3TransportResult out;
    double area=0.;
    for(std::size_t k=0;k<3;++k){
        if(!(vertexAreas_cm2[k]>=0.))throw std::invalid_argument("Negative transport averaging area");
        if(!(density_cm3[k].value>0.) || !(vertexMobility_cm2_V_s[k].value>=0.))
            throw std::domain_error("Invalid transport density or mobility");
        area+=vertexAreas_cm2[k];
        out.elementMobility_cm2_V_s=out.elementMobility_cm2_V_s
            +AvJet(vertexAreas_cm2[k])*vertexMobility_cm2_V_s[k];
    }
    if(!(area>0.))throw std::invalid_argument("Transport requires positive averaging area");
    out.elementMobility_cm2_V_s=out.elementMobility_cm2_V_s/AvJet(area);
    for(std::size_t k=0;k<3;++k){
        if(!(edgeCoefficients[k]>=0.) || !(einsteinFactor[k].value>0.))
            throw std::invalid_argument("Invalid transport coefficient or Einstein factor");
        const auto j=(k+1)%3;
        const auto thermal=AvJet(thermalVoltage_V)*einsteinFactor[k];
        AvJet flux;
        if(!directArguments){
            const auto delta=AvJet(electron?1.:-1.)*(quasiFermi_V[j]-quasiFermi_V[k])/thermal;
            const auto argument=avLogRatio(density_cm3[j],density_cm3[k])+delta;
            flux=avSgFlux(AvJet(edgeCoefficients[k])*out.elementMobility_cm2_V_s
                *thermal,density_cm3[k],density_cm3[j],delta,argument);
        }
        if(directArguments) {
            const bool reverse=nodeOrder[k]>nodeOrder[j];const auto u=reverse?j:k,v=reverse?k:j;
            flux=populationSgFlux(AvJet(edgeCoefficients[k])*out.elementMobility_cm2_V_s
                *thermal,density_cm3[u],density_cm3[v],(*directArguments)[k]);
            if(reverse)flux=-flux;
        }
        out.edgeFlux_per_cm_s[k]=flux;
        out.nodeFlux_per_cm_s[k]=out.nodeFlux_per_cm_s[k]+flux;
        out.nodeFlux_per_cm_s[j]=out.nodeFlux_per_cm_s[j]-flux;
    }
    return out;
}

} // namespace vela::experimental
