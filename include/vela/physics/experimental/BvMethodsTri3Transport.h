#pragma once

#include "vela/discretization/Tri3Operators.h"
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
    const std::array<std::size_t,3>& nodeOrder={0,1,2},
    const std::array<bool,3>& freezePopulationFeedback={false,false,false})
{
    if(!(thermalVoltage_V>0.))throw std::invalid_argument("Transport requires positive thermal voltage");
    Tri3TransportResult out;
    for (std::size_t k = 0; k < 3; ++k) {
        if(!(vertexAreas_cm2[k]>=0.))throw std::invalid_argument("Negative transport averaging area");
        if(!(density_cm3[k].value>0.) || !(vertexMobility_cm2_V_s[k].value>=0.))
            throw std::domain_error("Invalid transport density or mobility");
    }
    out.elementMobility_cm2_V_s = detail::tri3VertexAverage<AvJet>(
        vertexAreas_cm2, [&](std::size_t vertex) { return vertexMobility_cm2_V_s[vertex]; });
    for(std::size_t k=0;k<3;++k){
        if(!(edgeCoefficients[k]>=0.) || !(einsteinFactor[k].value>0.))
            throw std::invalid_argument("Invalid transport coefficient or Einstein factor");
        const auto j=(k+1)%3;
        const auto thermal=AvJet(thermalVoltage_V)*einsteinFactor[k];
        AvJet flux;
        auto prefactor=AvJet(edgeCoefficients[k])*out.elementMobility_cm2_V_s*thermal;
        // Population slots are 3..8, potential slots are 0..2. Freeze the
        // coefficient AND its argument input as a pair, retaining all values.
        if(freezePopulationFeedback[k])
            for(int z=3;z<9;++z)prefactor.derivative[z]=0.;
        if(!directArguments){
            auto delta=AvJet(electron?1.:-1.)*(quasiFermi_V[j]-quasiFermi_V[k])/thermal;
            const auto logRatio=avLogRatio(density_cm3[j],density_cm3[k]);
            auto argument=logRatio+delta;
            if(freezePopulationFeedback[k])for(int z=3;z<9;++z){
                argument.derivative[z]=0.;
                // avSgFlux's expm1 branch uses delta, with a=logRatio+delta.
                delta.derivative[z]=-logRatio.derivative[z];
            }
            flux=avSgFlux(prefactor,density_cm3[k],density_cm3[j],delta,argument);
        }
        if(directArguments) {
            const bool reverse=nodeOrder[k]>nodeOrder[j];const auto u=reverse?j:k,v=reverse?k:j;
            auto argument=(*directArguments)[k];
            if(freezePopulationFeedback[k])for(int z=3;z<9;++z)argument.derivative[z]=0.;
            flux=populationSgFlux(prefactor,density_cm3[u],density_cm3[v],argument);
            if(reverse)flux=-flux;
        }
        out.edgeFlux_per_cm_s[k]=flux;
        out.nodeFlux_per_cm_s[k]=out.nodeFlux_per_cm_s[k]+flux;
        out.nodeFlux_per_cm_s[j]=out.nodeFlux_per_cm_s[j]-flux;
    }
    return out;
}

// The population columns for different vertices have disjoint support.
// Reproduce the dual log-ratio branch without constructing either jet.
inline double populationTransportLogRatio(double to,double from,double toTangent,double fromTangent)
{
    if(!(to>0.) || !(from>0.))throw std::domain_error("AV log requires positive density");
    const double ratio=to/from;
    if(avRatioInOrdinaryRange(ratio) && std::isfinite(toTangent/from) &&
       std::isfinite((0.-ratio*fromTangent)/from))
        return std::log(ratio);
    return std::log(to)-std::log(from);
}

inline double tri3BernoulliValue(double x)
{
    if(std::abs(x)<1e-5) {
        const double z=x*x;
        return 1.-x/2.+z/12.-z*z/720.+z*z*z/30240.;
    }
    if(x>50.) {
        const double e=std::exp(-x);
        return x*e/(1.-e);
    }
    if(x<-50.)return -x/(1.-std::exp(x));
    return x/std::expm1(x);
}

struct Tri3TransportValues {
    std::array<double,3> edges{},nodes{};
};

inline Tri3TransportValues evaluateBvMethodsTri3TransportValues(
    const std::array<double,3>& areas,const std::array<double,3>& edgeCoefficients,
    const std::array<double,3>& density,const std::array<double,3>& qf,
    const std::array<double,3>& mobility,double thermalVoltage,bool electron,
    const std::array<double,3>& einstein,const std::array<double,3>& logRatio)
{
    if(!(thermalVoltage>0.))throw std::invalid_argument("Transport requires positive thermal voltage");
    for (std::size_t k = 0; k < 3; ++k) {
        if(!(areas[k]>=0.))throw std::invalid_argument("Negative transport averaging area");
        if(!(density[k]>0.) || !(mobility[k]>=0.))
            throw std::domain_error("Invalid transport density or mobility");
    }
    const double average = detail::tri3VertexAverage<double>(
        areas, [&](std::size_t vertex) { return mobility[vertex]; });
    Tri3TransportValues out;
    for(std::size_t k=0;k<3;++k) {
        if(!(edgeCoefficients[k]>=0.) || !(einstein[k]>0.))
            throw std::invalid_argument("Invalid transport coefficient or Einstein factor");
        const auto j=(k+1)%3;
        const double thermal=thermalVoltage*einstein[k];
        const double coefficient=edgeCoefficients[k]*average*thermal;
        const double delta=(electron?1.:-1.)*(qf[j]-qf[k])/thermal;
        const double argument=logRatio[k]+delta;
        const double flux = detail::sgBalancedEdgeFlux(coefficient, density[k], density[j],
            delta, delta, argument, tri3BernoulliValue, [](double x) { return std::expm1(x); });
        out.edges[k]=flux;
        out.nodes[k]=out.nodes[k]+flux;
        out.nodes[j]=out.nodes[j]-flux;
    }
    return out;
}

} // namespace vela::experimental
