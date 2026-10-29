#pragma once

#include "vela/physics/experimental/BvMethodsPopulationSg.h"

namespace vela::experimental {

// All energies below are explicitly normalized by VT; concentrations by N0.
// This reference is independent of nonlinear row/column equilibration.
struct BvMethodsPrimitiveReference {
    double thermalVoltage_V;
    double concentration_cm3;
    double affinityAtTemperature_eV;
    double energyReference_overVT;
    double baseBandgap_eV;
    bool stableStatistics = false; // Explicit arithmetic policy; preserve prior primitive default.
    void validate() const {
        if(!(thermalVoltage_V>0.) || !(concentration_cm3>0.) || !(baseBandgap_eV>0.) ||
           !std::isfinite(thermalVoltage_V) || !std::isfinite(concentration_cm3) ||
           !std::isfinite(affinityAtTemperature_eV) || !std::isfinite(energyReference_overVT) ||
           !std::isfinite(baseBandgap_eV))throw std::invalid_argument("Invalid BVMethods primitive reference");
    }
};
inline AvJet primitiveStatisticalGamma(const AvJet& population,double dos) {
    if(!(population.value>0.) || !(dos>0.) || !std::isfinite(population.value) || !std::isfinite(dos))
        throw std::domain_error("Invalid primitive population/DOS");
    const double ratio=population.value/dos;
    AvJet gamma(JoyceDixonStatistics::primitiveLogGamma(ratio));
    const double slope=JoyceDixonStatistics::primitiveLogGammaSlope(ratio)/dos;
    for(std::size_t k=0;k<gamma.derivative.size();++k)gamma.derivative[k]=slope*population.derivative[k];
    return gamma;
}
struct BvMethodsPrimitiveCarrier {
    AvJet potential_overVT,population_overN0,logGamma;
    double dos_overN0,conductionOffset_overVT,narrowing_overVT;
};
inline BvMethodsPrimitiveCarrier populationPrimitiveFromNormalized(
    const AvJet& potential_overVT,const AvJet& population_overN0,
    double dos_overN0,double narrowing_overVT,const BvMethodsPrimitiveReference& ref)
{
    ref.validate();
    if(!std::isfinite(potential_overVT.value) || !(population_overN0.value>0.) ||
       !std::isfinite(population_overN0.value) || !(dos_overN0>0.) || !std::isfinite(dos_overN0) ||
       !(narrowing_overVT>=0.) || !std::isfinite(narrowing_overVT))
        throw std::domain_error("Invalid normalized BVMethods primitive state");
    const double offset=ref.affinityAtTemperature_eV/ref.thermalVoltage_V+.5*narrowing_overVT;
    return {potential_overVT,population_overN0,(ref.stableStatistics?primitiveStatisticalGamma(population_overN0,dos_overN0):populationLogGamma(population_overN0,dos_overN0)),
            dos_overN0,offset,narrowing_overVT};
}
inline BvMethodsPrimitiveCarrier populationPrimitiveFromPhysical(
    const AvJet& potential_V,const AvJet& population_cm3,double dos_cm3,
    double narrowing_eV,const BvMethodsPrimitiveReference& ref)
{
    ref.validate();
    return populationPrimitiveFromNormalized(potential_V/AvJet(ref.thermalVoltage_V),
        population_cm3/AvJet(ref.concentration_cm3),dos_cm3/ref.concentration_cm3,
        narrowing_eV/ref.thermalVoltage_V,ref);
}
inline AvJet primitiveSgArgument(const BvMethodsPrimitiveCarrier& from,
    const BvMethodsPrimitiveCarrier& to,bool electron,const AvJet& einstein=AvJet(1))
{
    if(!(einstein.value>0.))throw std::invalid_argument("Invalid primitive Einstein factor");
    const double sign=electron?-1.:1.;
    auto difference=(from.potential_overVT-to.potential_overVT)
                   +AvJet(from.conductionOffset_overVT-to.conductionOffset_overVT);
    if(!electron)difference=difference-AvJet(from.narrowing_overVT-to.narrowing_overVT);
    difference=difference+AvJet(sign)*avLogRatio(AvJet(to.dos_overN0),AvJet(from.dos_overN0));
    const auto argument=AvJet(sign)*difference+(to.logGamma-from.logGamma);
    return (argument+(einstein-AvJet(1))*avLogRatio(to.population_overN0,from.population_overN0))/einstein;
}
inline AvJet primitiveSgParticleFlux(const BvMethodsPrimitiveCarrier& from,
    const BvMethodsPrimitiveCarrier& to,const BvMethodsPrimitiveReference& ref,
    double length_cm,bool electron)
{
    ref.validate();if(!(length_cm>0.) || !std::isfinite(length_cm))throw std::invalid_argument("Invalid primitive edge length");
    const auto difference=populationSgFlux(AvJet(electron?-1.:1.),from.population_overN0,
        to.population_overN0,primitiveSgArgument(from,to,electron));
    return difference*AvJet(ref.concentration_cm3)*AvJet(ref.thermalVoltage_V)/AvJet(length_cm);
}
// A mathematical band/QF relation with explicit origin. Native floating-order
// equivalence is a separate qualification from this reference-frame contract.
inline AvJet primitiveQuasiFermiPotential(const BvMethodsPrimitiveCarrier& node,
    const BvMethodsPrimitiveReference& ref,bool electron)
{
    ref.validate();
    auto band=node.potential_overVT+AvJet(node.conductionOffset_overVT)-AvJet(ref.energyReference_overVT);
    const auto eta=avLogRatio(node.population_overN0,AvJet(node.dos_overN0))-node.logGamma;
    if(electron)band=band-eta;
    else band=band+AvJet(ref.baseBandgap_eV/ref.thermalVoltage_V)-AvJet(node.narrowing_overVT)+eta;
    return band*AvJet(ref.thermalVoltage_V);
}
inline std::array<AvJet,2> primitiveAvCurrent(const AvTri3Geometry& g,
    const std::array<BvMethodsPrimitiveCarrier,3>& nodes,
    const BvMethodsPrimitiveReference& ref,bool electron)
{
    if(g.representative>=3)throw std::invalid_argument("Invalid primitive AV representative");
    std::array<std::array<double,2>,2> t{};std::array<AvJet,2> j{};
    std::size_t edge=0;const auto r=g.representative;
    for(std::size_t i=0;i<3;++i)if(i!=r) {
        const double dx=g.points_cm[i][0]-g.points_cm[r][0],dy=g.points_cm[i][1]-g.points_cm[r][1];
        const double length=std::hypot(dx,dy);
        if(!(length>0.))throw std::invalid_argument("Degenerate primitive AV edge");
        t[edge]={dx/length,dy/length};
        auto u=r,v=i;const bool reverse=g.nodeOrder[u]>g.nodeOrder[v];if(reverse)std::swap(u,v);
        j[edge]=primitiveSgParticleFlux(nodes[u],nodes[v],ref,length,electron);
        if(reverse)j[edge]=-j[edge];++edge;
    }
    const double determinant=t[0][0]*t[1][1]-t[0][1]*t[1][0];
    if(std::abs(determinant)<1e-14)throw std::invalid_argument("Degenerate primitive AV triangle");
    const auto x=(j[0]*AvJet(t[1][1])-AvJet(t[0][1])*j[1])/AvJet(determinant);
    const auto y=(AvJet(t[0][0])*j[1]-j[0]*AvJet(t[1][0]))/AvJet(determinant);
    return {AvJet(g.currentProjection[0])*x+AvJet(g.currentProjection[1])*y,
            AvJet(g.currentProjection[2])*x+AvJet(g.currentProjection[3])*y};
}
} // namespace vela::experimental
