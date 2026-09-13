#include "vela/physics/SiliconThermalPhysics.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/core/PhysicalConstants.h"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace vela {
namespace {
using D=ThermalQuantity;
D scalar(Real x) {D r;r.value=x;return r;}
D variable(Real x,int i) {D r=scalar(x);r.derivative[i]=1.;return r;}
D operator+(const D&a,const D&b) {D r=scalar(a.value+b.value);for(int i=0;i<4;++i)r.derivative[i]=a.derivative[i]+b.derivative[i];return r;}
D operator-(const D&a,const D&b) {D r=scalar(a.value-b.value);for(int i=0;i<4;++i)r.derivative[i]=a.derivative[i]-b.derivative[i];return r;}
D operator*(const D&a,const D&b) {D r=scalar(a.value*b.value);for(int i=0;i<4;++i)r.derivative[i]=a.derivative[i]*b.value+a.value*b.derivative[i];return r;}
D operator/(const D&a,const D&b) {D r=scalar(a.value/b.value);for(int i=0;i<4;++i)r.derivative[i]=(a.derivative[i]-r.value*b.derivative[i])/b.value;return r;}
D power(const D&a,Real b) {D r=scalar(std::pow(a.value,b));for(int i=0;i<4;++i)r.derivative[i]=b*std::pow(a.value,b-1.)*a.derivative[i];return r;}
D exp(const D&a) {D r=scalar(std::exp(a.value));for(int i=0;i<4;++i)r.derivative[i]=r.value*a.derivative[i];return r;}
D expm1(const D&a) {D r=scalar(std::expm1(a.value));for(int i=0;i<4;++i)r.derivative[i]=std::exp(a.value)*a.derivative[i];return r;}
D log(const D&a) {D r=scalar(std::log(a.value));for(int i=0;i<4;++i)r.derivative[i]=a.derivative[i]/a.value;return r;}
D fermi(const D&a) {D r=scalar(fermiDiracHalf(a.value));const Real slope=fermiDiracHalfDerivative(a.value);for(int i=0;i<4;++i)r.derivative[i]=slope*a.derivative[i];return r;}
template<std::size_t N>D polynomial(const std::array<Real,N>& c,const D&t) {D r=scalar(c.back());for(std::size_t i=N-1;i>0;--i)r=r*t+scalar(c[i-1]);return r;}
std::array<D,3> bands(const SiliconThermalParameters&p,const D&t) {
    const Real eg0=p.bandgap300_eV+p.bandgapAlpha_eV_per_K*90000./(p.bandgapBeta_K+300.);
    D eg=scalar(eg0)-scalar(p.bandgapAlpha_eV_per_K)*t*t/(t+scalar(p.bandgapBeta_K));
    D me=power(power(scalar(6.*p.electronTransverseMass*eg0)/eg,2.)*scalar(p.electronLongitudinalMass),1./3.)+scalar(p.electronMassOffset);
    D mh=power(polynomial(p.holeMassNumerator,t)/polynomial(p.holeMassDenominator,t),2./3.)+scalar(p.holeMassOffset);
    const Real c=2.*std::pow(2.*std::numbers::pi*constants::m0*constants::kb/(constants::h*constants::h),1.5);
    if (!(eg.value>0.)||!(me.value>0.)||!(mh.value>0.))throw std::invalid_argument("Invalid thermal band/DOS state");
    return {eg,scalar(c)*power(me*t,1.5),scalar(c)*power(mh*t,1.5)};
}
void positive(Real x){if(!std::isfinite(x)||x<=0.)throw std::invalid_argument("Thermal silicon requires finite positive scales");}
void validateState(const SiliconThermalState&s){
    for(Real x:{s.potential_V,s.electronQf_V,s.holeQf_V,s.temperature_K,s.donors_m3,s.acceptors_m3,s.electronQfReference_V,s.holeQfReference_V})
        if(!std::isfinite(x))throw std::invalid_argument("Thermal silicon state must be finite");
    if(s.temperature_K<50.||s.donors_m3<0.||s.acceptors_m3<0.)throw std::invalid_argument("Thermal silicon requires T>=50 K and nonnegative doping");
}
}
SiliconThermalPhysics::SiliconThermalPhysics(SiliconThermalParameters p):p_(p) {
    for(Real x:{p.bandgap300_eV,p.affinity300_eV,p.bandgapBeta_K,p.bgnReference_m3,p.electronTransverseMass,p.electronLongitudinalMass})positive(x);
    for(Real x:{p.bandgapAlpha_eV_per_K,p.bgnToAffinity,p.bgnCoefficient_eV,p.electronMassOffset,p.holeMassOffset})
        if(!std::isfinite(x)||x<0.)throw std::invalid_argument("Invalid thermal silicon coefficient");
    if(p.bgnToAffinity>1.)throw std::invalid_argument("BGN affinity fraction must be <=1");
    for(int i=0;i<2;++i){positive(p.srhTau300_s[i]);positive(p.srhReference_m3[i]);positive(p.augerReference_m3[i]);
        if(!std::isfinite(p.srhTemperatureExponent[i])||!std::isfinite(p.augerEnhancement[i])||p.augerEnhancement[i]<0.)throw std::invalid_argument("Invalid thermal recombination coefficient");}
    for(const auto& c:{p.holeMassNumerator,p.holeMassDenominator})for(Real x:c)if(!std::isfinite(x))throw std::invalid_argument("Invalid DOS coefficient");
    for(const auto& c:{p.augerElectron_m6_per_s,p.augerHole_m6_per_s})for(Real x:c)if(!std::isfinite(x))throw std::invalid_argument("Invalid Auger polynomial");
    const auto b=bands(p,scalar(300.));Nc300_=b[1].value;Nv300_=b[2].value;
    reference_=p.affinity300_eV+p.bandgap300_eV/2.+constants::Vt_300/2.*std::log(Nc300_/Nv300_);
}
SiliconThermalPhysics::DopingPreparation SiliconThermalPhysics::prepareDoping(Real donors,Real acceptors) const {
    if(!std::isfinite(donors)||!std::isfinite(acceptors)||donors<0.||acceptors<0.)
        throw std::invalid_argument("Invalid preparation doping");
    const Real doping=donors+acceptors;
    Real delta=0.;if(doping>0.){const Real x=std::log(doping/p_.bgnReference_m3);delta=p_.bgnCoefficient_eV*(x+std::hypot(x,std::sqrt(.5)));}
    if(p_.fermiBgnCorrection)delta+=fermiStatisticsBandgapCorrection(donors,acceptors,Nc300_,Nv300_,constants::Vt_300);
    if(!std::isfinite(delta))throw std::invalid_argument("Nonfinite prepared BGN");
    DopingPreparation prepared;prepared.donors_=donors;prepared.acceptors_=acceptors;prepared.delta_=delta;prepared.owner_=identity_;return prepared;
}
SiliconThermalPhysics::TemperaturePreparation SiliconThermalPhysics::prepareTemperature(Real temperature,const DopingPreparation& dopingState) const {
    if(dopingState.owner_!=identity_)throw std::invalid_argument("Foreign doping preparation");
    if(!std::isfinite(temperature)||temperature<50.)throw std::invalid_argument("Invalid preparation temperature");
    TemperaturePreparation prepared;prepared.owner_=identity_;prepared.temperature_=temperature;
    prepared.donors_=dopingState.donors_;prepared.acceptors_=dopingState.acceptors_;
    const D t=variable(temperature,3),vt=scalar(constants::kb/constants::q)*t,ratio=t/scalar(300.);
    const auto b=bands(p_,t);auto& r=prepared.base_;r.bandgap_eV=b[0];r.Nc_m3=b[1];r.Nv_m3=b[2];
    const Real delta=dopingState.delta_,doping=prepared.donors_+prepared.acceptors_;
    r.bandgapNarrowing_eV=delta;
    const D effectiveEg=b[0]-scalar(delta);
    prepared.vt_=vt;prepared.effectiveEg_=effectiveEg;
    r.affinity_eV=scalar(p_.affinity300_eV)+(scalar(p_.bandgap300_eV)-b[0])*scalar(.5)+scalar(p_.bgnToAffinity*delta);
    r.ni_m3=power(b[1]*b[2],.5)*exp(scalar(-.5)*b[0]/vt);
    r.effectiveNi_m3=power(b[1]*b[2],.5)*exp(scalar(-.5)*effectiveEg/vt);
    r.electronLifetime_s=scalar(p_.srhTau300_s[0]/(1.+doping/p_.srhReference_m3[0]))*power(ratio,p_.srhTemperatureExponent[0]);
    r.holeLifetime_s=scalar(p_.srhTau300_s[1]/(1.+doping/p_.srhReference_m3[1]))*power(ratio,p_.srhTemperatureExponent[1]);
    r.augerElectron_m6_per_s=polynomial(p_.augerElectron_m6_per_s,ratio);
    r.augerHole_m6_per_s=polynomial(p_.augerHole_m6_per_s,ratio);
    return prepared;
}
SiliconThermalPhysics::CarrierEvaluation SiliconThermalPhysics::carriers(const SiliconThermalState&s,const TemperaturePreparation& prepared) const {
    validateState(s);
    if(prepared.owner_!=identity_ || !prepared.matches(s.temperature_K,s.donors_m3,s.acceptors_m3))
        throw std::invalid_argument("Foreign or stale thermal preparation");
    CarrierEvaluation evaluation;evaluation.result=prepared.base_;auto& r=evaluation.result;
    const D psi=variable(s.potential_V,0),fn=variable(s.electronQf_V,1),fp=variable(s.holeQf_V,2);
    const D& vt=prepared.vt_;const D& effectiveEg=prepared.effectiveEg_;
    r.conductionBand_eV=scalar(reference_)-psi-r.affinity_eV;
    r.valenceBand_eV=r.conductionBand_eV-effectiveEg;
    const D en=s.electronQfReference_V==0.?(scalar(0.)-fn-r.conductionBand_eV)/vt:
        ((psi-scalar(s.electronQfReference_V))-fn-scalar(reference_)+r.affinity_eV)/vt;
    const D ep=s.holeQfReference_V==0.?(r.valenceBand_eV+fp)/vt:
        ((scalar(s.holeQfReference_V)-psi)+fp+scalar(reference_)-r.affinity_eV-effectiveEg)/vt;
    r.electronEta=en;r.holeEta=ep;
    const D fN=fermi(en),fP=fermi(ep);
    r.electrons_m3=r.Nc_m3*fN;r.holes_m3=r.Nv_m3*fP;
    evaluation.fN=fN;evaluation.fP=fP;
    for(const auto& q:{r.electrons_m3,r.holes_m3}){
        if(!std::isfinite(q.value))throw std::runtime_error("Nonfinite thermal density");
        for(Real d:q.derivative)if(!std::isfinite(d))throw std::runtime_error("Nonfinite thermal density derivative");
    }
    return evaluation;
}
std::array<ThermalQuantity,2> SiliconThermalPhysics::carrierDensities(const SiliconThermalState&s,const TemperaturePreparation& prepared) const {
    const auto result=carriers(s,prepared).result;return {result.electrons_m3,result.holes_m3};
}
SiliconThermalResult SiliconThermalPhysics::evaluate(const SiliconThermalState&s) const {
    validateState(s);return evaluate(s,prepareTemperature(s.temperature_K,prepareDoping(s.donors_m3,s.acceptors_m3)));
}
SiliconThermalResult SiliconThermalPhysics::evaluate(const SiliconThermalState&s,const TemperaturePreparation& prepared) const {
    auto evaluation=carriers(s,prepared);auto& r=evaluation.result;
    const D fn=variable(s.electronQf_V,1),fp=variable(s.holeQf_V,2);
    const auto& vt=prepared.vt_;const auto& en=r.electronEta;const auto& ep=r.holeEta;
    const auto& fN=evaluation.fN;const auto& fP=evaluation.fP;
    const D gammaN=exp(log(fN)-en),gammaP=exp(log(fP)-ep);
    const D splitting=(scalar(s.holeQfReference_V-s.electronQfReference_V)+fp-fn)/vt;
    const D excess=splitting.value>=0.?r.electrons_m3*r.holes_m3*(scalar(0.)-expm1(scalar(0.)-splitting)):
        gammaN*gammaP*r.effectiveNi_m3*r.effectiveNi_m3*expm1(splitting);
    r.srhRate_m3_per_s=excess/(r.holeLifetime_s*(r.electrons_m3+gammaN*r.effectiveNi_m3)+
        r.electronLifetime_s*(r.holes_m3+gammaP*r.effectiveNi_m3));
    r.augerElectron_m6_per_s=r.augerElectron_m6_per_s*(scalar(1.)+scalar(p_.augerEnhancement[0])*exp(scalar(0.)-r.electrons_m3/scalar(p_.augerReference_m3[0])));
    r.augerHole_m6_per_s=r.augerHole_m6_per_s*(scalar(1.)+scalar(p_.augerEnhancement[1])*exp(scalar(0.)-r.holes_m3/scalar(p_.augerReference_m3[1])));
    if(r.augerElectron_m6_per_s.value<0.||r.augerHole_m6_per_s.value<0.)throw std::invalid_argument("Auger polynomial negative at requested temperature");
    r.augerRate_m3_per_s=(r.augerElectron_m6_per_s*r.electrons_m3+r.augerHole_m6_per_s*r.holes_m3)*excess;
    // The inactive branch has zero derivatives; at equilibrium select its
    // one-sided derivative. WithGeneration retains the signed legacy law.
    if(!p_.augerWithGeneration && splitting.value<=0.)r.augerRate_m3_per_s=scalar(0.);
    for(const auto& q:{r.bandgap_eV,r.affinity_eV,r.Nc_m3,r.Nv_m3,r.ni_m3,r.effectiveNi_m3,r.electrons_m3,r.holes_m3,r.srhRate_m3_per_s,r.augerRate_m3_per_s}) {
        if(!std::isfinite(q.value))throw std::runtime_error("Nonfinite thermal silicon result");
        for(Real d:q.derivative)if(!std::isfinite(d))throw std::runtime_error("Nonfinite thermal silicon derivative");
    }
    return r;
}
} // namespace vela
