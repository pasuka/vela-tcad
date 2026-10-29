#pragma once
#include "vela/numerics/PopulationStepControl.h"
#include "vela/solver/SlotPopulationNumerics.h"
#include "vela/equation/CoupledDDAssembler.h"

namespace vela::population_step {
struct Frame {
    VectorXd n,p,etaN,etaP,tangentN,tangentP,relativeN,relativeP;
    std::vector<bool> active,nFree,pFree;
    Real scaledVt=1.;
    bool fermi=false;
};
inline Frame frame(const CoupledDDAssembler& a,const VectorXd& x,
                   const VectorXd& direction,const CoupledDDBoundaryConditions& bc) {
    const int N=a.numNodes();
    if(x.size()!=3*N||direction.size()!=x.size()||!direction.allFinite())
        throw std::invalid_argument("invalid population step coordinates");
    Frame f;f.n=a.electronDensity(x);f.p=a.holeDensity(x);
    f.fermi=usesFermiDirac(a.carrierStatistics());
    f.scaledVt=a.populationThermalVoltage_V()/(a.usesScaledState()?a.potentialScale():1.);
    f.active.resize(N);f.nFree.resize(N);f.pFree.resize(N);
    f.etaN=f.etaP=f.relativeN=f.relativeP=VectorXd::Zero(N);
    f.tangentN=(direction.head(N)-direction.segment(N,N))/f.scaledVt;
    f.tangentP=(direction.tail(N)-direction.head(N))/f.scaledVt;
    for(int i=0;i<N;++i) {
        const bool active=a.intrinsicDensity()[i]>0.;
        f.active[i]=active;f.nFree[i]=active&&!bc.phin.contains(i);f.pFree[i]=active&&!bc.phip.contains(i);
        if(!active)continue;
        if(!(f.n[i]>0.)||!(f.p[i]>0.)||!std::isfinite(f.n[i])||!std::isfinite(f.p[i]))
            throw std::domain_error("invalid population before global step");
        if(f.fermi) {
            f.etaN[i]=slot_population::populationEta(f.n[i],a.electronDensityOfStates()[i]);
            f.etaP[i]=slot_population::populationEta(f.p[i],a.holeDensityOfStates()[i]);
        }
        if(f.nFree[i])f.relativeN[i]=f.tangentN[i]*(f.fermi?slot_population::populationLogSlope(f.etaN[i]):1.);
        if(f.pFree[i])f.relativeP[i]=f.tangentP[i]*(f.fermi?slot_population::populationLogSlope(f.etaP[i]):1.);
    }
    return f;
}
inline VectorXd finiteStep(const CoupledDDAssembler& a,const Frame& f,
                           const VectorXd& direction,Real alpha) {
    if(!(alpha>=0.)||!std::isfinite(alpha))throw std::domain_error("invalid common population damping");
    const int N=a.numNodes();VectorXd mapped=alpha*direction;
    for(int b=0;b<2;++b)for(int i=0;i<N;++i) {
        if(!(b?f.pFree[i]:f.nFree[i]))continue;
        const double rel=alpha*(b?f.relativeP[i]:f.relativeN[i]);
        if(!(rel>-1.)||!std::isfinite(rel))throw std::domain_error("global population step lost positivity");
        const double t=alpha*(b?f.tangentP[i]:f.tangentN[i]);
        double remainder;
        if(!f.fermi)remainder=log1pRemainder(rel);
        else if(std::abs(t)<=1e-4)remainder=slot_population::inverseTangentRemainder(b?f.etaP[i]:f.etaN[i],t);
        else remainder=slot_population::relativeInversePopulation(b?f.p[i]:f.n[i],
            b?a.holeDensityOfStates()[i]:a.electronDensityOfStates()[i],rel,true)-t;
        mapped[(b+1)*N+i]+=(b?1.:-1.)*f.scaledVt*remainder;
    }
    return mapped;
}
}
