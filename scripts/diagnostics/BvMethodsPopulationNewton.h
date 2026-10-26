#pragma once
// Diagnostic coordinates and feasibility policy; no production dispatch.
#include "vela/core/Types.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace vela::diagnostics {
struct PopulationCoordinates {
    SparseMatrixd qfFromScaledPopulation;
    VectorXd psiScale_V, electronScale_cm3, holeScale_cm3;
    double voltageScale_V;
};
inline PopulationCoordinates populationCoordinates(
    const VectorXd& psi, const VectorXd& electrons, const VectorXd& holes,
    const VectorXd& electronSlope, const VectorXd& holeSlope,
    const std::vector<bool>& active, double potentialScale, double thermalVoltage,
    double voltage, double densityReference=1e10) {
    const int n=psi.size();
    if(n==0||electrons.size()!=n||holes.size()!=n||electronSlope.size()!=n||holeSlope.size()!=n||
       active.size()!=static_cast<std::size_t>(n)||!(potentialScale>0.)||!(thermalVoltage>0.)||!(densityReference>0.))
        throw std::invalid_argument("invalid population coordinate inputs");
    PopulationCoordinates out;
    out.psiScale_V=psi.cwiseAbs().array()+thermalVoltage;
    out.electronScale_cm3=electrons.cwiseAbs().array()+densityReference;
    out.holeScale_cm3=holes.cwiseAbs().array()+densityReference;
    out.voltageScale_V=std::abs(voltage)+thermalVoltage;
    std::vector<Eigen::Triplet<double>> t;
    for(int i=0;i<n;++i) {
        const double a=out.psiScale_V[i]/potentialScale;
        t.emplace_back(i,i,a);
        if(active[i]) {
            if(!(electrons[i]>0.)||!(holes[i]>0.)||!(electronSlope[i]>0.)||!(holeSlope[i]>0.))
                throw std::domain_error("population coordinates require positive densities and slopes");
            const double b=out.electronScale_cm3[i]/electronSlope[i]/potentialScale;
            const double c=out.holeScale_cm3[i]/holeSlope[i]/potentialScale;
            if(!std::isfinite(b)||!std::isfinite(c))throw std::domain_error("nonfinite population coordinate derivative");
            t.emplace_back(n+i,i,a);t.emplace_back(n+i,n+i,-b);
            t.emplace_back(2*n+i,i,a);t.emplace_back(2*n+i,2*n+i,c);
        } else {
            // Oxide/gauge slots remain QF coordinates; no fictitious population.
            t.emplace_back(n+i,n+i,1.);t.emplace_back(2*n+i,2*n+i,1.);
        }
    }
    t.emplace_back(3*n,3*n,out.voltageScale_V);
    out.qfFromScaledPopulation.resize(3*n+1,3*n+1);
    out.qfFromScaledPopulation.setFromTriplets(t.begin(),t.end());
    return out;
}
struct FractionToBoundary {
    double alpha=1.;
    int limitingNode=-1, limitingCarrier=-1;
};
inline FractionToBoundary populationFractionToBoundary(
    const VectorXd& n,const VectorXd& p,const VectorXd& dn,const VectorXd& dp,
    const std::vector<bool>& active,double tau) {
    if(!(tau>0.&&tau<1.)||n.size()!=p.size()||n.size()!=dn.size()||n.size()!=dp.size()||
       active.size()!=static_cast<std::size_t>(n.size()))throw std::invalid_argument("invalid fraction-to-boundary inputs");
    FractionToBoundary out;
    for(int i=0;i<n.size();++i)if(active[i])for(int b=0;b<2;++b) {
        const double c=b?p[i]:n[i],dc=b?dp[i]:dn[i];
        if(!(c>0.)||!std::isfinite(c)||!std::isfinite(dc))throw std::domain_error("invalid positive population direction");
        if(dc<0.) {const double bound=tau*(c/(-dc));if(bound<out.alpha){out.alpha=bound;out.limitingNode=i;out.limitingCarrier=b;}}
    }
    if(!(out.alpha>0.)||!std::isfinite(out.alpha))throw std::domain_error("no representable feasible population step");
    return out;
}
struct EquilibratedPopulationSystem {
    SparseMatrixd matrix;
    VectorXd rowWeights,columnWeights;
};
inline EquilibratedPopulationSystem equilibratePopulationSystem(const SparseMatrixd& a) {
    if(a.rows()==0||a.rows()!=a.cols())throw std::invalid_argument("square population system required");
    EquilibratedPopulationSystem out;out.matrix=a;
    out.rowWeights=VectorXd::Zero(a.rows());out.columnWeights=VectorXd::Zero(a.cols());
    for(int k=0;k<a.outerSize();++k)for(SparseMatrixd::InnerIterator v(a,k);v;++v) {
        if(!std::isfinite(v.value()))throw std::domain_error("nonfinite population matrix");
        out.rowWeights[v.row()]=std::max(out.rowWeights[v.row()],std::abs(v.value()));
    }
    for(int i=0;i<a.rows();++i){if(!(out.rowWeights[i]>0.))throw std::domain_error("zero population matrix row");out.rowWeights[i]=1./out.rowWeights[i];}
    for(int k=0;k<out.matrix.outerSize();++k)for(SparseMatrixd::InnerIterator v(out.matrix,k);v;++v){v.valueRef()*=out.rowWeights[v.row()];out.columnWeights[v.col()]=std::max(out.columnWeights[v.col()],std::abs(v.value()));}
    for(int i=0;i<a.cols();++i){if(!(out.columnWeights[i]>0.))throw std::domain_error("zero population matrix column");out.columnWeights[i]=1./out.columnWeights[i];}
    for(int k=0;k<out.matrix.outerSize();++k)for(SparseMatrixd::InnerIterator v(out.matrix,k);v;++v)v.valueRef()*=out.columnWeights[v.col()];
    return out;
}
inline bool restorationNoWorse(double candidateMerit,double candidateRms,double baselineMerit,double baselineRms) {
    return std::isfinite(candidateMerit)&&std::isfinite(candidateRms)&&candidateMerit<=baselineMerit&&candidateRms<=baselineRms;
}
inline bool restorationWithinTargets(double equationInf,double loadResidual,double rawRms,double equationTolerance,double loadTolerance) {
    return std::isfinite(equationInf)&&std::isfinite(loadResidual)&&std::isfinite(rawRms)&&
        equationInf<=equationTolerance&&std::abs(loadResidual)<=loadTolerance&&rawRms<1.;
}
inline bool derivativeBootstrapRequired(bool valid,double alpha,double merit,double beforeMerit,bool withinTargets) {
    return !valid||(!withinTargets&&(alpha<1.||!std::isfinite(merit)||merit>beforeMerit));
}
} // namespace vela::diagnostics
