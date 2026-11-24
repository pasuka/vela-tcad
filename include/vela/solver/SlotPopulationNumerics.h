#pragma once
#include "vela/physics/CarrierStatistics.h"
#include "vela/core/Types.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <limits>
namespace vela::slot_population {
enum class AcceptancePolicy { PrimaryOnly, PrimaryWithGuards };
// Diagnostic population constraint in internal potential coordinates. With the
// fixed DOS/bandgap used by SLOT, dp=dn=0 means dQF=dpsi, not dQF=0.
inline SparseMatrixd populationConstraintMatrix(const SparseMatrixd& a,
                                                const std::vector<int>& rows,
                                                int width) {
 if(width<=0 || a.rows()!=3*width || a.cols()!=3*width)
  throw std::invalid_argument("population constraint requires three nodal blocks");
 std::vector<bool> held(a.rows(),false);
 for(int i:rows){if(i<width || i>=3*width || held[i])
   throw std::invalid_argument("invalid or duplicate population constraint");held[i]=true;}
 std::vector<Eigen::Triplet<Real>> entries;entries.reserve(a.nonZeros()+2*rows.size());
 for(int j=0;j<a.outerSize();++j)for(SparseMatrixd::InnerIterator it(a,j);it;++it)
  if(!held[it.row()])entries.emplace_back(it.row(),it.col(),it.value());
 for(int i:rows){entries.emplace_back(i,i,1.);entries.emplace_back(i,i%width,-1.);}
 SparseMatrixd result(a.rows(),a.cols());result.setFromTriplets(entries.begin(),entries.end());
 return result;
}
inline bool retainPopulationHold(bool previouslyProjected,double unconstrainedRelative) {
 if(!std::isfinite(unconstrainedRelative))throw std::invalid_argument("nonfinite release direction");
 return previouslyProjected && unconstrainedRelative<=0.;
}
// Numerical validity, inverse mapping, linear quality and persistence guards
// remain mandatory at the caller. This selects only nonlinear acceptance.
inline bool accepts(double rho, double rawError, bool closure, bool boundary,
                    AcceptancePolicy policy) {
 if(!std::isfinite(rho) || !std::isfinite(rawError) || rho<0 || rawError<0)
  return false;
 return (rho<=1e-10 || rawError<1.) &&
        (policy==AcceptancePolicy::PrimaryOnly || (closure && boundary));
}
inline double fermiChange(double eta,double d) {
 if(!std::isfinite(eta) || !std::isfinite(d) || std::abs(d)>1024.)
  throw std::runtime_error("relative Fermi integration outside bounded range");
 static constexpr double x[]={.1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
 static constexpr double w[]={.3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
 const int parts=std::max(1,int(std::ceil(std::abs(d)/.25)));
 if(parts>4096)throw std::runtime_error("diagnostic integration range exceeded");
 long double sum=0;
 for(int p=0;p<parts;++p){double width=d/parts,mid=eta+(p+.5)*width;
  for(int j=0;j<4;++j)sum+=(long double)w[j]*(vela::fermiDiracHalfDerivative(mid+.5*width*x[j])+vela::fermiDiracHalfDerivative(mid-.5*width*x[j]))*.5*width;
 }
 return double(sum);
}
inline double relativeInverse(double z,double rel) {
 if(!(z>0) || !(rel>-1) || !std::isfinite(z) || !std::isfinite(rel))throw std::runtime_error("positive population required");
 if(rel==0)return 0; // Required identity property, not a density dead band.
 // Both endpoints must be dilute. In this range the inverse fugacity
 // correction to the log increment is below 4e-15; forming z*rel or an
 // absolute residual instead loses precision (and its tolerance) for
 // subnormal z. Restrict the repair to the extreme tail so ordinary accepted
 // states retain their existing arithmetic. This is not a population floor.
 const double logRatio=std::log1p(rel);
 if(z<1e-280 && std::log(z)+logRatio<=std::log(1e-14))return logRatio;
 double eta=vela::inverseFermiDiracHalf(z),target=z*rel;
 double d=std::log1p(rel)*z/vela::fermiDiracHalfDerivative(eta),lo=std::min(0.,2*d),hi=std::max(0.,2*d);
 for(int i=0;fermiChange(eta,lo)>target && i<60;++i)lo*=2;
 for(int i=0;fermiChange(eta,hi)<target && i<60;++i)hi*=2;
 if(fermiChange(eta,lo)>target || fermiChange(eta,hi)<target)throw std::runtime_error("unbracketed");
 for(int i=0;i<60;++i){double r=fermiChange(eta,d)-target;
  if(std::abs(r)<=2e-13*std::abs(target))return d;
  if(r>0)hi=d;else lo=d;
  double next=d-r/vela::fermiDiracHalfDerivative(eta+d);
  if(!(next>lo && next<hi) || !std::isfinite(next))next=lo+.5*(hi-lo);
  if(next==d)return d;d=next;
 }
 throw std::runtime_error("relative inverse iteration cap");
}

// Opt-in repair of the ill-conditioned population-increment inverse. Restrict
// both endpoints to the dilute region; this selects arithmetic, not a floor.
inline bool usesDiluteLogInverse(double z,double rel) {
 return z>0. && z<=1e-12 && rel>-1. && std::isfinite(rel) &&
        (long double)z*(1.L+(long double)rel)<=1e-12L;
}
inline double relativeInverseDiluteLog(double z,double rel) {
 if(!(z>0.) || !std::isfinite(z) || !(rel>-1.) || !std::isfinite(rel))
  throw std::runtime_error("positive population required");
 if(!usesDiluteLogInverse(z,rel))return relativeInverse(z,rel);
 if(rel==0.)return 0.;
 // log F(eta) = eta + log(1 + C(exp(eta))). Six fugacity terms
 // agree with the production dilute branch; evaluate C without forming F.
 auto correction=[](long double eta) {
  const long double t=std::exp(eta);long double c=0.;
  for(int k=6;k>=2;--k)
   c=(k%2 ? 1.L : -1.L)/(k*std::sqrt((long double)k))+t*c;
  return std::log1p(t*c);
 };
 const long double logz=std::log((long double)z);
 long double eta=logz;
 for(int k=0;k<4;++k)eta=logz-correction(eta);
 const long double target=std::log1p((long double)rel);
 const long double c0=correction(eta);long double d=target;
 for(int k=0;k<12;++k){
  const long double next=target+c0-correction(eta+d);
  if(next==d || std::abs(next-d)<=4*std::numeric_limits<long double>::epsilon()*std::abs(next))
   return (double)next;
  d=next;
 }
 throw std::runtime_error("dilute log inverse iteration cap");
}

// Population/DOS can underflow even when both physical inputs are positive and
// finite. Keep ordinary arithmetic unchanged; only the subnormal ratio tail
// uses logarithms. No population is clamped or excluded from any guard.
inline void validatePopulationCoordinates(double density,double dos) {
 if(!(density>0.) || !std::isfinite(density) || !(dos>0.) || !std::isfinite(dos))
  throw std::invalid_argument("positive finite population and DOS required");
}
inline double populationEta(double density,double dos) {
 validatePopulationCoordinates(density,dos);
 const double z=density/dos;
 if(z<std::numeric_limits<double>::min())return std::log(density)-std::log(dos);
 return inverseFermiDiracHalf(z);
}
inline double populationLogSlope(double eta) {
 if(!std::isfinite(eta))throw std::invalid_argument("finite population eta required");
 // In this tail the fugacity correction is below binary64 precision. Forming
 // F'/F instead can produce 0/0 or a ratio of inaccurately rounded subnormals.
 if(eta < -700.)return 1.;
 return fermiDiracHalfDerivative(eta)/fermiDiracHalf(eta);
}
// Population-primary coordinates require relative inversion accuracy. The
// historical inverse's absolute function tolerance is only a seed here.
inline double populationEtaAccurate(double density,double dos) {
 double eta=populationEta(density,dos);
 if(eta < -40.)return std::log(density)-std::log(dos);
 for(int k=0;k<12;++k) {
  const double value=fermiDiracHalf(eta),slope=populationLogSlope(eta);
  if(!(value>0.) || !(slope>0.) || !std::isfinite(slope))
   throw std::runtime_error("invalid relative population inverse");
  const long double mismatch=std::log(static_cast<long double>(value)*dos/density);
  const double next=eta-static_cast<double>(mismatch/slope);
  if(!std::isfinite(next))throw std::runtime_error("nonfinite relative population inverse");
  if(next==eta)return eta;
  eta=next;
 }
 const long double mismatch=std::log(static_cast<long double>(fermiDiracHalf(eta))*dos/density);
 if(std::abs(mismatch)>1e-12L)throw std::runtime_error("relative population inverse did not converge");
 return eta;
}
inline double populationEtaAccurate(double density, double dos, CarrierStatisticsModel model) {
    if (model != CarrierStatisticsModel::FermiDiracJoyceDixon)
        return populationEtaAccurate(density, dos);
    validatePopulationCoordinates(density, dos);
    const double ratio = density / dos;
    if (ratio < std::numeric_limits<double>::min())
        return std::log(density) - std::log(dos);
    return inverseFermiDiracHalf(ratio, model);
}
inline bool usesDiluteLogInversePopulation(double density,double dos,double rel) {
 validatePopulationCoordinates(density,dos);
 const double z=density/dos;
 if(z>=std::numeric_limits<double>::min())return usesDiluteLogInverse(z,rel);
 return rel>-1. && std::isfinite(rel) &&
        std::log(density)-std::log(dos)+std::log1p(rel)<=std::log(1e-12);
}
inline double relativeInversePopulation(double density,double dos,double rel,
                                       bool diluteLogRepair) {
 validatePopulationCoordinates(density,dos);
 if(!(rel>-1.) || !std::isfinite(rel))throw std::invalid_argument("valid population increment required");
 const double z=density/dos;
 if(z>=std::numeric_limits<double>::min())
  return diluteLogRepair?relativeInverseDiluteLog(z,rel):relativeInverse(z,rel);
 if(rel==0.)return 0.;
 const double logz=std::log(density)-std::log(dos),d=std::log1p(rel);
 const double logTarget=logz+d;
 // Both endpoints are dilute: inverse fugacity corrections are below 4e-15.
 if(logTarget<=std::log(1e-14))return d;
 // A very large positive increment can leave the tail. The endpoint is now
 // representable; use the full inverse rather than a Boltzmann extrapolation.
 return inverseFermiDiracHalf(std::exp(logTarget))-logz;
}

// For a small density tangent F'(eta)*t, return only the nonlinear part
// c of the inverse step d=t+c. Computing inverse(...)-t would lose c when
// the electrostatic and carrier steps nearly cancel. Taylor's integral
// remainder gives c=-d^2 integral_0^1 (1-s) F''(eta+s*d)/F'(eta) ds.
// This is an opt-in update arithmetic helper, not a density floor.
inline double inverseTangentRemainder(double eta,double t) {
 if(!std::isfinite(eta)||!std::isfinite(t)||std::abs(t)>1e-4)
  throw std::invalid_argument("small Fermi tangent required");
 if(t==0)return 0.;
 if(eta < -640.) {
  // log1p(t)-t, evaluated without subtracting two O(t) numbers.
  return t*t*(-.5+t*(1./3.+t*(-.25+t*(.2-t/6.))));
 }
 const double first=fermiDiracHalfDerivative(eta);
 if(!(first>0.)||!std::isfinite(first))throw std::runtime_error("invalid Fermi tangent slope");
 static constexpr double x[]={.1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
 static constexpr double w[]={.3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
 double c=0.;
 for(int k=0;k<12;++k){
  const double d=t+c;long double integral=0.;
  for(int j=0;j<4;++j)for(int side:{-1,1}){
   const double s=.5*(1.+side*x[j]);
   integral+=(long double)(.5*w[j])*(1.-s)*(fermiDiracHalfSecondDerivative(eta+s*d)/first);
  }
  const double next=(double)(-(long double)d*d*integral);
  if(next==c || std::abs(next-c)<=4*std::numeric_limits<double>::epsilon()*std::abs(next))return next;
  c=next;
 }
 throw std::runtime_error("small Fermi remainder did not converge");
}
inline double tangentQuasiFermiStep(double eta,double sign,double vt,
                                    double psiDirection,double qfDirection,
                                    double appliedPsi) {
 if((sign!=-1. && sign!=1.) || !(vt>0.) || !std::isfinite(vt))
  throw std::invalid_argument("invalid carrier sign or thermal scale");
 const double t=sign*(qfDirection-psiDirection);
 const double c=inverseTangentRemainder(eta,t);
 // Preserve the small QF component before adding damping/writeback effects.
 const double requestedPsi=vt*psiDirection;
 return (double)((long double)vt*qfDirection+
     ((long double)appliedPsi-requestedPsi)+(long double)sign*vt*c);
}

struct LinearQuality {
 bool accepted=true;
 double mixedRatio=0, componentBackward=0, constraintV=0;
 std::vector<double> blockNormwise;
};
// Same mixed backward guard as scripts/slot_linear_quality.py, accumulated
// directly from the assembled matrix in its column-major export order.
inline LinearQuality linearQuality(const SparseMatrixd& A,const VectorXd& b,
                                  const VectorXd& x,double vt,int width) {
 if(A.rows()!=A.cols() || A.rows()!=b.size() || x.size()!=b.size() ||
    width<=0 || b.size()%width || !(vt>0) || !std::isfinite(vt))
  throw std::invalid_argument("invalid population linear audit dimensions/scale");
 LinearQuality q;
 if(!b.allFinite() || !x.allFinite()) {q.accepted=false;return q;}
 const int n=b.size();std::vector<long double> ax(n),scale(n),sum(n),diag(n);
 std::vector<int> count(n);
 for(int i=0;i<n;++i)scale[i]=std::abs((long double)b[i]);
 for(int j=0;j<A.outerSize();++j)for(SparseMatrixd::InnerIterator it(A,j);it;++it){
  long double v=it.value();int i=it.row();
  if(!std::isfinite(v)){q.accepted=false;return q;}
  ax[i]+=v*(long double)x[it.col()];scale[i]+=std::abs(v)*std::abs((long double)x[it.col()]);sum[i]+=std::abs(v);
  if(v!=0)++count[i];if(i==it.col())diag[i]=v;
 }
 auto ratio=[](long double a,long double b){return b>0?a/b:(a==0?0.L:std::numeric_limits<long double>::infinity());};
 for(int k=0;k<n;k+=width){long double mr=0,ms=0;
  for(int i=k;i<k+width;++i){auto r=std::abs(ax[i]-(long double)b[i]);
   const auto m=ratio(r,1e-6L*scale[i]+sum[i]*1e-14L/vt);
   q.mixedRatio=std::max(q.mixedRatio,(double)m);q.componentBackward=std::max(q.componentBackward,(double)ratio(r,scale[i]));
   mr=std::max(mr,r);ms=std::max(ms,scale[i]);
   if(count[i]==1 && diag[i]!=0)q.constraintV=std::max(q.constraintV,(double)(std::abs((long double)x[i]-(long double)b[i]/diag[i])*vt));
  }
  q.blockNormwise.push_back((double)ratio(mr,ms));
 }
 q.accepted=q.mixedRatio<=1 && q.constraintV<=1e-14 && *std::max_element(q.blockNormwise.begin(),q.blockNormwise.end())<=1e-12;
 return q;
}
struct RowRefinement {
 LinearQuality before, after;
 int sweeps=0, correctedRows=0;
};
// Reuse the caller's factorized original system for correction RHSs. No
// coefficient or tolerance changes; unsuccessful trials never replace x.
template<class Solve>
inline RowRefinement refineWithFactors(const SparseMatrixd& A,const VectorXd& b,
    VectorXd& x,double vt,int width,int maxSteps,Solve&& solve) {
 if(maxSteps<0 || maxSteps>3)throw std::invalid_argument("factor refinement requires 0..3 steps");
 RowRefinement out;out.before=out.after=linearQuality(A,b,x,vt,width);
 if(out.before.accepted || maxSteps==0 || !b.allFinite() || !x.allFinite())return out;
 VectorXd trial=x;
 for(int k=0;k<maxSteps;++k){
  std::vector<long double> defect(b.size());for(int i=0;i<b.size();++i)defect[i]=b[i];
  for(int j=0;j<A.outerSize();++j)for(SparseMatrixd::InnerIterator it(A,j);it;++it)defect[it.row()]-=(long double)it.value()*trial[j];
  VectorXd rhs(b.size());for(int i=0;i<b.size();++i)rhs[i]=(double)defect[i];
  if(!rhs.allFinite())return out;
  VectorXd dx=solve(rhs);if(dx.size()!=x.size() || !dx.allFinite())return out;
  for(int i=0;i<x.size();++i)trial[i]=(double)((long double)trial[i]+dx[i]);
  ++out.sweeps;const auto q=linearQuality(A,b,trial,vt,width);
  if(q.accepted){x=trial;out.after=q;return out;}
 }
 return out;
}
// Bounded diagnostic coordinate relaxation, not a change to A, b or the guard.
// Each sweep selects failing mixed rows from the current defect. Commit only
// if ALL original guards pass; otherwise return the original direction intact.
inline RowRefinement refineMixedRows(const SparseMatrixd& A,const VectorXd& b,
                                    VectorXd& x,double vt,int width,int maxSweeps) {
 if(maxSweeps<0 || maxSweeps>4)throw std::invalid_argument("row refinement requires 0..4 sweeps");
 RowRefinement out;out.before=out.after=linearQuality(A,b,x,vt,width);
 if(out.before.accepted || maxSweeps==0)return out;
 if(!x.allFinite() || !b.allFinite())return out;
 VectorXd trial=x;const int n=b.size();
 for(int k=0;k<maxSweeps;++k){
  std::vector<long double> ax(n),scale(n),sum(n),diag(n);std::vector<int> counts(n);
  for(int i=0;i<n;++i)scale[i]=std::abs((long double)b[i]);
  for(int j=0;j<A.outerSize();++j)for(SparseMatrixd::InnerIterator it(A,j);it;++it){
   if(!std::isfinite(it.value()))return out;
   const int i=it.row();const long double v=it.value();
   ax[i]+=v*(long double)trial[j];scale[i]+=std::abs(v)*std::abs((long double)trial[j]);sum[i]+=std::abs(v);
   if(v!=0)++counts[i];if(i==j)diag[i]=v;
  }
  int changed=0;
  for(int i=0;i<n;++i){
   const auto defect=ax[i]-(long double)b[i];
   // Never repair boundary singleton rows by silently changing constraints.
   if(counts[i]<=1 || diag[i]==0 || std::abs(defect)<=1e-6L*scale[i]+sum[i]*1e-14L/vt)continue;
   const double next=(double)((long double)trial[i]-defect/diag[i]);
   if(!std::isfinite(next))return out;
   if(next!=trial[i]){trial[i]=next;++changed;}
  }
  ++out.sweeps;out.correctedRows+=changed;
  const auto q=linearQuality(A,b,trial,vt,width);
  if(q.accepted){x=trial;out.after=q;return out;}
  if(!changed)break;
 }
 return out;
}
} // namespace vela::slot_population
