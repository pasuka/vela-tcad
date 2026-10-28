#pragma once

#include "vela/equation/AssemblerUtils.h"
// Opt-in BVMethods Tri3 element contract. Physical fields are volts and
// caller-internal populations; local jets use physical-voltage columns.
// SG uses stable log-density/expm1 QF differences. No native/cache input.

#include "vela/physics/experimental/BvMethodsTri3Avalanche.h"
#include "vela/physics/ElementLombardi.h"
#include "vela/physics/experimental/BvMethodsTri3Transport.h"
#include "vela/physics/experimental/JoyceDixonStatistics.h"

namespace vela::detail {
struct BvMethodsQfState {
 VectorXd psi,phin,phip,n,p;
 bool populationCoordinates=false;
};
struct BvMethodsTri3Context {
 const DeviceMesh& mesh;
 const MaterialDatabase& matdb;
 const DopingModel& doping;
 struct Controls {
  MobilityModelConfig mobility;
  ImpactIonizationModelConfig impactIonization;
  UnitScalingConfig inputScaling;
  CarrierStatisticsModel statistics;
 } newton;
};
using experimental::AvJet;
struct BvMethodsTri3Evaluation {
 VectorXd absoluteSource;
 VectorXd electronGeneration,holeGeneration;
 VectorXd source;
 VectorXd allGeneration;
 SparseMatrixd jacobian;
 int zeroActiveElectron=0,zeroActiveHole=0;
};

class BvMethodsTri3Operator {
 const BvMethodsTri3Context& problem_;
 std::vector<Material> materials_;
 MobilityModelConfig cfg_,geo_;
 std::shared_ptr<MobilityModel> bulk_;
 std::vector<bool> electronActive_,holeActive_;
 std::vector<double> nc_,nv_;
 double potentialScale_,factor_,vt_;
 struct CellCache {
  Index id;
  std::array<Index,3> nodes;
  experimental::AvTri3Geometry geometry;
  std::array<Point2,3> shape;
  Point2 distanceGradient;
  std::array<double,3> bulkN,bulkP;
 };
 std::vector<CellCache> cells_;
 experimental::AvCoefficient electron_,hole_;
public:
 BvMethodsTri3Operator(const BvMethodsTri3Context& problem,const BvMethodsQfState& initial,
        const std::vector<double>& nc,const std::vector<double>& nv,
        double potentialScale,double physicalLineFactor,double thermalVoltage)
 :problem_(problem),materials_(detail::buildCellMaterials(problem.mesh,problem.matdb,300.)),
  cfg_(problem.newton.mobility),geo_(cfg_),potentialScale_(potentialScale),
  factor_(100./physicalLineFactor),vt_(thermalVoltage)
 {
  const auto& mesh=problem.mesh;const int n=mesh.numNodes();const auto units=problem.newton.inputScaling.unitSystem();
  electronActive_.resize(n);holeActive_.resize(n);nc_.resize(n);nv_.resize(n);
  for(int i=0;i<n;++i){nc_[i]=units.internalConcentrationToM3(nc[i])/1e6;nv_[i]=units.internalConcentrationToM3(nv[i])/1e6;}
  geo_.model="phumob_field_lombardi";geo_.edgeAveraging="element_box_phumob";
  geo_.surface.discretization="element_distance_gradient";geo_.highFieldGradientDiscretization="element_vertex_partial_layer";
  geo_.surface.cellNormalX.clear();geo_.surface.cellNormalY.clear();geo_.surface.cellDistances.clear();geo_.surface.nodeDistances.clear();geo_.highFieldCells.clear();
  auto ec=detail::buildEdgeCellMap(mesh);detail::updateSurfaceMobilityCellGeometry(geo_,mesh,ec,initial.psi,cfg_.surface.coordinateFieldFactor,&materials_);
  auto bc=cfg_;bc.model="masetti";bulk_=makeMobilityModel(bc);
  const auto& ion=problem.newton.impactIonization;
  auto invcm=[&](double a){return units.internalInverseLengthToMInv(a)/100.;};
  auto fieldcm=[&](double a){return units.internalElectricFieldToVPerM(a)/100.;};
  electron_={invcm(ion.electronALow),fieldcm(ion.electronBLow),invcm(ion.electronAHigh),fieldcm(ion.electronBHigh),fieldcm(ion.switchField)};
  hole_={invcm(ion.holeALow),fieldcm(ion.holeBLow),invcm(ion.holeAHigh),fieldcm(ion.holeBHigh),fieldcm(ion.switchField)};
  for(Index cid=0;cid<mesh.numCells();++cid){if(materials_[cid].ni<=0)continue;
   const auto& c=mesh.getCell(cid);CellCache a;a.id=cid;
   for(auto i:c.node_ids){electronActive_[i]=true;holeActive_[i]=true;}
   for(int k=0;k<3;++k){a.nodes[k]=c.node_ids[k];const auto& p=mesh.getNode(a.nodes[k]);a.geometry.points_cm[k]={units.internalLengthToMeters(p.x)*100.,units.internalLengthToMeters(p.y)*100.};}
   const auto& p0=mesh.getNode(a.nodes[0]);const auto& p1=mesh.getNode(a.nodes[1]);const auto& p2=mesh.getNode(a.nodes[2]);
   const double x1=p1.x-p0.x,y1=p1.y-p0.y,x2=p2.x-p0.x,y2=p2.y-p0.y,det=x1*y2-x2*y1;
   a.shape[1]={y2/det,-x2/det};a.shape[2]={-y1/det,x1/det};a.shape[0]=-a.shape[1]-a.shape[2];
   const auto weights=detail::tri3ElementVertexBoxMeasures(mesh,c);
   double total=0;for(double w:weights)total+=w;
   int rep=0;for(int k=1;k<3;++k)if(weights[k]>weights[rep]+128*std::numeric_limits<double>::epsilon()*total)rep=k;
   a.geometry.representative=rep;
   for(int k=0;k<3;++k){auto i=a.nodes[k];a.geometry.boxAreas_cm2[k]=weights[k]*units.areaM2PerInternal()*1e4;
    auto imp=problem.doping.totalImpurity(i);a.bulkN[k]=bulk_->electronMobility(materials_[cid],imp,0,0);a.bulkP[k]=bulk_->holeMobility(materials_[cid],imp,0,0);
   }
   a.distanceGradient=a.shape[1]*(geo_.surface.nodeDistances[a.nodes[1]]-geo_.surface.nodeDistances[a.nodes[0]])+a.shape[2]*(geo_.surface.nodeDistances[a.nodes[2]]-geo_.surface.nodeDistances[a.nodes[0]]);
   Eigen::Matrix2d G,S;G.col(0)=a.shape[1];G.col(1)=a.shape[2];const auto& z=geo_.highFieldCells[cid].stencil;
   S<<z[0],z[2],z[1],z[3];auto P=(S/cfg_.surface.coordinateFieldFactor*G.inverse()).eval();
   a.geometry.currentProjection={P(0,0),P(0,1),P(1,0),P(1,1)};cells_.push_back(a);
  }
 }

 BvMethodsTri3Evaluation evaluate(const BvMethodsQfState& state,bool derivatives) const {
  const auto units=problem_.newton.inputScaling.unitSystem();const int n=problem_.mesh.numNodes();BvMethodsTri3Evaluation out;
  out.electronGeneration=VectorXd::Zero(n);out.holeGeneration=VectorXd::Zero(n);out.absoluteSource=VectorXd::Zero(3*n);out.source=VectorXd::Zero(3*n);out.allGeneration=VectorXd::Zero(n);
  std::vector<Eigen::Triplet<double>> entries;if(derivatives)entries.reserve(cells_.size()*54);
  for(const auto&a:cells_){std::array<AvJet,3> psi,fn,fp,ns,ps;
   for(int k=0;k<3;++k){auto i=a.nodes[k];psi[k]=AvJet::variable(state.psi[i],k);fn[k]=AvJet::variable(state.phin[i],3+k);fp[k]=AvJet::variable(state.phip[i],6+k);
    const double nd=units.internalConcentrationToM3(state.n[i])/1e6,pd=units.internalConcentrationToM3(state.p[i])/1e6;
    ns[k]=AvJet(nd);ps[k]=AvJet(pd);
    const double dn=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(nd,nc_[i])/vt_;
    const double dp=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(pd,nv_[i])/vt_;
    ns[k].derivative[k]=dn;ns[k].derivative[3+k]=-dn;ps[k].derivative[k]=-dp;ps[k].derivative[6+k]=dp;
    if(state.populationCoordinates){
     // Jets retain the historical final V0 scatter; concentration columns
     // therefore carry their physical internal-unit factor divided by V0.
     const double dc=units.internalConcentrationToM3(1.)/1e6/potentialScale_;
     ns[k].derivative.fill(0.);ps[k].derivative.fill(0.);
     ns[k].derivative[3+k]=dc;ps[k].derivative[6+k]=dc;
     fn[k].derivative.fill(0.);fp[k].derivative.fill(0.);
     fn[k].derivative[k]=1.;fp[k].derivative[k]=1.;
     fn[k].derivative[3+k]=-dc/dn;fp[k].derivative[6+k]=dc/dp;
    }
   }
   const auto& s=geo_.highFieldCells[a.id].stencil;
   auto field=[&](const std::array<AvJet,3>& z){return experimental::avNorm({(z[1]-z[0])*AvJet(s[0])+(z[2]-z[0])*AvJet(s[2]),(z[1]-z[0])*AvJet(s[1])+(z[2]-z[0])*AvJet(s[3])});};
   const auto en=field(geo_.highFieldCells[a.id].contact?psi:fn),hp=field(geo_.highFieldCells[a.id].contact?psi:fp);
   const auto rawNormal=(psi[1]-psi[0])*AvJet(a.shape[1].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor)+(psi[2]-psi[0])*AvJet(a.shape[2].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor);
   const auto normal=detail::dualAbs(rawNormal);int rep=a.geometry.representative;auto node=a.nodes[rep];
   auto mobility=[&](bool e,const AvJet& drive){
    const double low=e?a.bulkN[rep]:a.bulkP[rep];const auto& lp=e?cfg_.electronLombardi:cfg_.holeLombardi;const auto& hf=e?cfg_.electronField:cfg_.holeField;
    const auto surface=evaluateElementLombardi(low,normal.value,geo_.surface.nodeDistances[node],problem_.doping.totalImpurity(node),300.,lp,cfg_);
    const auto high=evaluateElementCanali(surface.mobility,drive.value,hf,100./cfg_.internalFieldToVPerM);
    const double scale=units.internalMobilityToM2PerVS(1.)*1e4;AvJet mu(high.mobility*scale);
    for(int k=0;k<9;++k)mu.derivative[k]=(high.lowDerivative*surface.fieldDerivative*normal.derivative[k]+high.fieldDerivative*drive.derivative[k])*scale;
    return mu;
   };
   const auto mn=mobility(true,en),mp=mobility(false,hp);
   const auto result=experimental::evaluateBvMethodsTri3Avalanche(a.geometry,ns,ps,fn,fp,mn,mp,en,hp,vt_,electron_,hole_);
   if(result.electronZeroNorm&&experimental::avCoefficient(en,electron_).value>0)++out.zeroActiveElectron;
   if(result.holeZeroNorm&&experimental::avCoefficient(hp,hole_).value>0)++out.zeroActiveHole;
   const auto genN=experimental::avCoefficient(en,electron_)*mn*experimental::avNorm(result.electronParticleCurrentWithoutMobility);
   const auto genP=experimental::avCoefficient(hp,hole_)*mp*experimental::avNorm(result.holeParticleCurrentWithoutMobility);
   for(int k=0;k<3;++k){auto node=a.nodes[k];out.electronGeneration[node]+=a.geometry.boxAreas_cm2[k]*genN.value*factor_;out.holeGeneration[node]+=a.geometry.boxAreas_cm2[k]*genP.value*factor_;const auto& g=result.nodeGeneration_per_cm_s[k];if(!std::isfinite(g.value))throw std::runtime_error("Nonfinite candidate AV source");out.allGeneration[node]+=g.value;
    for(int block=1;block<=2;++block){if(!(block==1?electronActive_[node]:holeActive_[node]))continue;const int row=block*n+node;out.source[row]-=g.value*factor_;
     if(derivatives)for(int j=0;j<9;++j){const double d=-g.derivative[j]*factor_*potentialScale_;if(d!=0)entries.emplace_back(row,(j/3)*n+a.nodes[j%3],d);}
    }
   }
  }
  out.jacobian.resize(3*n,3*n);if(derivatives)out.jacobian.setFromTriplets(entries.begin(),entries.end());return out;
 }
 // Diagnostic approximate linearization only: retain the live mobility value,
 // but omit its tangent. Production calls retain the full derivative default.
 BvMethodsTri3Evaluation transport(const BvMethodsQfState& state,bool derivatives,bool generalized,bool coefficientAveraging=false,bool freezeMobilityDerivatives=false) const {
  const auto units=problem_.newton.inputScaling.unitSystem();const int n=problem_.mesh.numNodes();BvMethodsTri3Evaluation out;
  out.electronGeneration=VectorXd::Zero(n);out.holeGeneration=VectorXd::Zero(n);out.absoluteSource=VectorXd::Zero(3*n);out.source=VectorXd::Zero(3*n);out.allGeneration=VectorXd::Zero(n);
  std::vector<Eigen::Triplet<double>> entries;if(derivatives)entries.reserve(cells_.size()*54);
  for(const auto&a:cells_){std::array<AvJet,3> psi,fn,fp,ns,ps;
   for(int k=0;k<3;++k){auto i=a.nodes[k];psi[k]=AvJet::variable(state.psi[i],k);fn[k]=AvJet::variable(state.phin[i],3+k);fp[k]=AvJet::variable(state.phip[i],6+k);
    const double nd=units.internalConcentrationToM3(state.n[i])/1e6,pd=units.internalConcentrationToM3(state.p[i])/1e6;
    ns[k]=AvJet(nd);ps[k]=AvJet(pd);
    const double dn=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(nd,nc_[i])/vt_;
    const double dp=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(pd,nv_[i])/vt_;
    ns[k].derivative[k]=dn;ns[k].derivative[3+k]=-dn;ps[k].derivative[k]=-dp;ps[k].derivative[6+k]=dp;
    if(state.populationCoordinates){
     // Jets retain the historical final V0 scatter; concentration columns
     // therefore carry their physical internal-unit factor divided by V0.
     const double dc=units.internalConcentrationToM3(1.)/1e6/potentialScale_;
     ns[k].derivative.fill(0.);ps[k].derivative.fill(0.);
     ns[k].derivative[3+k]=dc;ps[k].derivative[6+k]=dc;
     fn[k].derivative.fill(0.);fp[k].derivative.fill(0.);
     fn[k].derivative[k]=1.;fp[k].derivative[k]=1.;
     fn[k].derivative[3+k]=-dc/dn;fp[k].derivative[6+k]=dc/dp;
    }
   }
   const auto& s=geo_.highFieldCells[a.id].stencil;
   auto field=[&](const std::array<AvJet,3>& z){return experimental::avNorm({(z[1]-z[0])*AvJet(s[0])+(z[2]-z[0])*AvJet(s[2]),(z[1]-z[0])*AvJet(s[1])+(z[2]-z[0])*AvJet(s[3])});};
   const auto en=field(geo_.highFieldCells[a.id].contact?psi:fn),hp=field(geo_.highFieldCells[a.id].contact?psi:fp);
   const auto rawNormal=(psi[1]-psi[0])*AvJet(a.shape[1].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor)+(psi[2]-psi[0])*AvJet(a.shape[2].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor);
   const auto normal=detail::dualAbs(rawNormal);int rep=a.geometry.representative;auto node=a.nodes[rep];
   auto mobility=[&](bool e,const AvJet& drive,int local){auto node=a.nodes[local];
    const double low=e?a.bulkN[local]:a.bulkP[local];const auto& lp=e?cfg_.electronLombardi:cfg_.holeLombardi;const auto& hf=e?cfg_.electronField:cfg_.holeField;
    const auto surface=evaluateElementLombardi(low,normal.value,geo_.surface.nodeDistances[node],problem_.doping.totalImpurity(node),300.,lp,cfg_);
    const auto high=evaluateElementCanali(surface.mobility,drive.value,hf,100./cfg_.internalFieldToVPerM);
    const double scale=units.internalMobilityToM2PerVS(1.)*1e4;AvJet mu(high.mobility*scale);
    for(int k=0;k<9;++k)mu.derivative[k]=(high.lowDerivative*surface.fieldDerivative*normal.derivative[k]+high.fieldDerivative*drive.derivative[k])*scale;
    return mu;
   };

   std::array<AvJet,3> mn,mp,etaN,etaP,gn{AvJet(1),AvJet(1),AvJet(1)},gp=gn;
   for(int k=0;k<3;++k){
    mn[k]=mobility(true,en,k);mp[k]=mobility(false,hp,k);
    if(freezeMobilityDerivatives){mn[k].derivative.fill(0.);mp[k].derivative.fill(0.);}
    const auto i=a.nodes[k];etaN[k]=AvJet(experimental::JoyceDixonStatistics::inverseDensity(ns[k].value,nc_[i]));etaP[k]=AvJet(experimental::JoyceDixonStatistics::inverseDensity(ps[k].value,nv_[i]));
    etaN[k].derivative[k]=1/vt_;etaN[k].derivative[3+k]=-1/vt_;
    etaP[k].derivative[k]=-1/vt_;etaP[k].derivative[6+k]=1/vt_;
    if(state.populationCoordinates){
     etaN[k].derivative.fill(0.);etaP[k].derivative.fill(0.);
     etaN[k].derivative[3+k]=-fn[k].derivative[3+k]/vt_;
     etaP[k].derivative[6+k]=fp[k].derivative[6+k]/vt_;
    }
   }
   auto einstein=[&](const std::array<AvJet,3>& population,const std::array<AvJet,3>& eta,int k){
    const int j=(k+1)%3;auto log=experimental::avLogRatio(population[j],population[k]);
    if(std::abs(log.value)>1e-10)return (eta[j]-eta[k])/log;
    const auto mid=(eta[j]+eta[k])/AvJet(2);const auto jet=experimental::JoyceDixonStatistics::evaluate(mid.value);
    AvJet value(std::max(1.,jet.value/jet.first));
    if(value.value>1.)for(int z=0;z<9;++z)value.derivative[z]=(1-jet.value*jet.second/(jet.first*jet.first))*mid.derivative[z];
    return value;
   };
   if(generalized)for(int k=0;k<3;++k){gn[k]=einstein(ns,etaN,k);gp[k]=einstein(ps,etaP,k);}
   const auto& coefficients=problem_.mesh.poissonCellEdgeCoefficients(a.id);
   std::array<double,3> edge{coefficients[0],coefficients[1],coefficients[2]};
   auto averagingAreas=a.geometry.boxAreas_cm2;
   if(coefficientAveraging){auto weights=detail::cellBoxNodeMeasures(problem_.mesh,a.id);for(int k=0;k<3;++k)averagingAreas[k]=weights[k]*units.areaM2PerInternal()*1e4;}
   auto tn=experimental::evaluateBvMethodsTri3Transport(averagingAreas,edge,ns,fn,mn,vt_,true,gn);
   auto tp=experimental::evaluateBvMethodsTri3Transport(averagingAreas,edge,ps,fp,mp,vt_,false,gp);
   for(int k=0;k<3;++k)for(int block=1;block<=2;++block){
    const int row=block*n+a.nodes[k];const auto& f=block==1?tn.nodeFlux_per_cm_s[k]:tp.nodeFlux_per_cm_s[k];out.source[row]+=f.value*factor_;out.absoluteSource[row]+=(std::abs((block==1?tn.edgeFlux_per_cm_s[k]:tp.edgeFlux_per_cm_s[k]).value)+std::abs((block==1?tn.edgeFlux_per_cm_s[(k+2)%3]:tp.edgeFlux_per_cm_s[(k+2)%3]).value))*factor_;
    if(derivatives)for(int j=0;j<9;++j){double value=f.derivative[j]*factor_*potentialScale_;if(value!=0)entries.emplace_back(row,(j/3)*n+a.nodes[j%3],value);}
   }
  }
  out.jacobian.resize(3*n,3*n);if(derivatives)out.jacobian.setFromTriplets(entries.begin(),entries.end());return out;
 }
};
} // namespace vela::detail
