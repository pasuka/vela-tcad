#pragma once
#include <optional>
#include <bit>
#include <cstdint>
#include "vela/core/ArithmeticPrecision.h"
#include "vela/physics/experimental/BvMethodsAvValues.h"

#include "vela/equation/AssemblerUtils.h"
#include <limits>
// Opt-in BVMethods Tri3 element contract. Physical fields are volts and
// caller-internal populations; local jets use physical-voltage columns.
// SG uses stable log-density/expm1 QF differences. No native/cache input.

#include "vela/physics/experimental/BvMethodsTri3Avalanche.h"
#include "vela/physics/PreparedElementLombardi.h"
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
  BvMethodsNumericalOptions numerics{};
 } newton;
};
using experimental::AvJet;

struct CanaliEvaluationOptions {
 bool avalancheValueOnly=false;
 bool sharedPower=false;
 bool audit=false;

 ElementCanaliValue full(Real low,Real field,const FieldMobilityParameters& p,Real cutoff) const {
  if(!sharedPower)return evaluateElementCanali(low,field,p,cutoff);
  const auto value=evaluateElementCanaliSharedPower(low,field,p,cutoff);
  if(audit) {
   const auto ref=evaluateElementCanali(low,field,p,cutoff);
   if(!same(value.mobility,ref.mobility) || !same(value.lowDerivative,ref.lowDerivative) ||
      !same(value.fieldDerivative,ref.fieldDerivative))
    throw std::runtime_error("Canali shared-power bit mismatch");
  }
  return value;
 }

 static bool same(Real a,Real b) {
  return std::bit_cast<std::uint64_t>(a)==std::bit_cast<std::uint64_t>(b);
 }
};

using BvMethodsRepresentativeMobilities=std::vector<std::array<double,2>>;

struct BvMethodsTri3Evaluation {
 VectorXd absoluteSource;
 VectorXd electronGeneration,holeGeneration;
 VectorXd source;
 VectorXd allGeneration;
 SparseMatrixd jacobian;
 int zeroActiveElectron=0,zeroActiveHole=0;
 int zeroNormElectron=0,zeroNormHole=0;
 int scalarCells=0,fallbackCells=0,reusedMobilityCells=0,auditedMobilityPairs=0;
};

class BvMethodsTri3Operator {
 const BvMethodsTri3Context& problem_;
 std::vector<Material> materials_;
 MobilityModelConfig cfg_,geo_;
 std::shared_ptr<MobilityModel> bulk_;
 std::vector<bool> electronActive_,holeActive_;
 std::vector<double> nc_,nv_,ni_;
 std::vector<PreparedElementLombardi> lombardiN_,lombardiP_;
 bool directPopulationSg_=false;
 double potentialScale_,factor_,vt_;
 struct CellCache {
  Index id;
  std::array<Index,3> nodes;
  experimental::AvTri3Geometry geometry;
  std::array<double,3> transportAreas_cm2;
  std::array<Point2,3> shape;
  Point2 distanceGradient;
  std::array<double,3> bulkN,bulkP;
 };
 std::vector<CellCache> cells_;
 experimental::AvCoefficient electron_,hole_;
 experimental::JoyceDixonJet statisticsJet(double eta) const {
  if (problem_.newton.statistics == CarrierStatisticsModel::FermiDiracJoyceDixon)
   return experimental::JoyceDixonStatistics::evaluate(eta);
  return {fermiDiracHalf(eta), fermiDiracHalfDerivative(eta), fermiDiracHalfSecondDerivative(eta)};
 }
 struct NodeStatistics {double n=0.,p=0.,dn=0.,dp=0.,etaN=0.,etaP=0.,tangentN=1.,tangentP=1.;};
 // Values depend on the node state, not on the cell's local AD slot. Keep
 // this scratch local to one assembly; no reuse across states or owners.
 std::vector<NodeStatistics> nodeStatistics(const BvMethodsQfState& state,bool inverse,const VectorXd* columns=nullptr,bool slopes=true) const {
  const auto units=problem_.newton.inputScaling.unitSystem();
  std::vector<NodeStatistics> result(problem_.mesh.numNodes());
  for(std::size_t i=0;i<result.size();++i)if(electronActive_[i]) {
   auto& a=result[i];a.n=units.internalConcentrationToCm3(state.n[i]);a.p=units.internalConcentrationToCm3(state.p[i]);
   if(problem_.newton.statistics == CarrierStatisticsModel::FermiDirac && (slopes || inverse)) {
    const auto n = integralFermiPopulation(a.n,nc_[i]);
    const auto p = integralFermiPopulation(a.p,nv_[i]);
    a.dn=n.densitySlope/vt_;a.dp=p.densitySlope/vt_;
    a.etaN=n.eta;a.etaP=p.eta;
   } else {
    if(slopes) {
     a.dn=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(a.n,nc_[i])/vt_;
     a.dp=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(a.p,nv_[i])/vt_;
    }
    if(inverse){a.etaN=experimental::JoyceDixonStatistics::inverseDensity(a.n,nc_[i]);a.etaP=experimental::JoyceDixonStatistics::inverseDensity(a.p,nv_[i]);}
   }
   // Evaluate extremely depleted columns in a relative local direction,
   // then undo that direction scale only after flux/source cancellation.
   // State values and global population coordinates remain unchanged.
   if(state.populationCoordinates) {
    if(a.n<std::sqrt(std::numeric_limits<double>::min()))a.tangentN=a.n;
    if(a.p<std::sqrt(std::numeric_limits<double>::min()))a.tangentP=a.p;
    if(columns){a.tangentN=(*columns)[result.size()+i];a.tangentP=(*columns)[2*result.size()+i];}
   }
  }
  return result;
 }
public:
 BvMethodsTri3Operator(const BvMethodsTri3Context& problem,const BvMethodsQfState& initial,
        const std::vector<double>& nc,const std::vector<double>& nv,
        double potentialScale,double physicalLineFactor,double thermalVoltage,
        const std::vector<double>& ni={},bool directPopulationSg=false)
 :problem_(problem),materials_(detail::buildCellMaterials(problem.mesh,problem.matdb,300.)),
  cfg_(problem.newton.mobility),geo_(cfg_),potentialScale_(potentialScale),
  factor_(100./physicalLineFactor),vt_(thermalVoltage)
 {
  const auto& mesh=problem.mesh;const int n=mesh.numNodes();const auto units=problem.newton.inputScaling.unitSystem();
  electronActive_.resize(n);holeActive_.resize(n);nc_.resize(n);nv_.resize(n);
  directPopulationSg_=directPopulationSg;
  if(directPopulationSg_ && ni.size()!=static_cast<std::size_t>(n))throw std::invalid_argument("Population SG requires intrinsic densities");
  ni_.resize(n);for(int i=0;i<n;++i){nc_[i]=units.internalConcentrationToCm3(nc[i]);nv_[i]=units.internalConcentrationToCm3(nv[i]);if(directPopulationSg_)ni_[i]=units.internalConcentrationToCm3(ni[i]);}
  geo_.model="phumob_field_lombardi";geo_.edgeAveraging="element_box_phumob";
  geo_.surface.discretization="element_distance_gradient";geo_.highFieldGradientDiscretization="element_vertex_partial_layer";
  if(cfg_.model=="constant"){
   // Build the shared stencils even though constant transport does not use the field.
   geo_.highFieldDrivingForce="quasi_fermi_gradient";
   geo_.jacobianFieldDerivatives=true;
  }
  geo_.surface.cellNormalX.clear();geo_.surface.cellNormalY.clear();geo_.surface.cellDistances.clear();geo_.surface.nodeDistances.clear();geo_.highFieldCells.clear();
  const auto& ec = detail::buildEdgeCellMap(mesh);detail::updateSurfaceMobilityCellGeometry(geo_,mesh,ec,initial.psi,cfg_.surface.coordinateFieldFactor,&materials_);
  auto bc=cfg_;if(bc.model!="constant")bc.model="masetti";bulk_=makeMobilityModel(bc);
  const auto& ion=problem.newton.impactIonization;
  auto invcm=[&](double a){return units.internalInverseLengthToMInv(a)/100.;};
  auto fieldcm=[&](double a){return units.internalElectricFieldToVPerM(a)/100.;};
  electron_={invcm(ion.electronALow),fieldcm(ion.electronBLow),invcm(ion.electronAHigh),fieldcm(ion.electronBHigh),fieldcm(ion.switchField)};
  hole_={invcm(ion.holeALow),fieldcm(ion.holeBLow),invcm(ion.holeAHigh),fieldcm(ion.holeBHigh),fieldcm(ion.switchField)};
  electron_.smoothSwitch=hole_.smoothSwitch=problem.newton.numerics.avalancheSwitch==BvMethodsAvalancheSwitch::Smooth;
  for(Index cid=0;cid<mesh.numCells();++cid){if(materials_[cid].ni<=0)continue;
   const auto& c=mesh.getCell(cid);CellCache a;a.id=cid;
   for(auto i:c.node_ids){electronActive_[i]=true;holeActive_[i]=true;}
   for(int k=0;k<3;++k){a.nodes[k]=c.node_ids[k];a.geometry.nodeOrder[k]=a.nodes[k];const auto& p=mesh.getNode(a.nodes[k]);a.geometry.points_cm[k]={units.internalLengthToMeters(p.x)*100.,units.internalLengthToMeters(p.y)*100.};}
   const auto& p0=mesh.getNode(a.nodes[0]);const auto& p1=mesh.getNode(a.nodes[1]);const auto& p2=mesh.getNode(a.nodes[2]);
   const double x1=p1.x-p0.x,y1=p1.y-p0.y,x2=p2.x-p0.x,y2=p2.y-p0.y,det=x1*y2-x2*y1;
   a.shape[1]={y2/det,-x2/det};a.shape[2]={-y1/det,x1/det};a.shape[0]=-a.shape[1]-a.shape[2];
   const auto weights=detail::tri3ElementVertexBoxMeasures(mesh,c);
   const auto avalancheWeights = problem.newton.numerics.avalancheWeights == BvMethodsAvalancheWeights::CoefficientBox
       ? conservativeCoefficientBoxMeasures(mesh,cid) : weights;
   double total=0;for(double w:weights)total+=w;
   int rep=0;for(int k=1;k<3;++k)if(weights[k]>weights[rep]+128*std::numeric_limits<double>::epsilon()*total)rep=k;
   a.geometry.representative=rep;
   if(problem.newton.numerics.avalancheCurrent==BvMethodsAvalancheCurrent::GssLaux) {
    a.geometry.usesGssLaux=true;
    const auto partial=tri3ElementEdgeBoxPartialVolumes(mesh,c);
    for(int k=0;k<3;++k)a.geometry.gssEdgePartialVolumes_cm2[k]=partial[k]*units.areaM2PerInternal()*1e4;
   }
   for(int k=0;k<3;++k){auto i=a.nodes[k];
    a.transportAreas_cm2[k]=weights[k]*units.areaM2PerInternal()*1e4;
    a.geometry.boxAreas_cm2[k]=avalancheWeights[k]*units.areaM2PerInternal()*1e4;
    auto imp=problem.doping.totalImpurity(i);a.bulkN[k]=bulk_->electronMobility(materials_[cid],imp,0,0);a.bulkP[k]=bulk_->holeMobility(materials_[cid],imp,0,0);
   }
   a.distanceGradient=a.shape[1]*(geo_.surface.nodeDistances[a.nodes[1]]-geo_.surface.nodeDistances[a.nodes[0]])+a.shape[2]*(geo_.surface.nodeDistances[a.nodes[2]]-geo_.surface.nodeDistances[a.nodes[0]]);
   Eigen::Matrix2d G,S;G.col(0)=a.shape[1];G.col(1)=a.shape[2];const auto& z=geo_.highFieldCells[cid].stencil;
   S<<z[0],z[2],z[1],z[3];auto P=(S/cfg_.surface.coordinateFieldFactor*G.inverse()).eval();
   a.geometry.currentProjection={P(0,0),P(0,1),P(1,0),P(1,1)};cells_.push_back(a);
  }
  lombardiN_.resize(n);lombardiP_.resize(n);
  for(int i=0;i<n;++i)if(electronActive_[i]) {
   const double distance=geo_.surface.nodeDistances[i],impurity=problem.doping.totalImpurity(i);
   lombardiN_[i]=prepareElementLombardi(distance,impurity,300.,cfg_.electronLombardi,cfg_);
   lombardiP_[i]=prepareElementLombardi(distance,impurity,300.,cfg_.holeLombardi,cfg_);
  }
 }

 BvMethodsTri3Evaluation evaluate(const BvMethodsQfState& state,bool derivatives,const VectorXd* columns=nullptr,bool reuseGeneration=true,CanaliEvaluationOptions canali={},bool scalar=false,
     const BvMethodsRepresentativeMobilities* mobilityValues=nullptr,bool auditMobility=false) const {
  const auto units=problem_.newton.inputScaling.unitSystem();const int n=problem_.mesh.numNodes();BvMethodsTri3Evaluation out;
  out.electronGeneration=VectorXd::Zero(n);out.holeGeneration=VectorXd::Zero(n);out.absoluteSource=VectorXd::Zero(3*n);out.source=VectorXd::Zero(3*n);out.allGeneration=VectorXd::Zero(n);
  std::vector<Eigen::Triplet<double>> entries;if(derivatives)entries.reserve(cells_.size()*54);
  const auto nodes=nodeStatistics(state,false,columns);
  const bool useScalar=scalar && !derivatives && supportsTransportValues(state) &&
      reuseGeneration && canali.avalancheValueOnly;
  if(mobilityValues && (!useScalar || mobilityValues->size()!=cells_.size()))
   throw std::invalid_argument("Invalid representative mobility buffer for AV");
  std::size_t cellIndex=0;
  for(const auto&a:cells_){
   const auto index=cellIndex++;
   auto accumulate=[&](double genN,double genP,const auto& nodeGeneration) {
    for(int k=0;k<3;++k) {
     const auto node=a.nodes[k];
     out.electronGeneration[node]+=a.geometry.boxAreas_cm2[k]*genN*factor_;
     out.holeGeneration[node]+=a.geometry.boxAreas_cm2[k]*genP*factor_;
     const double g=nodeGeneration(k);
     if(!std::isfinite(g))throw std::runtime_error("Nonfinite candidate AV source");
     out.allGeneration[node]+=g;
     if(electronActive_[node])out.source[n+node]-=g*factor_;
     if(holeActive_[node])out.source[2*n+node]-=g*factor_;
    }
   };
   if(useScalar) {
    std::array<double,3> psi,fn,fp,ns,ps,tn,tp;
    const double dc=units.internalConcentrationToCm3(1.)/potentialScale_;
    bool safe=experimental::avCoefficientsSafe(electron_)&&experimental::avCoefficientsSafe(hole_);
    for(int k=0;k<3;++k) {
     const auto i=a.nodes[k];
     const auto& node=nodes[i];
     psi[k]=state.psi[i];fn[k]=state.phin[i];fp[k]=state.phip[i];
     ns[k]=node.n;ps[k]=node.p;
     tn[k]=dc*node.tangentN;
     tp[k]=dc*node.tangentP;
     safe=safe&&std::abs(tn[k]/node.dn)<=1e30&&std::abs(tp[k]/node.dp)<=1e30;
    }
    const auto& cell=geo_.highFieldCells[a.id];
    const auto& s=cell.stencil;
    for(double v:s)safe=safe&&std::abs(v)<=1e30;
    const auto field = [&](const auto& values) {
     const auto [gx, gy] = tri3Gradient(values, s);
     return std::hypot(gx, gy);
    };
    const double en=field(cell.contact?psi:fn),hp=field(cell.contact?psi:fp);
    const double rawNormal=(psi[1]-psi[0])*(a.shape[1].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor)
                          +(psi[2]-psi[0])*(a.shape[2].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor);
    const double normal=rawNormal>0.?rawNormal:(rawNormal<0.?-rawNormal:0.);
    if(safe&&experimental::avFieldSafe(en)&&experimental::avFieldSafe(hp)) {
     const auto rep=a.geometry.representative;
     const auto node=a.nodes[rep];
     const double scale=units.internalMobilityToM2PerVS(1.)*1e4;
     const double cutoff=100./cfg_.internalFieldToVPerM;
     auto mobility=[&](bool e,double drive) {
      const auto& lp=e?cfg_.electronLombardi:cfg_.holeLombardi;
      const auto& hf=e?cfg_.electronField:cfg_.holeField;
      const auto surface=evaluatePreparedElementLombardi(e?a.bulkN[rep]:a.bulkP[rep],normal,
          e?lombardiN_[node]:lombardiP_[node],lp,cfg_);
      const auto value=static_cast<Real>(elementCanali(static_cast<vela::ArithmeticReal>(surface.mobility),
          static_cast<vela::ArithmeticReal>(drive),hf,cutoff));
      if(canali.audit && !CanaliEvaluationOptions::same(value,
          evaluateElementCanali(surface.mobility,drive,hf,cutoff).mobility))
       throw std::runtime_error("Canali avalanche value bit mismatch");
      return value*scale;
     };
     double mn,mp;
     if(mobilityValues) {
      const auto& stored=(*mobilityValues)[index];
      mn=stored[0];mp=stored[1];
      ++out.reusedMobilityCells;
      if(auditMobility || canali.audit) {
       if(!CanaliEvaluationOptions::same(mn,mobility(true,en)) ||
          !CanaliEvaluationOptions::same(mp,mobility(false,hp)))
        throw std::runtime_error("AV representative mobility bit mismatch");
       ++out.auditedMobilityPairs;
      }
     } else {
      mn=mobility(true,en);mp=mobility(false,hp);
     }
     const auto jn=experimental::avParticleCurrentValues(a.geometry,ns,fn,tn,vt_,true);
     const auto jp=experimental::avParticleCurrentValues(a.geometry,ps,fp,tp,vt_,false);
     const auto value=experimental::avGenerationValues(a.geometry,jn,jp,mn,mp,en,hp,electron_,hole_);
     out.zeroNormElectron+=value.electronZero;
     out.zeroNormHole+=value.holeZero;
     out.zeroActiveElectron+=value.electronActive;
     out.zeroActiveHole+=value.holeActive;
     accumulate(value.electron,value.hole,[&](int k){return value.node[k];});
     ++out.scalarCells;
     continue;
    }
    ++out.fallbackCells;
   }
   std::array<AvJet,3> psi,fn,fp,ns,ps;
   for(int k=0;k<3;++k){auto i=a.nodes[k];psi[k]=AvJet::variable(state.psi[i],k);fn[k]=AvJet::variable(state.phin[i],3+k);fp[k]=AvJet::variable(state.phip[i],6+k);
    const auto& node=nodes[i];ns[k]=AvJet(node.n);ps[k]=AvJet(node.p);
    const double dn=node.dn,dp=node.dp;
    ns[k].derivative[k]=dn;ns[k].derivative[3+k]=-dn;ps[k].derivative[k]=-dp;ps[k].derivative[6+k]=dp;
    if(state.populationCoordinates){
     // Jets retain the historical final V0 scatter; concentration columns
     // therefore carry their physical internal-unit factor divided by V0.
     const double dc=units.internalConcentrationToCm3(1.)/potentialScale_;
     ns[k].derivative.fill(0.);ps[k].derivative.fill(0.);
     ns[k].derivative[3+k]=dc*node.tangentN;ps[k].derivative[6+k]=dc*node.tangentP;
     fn[k].derivative.fill(0.);fp[k].derivative.fill(0.);
     fn[k].derivative[k]=1.;fp[k].derivative[k]=1.;
     fn[k].derivative[3+k]=-(dc*node.tangentN)/dn;fp[k].derivative[6+k]=(dc*node.tangentP)/dp;
    }
   }
   const auto& fieldCell = geo_.highFieldCells[a.id];
   const auto en = experimental::avNorm(tri3Gradient(fieldCell.contact ? psi : fn, fieldCell.stencil));
   const auto hp = experimental::avNorm(tri3Gradient(fieldCell.contact ? psi : fp, fieldCell.stencil));
   const auto rawNormal=(psi[1]-psi[0])*AvJet(a.shape[1].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor)+(psi[2]-psi[0])*AvJet(a.shape[2].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor);
   const auto normal=detail::dualAbs(rawNormal);int rep=a.geometry.representative;auto node=a.nodes[rep];
   auto mobility=[&](bool e,const AvJet& drive){
    const double low=e?a.bulkN[rep]:a.bulkP[rep];const auto& lp=e?cfg_.electronLombardi:cfg_.holeLombardi;const auto& hf=e?cfg_.electronField:cfg_.holeField;
    const auto surface=evaluatePreparedElementLombardi(low,normal.value,e?lombardiN_[node]:lombardiP_[node],lp,cfg_);
    const double cutoff=100./cfg_.internalFieldToVPerM;
    const double scale=units.internalMobilityToM2PerVS(1.)*1e4;
    if(!derivatives && canali.avalancheValueOnly) {
     const auto value=static_cast<Real>(elementCanali(static_cast<vela::ArithmeticReal>(surface.mobility),
        static_cast<vela::ArithmeticReal>(drive.value),hf,cutoff));
     if(canali.audit && !CanaliEvaluationOptions::same(value,evaluateElementCanali(surface.mobility,drive.value,hf,cutoff).mobility))
      throw std::runtime_error("Canali avalanche value bit mismatch");
     return AvJet(value*scale);
    }
    const auto high=canali.full(surface.mobility,drive.value,hf,cutoff);
    AvJet mu(high.mobility*scale);
    for(int k=0;k<9;++k)mu.derivative[k]=(high.lowDerivative*surface.fieldDerivative*normal.derivative[k]+high.fieldDerivative*drive.derivative[k])*scale;
    return mu;
   };
   const auto mn=mobility(true,en),mp=mobility(false,hp);
   std::array<AvJet,2> directN,directP;
   if(directPopulationSg_){
    if(!state.populationCoordinates)throw std::invalid_argument("Direct SG requires population coordinates");
    std::array<double,3> ni,nc,nv;for(int k=0;k<3;++k){ni[k]=ni_[a.nodes[k]];nc[k]=nc_[a.nodes[k]];nv[k]=nv_[a.nodes[k]];}
    directN=experimental::populationAvCurrent(a.geometry,ns,psi,ni,nc,vt_,true);
    directP=experimental::populationAvCurrent(a.geometry,ps,psi,ni,nv,vt_,false);
   }
   const auto result=experimental::evaluateBvMethodsTri3Avalanche(a.geometry,ns,ps,fn,fp,mn,mp,en,hp,vt_,electron_,hole_,
       directPopulationSg_?&directN:nullptr,directPopulationSg_?&directP:nullptr);
   out.zeroNormElectron+=result.electronZeroNorm;
   out.zeroNormHole+=result.holeZeroNorm;
   double genN=result.electronGeneration_cm3_s,genP=result.holeGeneration_cm3_s;
   if(reuseGeneration) {
    out.zeroActiveElectron+=result.electronActiveAtZero;
    out.zeroActiveHole+=result.holeActiveAtZero;
   } else {
    if(result.electronZeroNorm&&experimental::avCoefficient(en,electron_).value>0)++out.zeroActiveElectron;
    if(result.holeZeroNorm&&experimental::avCoefficient(hp,hole_).value>0)++out.zeroActiveHole;
    genN=(experimental::avCoefficient(en,electron_)*mn*experimental::avNorm(result.electronParticleCurrentWithoutMobility)).value;
    genP=(experimental::avCoefficient(hp,hole_)*mp*experimental::avNorm(result.holeParticleCurrentWithoutMobility)).value;
   }
   accumulate(genN,genP,[&](int k){return result.nodeGeneration_per_cm_s[k].value;});
   if(derivatives)for(int k=0;k<3;++k) {
    const auto node=a.nodes[k];
    const auto& g=result.nodeGeneration_per_cm_s[k];
    for(int block=1;block<=2;++block) {
     if(!(block==1?electronActive_[node]:holeActive_[node]))continue;
     const int row=block*n+node;
     for(int j=0;j<9;++j){const double tangent=j<3?1.:(j<6?nodes[a.nodes[j%3]].tangentN:nodes[a.nodes[j%3]].tangentP);const int column=(j/3)*n+a.nodes[j%3];const double d=columns?(-g.derivative[j]*factor_*potentialScale_)*((*columns)[column]/tangent):(-g.derivative[j]*factor_*potentialScale_)/tangent;if(d!=0)entries.emplace_back(row,(j/3)*n+a.nodes[j%3],d);}
    }
   }
  }
  out.jacobian.resize(3*n,3*n);if(derivatives)out.jacobian.setFromTriplets(entries.begin(),entries.end());return out;
 }
 // Diagnostic approximate linearization only: retain the live mobility value,
 // but omit its tangent. Production calls retain the full derivative default.
 BvMethodsTri3Evaluation transport(const BvMethodsQfState& state,bool derivatives,bool generalized,bool coefficientAveraging=false,bool freezeMobilityDerivatives=false,const VectorXd* columns=nullptr,double statisticalFloor_cm3=0.,
    bool reuseLombardiPowers=false,CanaliEvaluationOptions canali={}) const {
  const auto units=problem_.newton.inputScaling.unitSystem();const int n=problem_.mesh.numNodes();BvMethodsTri3Evaluation out;
  out.electronGeneration=VectorXd::Zero(n);out.holeGeneration=VectorXd::Zero(n);out.absoluteSource=VectorXd::Zero(3*n);out.source=VectorXd::Zero(3*n);out.allGeneration=VectorXd::Zero(n);
  std::vector<Eigen::Triplet<double>> entries;if(derivatives)entries.reserve(cells_.size()*54);
  const auto nodes=nodeStatistics(state,true,columns);
  for(const auto&a:cells_){std::array<AvJet,3> psi,fn,fp,ns,ps;
   for(int k=0;k<3;++k){auto i=a.nodes[k];psi[k]=AvJet::variable(state.psi[i],k);fn[k]=AvJet::variable(state.phin[i],3+k);fp[k]=AvJet::variable(state.phip[i],6+k);
    const auto& node=nodes[i];ns[k]=AvJet(node.n);ps[k]=AvJet(node.p);
    const double dn=node.dn,dp=node.dp;
    ns[k].derivative[k]=dn;ns[k].derivative[3+k]=-dn;ps[k].derivative[k]=-dp;ps[k].derivative[6+k]=dp;
    if(state.populationCoordinates){
     // Jets retain the historical final V0 scatter; concentration columns
     // therefore carry their physical internal-unit factor divided by V0.
     const double dc=units.internalConcentrationToCm3(1.)/potentialScale_;
     ns[k].derivative.fill(0.);ps[k].derivative.fill(0.);
     ns[k].derivative[3+k]=dc*node.tangentN;ps[k].derivative[6+k]=dc*node.tangentP;
     fn[k].derivative.fill(0.);fp[k].derivative.fill(0.);
     fn[k].derivative[k]=1.;fp[k].derivative[k]=1.;
     fn[k].derivative[3+k]=-(dc*node.tangentN)/dn;fp[k].derivative[6+k]=(dc*node.tangentP)/dp;
    }
   }
   const auto& fieldCell = geo_.highFieldCells[a.id];
   const auto en = experimental::avNorm(tri3Gradient(fieldCell.contact ? psi : fn, fieldCell.stencil));
   const auto hp = experimental::avNorm(tri3Gradient(fieldCell.contact ? psi : fp, fieldCell.stencil));
   const auto rawNormal=(psi[1]-psi[0])*AvJet(a.shape[1].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor)+(psi[2]-psi[0])*AvJet(a.shape[2].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor);
   const auto normal=detail::dualAbs(rawNormal);int rep=a.geometry.representative;auto node=a.nodes[rep];
   std::optional<LombardiCellPowers> powers;
   if (reuseLombardiPowers)
    powers.emplace();
   auto mobility=[&](bool e,const AvJet& drive,int local){auto node=a.nodes[local];
    const double low=e?a.bulkN[local]:a.bulkP[local];const auto& lp=e?cfg_.electronLombardi:cfg_.holeLombardi;const auto& hf=e?cfg_.electronField:cfg_.holeField;
    const double scale=units.internalMobilityToM2PerVS(1.)*1e4;
    if(cfg_.model=="constant")return AvJet(low*scale);
    const auto surface=evaluatePreparedElementLombardi(
        low,normal.value,e?lombardiN_[node]:lombardiP_[node],lp,cfg_,
        powers ? &*powers : nullptr,e);
    const auto high=canali.full(surface.mobility,drive.value,hf,100./cfg_.internalFieldToVPerM);
    AvJet mu(high.mobility*scale);
    for(int k=0;k<9;++k)mu.derivative[k]=(high.lowDerivative*surface.fieldDerivative*normal.derivative[k]+high.fieldDerivative*drive.derivative[k])*scale;
    return mu;
   };

   std::array<AvJet,3> mn,mp,etaN,etaP,gn{AvJet(1),AvJet(1),AvJet(1)},gp=gn;
   for(int k=0;k<3;++k){
    mn[k]=mobility(true,en,k);mp[k]=mobility(false,hp,k);
    if(freezeMobilityDerivatives){mn[k].derivative.fill(0.);mp[k].derivative.fill(0.);}
    const auto i=a.nodes[k];etaN[k]=AvJet(nodes[i].etaN);etaP[k]=AvJet(nodes[i].etaP);
    etaN[k].derivative[k]=1/vt_;etaN[k].derivative[3+k]=-1/vt_;
    etaP[k].derivative[k]=-1/vt_;etaP[k].derivative[6+k]=1/vt_;
    if(state.populationCoordinates){
     etaN[k].derivative.fill(0.);etaP[k].derivative.fill(0.);
     etaN[k].derivative[3+k]=-fn[k].derivative[3+k]/vt_;
     etaP[k].derivative[6+k]=fp[k].derivative[6+k]/vt_;
    }
   }
   auto einstein=[&](const std::array<AvJet,3>& population,const std::array<AvJet,3>& eta,int k){
    int j=(k+1)%3;if(directPopulationSg_ && a.geometry.nodeOrder[k]>a.geometry.nodeOrder[j])std::swap(k,j);
    auto log=experimental::avLogRatio(population[j],population[k]);
    if(std::abs(log.value)>1e-10)return (eta[j]-eta[k])/log;
    const auto mid=(eta[j]+eta[k])/AvJet(2);
    const auto stats=statisticsJet(mid.value);
    AvJet value(std::max(1.,stats.value/stats.first));
    if(value.value>1.)for(int z=0;z<9;++z)value.derivative[z]=(1-stats.value*stats.second/(stats.first*stats.first))*mid.derivative[z];
    return value;
   };
   if(generalized)for(int k=0;k<3;++k){gn[k]=einstein(ns,etaN,k);gp[k]=einstein(ps,etaP,k);}
   const auto& coefficients=problem_.mesh.poissonCellEdgeCoefficients(a.id);
   std::array<double,3> edge{coefficients[0],coefficients[1],coefficients[2]};
   auto averagingAreas=a.transportAreas_cm2;
   if(coefficientAveraging){auto weights=detail::cellBoxNodeMeasures(problem_.mesh,a.id);for(int k=0;k<3;++k)averagingAreas[k]=weights[k]*units.areaM2PerInternal()*1e4;}
   std::array<AvJet,3> an,ap;
   if(directPopulationSg_){
    if(!state.populationCoordinates)throw std::invalid_argument("Direct SG requires population coordinates");
    for(int k=0;k<3;++k){int u=k,v=(k+1)%3;if(a.geometry.nodeOrder[u]>a.geometry.nodeOrder[v])std::swap(u,v);
     an[k]=experimental::populationSgArgument(ns[u],ns[v],psi[u],psi[v],ni_[a.nodes[u]],ni_[a.nodes[v]],nc_[a.nodes[u]],nc_[a.nodes[v]],vt_,true,gn[k]);
     ap[k]=experimental::populationSgArgument(ps[u],ps[v],psi[u],psi[v],ni_[a.nodes[u]],ni_[a.nodes[v]],nv_[a.nodes[u]],nv_[a.nodes[v]],vt_,false,gp[k]);}
   }
   std::array<bool,3> freezeN{},freezeP{};
   if(derivatives && state.populationCoordinates && statisticalFloor_cm3>0.)
    for(int k=0;k<3;++k){const int j=(k+1)%3;
     // Permit only rounding incurred by the internal <-> cm^-3 unit mapping.
     const double floor=statisticalFloor_cm3*(1.+16.*std::numeric_limits<double>::epsilon());
     freezeN[k]=ns[k].value<=floor || ns[j].value<=floor;
     freezeP[k]=ps[k].value<=floor || ps[j].value<=floor;
    }
   auto tn=experimental::evaluateBvMethodsTri3Transport(averagingAreas,edge,ns,fn,mn,vt_,true,gn,directPopulationSg_?&an:nullptr,a.geometry.nodeOrder,freezeN);
   auto tp=experimental::evaluateBvMethodsTri3Transport(averagingAreas,edge,ps,fp,mp,vt_,false,gp,directPopulationSg_?&ap:nullptr,a.geometry.nodeOrder,freezeP);
   for(int k=0;k<3;++k)for(int block=1;block<=2;++block){
    const int row=block*n+a.nodes[k];const auto& f=block==1?tn.nodeFlux_per_cm_s[k]:tp.nodeFlux_per_cm_s[k];out.source[row]+=f.value*factor_;out.absoluteSource[row]+=(std::abs((block==1?tn.edgeFlux_per_cm_s[k]:tp.edgeFlux_per_cm_s[k]).value)+std::abs((block==1?tn.edgeFlux_per_cm_s[(k+2)%3]:tp.edgeFlux_per_cm_s[(k+2)%3]).value))*factor_;
    if(derivatives)for(int j=0;j<9;++j){const double tangent=j<3?1.:(j<6?nodes[a.nodes[j%3]].tangentN:nodes[a.nodes[j%3]].tangentP);const int column=(j/3)*n+a.nodes[j%3];double value=columns?(f.derivative[j]*factor_*potentialScale_)*((*columns)[column]/tangent):(f.derivative[j]*factor_*potentialScale_)/tangent;if(value!=0)entries.emplace_back(row,(j/3)*n+a.nodes[j%3],value);}
   }
  }
  out.jacobian.resize(3*n,3*n);if(derivatives)out.jacobian.setFromTriplets(entries.begin(),entries.end());return out;
 }
 bool supportsTransportValues(const BvMethodsQfState& state) const {
  return state.populationCoordinates && !directPopulationSg_;
 }

 BvMethodsTri3Evaluation transportValues(const BvMethodsQfState& state,bool generalized,
    bool coefficientAveraging,bool reuseLombardiPowers,bool valueOnlyCanali,bool auditCanali,
    BvMethodsRepresentativeMobilities* mobilityValues=nullptr) const {
  if(!supportsTransportValues(state))throw std::invalid_argument("Scalar Tri3 transport requires population QF SG");
  const auto units=problem_.newton.inputScaling.unitSystem();
  const int n=problem_.mesh.numNodes();
  BvMethodsTri3Evaluation out;
  out.electronGeneration=VectorXd::Zero(n);out.holeGeneration=VectorXd::Zero(n);
  out.absoluteSource=VectorXd::Zero(3*n);out.source=VectorXd::Zero(3*n);out.allGeneration=VectorXd::Zero(n);
  out.jacobian.resize(3*n,3*n);
  const auto nodes=nodeStatistics(state,true,nullptr,false);
  const double dc=units.internalConcentrationToCm3(1.)/potentialScale_;
  const double mobilityScale=units.internalMobilityToM2PerVS(1.)*1e4;
  if(mobilityValues){mobilityValues->clear();mobilityValues->reserve(cells_.size());}
  for(const auto& a:cells_) {
   if(mobilityValues)mobilityValues->emplace_back();
   std::array<double,3> psi,fn,fp;
   for(int k=0;k<3;++k) {
    const auto i=a.nodes[k];psi[k]=state.psi[i];fn[k]=state.phin[i];fp[k]=state.phip[i];
   }
   const auto& cell=geo_.highFieldCells[a.id];
   const auto& s=cell.stencil;
   const auto field = [&](const auto& values) {
    const auto [gx, gy] = tri3Gradient(values, s);
    return std::hypot(gx, gy);
   };
   const double en=field(cell.contact?psi:fn),hp=field(cell.contact?psi:fp);
   const double rawNormal=(psi[1]-psi[0])*(a.shape[1].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor)
                         +(psi[2]-psi[0])*(a.shape[2].dot(a.distanceGradient)*cfg_.surface.coordinateFieldFactor);
   const double normal=rawNormal>0.?rawNormal:(rawNormal<0.?-rawNormal:0.);
   std::optional<LombardiCellPowers> powers;
   if(reuseLombardiPowers)powers.emplace();
   const auto& coefficients=problem_.mesh.poissonCellEdgeCoefficients(a.id);
   const std::array<double,3> edge{coefficients[0],coefficients[1],coefficients[2]};
   auto areas=a.transportAreas_cm2;
   if(coefficientAveraging) {
    const auto weights=detail::cellBoxNodeMeasures(problem_.mesh,a.id);
    for(int k=0;k<3;++k)areas[k]=weights[k]*units.areaM2PerInternal()*1e4;
   }
   for(int block=1;block<=2;++block) {
    const bool electron=block==1;
    const auto& lp=electron?cfg_.electronLombardi:cfg_.holeLombardi;
    const auto& hf=electron?cfg_.electronField:cfg_.holeField;
    const auto& prepared=electron?lombardiN_:lombardiP_;
    const auto& bulk=electron?a.bulkN:a.bulkP;
    const double drive=electron?en:hp;
    std::array<double,3> density,tangent,eta,mobility,logRatio,einstein{1.,1.,1.};
    for(int k=0;k<3;++k) {
     const auto i=a.nodes[k];const auto& node=nodes[i];
     density[k]=electron?node.n:node.p;
     tangent[k]=dc*(electron?node.tangentN:node.tangentP);
     eta[k]=electron?node.etaN:node.etaP;
     if(cfg_.model=="constant"){
      mobility[k]=bulk[k]*mobilityScale;
      continue;
     }
     const auto surface=evaluatePreparedElementLombardi(bulk[k],normal,prepared[i],lp,cfg_,powers?&*powers:nullptr,electron);
     const double cutoff=100./cfg_.internalFieldToVPerM;
     const auto valueOnly=[&] {
      return static_cast<Real>(elementCanali(static_cast<vela::ArithmeticReal>(surface.mobility),
          static_cast<vela::ArithmeticReal>(drive),hf,cutoff));
     };
     const double mobilityValue=valueOnlyCanali ? valueOnly()
        : evaluateElementCanali(surface.mobility,drive,hf,cutoff).mobility;
     if(auditCanali) {
      const double reference=evaluateElementCanali(surface.mobility,drive,hf,cutoff).mobility;
      if(std::bit_cast<std::uint64_t>(mobilityValue)!=std::bit_cast<std::uint64_t>(reference))
       throw std::runtime_error("BVMethods value-only Canali mobility mismatch");
     }
     mobility[k]=mobilityValue*mobilityScale;
     if(mobilityValues && k==a.geometry.representative)
      mobilityValues->back()[electron?0:1]=mobility[k];
    }
    for(int k=0;k<3;++k) {
     const int j=(k+1)%3;
     logRatio[k]=experimental::populationTransportLogRatio(density[j],density[k],tangent[j],tangent[k]);
     if(!generalized)continue;
     if(std::abs(logRatio[k])>1e-10)einstein[k]=(eta[j]-eta[k])/logRatio[k];
     else {
      const double mid=(eta[j]+eta[k])/2.;
      const auto stats=statisticsJet(mid);
      einstein[k]=std::max(1.,stats.value/stats.first);
     }
    }
    const auto flux=experimental::evaluateBvMethodsTri3TransportValues(areas,edge,density,electron?fn:fp,
       mobility,vt_,electron,einstein,logRatio);
    for(int k=0;k<3;++k) {
     const int row=block*n+a.nodes[k];
     out.source[row]+=flux.nodes[k]*factor_;
     out.absoluteSource[row]+=(std::abs(flux.edges[k])+std::abs(flux.edges[(k+2)%3]))*factor_;
    }
   }
  }
  return out;
 }
};
} // namespace vela::detail
