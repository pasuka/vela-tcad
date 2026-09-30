#include "vela/equation/SplitDDRuntime.h"
#include "vela/equation/SplitDDOperator.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/solver/GummelSolver.h"
#include "vela/post/ContactCurrent.h"
#include "vela/core/PhysicalConstants.h"
#include <boost/uuid/detail/sha1.hpp>
#include <iomanip>
#include <sstream>
#include <algorithm>

namespace vela {
using split_dd::Wide;
namespace {
VectorXd project(const std::vector<Wide>& values) {
    VectorXd out(values.size());for(std::size_t i=0;i<values.size();++i)out(i)=static_cast<Real>(values[i]);return out;
}
split_dd::State savedState(const DDSolution& s) {
    const int N=s.psi.size();
    if(s.packedState.size()!=3*N || s.packedLow.size()!=3*N ||
       s.electronQfReference.size()!=N || s.holeQfReference.size()!=N)
        throw std::invalid_argument("Incomplete split DD checkpoint");
    std::vector<numerics::SplitCoordinate> pairs;std::vector<double> en,hp;
    for(int k=0;k<3*N;++k)pairs.push_back({s.packedState(k),s.packedLow(k)});
    for(int i=0;i<N;++i){en.push_back(s.electronQfReference(i));hp.push_back(s.holeQfReference(i));}
    return split_dd::State(pairs,en,hp,s.packedPotentialScale_V,s.packedMeshFingerprint);
}
}
bool DDSolution::hasConsistentSplitPackedState() const {
    try {
        auto state=savedState(*this);const int N=state.size();
        if(N==0 || phin.size()!=N || phip.size()!=N || phinIncrement.size()!=N || phipIncrement.size()!=N)return false;
        for(int i=0;i<N;++i) {
            if(psi(i)!=static_cast<Real>(state.potential(0,i)) ||
               phin(i)!=static_cast<Real>(state.potential(1,i)) ||
               phip(i)!=static_cast<Real>(state.potential(2,i)) ||
               phinIncrement(i)!=static_cast<Real>(state.coordinate(N+i)*Wide(packedPotentialScale_V)) ||
               phipIncrement(i)!=static_cast<Real>(state.coordinate(2*N+i)*Wide(packedPotentialScale_V)))return false;
        }
        return true;
    } catch(const std::exception&) {return false;}
}
namespace {
Wide physicalAt(const DDSolution& s,int block,int node) {
    const int N=s.psi.size();
    if(block<0 || block>2 || node<0 || node>=N || s.packedState.size()!=3*N || s.packedLow.size()!=3*N)
        throw std::invalid_argument("Invalid split field index");
    const int k=block*N+node;
    Wide result=(Wide(s.packedState(k))+Wide(s.packedLow(k)))*Wide(s.packedPotentialScale_V);
    if(block==1)result+=Wide(s.electronQfReference(node));
    if(block==2)result+=Wide(s.holeQfReference(node));
    return result;
}
}
Real DDSolution::splitPotentialAt(int block,int node) const {return static_cast<Real>(physicalAt(*this,block,node));}
Real DDSolution::splitQuasiFermiDifferenceAt(int node) const {return static_cast<Real>(physicalAt(*this,2,node)-physicalAt(*this,1,node));}

struct SplitDDRuntime::Impl {
    nlohmann::json model;split_dd::Operator op;VectorXd low;
    mutable VectorXd cachedX,cachedLow;mutable split_dd::Evaluation cached;
    std::string key;
    explicit Impl(nlohmann::json j):model(std::move(j)),op(model),low(VectorXd::Zero(3*model["nodes"].size())),key(model["mesh_fingerprint"]) {}
    split_dd::State state(const VectorXd& x) const {
        const int N=model["nodes"].size();
        if(x.size()!=3*N || low.size()!=x.size())throw std::invalid_argument("Split runtime vector size mismatch");
        std::vector<numerics::SplitCoordinate> pairs;std::vector<double> en,hp;
        for(int k=0;k<x.size();++k)pairs.push_back({x(k),low(k)});
        for(int i=0;i<N;++i){en.push_back(model["nodes"][i][7]);hp.push_back(model["nodes"][i][8]);}
        return split_dd::State(pairs,en,hp,model["potential_scale"],key);
    }
    const split_dd::Evaluation& evaluate(const VectorXd& x) const {
        if(cachedX.size()!=x.size() || cachedLow.size()!=low.size() ||
           !(cachedX.array()==x.array()).all() || !(cachedLow.array()==low.array()).all()) {
            cached=op.evaluate(state(x));cachedX=x;cachedLow=low;
        }
        return cached;
    }
    std::vector<Wide> residual(const VectorXd& x,const CoupledDDBoundaryConditions& bcs) const {
        if(!bcs.thermionic.empty())throw std::invalid_argument("Split runtime does not support thermionic boundaries");
        auto r=evaluate(x).residual;auto s=state(x);const int N=s.size();const double scale=model["potential_scale"];
        for(const auto& [i,z]:bcs.psi)r[i]=s.coordinate(i)-Wide(z);
        for(const auto& [i,z]:bcs.phin)r[N+i]=s.coordinate(N+i)-Wide(z-model["nodes"][i][7].get<double>()/scale);
        for(const auto& [i,z]:bcs.phip)r[2*N+i]=s.coordinate(2*N+i)-Wide(z-model["nodes"][i][8].get<double>()/scale);
        return r;
    }
};
SplitDDRuntime::SplitDDRuntime(nlohmann::json model):impl_(std::make_unique<Impl>(std::move(model))) {}
SplitDDRuntime::~SplitDDRuntime()=default;
void SplitDDRuntime::setLow(const VectorXd& low) {if(low.size()!=impl_->low.size() || !low.allFinite())throw std::invalid_argument("Invalid split coordinate tail");impl_->low=low;}
VectorXd SplitDDRuntime::low() const {return impl_->low;}
const std::string& SplitDDRuntime::fingerprint() const {return impl_->key;}
VectorXd SplitDDRuntime::candidate(const VectorXd& x,const VectorXd& step,Real alpha,const VectorXd& baseLow) {
    setLow(baseLow);auto s=impl_->state(x);std::vector<double> dx(step.data(),step.data()+step.size());auto next=s.shifted(dx,alpha).checkpoint();VectorXd hi(x.size()),lo(x.size());
    for(int i=0;i<x.size();++i){hi(i)=next["coordinates"][i][0];lo(i)=next["coordinates"][i][1];}setLow(lo);return hi;
}
VectorXd SplitDDRuntime::residual(const VectorXd& x,const CoupledDDBoundaryConditions& bcs) const {return project(impl_->residual(x,bcs));}
VectorXd SplitDDRuntime::symmetricDifference(const VectorXd& x,const VectorXd& step,
    const CoupledDDBoundaryConditions& bcs,VectorXd& forward,VectorXd& backward) {
    const VectorXd base=low();
    try {
        auto plus=candidate(x,step,1.,base);auto a=impl_->residual(plus,bcs);
        auto minus=candidate(x,step,-1.,base);auto b=impl_->residual(minus,bcs);
        forward=project(a);backward=project(b);
        for(std::size_t i=0;i<a.size();++i)a[i]=(a[i]-b[i])/2;
        setLow(base);return project(a);
    } catch(...) {setLow(base);throw;}
}
VectorXd SplitDDRuntime::density(const VectorXd& x,bool electron) const {const auto& o=impl_->evaluate(x);return project(electron?o.n:o.p);}
Real SplitDDRuntime::current(const VectorXd& x,const std::string& name) const {return static_cast<Real>(impl_->evaluate(x).currents.at(name));}
std::vector<CoupledDDCarrierTermDiagnostic> SplitDDRuntime::terms(const VectorXd& x,const CoupledDDBoundaryConditions& bcs) const {
    const auto& o=impl_->evaluate(x);const int N=o.n.size();auto rr=impl_->residual(x,bcs);std::vector<CoupledDDCarrierTermDiagnostic> out(N);
    std::vector<Wide> ef(N),hf(N),ea(N),ha(N);std::vector<bool> active(N,false);
    for(std::size_t k=0;k<o.electronFlux.size();++k){const auto& e=impl_->model["transport_edges"][k];int i=e["i"],j=e["j"];auto en=o.electronFlux[k],hp=o.holeFlux[k];ef[i]+=en;ef[j]-=en;hf[i]+=hp;hf[j]-=hp;ea[i]+=abs(en);ea[j]+=abs(en);ha[i]+=abs(hp);ha[j]+=abs(hp);if(!e["weights"].empty())active[i]=active[j]=true;}
    for(int i=0;i<N;++i){auto& t=out[i];t.nodeId=i;t.electronDensity_m3=static_cast<Real>(o.n[i]);t.holeDensity_m3=static_cast<Real>(o.p[i]);
        t.electronContinuityActive=active[i]&&!bcs.phin.contains(i);t.holeContinuityActive=active[i]&&!bcs.phip.contains(i);
        t.electronFlux=static_cast<Real>(ef[i]);t.holeFlux=static_cast<Real>(hf[i]);t.electronFluxAbsSum=static_cast<Real>(ea[i]);t.holeFluxAbsSum=static_cast<Real>(ha[i]);
        Wide source=o.srh[i]*Wide(impl_->model["node_physics"][i][4].get<double>());t.electronRecombination=t.holeRecombination=static_cast<Real>(source);
        t.electronResidual=static_cast<Real>(rr[N+i]);t.holeResidual=static_cast<Real>(rr[2*N+i]);
        if(!active[i]){t.electronGauge=t.electronResidual;t.holeGauge=t.holeResidual;}
        if(bcs.phin.contains(i)){t.electronBoundary=t.electronResidual;t.electronFlux=t.electronFluxAbsSum=t.electronRecombination=t.electronGauge=0.;}
        if(bcs.phip.contains(i)){t.holeBoundary=t.holeResidual;t.holeFlux=t.holeFluxAbsSum=t.holeRecombination=t.holeGauge=0.;}
    }return out;
}
std::vector<CoupledDDEdgeFluxDiagnostic> SplitDDRuntime::edges(const VectorXd& x) const {
    const auto& o=impl_->evaluate(x);const auto& model=impl_->model;
    std::vector<CoupledDDEdgeFluxDiagnostic> out;
    auto cast=[](const Wide& v){return static_cast<Real>(v);};
    auto bernoulli=[](const Wide& u)->Wide{return u==0?Wide(1):u/expm1(u);};
    const Wide lineFactor=Wide(model["current_factor"].get<double>())*Wide(1e6)/Wide(constants::q);
    for(std::size_t k=0;k<o.electronFlux.size();++k) {
        const auto& e=model["transport_edges"][k];int i=e["i"],j=e["j"];
        CoupledDDEdgeFluxDiagnostic d;d.edgeId=e["id"];d.node0=i;d.node1=j;d.length_m=e["length"];d.couple_m=e["couple"];
        d.psi0_V=cast(o.psi[i]);d.psi1_V=cast(o.psi[j]);d.phin0_V=cast(o.phin[i]);d.phin1_V=cast(o.phin[j]);d.phip0_V=cast(o.phip[i]);d.phip1_V=cast(o.phip[j]);
        d.electronDensity0_m3=cast(o.n[i]);d.electronDensity1_m3=cast(o.n[j]);d.holeDensity0_m3=cast(o.p[i]);d.holeDensity1_m3=cast(o.p[j]);d.ni0_m3=model["nodes"][i][9];d.ni1_m3=model["nodes"][j][9];
        const Wide field=(o.psi[j]-o.psi[i])*Wide(model["field_factor"].get<double>())/Wide(e["length"].get<double>());
        d.electricField_V_m=cast(abs(field));d.electronMobilityDriveInternal=model.contains("hfs")?std::numeric_limits<Real>::quiet_NaN():cast(abs(field));
        const Wide mn=o.edgeMun[k],mp=o.edgeMup[k];
        d.electronMobility_m2_V_s=cast(mn);d.holeMobility_m2_V_s=cast(mp);
        d.electronFlux=cast(o.electronFlux[k]);d.holeFlux=cast(o.holeFlux[k]);
        d.electronParticleLineFlux_per_m_s=cast(o.electronFlux[k]*lineFactor);d.holeParticleLineFlux_per_m_s=cast(o.holeFlux[k]*lineFactor);
        if(!e["weights"].empty()) {
            Wide vt=model["nodes"][i][10].get<double>(),eta=(o.psi[j]-o.psi[i])/vt+log(Wide(d.ni1_m3)/Wide(d.ni0_m3));
            Wide left=o.n[i]*bernoulli(-eta),right=o.n[j]*bernoulli(eta),coef=mn*Wide(e["coefficient"].get<double>());
            d.electronSgBoltzmannDecompositionAvailable=true;d.electronSgEta=cast(eta);d.electronSgBernoulliMinusEta=cast(bernoulli(-eta));d.electronSgBernoulliEta=cast(bernoulli(eta));
            d.electronSgLeftTermInternal=cast(left);d.electronSgRightTermInternal=cast(right);d.electronSgSignedDifferenceInternal=cast(left-right);
            d.electronSgReconstructedFluxScaled=cast(coef*(left-right));d.electronSgStableFluxScaled=d.electronSgHighPrecisionReferenceFluxScaled=d.electronFlux;
            d.electronSgHighPrecisionReferenceTermScaleScaled=cast(abs(coef)*(abs(left)+abs(right)));
            d.electronSgLogLeftOverRight=cast((o.phin[j]-o.phin[i])/vt);d.electronSgRightFactorFluxScaled=cast(coef*right);
            d.electronQfReference0_V=model["nodes"][i][7];d.electronQfReference1_V=model["nodes"][j][7];
            d.electronSgPhin0Relative_V=cast(o.phin[i]-Wide(d.electronQfReference0_V));d.electronSgPhin1Relative_V=cast(o.phin[j]-Wide(d.electronQfReference0_V));
            d.electronSgCancellationCondition=cast((abs(left)+abs(right))/std::max(Wide(abs(left-right)),Wide("1e-100")));
        }
        out.push_back(d);
    }return out;
}
ContactCurrentDetailedResult SplitDDRuntime::contact(const VectorXd& x,const std::string& name) const {
    const auto& o=impl_->evaluate(x);const auto& model=impl_->model;
    const Wide factor=Wide(model["current_factor"].get<double>())*Wide(1e6);
    const Wide cs=model["continuity_scale"].get<double>();
    ContactCurrentDetailedResult result;Wide en=0,hp=0,ed=0,hd=0,ef=0,hf=0;
    auto cast=[](const Wide& value){return static_cast<Real>(value);};
    auto bernoulli=[](const Wide& u)->Wide{return u==0?Wide(1):u/expm1(u);};
    for(std::size_t k=0;k<o.electronFlux.size();++k) {
        const auto& edge=model["transport_edges"][k];
        if(!edge["ports"].contains(name))continue;
        const int i=edge["i"],j=edge["j"];const Wide sign=edge["ports"][name].get<double>();
        const Wide ec=-sign*factor*o.electronFlux[k],hc=-sign*factor*o.holeFlux[k];
        const Wide mn=o.edgeMun[k],mp=o.edgeMup[k];
        const Wide vt=model["nodes"][i][10].get<double>();
        const Wide u=(o.psi[j]-o.psi[i])/vt;
        // Electric drift uses the physical electric field. The complementary
        // diffusion component includes the intrinsic-density (BGN) gradient.
        const Wide eDrift=-sign*factor*mn*Wide(edge["coefficient"].get<double>())*u*(o.n[i]+o.n[j])/2;
        const Wide hDrift=sign*factor*mp*Wide(edge["coefficient"].get<double>())*u*(o.p[i]+o.p[j])/2;
        en+=ec;hp+=hc;ed+=eDrift;hd+=hDrift;ef+=ec-eDrift;hf+=hc-hDrift;
        ContactCurrentEdgeDiagnostic d;d.edgeId=edge["id"];d.node0=i;d.node1=j;
        d.edgeLength_m=edge["length_m"];d.edgeCouple_m=edge["couple_m"];d.outwardSign=cast(sign);
        d.bernoulliU=cast(u);d.bernoulliBplus=cast(bernoulli(u));d.bernoulliBminus=cast(bernoulli(-u));
        d.electronUsedQuasiFermi=d.holeUsedQuasiFermi=!edge["weights"].empty();
        d.psi0=cast(o.psi[i]);d.psi1=cast(o.psi[j]);d.phin0=cast(o.phin[i]);d.phin1=cast(o.phin[j]);d.phip0=cast(o.phip[i]);d.phip1=cast(o.phip[j]);
        d.electronQfReference0=model["nodes"][i][7];d.electronQfReference1=model["nodes"][j][7];
        d.electronSgPhin0Relative=cast(o.phin[i]-Wide(d.electronQfReference0));d.electronSgPhin1Relative=cast(o.phin[j]-Wide(d.electronQfReference0));
        d.electronSgPsi0Relative=cast(o.psi[i]-Wide(d.electronQfReference0));d.electronSgPsi1Relative=cast(o.psi[j]-Wide(d.electronQfReference0));
        d.n0=cast(o.n[i]);d.n1=cast(o.n[j]);d.p0=cast(o.p[i]);d.p1=cast(o.p[j]);d.ni0=model["nodes"][i][9];d.ni1=model["nodes"][j][9];
        d.mun=cast(mn);d.mup=cast(mp);
        const Wide drive=abs(o.psi[j]-o.psi[i])*Wide(model["field_factor"].get<double>())/Wide(edge["length"].get<double>());
        d.electronMobilityDriveInternal=d.holeMobilityDriveInternal=model.contains("hfs")?std::numeric_limits<Real>::quiet_NaN():cast(drive);
        const Wide couple=edge["couple"].get<double>();
        if(couple!=0){d.electronContinuityFlux=cast(o.electronFlux[k]*cs/couple);d.holeContinuityFlux=cast(o.holeFlux[k]*cs/couple);}
        d.electronCurrent=d.electronCurrentLongDoubleReference=cast(ec);d.holeCurrent=d.holeCurrentLongDoubleReference=cast(hc);
        d.electronDriftCurrent=cast(eDrift);d.electronDiffusionCurrent=cast(ec-eDrift);d.holeDriftCurrent=cast(hDrift);d.holeDiffusionCurrent=cast(hc-hDrift);d.totalCurrent=cast(ec-hc);result.edges.push_back(d);
    }
    result.totals.electronCurrent=cast(en);result.totals.holeCurrent=cast(hp);result.totals.totalCurrent=cast(en-hp);
    result.totals.electronDriftCurrent=cast(ed);result.totals.electronDiffusionCurrent=cast(ef);result.totals.holeDriftCurrent=cast(hd);result.totals.holeDiffusionCurrent=cast(hf);
    result.precision.electronCurrentCompensated=result.precision.electronCurrentLongDoubleReference=cast(en);
    result.precision.holeCurrentCompensated=result.precision.holeCurrentLongDoubleReference=cast(hp);
    result.precision.totalCurrentCompensated=result.precision.totalCurrentLongDoubleReference=cast(en-hp);
    return result;
}
void SplitDDRuntime::save(const VectorXd& x,DDSolution& out) const {
    const auto& o=impl_->evaluate(x);auto s=impl_->state(x);int N=s.size();
    out.packedState=x;out.packedLow=impl_->low;out.packedPotentialScale_V=s.potentialScale();out.packedMeshFingerprint=impl_->key;
    out.psi=project(o.psi);out.phin=project(o.phin);out.phip=project(o.phip);out.n=project(o.n);out.p=project(o.p);
    out.phinIncrement.resize(N);out.phipIncrement.resize(N);out.electronQfReference.resize(N);out.holeQfReference.resize(N);
    for(int i=0;i<N;++i){out.phinIncrement(i)=static_cast<Real>(s.coordinate(N+i)*Wide(s.potentialScale()));out.phipIncrement(i)=static_cast<Real>(s.coordinate(2*N+i)*Wide(s.potentialScale()));out.electronQfReference(i)=s.reference(1,i);out.holeQfReference(i)=s.reference(2,i);}
}
VectorXd SplitDDRuntime::restore(const DDSolution& s) {
    if(!s.hasConsistentSplitPackedState() || s.packedMeshFingerprint!=impl_->key)throw std::invalid_argument("Split checkpoint field or mesh mismatch");
    auto old=savedState(s);auto frame=impl_->state(s.packedState);
    if(old.potentialScale()!=frame.potentialScale())throw std::invalid_argument("Split checkpoint scale mismatch");
    for(std::size_t i=0;i<old.size();++i)for(int b=1;b<=2;++b)if(old.reference(b,i)!=frame.reference(b,i))throw std::invalid_argument("Split checkpoint QF reference mismatch");
    setLow(s.packedLow);return s.packedState;
}

VectorXd SplitDDRuntime::restoreForContinuation(const DDSolution& s) {
    if(!s.hasConsistentSplitPackedState() || s.packedMeshFingerprint!=impl_->key)
        throw std::invalid_argument("Split continuation field or mesh mismatch");
    const auto old=savedState(s),frame=impl_->state(s.packedState);
    if(old.potentialScale()!=frame.potentialScale())
        throw std::invalid_argument("Split continuation scale mismatch");
    VectorXd hi=s.packedState,lo=s.packedLow;
    const int N=old.size();
    for(int i=0;i<N;++i)for(int b=1;b<=2;++b) {
        if(old.reference(b,i)==frame.reference(b,i))continue;
        const int k=b*N+i;
        const Wide value=old.coordinate(k)+
            (Wide(old.reference(b,i))-Wide(frame.reference(b,i)))/Wide(old.potentialScale());
        hi(k)=static_cast<Real>(value);
        lo(k)=static_cast<Real>(value-Wide(hi(k)));
    }
    if(!hi.allFinite() || !lo.allFinite())throw std::invalid_argument("Nonfinite split continuation");
    setLow(lo);return hi;
}

void SplitDDRuntime::applyBoundary(VectorXd& x,const CoupledDDBoundaryConditions& bcs) {
    // Use precisely the same normalized constraints as residual().
    const auto frame=impl_->state(x);const int N=frame.size();auto lo=low();
    for(const auto& [i,z]:bcs.psi){x(i)=z;lo(i)=0.;}
    for(const auto& [i,z]:bcs.phin){x(N+i)=z-frame.reference(1,i)/frame.potentialScale();lo(N+i)=0.;}
    for(const auto& [i,z]:bcs.phip){x(2*N+i)=z-frame.reference(2,i)/frame.potentialScale();lo(2*N+i)=0.;}
    setLow(lo);
}

DDSolution predictSplitDDState(const DDSolution& previous,const DDSolution& current,
    Real ratio,const std::array<bool,3>& fields) {
    if(!previous.hasConsistentSplitPackedState() || !current.hasConsistentSplitPackedState() ||
       previous.packedMeshFingerprint!=current.packedMeshFingerprint ||
       previous.packedPotentialScale_V!=current.packedPotentialScale_V ||
       previous.psi.size()!=current.psi.size() || !std::isfinite(ratio))
        throw std::invalid_argument("Split predictor incompatible states");
    DDSolution out=current;
    if(ratio==0. || std::none_of(fields.begin(),fields.end(),[](bool b){return b;}))return out;
    const auto old=savedState(previous),now=savedState(current);const int N=now.size();
    if(current.n.size()!=N || current.p.size()!=N || !current.n.allFinite() || !current.p.allFinite() ||
       (current.n.array()<0.).any() || (current.p.array()<0.).any())
        throw std::invalid_argument("Split predictor invalid densities");
    for(int b=0;b<3;++b)if(fields[b])for(int i=0;i<N;++i) {
        const int k=b*N+i;
        const Wide value=now.coordinate(k)+Wide(ratio)*(now.potential(b,i)-old.potential(b,i))/Wide(now.potentialScale());
        out.packedState(k)=static_cast<Real>(value);
        out.packedLow(k)=static_cast<Real>(value-Wide(out.packedState(k)));
    }
    const auto predicted=savedState(out);
    for(int i=0;i<N;++i) {
        const Wide psi=predicted.potential(0,i),en=predicted.potential(1,i),hp=predicted.potential(2,i);
        out.psi(i)=static_cast<Real>(psi);out.phin(i)=static_cast<Real>(en);out.phip(i)=static_cast<Real>(hp);
        out.phinIncrement(i)=static_cast<Real>(predicted.coordinate(N+i)*Wide(now.potentialScale()));
        out.phipIncrement(i)=static_cast<Real>(predicted.coordinate(2*N+i)*Wide(now.potentialScale()));
        // The qualified split runtime is classical Boltzmann at 300 K.
        // These projected densities are initial-guess fields; the assembler
        // independently reconstructs density from its material and packed state.
        out.n(i)=static_cast<Real>(Wide(current.n(i))*exp((psi-en-now.potential(0,i)+now.potential(1,i))/Wide(constants::Vt_300)));
        out.p(i)=static_cast<Real>(Wide(current.p(i))*exp((hp-psi-now.potential(2,i)+now.potential(0,i))/Wide(constants::Vt_300)));
    }
    if(!out.hasConsistentSplitPackedState() || !out.n.allFinite() || !out.p.allFinite())
        throw std::invalid_argument("Nonfinite split predictor");
    out.converged=false;out.iters=0;return out;
}

void CoupledDDAssembler::enableSplitDDState(bool enabled) {splitDDEnabled_=enabled;splitRuntime_.reset();}
std::shared_ptr<SplitDDRuntime> CoupledDDAssembler::splitRuntime() const {
    if(!splitDDEnabled_)return {};
    if(!splitRuntime_)splitRuntime_=buildSplitRuntime();return splitRuntime_;
}
std::shared_ptr<SplitDDRuntime> CoupledDDAssembler::buildSplitRuntime() const {
    const bool enormal=usesElementDistanceLombardi(mobilityConfig_);
    if((mobilityConfig_.model!="phumob" && !enormal) || mobilityConfig_.edgeAveraging!="element_box_phumob" || usesFermiDirac_ || !recombination_.srhEnabled() || recombination_.augerEnabled() || recombination_.bandToBandEnabled() || impactIonizationCoupled_ || electronQuantumPotentialConfig_.enabled || (electronQuantumPotential_V_.size()&&!electronQuantumPotential_V_.isZero(0.)))throw std::invalid_argument("Split DD requires qualified element-box PhuMob, Boltzmann and SRH without quantum/avalanche");
    const int N=mesh_.numNodes();const double scale=scaling_.enabled?scaling_.V0:1.;
    const Real cs=scaling_.enabled?scaling_.C0*scaling_.D0:1.;const Real fieldFactor=scaling_.enabled?scaling_.fieldFromCoordinateDeltaFactor:1.;
    const Real sourceFactor=scaling_.enabled?scaling_.unitSystem.continuitySourceIntegralFactor():1.;const Real area=scaling_.enabled?scaling_.chargeAreaFactor:1.;
    nlohmann::json j={{"schema","vela.split-dd-operator.v1"},{"temperature_K",300.},{"potential_scale",scale},{"gauge_policy","encoded_reference"}};
    nlohmann::json geometry=nlohmann::json::array();
    for(int i=0;i<N;++i){const auto& node=mesh_.getNode(i);geometry.push_back({node.x,node.y});
        j["nodes"].push_back({0.,0.,0.,0.,0.,0.,scale,electronQuasiFermiReferenceAt(i),holeQuasiFermiReferenceAt(i),ni_[i],Vt_,doping_.netDoping(i),poissonElectronVol_[i],poissonHoleVol_[i],poissonDopantVol_[i],constants::q,area,scaling_.enabled?scaling_.permittivityReference_F_per_m*scaling_.V0:1.,fixedInterfaceChargeRhs_(i)});
        Real doping=recombination_.srhDopingConcentration(doping_.donors(i),doping_.acceptors(i));j["node_physics"].push_back({doping_.donors(i),doping_.acceptors(i),recombination_.electronLifetime(doping),recombination_.holeLifetime(doping),srhVol_[i]*sourceFactor/cs});}
    for(Index c=0;c<mesh_.numCells();++c){const auto& cell=mesh_.getCell(c);geometry.push_back({cell.node_ids[0],cell.node_ids[1],cell.node_ids[2],cell.region_id});}
    for(const auto& contact:mesh_.contacts())geometry.push_back({contact.name,contact.node_ids});
    boost::uuids::detail::sha1 hasher;auto serialized=geometry.dump();hasher.process_bytes(serialized.data(),serialized.size());boost::uuids::detail::sha1::digest_type digest;hasher.get_digest(digest);std::ostringstream id;id<<"sha1:";for(auto z:digest)id<<std::hex<<std::setfill('0')<<std::setw(2)<<unsigned(z);j["mesh_fingerprint"]=id.str();
    const auto& pm=mobilityConfig_.phuMob;const auto& en=pm.donorSpecies==PhuMobDonorSpecies::Arsenic?pm.electronArsenic:pm.electronPhosphorus;const auto& hp=pm.holeBoron;
    j["phumob"]={{"electron",{en.muMax,en.muMin,en.theta,en.nRef,en.alpha}},{"hole",{hp.muMax,hp.muMin,hp.theta,hp.nRef,hp.alpha}},
      {"common",{pm.donorClusterReference,pm.acceptorClusterReference,pm.donorClusterCoefficient,pm.acceptorClusterCoefficient,pm.electronMassRatio,pm.holeMassRatio,pm.conwellWeisskopfFactor,pm.brooksHerringFactor,pm.electronHoleScatteringFactor,pm.holeElectronScatteringFactor,pm.gA,pm.gB,pm.gC,pm.gAlpha,pm.gAlphaPrime,pm.gBeta,pm.gGamma,pm.internalConcentrationToCm3,pm.internalMobilityToCm2PerVS}}};
    j["boundary"]=nlohmann::json::array({nlohmann::json::array(),nlohmann::json::array(),nlohmann::json::array()});
    std::vector<int> surfaceCell(mesh_.numCells(),-1);
    if(enormal) {
        if(mobilityConfig_.electronLombardi.alpha!=0 || mobilityConfig_.holeLombardi.alpha!=0)
            throw std::invalid_argument("Split element Lombardi requires alpha=0");
        detail::updateSurfaceMobilityCellGeometry(mobilityConfig_,mesh_,edgeCells_,VectorXd::Zero(N),mobilityConfig_.surface.coordinateFieldFactor,&cellMaterials_);
        if(usesElementHighField(mobilityConfig_)) {
            j["hfs"]={{"cutoff",100./mobilityConfig_.internalFieldToVPerM},
                {"parameters",{{mobilityConfig_.electronField.saturationVelocity,mobilityConfig_.electronField.beta},{mobilityConfig_.holeField.saturationVelocity,mobilityConfig_.holeField.beta}}}};
        }
        auto& s=j["enormal"];s["cells"]=nlohmann::json::array();s["parameters"]=nlohmann::json::array();
        s["mobility_to_SI"]=mobilityConfig_.internalMobilityToM2PerVS;s["concentration_to_SI"]=mobilityConfig_.internalConcentrationToM3;
        for(const auto& p:{mobilityConfig_.electronLombardi,mobilityConfig_.holeLombardi})s["parameters"].push_back({p.B,p.C,p.N0,p.N2,p.lambda,p.k,p.delta,p.A,p.eta,p.criticalLength,p.acousticFactor,p.roughnessFactor});
        for(Index cid=0;cid<mesh_.numCells();++cid) {
            const auto& mat=cellMaterials_[cid];if(mat.ni<=0&&mat.mun<=0&&mat.mup<=0)continue;
            if(mat.temperature_K.value_or(300.)!=300.||mat.mun<=0||mat.mup<=0)throw std::invalid_argument("Split Lombardi requires 300 K silicon transport");
            const auto& cell=mesh_.getCell(cid);const auto& p0=mesh_.getNode(cell.node_ids[0]);const auto& p1=mesh_.getNode(cell.node_ids[1]);const auto& p2=mesh_.getNode(cell.node_ids[2]);
            Real x1=p1.x-p0.x,y1=p1.y-p0.y,x2=p2.x-p0.x,y2=p2.y-p0.y,det=x1*y2-x2*y1;
            Real nx=mobilityConfig_.surface.cellNormalX[cid],ny=mobilityConfig_.surface.cellNormalY[cid];
            Real f=mobilityConfig_.surface.coordinateFieldFactor*mobilityConfig_.internalFieldToVPerM;
            const auto measures=detail::cellBoxNodeMeasures(mesh_,cid);long double volume=(long double)measures[0]+measures[1]+measures[2];
            nlohmann::json distances=nlohmann::json::array(),weights=nlohmann::json::array();
            for(int k=0;k<3;++k){distances.push_back(mobilityConfig_.surface.nodeDistances[cell.node_ids[k]]*mobilityConfig_.internalLengthToM);weights.push_back(static_cast<Real>((long double)measures[k]/volume));}
            surfaceCell[cid]=static_cast<int>(s["cells"].size());
            s["cells"].push_back({{"nodes",cell.node_ids},{"distances_m",distances},{"weights",weights},{"stencil",{(y2*nx-x2*ny)/det*f,(-y1*nx+x1*ny)/det*f}}});
            if(usesElementHighField(mobilityConfig_)) {
                const auto& h=mobilityConfig_.highFieldCells.at(cid);
                s["cells"].back()["hfs_stencil"]=h.stencil;s["cells"].back()["hfs_contact"]=h.contact;
            }
        }
    }
    j["continuity_scale"]=cs;j["field_factor"]=fieldFactor;j["current_factor"]=constants::q*cs*(scaling_.enabled?scaling_.currentDensityLineIntegralFactor:1.)*1e-6;
    for(Index e=0;e<mesh_.numEdges();++e){const auto& edge=mesh_.getEdge(e);if(edge.length<1e-30)continue;j["poisson_edges"].push_back({edge.n0,edge.n1,edgeAssemblyKernels_[e].poissonCoupling});
        nlohmann::json weights=nlohmann::json::array(),ports=nlohmann::json::object();long double total=0;
        for(Index cid:edgeCells_[e]){const auto& m=cellMaterials_[cid];if(m.ni<=0&&m.mun<=0&&m.mup<=0)continue;total+=detail::cellBoxEdgeCoefficient(mesh_,cid,e);}
        if(total>0)for(Index cid:edgeCells_[e]){const auto& m=cellMaterials_[cid];if(m.ni<=0&&m.mun<=0&&m.mup<=0)continue;if(m.temperature_K.value_or(300.)!=300.||m.mun<=0||m.mup<=0)throw std::invalid_argument("Split DD requires 300 K silicon transport");Real g=detail::cellBoxEdgeCoefficient(mesh_,cid,e);if(g==0)continue;const auto measures=detail::cellBoxNodeMeasures(mesh_,cid);long double volume=(long double)measures[0]+measures[1]+measures[2];for(int k=0;k<3;++k)weights.push_back({mesh_.getCell(cid).node_ids[k],static_cast<Real>((long double)g*measures[k]/volume/total)});}
        if(couple_[e]>0)for(const auto& c:mesh_.contacts()){bool i=std::find(c.node_ids.begin(),c.node_ids.end(),edge.n0)!=c.node_ids.end();bool h=std::find(c.node_ids.begin(),c.node_ids.end(),edge.n1)!=c.node_ids.end();if(i!=h)ports[c.name]=i?1.:-1.;}
        j["transport_edges"].push_back({{"id",e},{"i",edge.n0},{"j",edge.n1},{"coefficient",Vt_*fieldFactor*couple_[e]/edge.length/cs},{"weights",weights},{"ports",ports},{"length",edge.length},{"couple",couple_[e]},
            {"length_m",scaling_.unitSystem.internalLengthToMeters(edge.length)},{"couple_m",scaling_.unitSystem.internalLengthToMeters(couple_[e])}});
        if(enormal) {
            auto& ws=j["transport_edges"].back()["cell_weights"];ws=nlohmann::json::array();
            if(total>0)for(Index cid:edgeCells_[e])if(surfaceCell[cid]>=0) {
                const Real g=detail::cellBoxEdgeCoefficient(mesh_,cid,e);if(g!=0)ws.push_back({surfaceCell[cid],static_cast<Real>((long double)g/total)});
            }
        }
    }
    return std::make_shared<SplitDDRuntime>(std::move(j));
}
} // namespace vela
