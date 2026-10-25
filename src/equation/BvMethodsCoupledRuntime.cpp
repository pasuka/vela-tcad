#include "vela/equation/BvMethodsCoupledRuntime.h"
#include "vela/equation/BvMethodsTri3Operator.h"
#include <algorithm>
#include <stdexcept>

namespace vela {
struct BvMethodsCoupledRuntime::Impl {
    const CoupledDDAssembler& owner;
    std::unique_ptr<CoupledDDAssembler> poisson;
    detail::BvMethodsTri3Context context;
    std::unique_ptr<detail::BvMethodsTri3Operator> element;
    std::vector<double> volumes;
    std::vector<bool> active;
    double scale,continuity;
    struct Evaluation {
        detail::BvMethodsTri3Evaluation transport,av;
        VectorXd srh;
        SparseMatrixd srhJacobian;
        detail::BvMethodsQfState state;
    };
    explicit Impl(const CoupledDDAssembler& a)
    :owner(a),context{a.mesh_,a.matdb_,a.doping_,{a.mobilityConfig_,a.impactIonizationConfig_,{},a.carrierStatisticsModel_}},
     scale(a.scaling_.enabled?a.scaling_.V0:1.),continuity(a.scaling_.enabled?a.scaling_.C0*a.scaling_.D0:1.) {
        const auto fail=[](const char* s){throw std::invalid_argument(std::string("bvmethods_tri3_qf: ")+s);};
        if(a.carrierStatisticsModel_!=CarrierStatisticsModel::FermiDiracJoyceDixon)
            fail("requires explicit fermi_dirac_joyce_dixon statistics");
        if(std::abs(a.Vt_*constants::q/constants::kb-300.)>1e-8 ||
           a.electronQuantumPotentialConfig_.enabled || a.carrierDiagonalFloor_.enabled)
            fail("requires 300 K classical DD without carrier diagonal floors");
        if(a.mobilityConfig_.model!="masetti_field_lombardi")
            fail("requires masetti_field_lombardi mobility parameters");
        if(a.mobilityConfig_.highFieldDrivingForce!="quasi_fermi_gradient" || !a.mobilityConfig_.jacobianFieldDerivatives)
            fail("requires quasi_fermi_gradient mobility and live field derivatives");
        if(a.mobilityConfig_.dopingConcentrationBasis!="total_impurity")
            fail("requires total_impurity mobility concentration basis");
        if(a.recombination_.augerEnabled() || a.recombination_.bandToBandEnabled())
            fail("only SRH recombination is qualified");
        if(a.mesh_.poissonPermittivityPolicy()!=BoxGeometryBuilder::PoissonPermittivityPolicy::CellMaterial)
            fail("requires cell_material Poisson geometry");
        const auto& av=a.impactIonizationConfig_;
        if(av.model!="none" && (av.model!="van_overstraeten" ||
           av.couplingMode!="self_consistent" || av.sourceJacobianMode!="local_ad"))
            fail("requires self-consistent Van Overstraeten with local_ad");
        if(av.model!="none" && (av.drivingForce!="quasi_fermi_gradient" || av.generation!="current_density" ||
           av.minimumField!=0. || av.quasiFermiCarrierTruncation!=0. || av.electronDrivingForceRefDensity!=0. ||
           av.holeDrivingForceRefDensity!=0. || av.sourceGeometryScale!=1. || av.sourceVolumeFactor!=0. ||
           av.aScale!=1. || av.bScale!=1.))
            fail("profile owns GradQF current-source support; cutoffs, interpolation and source overrides unqualified");
        for(Index cid=0;cid<a.mesh_.numCells();++cid)
            if(a.mesh_.getCell(cid).node_ids.size()!=3)fail("requires Tri3 cells");
        auto legacy=a.scaling_;legacy.coupledOperator="legacy";
        auto noAv=av;noAv.model="none";
        poisson=std::make_unique<CoupledDDAssembler>(a.mesh_,a.matdb_,a.doping_,a.Vt_,
            a.mobilityConfig_,a.recombinationConfig_,a.bandgapNarrowingConfig_,noAv,
            std::vector<RegionFixedChargeSpec>{},std::vector<InterfaceSheetChargeSpec>{},legacy,
            a.carrierDiagonalFloor_,a.carrierStatistics_,a.electronQuantumPotentialConfig_);
        // The owner already integrated fixed/interface charge with its unit policy.
        poisson->fixedInterfaceChargeRhs_=a.fixedInterfaceChargeRhs_;
        context.newton.inputScaling.physicalUnitSystem=a.scaling_.enabled
            ?a.scaling_.unitSystem:PhysicalUnitSystem::legacySI();
        const int n=a.mesh_.numNodes();volumes.resize(n,0.);active.resize(n,false);
        for(Index cid=0;cid<a.mesh_.numCells();++cid)if(a.cellMaterials_[cid].ni>0.) {
            const auto w=detail::cellBoxNodeMeasures(a.mesh_,cid);const auto& cell=a.mesh_.getCell(cid);
            for(int k=0;k<3;++k){volumes[cell.node_ids[k]]+=w[k];active[cell.node_ids[k]]=true;}
        }
        auto state=makeState(VectorXd::Zero(3*n));
        const auto& u=context.newton.inputScaling.unitSystem();
        const double line=u.currentDensityAM2PerInternal()*u.lengthMPerInternal()*continuity;
        element=std::make_unique<detail::BvMethodsTri3Operator>(context,state,a.Nc_,a.Nv_,scale,line,a.Vt_);
    }
    void sync()const {
        if(!owner.electronQuasiFermiReferenceField().isZero(0.) || !owner.holeQuasiFermiReferenceField().isZero(0.))
            throw std::invalid_argument("bvmethods_tri3_qf v1 requires zero-reference ordinary QF coordinates");
        poisson->setQuasiFermiReferenceFields(owner.electronQuasiFermiReferenceField(),owner.holeQuasiFermiReferenceField());
    }
    detail::BvMethodsQfState makeState(const VectorXd& x) const {
        const int n=owner.mesh_.numNodes();if(x.size()!=3*n || !x.allFinite())
            throw std::invalid_argument("bvmethods_tri3_qf: invalid QF state");
        detail::BvMethodsQfState s;s.psi=x.head(n)*scale;
        s.phin=x.segment(n,n)*scale+owner.electronQuasiFermiReferenceField();
        s.phip=x.tail(n)*scale+owner.holeQuasiFermiReferenceField();
        s.n=owner.electronDensity(x);s.p=owner.holeDensity(x);
        for(int i=0;i<n;++i)if(active[i] && (!(s.n[i]/owner.Nc_[i]>0.) || !(s.p[i]/owner.Nv_[i]>0.) || !std::isfinite(s.n[i]) || !std::isfinite(s.p[i])))
            throw std::domain_error("BVMethods trial state outside positive finite carrier domain");
        return s;
    }
    Evaluation evaluate(const VectorXd& x,bool jac)const {
        Evaluation out;out.state=makeState(x);const auto& s=out.state;const int n=s.n.size();
        out.transport=element->transport(s,jac,true,true);
        if(owner.impactIonizationCoupled_)out.av=element->evaluate(s,jac);
        else {out.av.source=VectorXd::Zero(3*n);out.av.allGeneration=VectorXd::Zero(n);out.av.jacobian.resize(3*n,3*n);out.av.electronGeneration=VectorXd::Zero(n);out.av.holeGeneration=VectorXd::Zero(n);}
        out.srh=VectorXd::Zero(3*n);std::vector<Eigen::Triplet<double>> entries;
        const auto& u=context.newton.inputScaling.unitSystem();
        const double sourceFactor=u.continuitySourceIntegralFactor()/continuity;
        for(int i=0;i<n;++i)if(active[i]&&owner.recombination_.srhEnabled()) {
            using detail::BvMethodsQfState;using experimental::AvJet;
            const auto model=owner.carrierStatisticsModel_;const double vt=owner.Vt_;
            auto population=[&](bool electron){
                const double v=electron?s.n[i]:s.p[i],dos=electron?owner.Nc_[i]:owner.Nv_[i];
                AvJet a(v);const double eta=inverseFermiDiracHalf(v/dos,model),d=dos*fermiDiracHalfDerivative(eta,model)/vt;
                a.derivative[0]=electron?d:-d;a.derivative[electron?1:2]=electron?-d:d;return a;
            };
            const auto nj=population(true),pj=population(false);
            auto gamma=[&](const AvJet& p,double dos){
                const double eta=inverseFermiDiracHalf(p.value/dos,model),slope=dos*fermiDiracHalfDerivative(eta,model);
                AvJet g(std::exp(std::log(p.value/dos)-eta));
                for(int k=0;k<3;++k)g.derivative[k]=g.value*(1/p.value-1/slope)*p.derivative[k];return g;
            };
            const auto gn=gamma(nj,owner.Nc_[i]),gp=gamma(pj,owner.Nv_[i]);
            const auto splitting=(AvJet::variable(s.phip[i],2)-AvJet::variable(s.phin[i],1))/AvJet(vt);
            const auto excess=gn*gp*AvJet(owner.ni_[i]*owner.ni_[i])*detail::dualExpm1(splitting);
            const double dop=owner.recombination_.srhDopingConcentration(owner.doping_.donors(i),owner.doping_.acceptors(i));
            const double rate=owner.recombination_.srhRateGeneralizedFromExcessProduct(excess.value,nj.value,pj.value,
                owner.ni_[i],owner.ni_[i],gn.value,gp.value,dop);
            const double measure=volumes[i]*sourceFactor;out.srh[n+i]=out.srh[2*n+i]=rate*measure;
            if(jac){const auto d=owner.recombination_.srhRateGeneralizedDerivativesFromExcessProduct(excess.value,
                nj.value,pj.value,owner.ni_[i],owner.ni_[i],gn.value,gp.value,dop);
                for(int k=0;k<3;++k){const double v=(d.dRateDn*nj.derivative[k]+d.dRateDp*pj.derivative[k]+
                    d.dRateDExcess*excess.derivative[k]+d.dRateDElectronDegeneracy*gn.derivative[k]+
                    d.dRateDHoleDegeneracy*gp.derivative[k])*measure*scale;
                    for(int b:{1,2})if(v!=0)entries.emplace_back(b*n+i,k*n+i,v);
                }
            }
        }
        out.srhJacobian.resize(3*n,3*n);out.srhJacobian.setFromTriplets(entries.begin(),entries.end());return out;
    }
    void validateBoundary(const CoupledDDBoundaryConditions& b)const {
        if(!b.thermionic.empty())throw std::invalid_argument("bvmethods_tri3_qf: thermionic boundaries unqualified");
    }
    bool constrained(int row,const CoupledDDBoundaryConditions& b)const {
        const int n=owner.mesh_.numNodes();return row<n?b.psi.contains(row):row<2*n?b.phin.contains(row-n):b.phip.contains(row-2*n);
    }
    SparseMatrixd project(const SparseMatrixd& M,const CoupledDDBoundaryConditions& b)const {
        std::vector<Eigen::Triplet<double>> t;
        for(int k=0;k<M.outerSize();++k)for(SparseMatrixd::InnerIterator v(M,k);v;++v)
            if(!constrained(v.row(),b))t.emplace_back(v.row(),v.col(),v.value());
        SparseMatrixd a(M.rows(),M.cols());a.setFromTriplets(t.begin(),t.end());return a;
    }
};
BvMethodsCoupledRuntime::BvMethodsCoupledRuntime(const CoupledDDAssembler& a):impl_(std::make_unique<Impl>(a)){}
BvMethodsCoupledRuntime::~BvMethodsCoupledRuntime()=default;
VectorXd BvMethodsCoupledRuntime::residual(const VectorXd& x,const CoupledDDBoundaryConditions& b)const {
    auto& p=*impl_;p.validateBoundary(b);p.sync();
    Impl::Evaluation e;try{e=p.evaluate(x,false);}catch(const std::domain_error&){return VectorXd::Constant(x.size(),std::numeric_limits<double>::quiet_NaN());}
    const int n=p.active.size();
    VectorXd r=p.poisson->residual(x,b);const auto carrier=e.transport.source+e.av.source+e.srh;
    for(int i=0;i<n;++i)for(int block:{1,2})if(p.active[i]&&!p.constrained(block*n+i,b))r[block*n+i]=carrier[block*n+i];
    return r;
}
SparseMatrixd BvMethodsCoupledRuntime::jacobian(const VectorXd& x,const CoupledDDBoundaryConditions& b)const {
    auto& p=*impl_;p.validateBoundary(b);p.sync();auto e=p.evaluate(x,true);const int n=p.active.size();
    const auto base=p.poisson->assembleJacobian(x,b);const SparseMatrixd carrier=e.transport.jacobian+e.av.jacobian+e.srhJacobian;
    std::vector<Eigen::Triplet<double>> t;
    for(int k=0;k<base.outerSize();++k)for(SparseMatrixd::InnerIterator v(base,k);v;++v)
        if(v.row()<n || !p.active[v.row()%n] || p.constrained(v.row(),b))t.emplace_back(v.row(),v.col(),v.value());
    for(int k=0;k<carrier.outerSize();++k)for(SparseMatrixd::InnerIterator v(carrier,k);v;++v)
        if(!p.constrained(v.row(),b))t.emplace_back(v.row(),v.col(),v.value());
    SparseMatrixd result(3*n,3*n);result.setFromTriplets(t.begin(),t.end());return result;
}
std::vector<CoupledDDCarrierTermDiagnostic> BvMethodsCoupledRuntime::terms(const VectorXd& x,const CoupledDDBoundaryConditions& b,bool includeImpact)const {
    auto& p=*impl_;p.validateBoundary(b);auto e=p.evaluate(x,false);const int n=p.active.size();
    auto r=residual(x,b);const auto& u=p.context.newton.inputScaling.unitSystem();std::vector<CoupledDDCarrierTermDiagnostic> rows(n);
    for(int i=0;i<n;++i){auto& a=rows[i];a.nodeId=i;a.electronDensity_m3=u.internalConcentrationToM3(e.state.n[i]);a.holeDensity_m3=u.internalConcentrationToM3(e.state.p[i]);
        a.electronFlux=e.transport.source[n+i];a.holeFlux=e.transport.source[2*n+i];a.electronFluxAbsSum=e.transport.absoluteSource[n+i];a.holeFluxAbsSum=e.transport.absoluteSource[2*n+i];
        a.electronRecombination=a.holeRecombination=e.srh[n+i];a.electronImpact=a.holeImpact=includeImpact?e.av.source[n+i]:0.;
        a.impactElectronSource=includeImpact?e.av.electronGeneration[i]:0.;a.impactHoleSource=includeImpact?e.av.holeGeneration[i]:0.;
        a.impactCombinedSource=-a.electronImpact;a.electronResidual=r[n+i];a.holeResidual=r[2*n+i];
        a.electronContinuityActive=p.active[i]&&!b.phin.contains(i);a.holeContinuityActive=p.active[i]&&!b.phip.contains(i);
        if(!includeImpact){if(a.electronContinuityActive)a.electronResidual-=e.av.source[n+i];if(a.holeContinuityActive)a.holeResidual-=e.av.source[2*n+i];}
        if(!p.active[i]){a.electronGauge=r[n+i];a.holeGauge=r[2*n+i];}
        if(b.phin.contains(i))a.electronBoundary=r[n+i]-a.electronFlux-a.electronRecombination-a.electronImpact-a.electronGauge;
        if(b.phip.contains(i))a.holeBoundary=r[2*n+i]-a.holeFlux-a.holeRecombination-a.holeImpact-a.holeGauge;
    }return rows;
}
VectorXd BvMethodsCoupledRuntime::avalancheResidual(const VectorXd& x,const CoupledDDBoundaryConditions& b)const {
    auto& p=*impl_;auto e=p.evaluate(x,false);auto r=e.av.source;
    for(int i=0;i<r.size();++i)if(p.constrained(i,b))r[i]=0;return r;
}
SparseMatrixd BvMethodsCoupledRuntime::avalancheJacobian(const VectorXd& x,const CoupledDDBoundaryConditions& b)const {
    return impl_->project(impl_->evaluate(x,true).av.jacobian,b);
}
VectorXd BvMethodsCoupledRuntime::transportResidual(const VectorXd& x)const {
    auto& p=*impl_;p.sync();
    try{return p.element->transport(p.makeState(x),false,true,true).source;}
    catch(const std::domain_error&){return VectorXd::Constant(x.size(),std::numeric_limits<double>::quiet_NaN());}
}
SparseMatrixd BvMethodsCoupledRuntime::transportJacobian(const VectorXd& x)const {
    auto& p=*impl_;p.sync();return p.element->transport(p.makeState(x),true,true,true).jacobian;
}
VectorXd CoupledDDAssembler::bvMethodsTransportResidual(const VectorXd& x)const {
    auto r=bvMethodsRuntime();if(!r)throw std::invalid_argument("BVMethods transport requested for legacy operator");return r->transportResidual(x);
}
SparseMatrixd CoupledDDAssembler::bvMethodsTransportJacobian(const VectorXd& x)const {
    auto r=bvMethodsRuntime();if(!r)throw std::invalid_argument("BVMethods transport requested for legacy operator");return r->transportJacobian(x);
}
std::shared_ptr<BvMethodsCoupledRuntime> CoupledDDAssembler::bvMethodsRuntime()const {
    if(scaling_.coupledOperator=="legacy")return {};
    if(scaling_.coupledOperator!="bvmethods_tri3_qf")throw std::invalid_argument("Unknown coupled operator");
    if(!bvMethodsRuntime_)bvMethodsRuntime_=std::make_shared<BvMethodsCoupledRuntime>(*this);return bvMethodsRuntime_;
}
}
