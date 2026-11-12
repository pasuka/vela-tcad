#include "vela/equation/BvMethodsCoupledRuntime.h"
#include "vela/core/PerformanceProfiler.h"
#include "vela/equation/BvMethodsTri3Operator.h"
#include "vela/physics/experimental/PopulationRecombination.h"
#include <algorithm>
#include <cstdlib>
#include <string_view>
#include <stdexcept>
#include <cstring>
#include <optional>
#include <tuple>

namespace vela {
namespace {
bool enabledEnvironmentFlag(const char* name) {
    const char* value=std::getenv(name);
    return value && std::string_view(value)=="1";
}

enum class TransportValueMode { dual, shadow, scalar };
TransportValueMode transportValueMode() {
    static const auto mode=[] {
        const char* value=std::getenv("VELA_BVMETHODS_TRI3_VALUE_MODE");
        if(value && std::string_view(value)=="shadow")return TransportValueMode::shadow;
        if(value && std::string_view(value)=="scalar")return TransportValueMode::scalar;
        return TransportValueMode::dual;
    }();
    return mode;
}

bool equalStateBits(const VectorXd& a,const VectorXd& b) {
    return a.size()==b.size() && (a.size()==0 || std::memcmp(a.data(),b.data(),sizeof(Real)*a.size())==0);
}

std::string_view portTransportStage(BvMethodsTransportPurpose purpose) {
    switch(purpose) {
    case BvMethodsTransportPurpose::terminalValue:return "bvmethods.tri3.port_value";
    case BvMethodsTransportPurpose::terminalLinearization:return "bvmethods.tri3.port_linearization_value";
    case BvMethodsTransportPurpose::terminalValidation:return "bvmethods.tri3.port_validation";
    case BvMethodsTransportPurpose::contactOutput:return "bvmethods.tri3.contact_output";
    case BvMethodsTransportPurpose::other:return "bvmethods.tri3.transport_other";
    }
    return "bvmethods.tri3.transport_other";
}
} // namespace
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
    template<class T> struct Entry {
        VectorXd state;
        std::optional<T> value;
        bool matches(const VectorXd& x) const {
            return value.has_value() && state.size()==x.size() &&
                (x.size()==0 || std::memcmp(state.data(),x.data(),sizeof(Real)*x.size())==0);
        }
    };
    mutable std::array<Entry<Evaluation>,2> evaluations;
    mutable std::array<Entry<detail::BvMethodsTri3Evaluation>,2> transports;
    mutable VectorXd electronReference,holeReference;
    // A single runtime owns immutable geometry, coefficients and options.
    // Audit one previous source-free value state without changing dispatch.
    mutable VectorXd lastPortValueState;
    mutable BvMethodsTransportPurpose lastPortPurpose=BvMethodsTransportPurpose::other;
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
        if((a.scaling_.populationSgCurrent!="quasi_fermi" && a.scaling_.populationSgCurrent!="band_population") ||
           (a.scaling_.populationSgCurrent!="quasi_fermi"&&!a.usesPopulationState()))fail("invalid population SG arithmetic");
        if(!std::isfinite(a.scaling_.populationSgBootstrapFloor_cm3) || a.scaling_.populationSgBootstrapFloor_cm3<0. ||
           (a.scaling_.populationSgBootstrapFloor_cm3>0. &&
            (!a.usesPopulationState() || a.scaling_.populationSgCurrent!="quasi_fermi")))
            fail("statistical bootstrap requires population quasi_fermi SG and a finite nonnegative floor");
        auto legacy=a.scaling_;legacy.coupledOperator="legacy";legacy.stateVariables="potential_quasi_fermi";legacy.populationSgCurrent="quasi_fermi";legacy.populationSgBootstrapFloor_cm3=0.;
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
        VectorXd seed=VectorXd::Zero(3*n);
        if(a.usesPopulationState())for(int i=0;i<n;++i)if(active[i])seed[n+i]=seed[2*n+i]=a.ni_[i];
        auto state=makeState(seed);
        const auto& u=context.newton.inputScaling.unitSystem();
        const double line=u.currentDensityAM2PerInternal()*u.lengthMPerInternal()*continuity;
        element=std::make_unique<detail::BvMethodsTri3Operator>(context,state,a.Nc_,a.Nv_,scale,line,a.Vt_,a.ni_,a.scaling_.populationSgCurrent=="band_population");
    }
    void sync()const {
        const auto en=owner.electronQuasiFermiReferenceField(),hp=owner.holeQuasiFermiReferenceField();
        if(!en.isZero(0.) || !hp.isZero(0.))
            throw std::invalid_argument("bvmethods_tri3_qf v1 requires zero-reference ordinary QF coordinates");
        if(owner.bvMethodsEvaluationReuse_) {
            const auto same=[](const VectorXd& a,const VectorXd& b){return a.size()==b.size() &&
                (a.size()==0 || std::memcmp(a.data(),b.data(),sizeof(Real)*a.size())==0);};
            if(!same(en,electronReference)||!same(hp,holeReference)) {
                for(auto& e:evaluations)e.value.reset();for(auto& t:transports)t.value.reset();
                electronReference=en;holeReference=hp;
            }
        }
        poisson->setQuasiFermiReferenceFields(en,hp);
    }
    detail::BvMethodsQfState makeState(const VectorXd& x) const {
        const int n=owner.mesh_.numNodes();if(x.size()!=3*n || !x.allFinite())
            throw std::invalid_argument("bvmethods_tri3_qf: invalid QF state");
        detail::BvMethodsQfState s;s.psi=x.head(n)*scale;
        if(owner.usesPopulationState()) {
            const auto state=owner.populationPhysicalFields(x);
            s.phin=state.phin_V;s.phip=state.phip_V;
            s.n=state.n;s.p=state.p;s.populationCoordinates=true;
            return s;
        }
        s.phin=x.segment(n,n)*scale+owner.electronQuasiFermiReferenceField();
        s.phip=x.tail(n)*scale+owner.holeQuasiFermiReferenceField();
        s.n=owner.electronDensity(x);s.p=owner.holeDensity(x);
        for(int i=0;i<n;++i)if(active[i] && (!(s.n[i]>0.) || !(s.p[i]>0.) || !std::isfinite(s.n[i]) || !std::isfinite(s.p[i])))
            throw std::domain_error("BVMethods trial state outside positive finite carrier domain");
        return s;
    }
    detail::BvMethodsTri3Evaluation assembleTransport(const detail::BvMethodsQfState& state,
        bool jac,bool approximate,const VectorXd* columns,double floor,bool reusePowers,
        std::string_view stage) const {
        const auto mode=jac || !element->supportsTransportValues(state)
            ? TransportValueMode::dual : transportValueMode();
        if(mode==TransportValueMode::scalar) {
            ScopedPerformanceTimer timer(stage);
            return element->transportValues(state,true,true,reusePowers);
        }
        auto reference=[&] {
            ScopedPerformanceTimer timer(stage);
            return element->transport(state,jac,true,true,approximate,columns,floor,reusePowers);
        }();
        if(mode==TransportValueMode::shadow) {
            ScopedPerformanceTimer timer("bvmethods.tri3.scalar_shadow");
            const auto candidate=element->transportValues(state,true,true,reusePowers);
            if(!equalStateBits(reference.source,candidate.source) ||
               !equalStateBits(reference.absoluteSource,candidate.absoluteSource)) {
                incrementPerformanceCounter("bvmethods.tri3.scalar_shadow_mismatches");
                throw std::runtime_error("BVMethods scalar transport shadow mismatch");
            }
        }
        return reference;
    }

    detail::BvMethodsTri3Evaluation transport(const VectorXd& x,bool jac,
        BvMethodsTransportPurpose purpose=BvMethodsTransportPurpose::other) const {
        auto& cached=transports[jac?1:0];
        if(owner.bvMethodsEvaluationReuse_&&cached.matches(x))return *cached.value;
        static const bool audit=enabledEnvironmentFlag("VELA_BVMETHODS_TRI3_CONTACT_AUDIT");
        std::optional<ScopedPerformanceTimer> repeatTimer;
        if(audit && activePerformanceProfiler() && !jac && purpose!=BvMethodsTransportPurpose::other) {
            incrementPerformanceCounter("bvmethods.tri3.port_audited_values");
            if(lastPortPurpose!=BvMethodsTransportPurpose::other && equalStateBits(lastPortValueState,x)) {
                incrementPerformanceCounter("bvmethods.tri3.port_same_state");
                if(lastPortPurpose!=purpose)incrementPerformanceCounter("bvmethods.tri3.port_same_state_cross_purpose");
                repeatTimer.emplace("bvmethods.tri3.port_same_state_recompute");
            } else lastPortValueState=x;
            lastPortPurpose=purpose;
        }
        auto result=assembleTransport(makeState(x),jac,false,nullptr,0.,false,
            jac?"bvmethods.tri3.terminal_jacobian":portTransportStage(purpose));
        if(owner.bvMethodsEvaluationReuse_){cached.state=x;cached.value=result;}
        return result;
    }
    Evaluation evaluate(const VectorXd& x,bool jac)const {
        if(!owner.bvMethodsEvaluationReuse_)return compute(x,jac);
        sync();auto& cached=evaluations[jac?1:0];
        if(owner.bvMethodsEvaluationReuse_&&cached.matches(x))return *cached.value;
        auto result=compute(x,jac);
        if(owner.bvMethodsEvaluationReuse_){cached.state=x;cached.value=result;}
        return result;
    }
    Evaluation compute(const VectorXd& x,bool jac,bool approximate=false,const VectorXd* columns=nullptr,bool freezeStatisticalFloor=true)const {
        Evaluation out;out.state=makeState(x);const auto& s=out.state;const int n=s.n.size();
        static const bool enableLombardiPowReuse = [] {
            const char* flag = std::getenv("VELA_BVMETHODS_DISABLE_LOMBARDI_POW_REUSE");
            return !flag || std::string_view(flag) != "1";
        }();
        const bool reuseLombardiPowers = enableLombardiPowReuse && !jac &&
            s.populationCoordinates && owner.scaling_.populationSgCurrent == "quasi_fermi";
        out.transport=assembleTransport(s,jac,approximate,columns,
            approximate&&freezeStatisticalFloor?owner.scaling_.populationSgBootstrapFloor_cm3:0.,
            reuseLombardiPowers,jac?"bvmethods.tri3.device_jacobian":"bvmethods.tri3.device_residual");
        if(owner.impactIonizationCoupled_) {
            static const bool reuseGeneration=!enabledEnvironmentFlag("VELA_BVMETHODS_DISABLE_AV_GENERATION_REUSE");
            ScopedPerformanceTimer timer(jac&&!approximate?"bvmethods.tri3.avalanche_jacobian":"bvmethods.tri3.avalanche_residual");
            out.av=element->evaluate(s,jac&&!approximate,columns,reuseGeneration);
        }
        else {out.av.source=VectorXd::Zero(3*n);out.av.allGeneration=VectorXd::Zero(n);out.av.jacobian.resize(3*n,3*n);out.av.electronGeneration=VectorXd::Zero(n);out.av.holeGeneration=VectorXd::Zero(n);}
        out.srh=VectorXd::Zero(3*n);std::vector<Eigen::Triplet<double>> entries;
        const auto& u=context.newton.inputScaling.unitSystem();
        const double sourceFactor=u.continuitySourceIntegralFactor()/continuity;
        for(int i=0;i<n;++i)if(active[i]&&owner.recombination_.srhEnabled()) {
            using detail::BvMethodsQfState;using experimental::AvJet;
            const double vt=owner.Vt_;
            const double tangentN=columns?(*columns)[n+i]:(owner.usesPopulationState() && s.n[i]<std::sqrt(std::numeric_limits<double>::min())?s.n[i]:1.);
            const double tangentP=columns?(*columns)[2*n+i]:(owner.usesPopulationState() && s.p[i]<std::sqrt(std::numeric_limits<double>::min())?s.p[i]:1.);
            auto population=[&](bool electron){
                const double v=electron?s.n[i]:s.p[i],dos=electron?owner.Nc_[i]:owner.Nv_[i];
                AvJet a(v);const double d=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(v,dos)/vt;
                if(owner.usesPopulationState())a.derivative[electron?1:2]=(electron?tangentN:tangentP)/scale;
                else {a.derivative[0]=electron?d:-d;a.derivative[electron?1:2]=electron?-d:d;}
                return a;
            };
            const auto nj=population(true),pj=population(false);
            auto gamma=[&](const AvJet& p,double dos){
                if(p.value/dos<std::numeric_limits<double>::min())return AvJet(1.);
                const double logGamma=experimental::JoyceDixonStatistics::primitiveLogGamma(p.value/dos);
                AvJet g(std::exp(logGamma));
                // Preserve the qualified value path, but differentiate the JD
                // log-gamma primitive directly. In the dilute limit the old
                // 1/n - 1/(Nc*F'(eta)) loses the entire tangent to cancellation.
                const double logGammaSlope=experimental::JoyceDixonStatistics::primitiveLogGammaSlope(p.value/dos);
                for(int k=0;k<3;++k)
                    g.derivative[k]=g.value*logGammaSlope*(p.derivative[k]/dos);
                return g;
            };
            const auto gn=gamma(nj,owner.Nc_[i]),gp=gamma(pj,owner.Nv_[i]);
            AvJet excess;
            if(owner.usesPopulationState()) {
                excess=experimental::populationRecombinationExcess(nj,pj,gn,gp,owner.ni_[i]);
            } else {
                const auto fn=AvJet::variable(s.phin[i],1),fp=AvJet::variable(s.phip[i],2);
                const auto splitting=(fp-fn)/AvJet(vt);
                // Positive splitting must not materialize exp(+splitting).
                excess=splitting.value>=0.
                    ? nj*pj*(-detail::dualExpm1(-splitting))
                    : gn*gp*AvJet(owner.ni_[i]*owner.ni_[i])*detail::dualExpm1(splitting);
            }
            const double dop=owner.recombination_.srhDopingConcentration(owner.doping_.donors(i),owner.doping_.acceptors(i));
            const double rate=owner.recombination_.srhRateGeneralizedFromExcessProduct(excess.value,nj.value,pj.value,
                owner.ni_[i],owner.ni_[i],gn.value,gp.value,dop);
            const double measure=volumes[i]*sourceFactor;out.srh[n+i]=out.srh[2*n+i]=rate*measure;
            if(jac){const auto d=owner.recombination_.srhRateGeneralizedDerivativesFromExcessProduct(excess.value,
                nj.value,pj.value,owner.ni_[i],owner.ni_[i],gn.value,gp.value,dop);
                for(int k=0;k<3;++k){const double value=((d.dRateDn*nj.derivative[k]+d.dRateDp*pj.derivative[k]+
                    (d.dRateDExcess==0.?0.:d.dRateDExcess*excess.derivative[k])+d.dRateDElectronDegeneracy*gn.derivative[k]+
                    d.dRateDHoleDegeneracy*gp.derivative[k])*measure*scale);
                    const double tangent=k==1?tangentN:(k==2?tangentP:1.);
                    const double v=columns?value*((*columns)[k*n+i]/tangent):value/tangent;
                    for(int b:{1,2})if(v!=0)entries.emplace_back(b*n+i,k*n+i,v);
                }
            }
        }
        out.srhJacobian.resize(3*n,3*n);out.srhJacobian.setFromTriplets(entries.begin(),entries.end());return out;
    }
    // Direct fixed-population Poisson assembly. No QF forward reconstruction
    // or subtraction of Debye terms from a QF Jacobian is involved.
    std::pair<VectorXd,SparseMatrixd> populationBase(const VectorXd& x,
        const CoupledDDBoundaryConditions& b,bool derivatives,const VectorXd* columns=nullptr)const {
        const auto s=makeState(x);const int n=active.size();
        VectorXd r=VectorXd::Zero(3*n);
        std::vector<Eigen::Triplet<double>> entries;
        const double area=owner.scaling_.enabled?owner.scaling_.chargeAreaFactor:1.;
        const double epsRef=owner.scaling_.enabled?owner.scaling_.permittivityReference_F_per_m:1.;
        const double rowScale=epsRef*scale;
        for(Index e=0;e<owner.mesh_.numEdges();++e) {
            const auto& edge=owner.mesh_.getEdge(e);if(edge.length<1e-30)continue;
            const int i=edge.n0,j=edge.n1;const double g=owner.edgeAssemblyKernels_[e].poissonCoupling;
            const double flux=(g/epsRef)*(x[i]-x[j]);r[i]+=flux;r[j]-=flux;
            if(derivatives)for(auto [a,c,v]:{std::tuple{i,i,g},std::tuple{j,j,g},std::tuple{i,j,-g},std::tuple{j,i,-g}})
                if(!b.psi.contains(a))entries.emplace_back(a,c,(v/epsRef)*(columns?(*columns)[c]:1.));
        }
        for(int i=0;i<n;++i) {
            double charge=-constants::q*s.p[i]*owner.poissonHoleVol_[i]*area;
            charge-=constants::q*owner.doping_.netDoping(i)*owner.poissonDopantVol_[i]*area;
            charge+=constants::q*s.n[i]*owner.poissonElectronVol_[i]*area;
            charge-=owner.fixedInterfaceChargeRhs_[i];r[i]+=charge/rowScale;
            if(b.psi.contains(i)) {r[i]=x[i]-b.psi.at(i);if(derivatives)entries.emplace_back(i,i,columns?(*columns)[i]:1.);}
            else if(derivatives && active[i]) {
                entries.emplace_back(i,n+i,(constants::q*owner.poissonElectronVol_[i]*area/rowScale)*(columns?(*columns)[n+i]:1.));
                entries.emplace_back(i,2*n+i,(-constants::q*owner.poissonHoleVol_[i]*area/rowScale)*(columns?(*columns)[2*n+i]:1.));
            }
            for(int block:{1,2}) {
                const auto& bc=block==1?b.phin:b.phip;const int row=block*n+i;
                if(!active[i]) {r[row]=x[row];if(derivatives)entries.emplace_back(row,row,columns?(*columns)[row]:1.);continue;}
                if(!bc.contains(i))continue;
                r[row]=(block==1?s.phin[i]:s.phip[i])/scale-bc.at(i);
                if(derivatives) {
                    const double c=block==1?s.n[i]:s.p[i],dos=block==1?owner.Nc_[i]:owner.Nv_[i];
                    const double slope=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(c,dos);
                    entries.emplace_back(row,i,columns?(*columns)[i]:1.);
                    entries.emplace_back(row,row,columns?((block==1?-1.:1.)*((*columns)[row]/slope))*owner.Vt_/scale:(block==1?-1.:1.)*owner.Vt_/slope/scale);
                }
            }
        }
        SparseMatrixd j(3*n,3*n);if(derivatives)j.setFromTriplets(entries.begin(),entries.end());
        return {r,j};
    }
    VectorXd residualFrom(const VectorXd& x,const CoupledDDBoundaryConditions& b,const Evaluation& e)const {
        const int n=active.size();
        VectorXd r=owner.usesPopulationState()?populationBase(x,b,false).first:poisson->residual(x,b);
        const auto carrier=e.transport.source+e.av.source+e.srh;
        for(int i=0;i<n;++i)for(int block:{1,2})if(active[i]&&!constrained(block*n+i,b))r[block*n+i]=carrier[block*n+i];
        return r;
    }
    std::vector<CoupledDDCarrierTermDiagnostic> termsFrom(const Evaluation& e,const VectorXd& r,
        const CoupledDDBoundaryConditions& b,bool includeImpact)const {
        const int n=active.size();
        const auto& u=context.newton.inputScaling.unitSystem();std::vector<CoupledDDCarrierTermDiagnostic> rows(n);
    for(int i=0;i<n;++i){auto& a=rows[i];a.nodeId=i;a.electronDensity_m3=u.internalConcentrationToM3(e.state.n[i]);a.holeDensity_m3=u.internalConcentrationToM3(e.state.p[i]);
        a.electronFlux=e.transport.source[n+i];a.holeFlux=e.transport.source[2*n+i];a.electronFluxAbsSum=e.transport.absoluteSource[n+i];a.holeFluxAbsSum=e.transport.absoluteSource[2*n+i];
        a.electronRecombination=a.holeRecombination=e.srh[n+i];a.electronImpact=a.holeImpact=includeImpact?e.av.source[n+i]:0.;
        a.impactElectronSource=includeImpact?e.av.electronGeneration[i]:0.;a.impactHoleSource=includeImpact?e.av.holeGeneration[i]:0.;
        a.impactCombinedSource=-a.electronImpact;a.electronResidual=r[n+i];a.holeResidual=r[2*n+i];
        a.electronContinuityActive=active[i]&&!b.phin.contains(i);a.holeContinuityActive=active[i]&&!b.phip.contains(i);
        if(!includeImpact){if(a.electronContinuityActive)a.electronResidual-=e.av.source[n+i];if(a.holeContinuityActive)a.holeResidual-=e.av.source[2*n+i];}
        if(!active[i]){a.electronGauge=r[n+i];a.holeGauge=r[2*n+i];}
        if(b.phin.contains(i))a.electronBoundary=r[n+i]-a.electronFlux-a.electronRecombination-a.electronImpact-a.electronGauge;
        if(b.phip.contains(i))a.holeBoundary=r[2*n+i]-a.holeFlux-a.holeRecombination-a.holeImpact-a.holeGauge;
    }return rows;
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
    return p.residualFrom(x,b,e);
}
SparseMatrixd BvMethodsCoupledRuntime::jacobian(const VectorXd& x,const CoupledDDBoundaryConditions& b,bool approximate,const VectorXd* columns,bool freezeStatisticalFloor)const {
    auto& p=*impl_;p.validateBoundary(b);p.sync();if(columns && (!p.owner.usesPopulationState()||columns->size()!=x.size()||!columns->allFinite()||(columns->array()<=0.).any()))throw std::invalid_argument("Invalid population linearization columns");
    auto e=columns?p.compute(x,true,approximate,columns,freezeStatisticalFloor):(approximate?p.compute(x,true,true,nullptr,freezeStatisticalFloor):p.evaluate(x,true));const int n=p.active.size();
    const auto base=p.owner.usesPopulationState()?p.populationBase(x,b,true,columns).second:p.poisson->assembleJacobian(x,b);const SparseMatrixd carrier=e.transport.jacobian+e.av.jacobian+e.srhJacobian;
    std::vector<Eigen::Triplet<double>> t;
    for(int k=0;k<base.outerSize();++k)for(SparseMatrixd::InnerIterator v(base,k);v;++v)
        if(v.row()<n || !p.active[v.row()%n] || p.constrained(v.row(),b))t.emplace_back(v.row(),v.col(),v.value());
    for(int k=0;k<carrier.outerSize();++k)for(SparseMatrixd::InnerIterator v(carrier,k);v;++v)
        if(!p.constrained(v.row(),b))t.emplace_back(v.row(),v.col(),v.value());
    SparseMatrixd result(3*n,3*n);result.setFromTriplets(t.begin(),t.end());return result;
}
std::vector<CoupledDDCarrierTermDiagnostic> BvMethodsCoupledRuntime::terms(const VectorXd& x,const CoupledDDBoundaryConditions& b,bool includeImpact)const {
    auto& p=*impl_;p.validateBoundary(b);p.sync();const auto e=p.evaluate(x,false);
    return p.termsFrom(e,p.residualFrom(x,b,e),b,includeImpact);
}
CoupledDDPhysicalDiagnostics BvMethodsCoupledRuntime::physicalDiagnostics(
    const VectorXd& x,const CoupledDDBoundaryConditions& b)const {
    auto& p=*impl_;p.validateBoundary(b);p.sync();const auto e=p.evaluate(x,false);
    CoupledDDPhysicalDiagnostics out;out.residual=p.residualFrom(x,b,e);
    const bool impact=p.owner.impactIonizationCoupled_;
    out.constrainedTerms=p.termsFrom(e,out.residual,b,impact);
    out.unconstrainedTerms=p.termsFrom(e,p.residualFrom(x,{},e),{},impact);
    return out;
}

VectorXd BvMethodsCoupledRuntime::avalancheResidual(const VectorXd& x,const CoupledDDBoundaryConditions& b)const {
    auto& p=*impl_;auto e=p.evaluate(x,false);auto r=e.av.source;
    for(int i=0;i<r.size();++i)if(p.constrained(i,b))r[i]=0;return r;
}
SparseMatrixd BvMethodsCoupledRuntime::avalancheJacobian(const VectorXd& x,const CoupledDDBoundaryConditions& b)const {
    return impl_->project(impl_->evaluate(x,true).av.jacobian,b);
}
VectorXd BvMethodsCoupledRuntime::transportResidual(const VectorXd& x,BvMethodsTransportPurpose purpose)const {
    auto& p=*impl_;p.sync();
    try{return p.transport(x,false,purpose).source;}
    catch(const std::domain_error&){return VectorXd::Constant(x.size(),std::numeric_limits<double>::quiet_NaN());}
}
SparseMatrixd BvMethodsCoupledRuntime::transportJacobian(const VectorXd& x,const VectorXd* columns)const {
    auto& p=*impl_;p.sync();
    if(columns){if(!p.owner.usesPopulationState()||columns->size()!=x.size()||!columns->allFinite()||(columns->array()<=0.).any())throw std::invalid_argument("Invalid population transport columns");return p.assembleTransport(p.makeState(x),true,false,columns,0.,false,"bvmethods.tri3.terminal_jacobian").jacobian;}
    return p.transport(x,true).jacobian;
}
std::vector<CoupledDDPoissonTermDiagnostic> BvMethodsCoupledRuntime::populationPoissonTerms(
    const VectorXd& x,const CoupledDDBoundaryConditions& b)const {
    auto& p=*impl_;const auto& a=p.owner;
    if(!a.usesPopulationState())throw std::invalid_argument("population Poisson diagnostics require population state");
    const auto s=p.makeState(x);const int n=p.active.size();
    const auto actual=p.populationBase(x,b,false).first;
    const auto free=p.populationBase(x,{},false).first;
    const double rowScale=a.scaling_.enabled?a.scaling_.permittivityReference_F_per_m*p.scale:1.;
    const double area=a.scaling_.enabled?a.scaling_.chargeAreaFactor:1.;
    std::vector<CoupledDDPoissonTermDiagnostic> rows(n);
    for(int i=0;i<n;++i) {
        auto& t=rows[i];t.nodeId=i;t.intrinsicDensity=a.ni_[i];
        t.electronDensityOfStates=a.Nc_[i];t.holeDensityOfStates=a.Nv_[i];
        t.suppliedElectronDensity=t.reconstructedElectronDensity=s.n[i];
        t.suppliedHoleDensity=t.reconstructedHoleDensity=s.p[i];t.hasSuppliedCarrierState=true;
        t.electronRequiredQuasiFermiPotential=s.phin[i];t.holeRequiredQuasiFermiPotential=s.phip[i];
        t.electronCharge=constants::q*s.n[i]*a.poissonElectronVol_[i]*area/rowScale;
        t.holeCharge=-constants::q*s.p[i]*a.poissonHoleVol_[i]*area/rowScale;
        t.dopingCharge=-constants::q*a.doping_.netDoping(i)*a.poissonDopantVol_[i]*area/rowScale;
        t.fixedInterfaceCharge=-a.fixedInterfaceChargeRhs_[i]/rowScale;
        t.unconstrainedResidual=free[i];t.productionResidual=actual[i];t.constrained=b.psi.contains(i);
        t.boundaryReplacement=t.productionResidual-t.unconstrainedResidual;
    }
    for(Index e=0;e<a.mesh_.numEdges();++e) {
        const auto& edge=a.mesh_.getEdge(e);if(edge.length<1e-30)continue;
        const double f=a.edgeAssemblyKernels_[e].poissonCoupling*(s.psi[edge.n0]-s.psi[edge.n1])/rowScale;
        rows[edge.n0].dielectricFlux+=f;rows[edge.n1].dielectricFlux-=f;
        rows[edge.n0].dielectricFluxAbsSum+=std::abs(f);rows[edge.n1].dielectricFluxAbsSum+=std::abs(f);
    }
    for(auto& t:rows)t.closureError=t.productionResidual-(t.dielectricFlux+t.electronCharge+t.holeCharge+
        t.dopingCharge+t.fixedInterfaceCharge+t.boundaryReplacement);
    return rows;
}
CoupledDDPhysicalDiagnostics CoupledDDAssembler::bvMethodsPhysicalDiagnostics(
    const VectorXd& x,const CoupledDDBoundaryConditions& b)const {
    auto r=bvMethodsRuntime();if(!r)throw std::invalid_argument("BVMethods diagnostics requested for legacy operator");
    return r->physicalDiagnostics(x,b);
}
VectorXd CoupledDDAssembler::bvMethodsTransportResidual(const VectorXd& x,BvMethodsTransportPurpose purpose)const {
    auto r=bvMethodsRuntime();if(!r)throw std::invalid_argument("BVMethods transport requested for legacy operator");return r->transportResidual(x,purpose);
}
SparseMatrixd CoupledDDAssembler::bvMethodsTransportJacobian(const VectorXd& x)const {
    auto r=bvMethodsRuntime();if(!r)throw std::invalid_argument("BVMethods transport requested for legacy operator");return r->transportJacobian(x);
}
SparseMatrixd CoupledDDAssembler::populationProtectedJacobian(const VectorXd& x,const CoupledDDBoundaryConditions& b,const VectorXd* columns)const {
    auto r=bvMethodsRuntime();
    if(!r||!usesPopulationState())throw std::invalid_argument("Protected tangents require BVMethods population state");
    // Keep the current-state physical operator and statistical edge tangents;
    // only mobility and avalanche derivatives are omitted.
    return r->jacobian(x,b,true,columns,false);
}
SparseMatrixd CoupledDDAssembler::populationJacobianInColumns(const VectorXd& x,const CoupledDDBoundaryConditions& b,const VectorXd& columns,bool approximate)const {
    auto r=bvMethodsRuntime();if(!r||!usesPopulationState())throw std::invalid_argument("Direct population columns require BVMethods population state");return r->jacobian(x,b,approximate,&columns);
}
SparseMatrixd CoupledDDAssembler::bvMethodsTransportJacobianInColumns(const VectorXd& x,const VectorXd& columns)const {
    auto r=bvMethodsRuntime();if(!r||!usesPopulationState())throw std::invalid_argument("Direct population transport columns require BVMethods population state");return r->transportJacobian(x,&columns);
}
std::shared_ptr<BvMethodsCoupledRuntime> CoupledDDAssembler::bvMethodsRuntime()const {
    if(scaling_.coupledOperator=="legacy")return {};
    if(scaling_.coupledOperator!="bvmethods_tri3_qf")throw std::invalid_argument("Unknown coupled operator");
    if(!bvMethodsRuntime_)bvMethodsRuntime_=std::make_shared<BvMethodsCoupledRuntime>(*this);return bvMethodsRuntime_;
}
}
