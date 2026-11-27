#include "vela/equation/LegacyPopulationRuntime.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/discretization/Bernoulli.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/solver/SlotPopulationNumerics.h"
#include "vela/core/PhysicalConstants.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vela {
struct LegacyPopulationRuntime::Impl {
    const CoupledDDAssembler& owner;
    std::unique_ptr<CoupledDDAssembler> qf;
    explicit Impl(const CoupledDDAssembler& a):owner(a) {
        if(!a.usesPopulationState() || a.usesBvMethodsOperator() ||
           a.splitDDEnabled_ || a.electronQuantumPotentialConfig_.enabled || a.extendedPoissonResidual_ ||
           a.scaling_.populationSgCurrent!="quasi_fermi" ||
           (a.carrierStatisticsModel_==CarrierStatisticsModel::FermiDiracJoyceDixon &&
            !a.effectiveBandSg_))
            throw std::invalid_argument("legacy population runtime requires classical Boltzmann/Fermi SLOT physics without split or quantum coordinates");
        auto scaling=a.scaling_;scaling.stateVariables="potential_quasi_fermi";
        qf=std::make_unique<CoupledDDAssembler>(a.mesh_,a.matdb_,a.doping_,a.Vt_,
            a.mobilityConfig_,a.recombinationConfig_,a.bandgapNarrowingConfig_,
            a.impactIonizationConfig_,std::vector<RegionFixedChargeSpec>{},
            std::vector<InterfaceSheetChargeSpec>{},scaling,a.carrierDiagonalFloor_,
            a.carrierStatistics_,a.electronQuantumPotentialConfig_);
        // Preserve the already resolved physical charge and statistical inputs.
        qf->fixedInterfaceChargeRhs_=a.fixedInterfaceChargeRhs_;
        qf->ni_=a.ni_;qf->niElectron_=a.niElectron_;qf->niHole_=a.niHole_;
        qf->intrinsicReference_=a.intrinsicReference_;qf->Nc_=a.Nc_;qf->Nv_=a.Nv_;
        qf->setLegacyEdgeJacobianSteps(a.legacyTransportJacobianStep_,a.legacyAvalancheJacobianStep_);
    }
    VectorXd bind(const VectorXd& y,const CoupledDDBoundaryConditions& b) const {
        if(!b.thermionic.empty() || !b.psiLow.empty() || !b.phinLow.empty() || !b.phipLow.empty())
            throw std::invalid_argument("legacy population contacts require classical Ohmic constraints");
        const auto s=owner.unpack(y);
        const Real scale=owner.scaling_.enabled?owner.scaling_.V0:1.;
        for(int i=0;i<s.n.size();++i)if(owner.ni_[i]>0.) {
            const Real n=vela::electronDensity(owner.niElectron_[i],owner.Nc_[i],
                s.psi[i]*scale,s.phin[i]*scale,owner.Vt_,owner.carrierStatisticsModel_);
            const Real p=vela::holeDensity(owner.niHole_[i],owner.Nv_[i],
                s.psi[i]*scale,s.phip[i]*scale,owner.Vt_,owner.carrierStatisticsModel_);
            if(!std::isfinite(n) || !std::isfinite(p) || n<=0. || p<=0. ||
               std::abs((static_cast<long double>(n)-s.n[i])/s.n[i])>=1e-7L ||
               std::abs((static_cast<long double>(p)-s.p[i])/s.p[i])>=1e-7L)
                throw std::runtime_error("legacy population inverse guard: QF consumer is inconsistent with authoritative density");
        }
        qf->authoritativeElectronDensity_=s.n;
        qf->authoritativeHoleDensity_=s.p;
        return qf->pack(s);
    }
    bool usesDirectIalJacobian() const {
        return owner.mobilityConfig_.model=="ialmob" && !owner.usesFermiDirac_ &&
            owner.mobilityConfig_.carrierCurrentDiscretization=="scharfetter_gummel_edge" &&
            !owner.carrierDiagonalFloor_.enabled && !owner.recombination_.bandToBandEnabled() &&
            (!owner.impactIonizationEnabled_ ||
             detail::usesIalElementAvalancheSource(owner.impactIonizationConfig_));
    }
    SparseMatrixd directIalJacobian(
        const VectorXd& y, const CoupledDDBoundaryConditions& b, bool sourceOnly=false) const {
        using Triplet = Eigen::Triplet<Real>;
        const VectorXd x=bind(y,b);
        const auto s=owner.unpack(y);
        const int count=owner.mesh_.numNodes();
        const Real scale=owner.scaling_.enabled?owner.scaling_.V0:1.;
        const Real rowScale=qf->continuityResidualScale();
        const Real fieldFactor=owner.scaling_.enabled?owner.scaling_.fieldFromCoordinateDeltaFactor:1.;
        const Real sourceFactor=owner.scaling_.enabled?
            owner.scaling_.unitSystem.continuitySourceIntegralFactor()/rowScale:1.;
        const VectorXd psi=s.psi*scale,fn=s.phin*scale,fp=s.phip*scale;
        qf->refreshIalMobility(x);
        const auto& cfg=qf->mobilityConfig_;
        const auto& geometry=*cfg.ialmobGeometry;
        const auto& transport=*cfg.ialmobState;
        const Real concentration=cfg.internalConcentrationToM3;
        const auto vertexStates=[&](const Cell& cell) {
            std::array<IalElementVertexState,3> states;
            for (int j=0;j<3;++j) {
                const auto i=cell.node_ids[j];
                states[j]={psi[i],fn[i],fp[i],
                    owner.doping_.donors(i)*concentration,owner.doping_.acceptors(i)*concentration,
                    s.n[i]*concentration,s.p[i]*concentration};
            }
            return states;
        };
        std::vector<Triplet> entries;
        entries.reserve(owner.mesh_.numEdges()*24);
        std::vector<bool> activeN(count),activeP(count);
        const auto add=[&](int row, int col, Real value) {
            if (value==0.) return;
            if (row<count ? b.psi.contains(row) : row<2*count ?
                b.phin.contains(row-count) : b.phip.contains(row-2*count)) return;
            entries.emplace_back(row,col,value);
        };
        const auto convert=[&](const detail::Tri3LocalForwardDual& value,
                                const Cell& cell,int k) {
            const auto node=cell.node_ids[k/3];
            const int component=k%3;
            const long double d=value.derivative[k];
            if (component==0) return static_cast<Real>(d*scale);
            return d==0.L ? Real{0.} : static_cast<Real>(
                d/(component==1?s.n[node]:s.p[node]));
        };
        if (!sourceOnly) {
            const bool feedback=cfg.ialmobTransportFeedbackDerivatives &&
                qf->transportMobilityDerivativeEnabled_;
            std::vector<IalElementMobilityResult> derivatives;
            if (feedback) {
                derivatives.reserve(geometry.cells.size());
                for (const auto& support:geometry.cells)
                    derivatives.push_back(evaluateIalTransportCell(cfg,support.cellId,
                        vertexStates(owner.mesh_.getCell(support.cellId)),true,owner.Vt_));
            }
            for (Index e=0;e<owner.mesh_.numEdges();++e) {
                const auto& edge=owner.edgeAssemblyKernels_[e];
                const int i=edge.n0,j=edge.n1;
                if (edge.length<1e-30) continue;
                const auto fluxColumns=[&](bool electron, Real mu, const VectorXd& density,
                        const std::vector<Real>& ni, std::vector<bool>& active) {
                    if (mu<=0.) return;
                    active[i]=active[j]=true;
                    const int offset=electron?count:2*count;
                    const Real sign=electron?1.:-1.;
                    const Real coef=mu*owner.Vt_*fieldFactor*owner.couple_[e]/edge.length;
                    const auto drift = sgPopulationEdgeDrift(
                        s.psi[i], s.psi[j], sign * scale, owner.Vt_, ni[i], ni[j], owner.bgnEnabled_);
                    const Real a = drift.hi;
                    const auto scatter=[&](int col,Real v) {
                        add(offset+i,col,v);
                        add(offset+j,col,-v);
                    };
                    const Real potential=-sign*coef/owner.Vt_*scale/rowScale*
                        (bernoulliDerivative(-a)*density[i]+bernoulliDerivative(a)*density[j]);
                    scatter(i,-potential);
                    scatter(j,potential);
                    if (ni[i]>0.) scatter(offset+i,coef/rowScale*bernoulli(-a));
                    if (ni[j]>0.) scatter(offset+j,-coef/rowScale*bernoulli(a));
                    if (!feedback) return;
                    const Real flux = sgBoltzmannPopulationContinuityFluxFromDrift(
                        density[i], density[j], drift, owner.Vt_, coef, 0., 0., 0.);
                    for (const auto& contribution:geometry.edges[e]) {
                        const auto& cell=owner.mesh_.getCell(geometry.cells[contribution.support].cellId);
                        const auto& local=derivatives[contribution.support];
                        const auto& m=electron?local.electron:local.hole;
                        for (int k=0;k<9;++k)
                            scatter((k%3)*count+cell.node_ids[k/3],static_cast<Real>(
                                static_cast<long double>(flux)/mu*contribution.weight/rowScale*
                                convert(m,cell,k)));
                    }
                };
                fluxColumns(true,transport.electronEdges[e],s.n,owner.niElectron_,activeN);
                fluxColumns(false,transport.holeEdges[e],s.p,owner.niHole_,activeP);
            }
        }
        if (owner.impactIonizationCoupled_) {
            auto options=detail::ialElementAvalancheOptionsFromConfig(owner.impactIonizationConfig_);
            options.populationThermalVoltage_V=owner.Vt_;
            const Real measure=concentration*cfg.internalLengthToM*cfg.internalLengthToM;
            for (const auto& cell:owner.mesh_.cells()) {
                if (owner.mesh_.getRegion(cell.region_id).material!="Si") continue;
                const auto local=evaluateIalTransportAvalancheCell(
                    cfg,cell.id,vertexStates(cell),options);
                for (int j=0;j<3;++j) {
                    const int row=cell.node_ids[j];
                    const auto source=local.electron.source_per_m_s[j]+local.hole.source_per_m_s[j];
                    if (source.value!=0.) activeN[row]=activeP[row]=true;
                    for (int k=0;k<9;++k) {
                        const Real d=-convert(source,cell,k)/measure*sourceFactor;
                        const int col=(k%3)*count+cell.node_ids[k/3];
                        add(count+row,col,d);
                        add(2*count+row,col,d);
                    }
                }
            }
        }
        if (!sourceOnly) {
            for (int i=0;i<count;++i) {
                if (owner.ni_[i]>0. &&
                    (owner.recombination_.srhEnabled() || owner.recombination_.augerEnabled())) {
                    const auto r=qf->nodeRecombinationDifferential(i,x,s.n[i],s.p[i],true);
                    if (r.rate!=0. || activeN[i] || activeP[i]) {
                        for (int k=1;k<3;++k) {
                            const Real d=r.derivative[k]*owner.recombinationVol_[i]*sourceFactor;
                            add(count+i,k*count+i,d);
                            add(2*count+i,k*count+i,d);
                        }
                        activeN[i]=activeP[i]=true;
                    }
                }
                const auto constraint=[&](int offset, bool electron, bool active, bool pinned) {
                    if (active && !pinned) return;
                    if (owner.ni_[i]>0.) {
                        entries.emplace_back(offset+i,i,1.);
                        entries.emplace_back(offset+i,offset+i,
                            (electron?-1.:1.)*owner.Vt_/scale/(electron?s.n[i]:s.p[i]));
                    } else entries.emplace_back(offset+i,offset+i,1.);
                };
                constraint(count,true,activeN[i],b.phin.contains(i));
                constraint(2*count,false,activeP[i],b.phip.contains(i));
            }
        }
        if (!sourceOnly)
            appendPoisson(entries,b);
        SparseMatrixd result(3*count,3*count);
        result.setFromTriplets(entries.begin(),entries.end());
        return result;
    }
    void appendPoisson(std::vector<Eigen::Triplet<Real>>& entries,
                       const CoupledDDBoundaryConditions& b) const {
        const int count=owner.mesh_.numNodes();
        const Real scale=owner.scaling_.enabled?owner.scaling_.V0:1.;
        const Real rowScale=owner.scaling_.enabled?
            owner.scaling_.permittivityReference_F_per_m*scale:1.;
        const Real area=owner.scaling_.enabled?owner.scaling_.chargeAreaFactor:1.;
        for (Index e=0;e<owner.mesh_.numEdges();++e) {
            const auto& edge=owner.mesh_.getEdge(e);
            if (edge.length<1e-30) continue;
            const Real g=owner.edgeAssemblyKernels_[e].poissonCoupling*scale/rowScale;
            if (!b.psi.contains(edge.n0)) {
                entries.emplace_back(edge.n0,edge.n0,g);
                entries.emplace_back(edge.n0,edge.n1,-g);
            }
            if (!b.psi.contains(edge.n1)) {
                entries.emplace_back(edge.n1,edge.n1,g);
                entries.emplace_back(edge.n1,edge.n0,-g);
            }
        }
        for (int i=0;i<count;++i) {
            if (b.psi.contains(i)) entries.emplace_back(i,i,1.);
            else if (owner.ni_[i]>0.) {
                entries.emplace_back(i,count+i,constants::q*owner.poissonElectronVol_[i]*area/rowScale);
                entries.emplace_back(i,2*count+i,-constants::q*owner.poissonHoleVol_[i]*area/rowScale);
            }
        }
    }
    SparseMatrixd transform(const SparseMatrixd& a,const VectorXd& y,
                             const CoupledDDBoundaryConditions* poisson=nullptr) const {
        const int n=owner.mesh_.numNodes();
        const long double scale=owner.scaling_.enabled?owner.scaling_.V0:1.;
        std::vector<Eigen::Triplet<Real>> t;t.reserve(2*a.nonZeros());
        for(int j=0;j<a.outerSize();++j)for(SparseMatrixd::InnerIterator v(a,j);v;++v) {
            if(v.value()==0.)continue;
            if(poisson && v.row()<n)continue;
            if(v.col()<n || owner.ni_[v.col()%n]<=0.) {
                t.emplace_back(v.row(),v.col(),v.value());continue;
            }
            const int i=v.col()%n;const bool electron=v.col()<2*n;
            const Real density=y[v.col()],dos=electron?owner.Nc_[i]:owner.Nv_[i];
            long double logarithmicSlope=1.;
            if(owner.usesFermiDirac_) {
                const auto model=owner.carrierStatisticsModel_;
                const Real eta=slot_population::populationEtaAccurate(density,dos,model);
                if(eta>=-40.) {
                    const Real value=fermiDiracHalf(eta,model);
                    const Real derivative=fermiDiracHalfDerivative(eta,model);
                    if(!(value>0.) || !(derivative>0.))throw std::domain_error("invalid population statistical slope");
                    logarithmicSlope=static_cast<long double>(derivative)/value;
                }
            }
            // Do not materialize 1/density or an overflowing QF conversion
            // matrix in binary64. Multiply each nonzero before division.
            const long double c=static_cast<long double>(v.value())*(electron?-1.L:1.L)*
                owner.Vt_/scale/static_cast<long double>(density)/logarithmicSlope;
            if(!std::isfinite(c) || std::abs(c)>std::numeric_limits<Real>::max())
                throw std::overflow_error("unrepresentable legacy population Jacobian coefficient");
            t.emplace_back(v.row(),v.col(),static_cast<Real>(c));
            t.emplace_back(v.row(),i,v.value()); // d(phi_n,p)/d(psi) = 1
        }
        if(poisson)
            appendPoisson(t,*poisson);
        SparseMatrixd result(a.rows(),a.cols());result.setFromTriplets(t.begin(),t.end());
        return result;
    }
};
LegacyPopulationRuntime::LegacyPopulationRuntime(const CoupledDDAssembler& a):impl_(std::make_unique<Impl>(a)){}
LegacyPopulationRuntime::~LegacyPopulationRuntime()=default;
VectorXd LegacyPopulationRuntime::residual(const VectorXd& y,const CoupledDDBoundaryConditions& b)const {
    const auto x=impl_->bind(y,b);return impl_->qf->residual(x,b);
}
SparseMatrixd LegacyPopulationRuntime::jacobian(const VectorXd& y,const CoupledDDBoundaryConditions& b,bool approximate)const {
    if(approximate)throw std::invalid_argument("legacy SLOT population derivatives require the full physical Jacobian");
    if(impl_->usesDirectIalJacobian())return impl_->directIalJacobian(y,b);
    const auto x=impl_->bind(y,b);return impl_->transform(impl_->qf->assembleJacobian(x,b),y,&b);
}
std::vector<CoupledDDCarrierTermDiagnostic> LegacyPopulationRuntime::terms(const VectorXd& y,const CoupledDDBoundaryConditions& b,bool impact)const {
    const auto x=impl_->bind(y,b);return impl_->qf->carrierContinuityTermDiagnosticsImpl(x,b,nullptr,impact,nullptr);
}
CoupledDDPhysicalDiagnostics LegacyPopulationRuntime::physicalDiagnostics(const VectorXd& y,const CoupledDDBoundaryConditions& b)const {
    return {residual(y,b),terms(y,b),terms(y,{})};
}
VectorXd LegacyPopulationRuntime::avalancheResidual(const VectorXd& y,const CoupledDDBoundaryConditions& b)const {
    const auto x=impl_->bind(y,b);return impl_->qf->impactIonizationSourceResidual(x,b);
}
SparseMatrixd LegacyPopulationRuntime::avalancheJacobian(const VectorXd& y,const CoupledDDBoundaryConditions& b)const {
    if(impl_->usesDirectIalJacobian())return impl_->directIalJacobian(y,b,true);
    const auto x=impl_->bind(y,b);return impl_->transform(impl_->qf->impactIonizationSourceJacobian(x,b),y);
}
std::vector<CoupledDDPoissonTermDiagnostic> LegacyPopulationRuntime::populationPoissonTerms(const VectorXd& y,const CoupledDDBoundaryConditions& b)const {
    const auto x=impl_->bind(y,b);return impl_->qf->poissonTermDiagnostics(x,b);
}
VectorXd LegacyPopulationRuntime::transportResidual(const VectorXd&)const {
    throw std::invalid_argument("source-free BVMethods transport is not a legacy SLOT operator");
}
SparseMatrixd LegacyPopulationRuntime::transportJacobian(const VectorXd&)const {
    throw std::invalid_argument("source-free BVMethods transport is not a legacy SLOT operator");
}
}
