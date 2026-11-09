#include "vela/equation/LegacyPopulationRuntime.h"
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
           a.carrierStatisticsModel_==CarrierStatisticsModel::FermiDiracJoyceDixon)
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
                const Real eta=slot_population::populationEtaAccurate(density,dos);
                if(eta>=-40.) {
                    const Real value=fermiDiracHalf(eta),derivative=fermiDiracHalfDerivative(eta);
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
        if(poisson) {
            // Charge is linear in authoritative n,p. Assemble these columns
            // directly, rather than divide an already underflowed QF charge
            // derivative or cancel its potential columns in binary64.
            const Real rowScale=owner.scaling_.enabled?
                owner.scaling_.permittivityReference_F_per_m*owner.scaling_.V0:1.;
            const Real area=owner.scaling_.enabled?owner.scaling_.chargeAreaFactor:1.;
            for(Index e=0;e<owner.mesh_.numEdges();++e) {
                const auto& edge=owner.mesh_.getEdge(e);if(edge.length<1e-30)continue;
                const Real g=owner.edgeAssemblyKernels_[e].poissonCoupling*static_cast<Real>(scale)/rowScale;
                if(!poisson->psi.contains(edge.n0)) {t.emplace_back(edge.n0,edge.n0,g);t.emplace_back(edge.n0,edge.n1,-g);}
                if(!poisson->psi.contains(edge.n1)) {t.emplace_back(edge.n1,edge.n1,g);t.emplace_back(edge.n1,edge.n0,-g);}
            }
            for(int i=0;i<n;++i) {
                if(poisson->psi.contains(i))t.emplace_back(i,i,1.);
                else if(owner.ni_[i]>0.) {
                    t.emplace_back(i,n+i,constants::q*owner.poissonElectronVol_[i]*area/rowScale);
                    t.emplace_back(i,2*n+i,-constants::q*owner.poissonHoleVol_[i]*area/rowScale);
                }
            }
        }
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
