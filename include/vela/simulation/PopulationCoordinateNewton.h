#pragma once

#include "vela/simulation/PrimaryVariableCoordinates.h"

namespace vela {

/// Production population-coordinate entry. The Newton algorithm remains the
/// single CoupledLoadLineNewton implementation. BVMethods can select an
/// equivalent defect calculation for the PARDISO bordered path.
class PopulationCoordinateNewton {
public:
    PopulationCoordinateNewton(ArclengthSystem system,
        CoupledLoadLineNewton::CurrentFunctional current,CoupledLoadLineConfig config,
        ArclengthScalarFunctional port={},DensityUpdateErrorEvaluator metric={},
        PrimaryVariableScales references={},bool unit=false,bool approximateNativeStop=false,
        ApproximateNativeRhsScales mapping={})
        :system_(std::move(system)),port_(std::move(port)),config_(std::move(config)),
         metric_(std::move(metric)),references_(std::move(references)),unit_(unit),
         approximateNativeStop_(approximateNativeStop),mapping_(mapping)
    {
        if(!port_.value)port_.value=std::move(current);
        if(unit_ && (!system_.makeTrial || !system_.columnScales || !system_.prescribedRows ||
                     config_.linearSolver!="direct_bordered"))
            throw std::invalid_argument("Unit primary coordinates require a directly bordered population system");
        if(approximateNativeStop_ && !metric_)
            throw std::invalid_argument("approximate native stopping requires a physical update metric");
        if(approximateNativeStop_ && config_.updateConvergence.mode=="off")
            config_.updateConvergence.mode="report";
    }

    ArclengthState secantPredict(const ArclengthState& previous,const ArclengthState& current,
                                Real target,Real resistance) const {
        return CoupledLoadLineNewton(system_,port_.value,config_,port_,metric_)
            .secantPredict(previous,current,target,resistance);
    }

    CoupledLoadLineResult solve(const ArclengthState& initial,Real target,Real resistance,
        const derivative_recovery::DerivativeSwitch* history=nullptr) const {
        auto references=references_;
        if(!unit_)references={VectorXd::Ones(initial.x.size()),1.,1.};
        references.validate(initial.x.size());
        const PrimaryVariableCoordinates coordinates(references);
        const auto device=unit_?coordinates.system(system_):system_;
        const auto port=unit_?coordinates.port(port_):port_;
        const auto metric=unit_?coordinates.metric(metric_):metric_;
        auto config=config_;
        config.loadLineTolerance_V/=references.voltage;
        config.maxVoltageUpdate_V/=references.voltage;

        // Only the latest physical state is cached. A new target owns a fresh
        // cache, preventing stale accepted-state histories from authorizing stop.
        ArclengthState loadState;bool loadValid=false;LoadLineResidual load{0.};
        const auto physicalLoad=[&](const ArclengthState& state)->LoadLineResidual {
            const bool same=loadValid && state.x.size()==loadState.x.size() &&
                std::memcmp(state.x.data(),loadState.x.data(),state.x.size()*sizeof(Real))==0 &&
                std::memcmp(&state.lambda,&loadState.lambda,sizeof(Real))==0;
            if(!same) {
                loadState=state;
                load=LoadLineResidual::evaluate(config_.voltageCoefficient,state.lambda,
                    resistance,port_.value(state.x,state.lambda),target);
                loadValid=true;
            }
            return load;
        };
        // The mapped stopping report and the original-system audit request the
        // same physical border at a fixed state. Reuse its exact coefficients,
        // not a reconstructed/rounded Unit matrix. Entries, raw residual and
        // derivative choice remain separate; no cache survives this target.
        ArclengthState assemblyState;
        bool assemblyValid=false,fullValid=false,approximateValid=false,residualValid=false;
        SparseMatrixd fullBorder,approximateBorder;
        VectorXd physicalResidual;
        const auto atAssemblyState=[&](const ArclengthState& state) {
            const bool same=assemblyValid && state.x.size()==assemblyState.x.size() &&
                std::memcmp(state.x.data(),assemblyState.x.data(),state.x.size()*sizeof(Real))==0 &&
                std::memcmp(&state.lambda,&assemblyState.lambda,sizeof(Real))==0;
            if(!same) {
                assemblyState=state;assemblyValid=true;
                fullValid=approximateValid=residualValid=false;
            }
        };
        const auto physicalBorder=[&](const ArclengthState& state,bool full)->const SparseMatrixd& {
            atAssemblyState(state);
            auto& valid=full?fullValid:approximateValid;
            auto& matrix=full?fullBorder:approximateBorder;
            if(!valid) {matrix=border(state,resistance,full);valid=true;}
            return matrix;
        };
        CoupledLoadLineCallbacks hooks;
        hooks.canonicalLoadResidual=[&](const ArclengthState& state) {
            return physicalLoad(coordinates.physical(state)).inVoltageUnits(references.voltage);
        };
        hooks.physicalStepNorm=[&](const VectorXd& step) {
            return coordinates.physical(step).cwiseAbs().maxCoeff();
        };
        if(unit_)hooks.validatePhysicalDirection=[&](const ArclengthState& state,
            const VectorXd& dx,Real dv,bool full,CoupledLoadLineResult::IterationTrace& trace,
            std::string& reason) {
            const auto physical=coordinates.physical(state);
            const auto& A=physicalBorder(physical,full);
            if(!residualValid) {
                physicalResidual=system_.residual(physical.x,physical.lambda);residualValid=true;
            }
            VectorXd b(A.rows());b.head(physical.x.size())=-physicalResidual;
            b[b.size()-1]=physicalLoad(physical).physicalRhs();
            VectorXd direction(b.size());direction.head(dx.size())=coordinates.physical(dx);
            direction[dx.size()]=dv*references.voltage;
            const auto calculateDefect = config.bvMethodsEftDefect &&
                config.linearSolver != "direct_bordered_qr" &&
                LinearSolver::selectedBackend() == "pardiso_mkl"
                ? linear_refinement::eftDefect : linear_refinement::highPrecisionDefect;
            const auto defect=calculateDefect(A,-b,direction);
            const auto rows=auditRows(A,b,system_.columnScales(physical));
            trace.physicalLinearAbsolute=defect.values.cwiseProduct(rows).cwiseAbs().maxCoeff();
            trace.physicalLinearBackward=defect.backwardMax;
            const bool pass=std::isfinite(trace.physicalLinearAbsolute) &&
                trace.physicalLinearAbsolute<=1e-8 && trace.physicalLinearBackward<=1e-10;
            if(!pass)reason="original physical linear precision gate failed";
            return pass;
        };
        if(approximateNativeStop_)hooks.stoppingRule=[&](const ArclengthState& state,
            const VectorXd& f,CoupledLoadLineResult& result) {
            const auto physical=coordinates.physical(state);
            atAssemblyState(physical);
            physicalResidual=f;residualValid=true;
            const Real L=physicalLoad(physical).voltage_V;
            if(!physical.x.allFinite() || !std::isfinite(physical.lambda) || !f.allFinite() || !std::isfinite(L))
                return std::string{};
            const bool update=result.updateError.valid && std::isfinite(result.updateError.rms) &&
                result.updateError.rms>=0. && result.updateError.rms<1.;
            Real mapped=std::numeric_limits<Real>::quiet_NaN();
            if(mapping_.available()) {
                // Mapping is diagnostic, not a qualified native RHSMin branch.
                // A missing/invalid map never authorizes residual-based stop.
                try {
                    const auto& A=physicalBorder(physical,true);
                    VectorXd rhs(f.size()+1);rhs.head(f.size())=-f;rhs[f.size()]=-L;
                    const auto reduced=eliminateLinearConstraints(A,rhs,
                        system_.prescribedRows(physical.x,physical.lambda));
                    mapped=0.;
                    for(int i=0;i<reduced.rhs.size();++i) {
                        const int row=reduced.kept[i];
                        const Real scale=row==f.size()?mapping_.load:
                            (row<f.size()/3?mapping_.poisson:mapping_.carrier);
                        mapped=std::max(mapped,std::abs(reduced.rhs[i])*scale);
                    }
                } catch(const std::exception&) {mapped=std::numeric_limits<Real>::quiet_NaN();}
            }
            result.trace.back().mappedRhsNorm=mapped;
            const bool residual=std::isfinite(mapped) && mapped>=0. && mapped<1e-5;
            return std::string(update?(residual?"native_approx_both":"native_approx_update"):
                (residual?"native_approx_mapped_rhs":""));
        };
        auto result=CoupledLoadLineNewton(device,port.value,config,port,metric,std::move(hooks)).solve(
            unit_?coordinates.scaled(initial):initial,target/references.voltage,
            resistance*references.current/references.voltage,history);
        if(unit_)result=coordinates.physical(std::move(result));
        result.loadLineResidual_V=physicalLoad(result.state).voltage_V;
        return result;
    }

private:
    SparseMatrixd border(const ArclengthState& state,Real R,bool full) const {
        const int n=state.x.size();
        const auto J=full?system_.jacobian(state.x,state.lambda):
            (system_.postProtectionDerivativeRecovery?system_.protectedJacobian(state.x,state.lambda):
                system_.bootstrapJacobian(state.x,state.lambda));
        const auto f=system_.parameterDerivative(state.x,state.lambda);
        const auto port=port_.linearize(state.x,state.lambda);
        std::vector<Eigen::Triplet<Real>> entries;
        for(int k=0;k<J.outerSize();++k)for(SparseMatrixd::InnerIterator t(J,k);t;++t)entries.emplace_back(t.row(),t.col(),t.value());
        for(int i=0;i<n;++i) {
            if(f[i])entries.emplace_back(i,n,f[i]);
            if(port.stateDerivative[i])entries.emplace_back(n,i,R*port.stateDerivative[i]);
        }
        entries.emplace_back(n,n,config_.voltageCoefficient+R*port.parameterDerivative);
        SparseMatrixd A(n+1,n+1);A.setFromTriplets(entries.begin(),entries.end());return A;
    }
    static VectorXd auditRows(const SparseMatrixd& A,const VectorXd& b,const VectorXd& columns) {
        if(columns.size()!=b.size() || !columns.allFinite() || (columns.array()<=0.).any())
            throw std::invalid_argument("invalid physical audit columns");
        VectorXd scale=b.cwiseAbs();
        for(int k=0;k<A.outerSize();++k)for(SparseMatrixd::InnerIterator t(A,k);t;++t)
            scale[t.row()]=std::max(scale[t.row()],std::abs(t.value()*columns[t.col()]));
        if(!scale.allFinite() || (scale.array()<=0.).any())throw std::runtime_error("invalid physical audit row scales");
        return scale.cwiseInverse();
    }
    ArclengthSystem system_;
    ArclengthScalarFunctional port_;
    CoupledLoadLineConfig config_;
    DensityUpdateErrorEvaluator metric_;
    PrimaryVariableScales references_;
    bool unit_,approximateNativeStop_;
    ApproximateNativeRhsScales mapping_;
};

} // namespace vela
