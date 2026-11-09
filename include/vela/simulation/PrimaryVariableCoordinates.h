#pragma once

#include "vela/simulation/CoupledLoadLine.h"
#include "vela/solver/ConstraintElimination.h"
#include <cstring>

namespace vela {

/// Numerical coordinates only. Public states, protection and error references
/// retain their established physical/packed contract.
struct PrimaryVariableScales {
    VectorXd device;
    Real voltage = 1.;
    Real current = 1.;

    void validate(Eigen::Index size) const {
        if(device.size()!=size || !device.allFinite() || (device.array()<=0.).any() ||
           !std::isfinite(voltage) || voltage<=0. || !std::isfinite(current) || current<=0.)
            throw std::invalid_argument("invalid primary-variable Unit scales");
    }
};

struct ApproximateNativeRhsScales {
    Real poisson=0., carrier=0., load=0.;
    bool available() const {
        return std::isfinite(poisson)&&poisson>0. && std::isfinite(carrier)&&carrier>0. &&
            std::isfinite(load)&&load>0.;
    }
};

class PrimaryVariableCoordinates {
public:
    explicit PrimaryVariableCoordinates(PrimaryVariableScales scales):scales_(std::move(scales)) {
        scales_.validate(scales_.device.size());
    }
    VectorXd physical(const VectorXd& x) const { return x.cwiseProduct(scales_.device); }
    VectorXd scaled(const VectorXd& x) const { return x.cwiseQuotient(scales_.device); }
    ArclengthState physical(const ArclengthState& x) const { return {physical(x.x),x.lambda*scales_.voltage}; }
    ArclengthState scaled(const ArclengthState& x) const { return {scaled(x.x),x.lambda/scales_.voltage}; }
    Real voltage() const { return scales_.voltage; }
    Real current() const { return scales_.current; }
    SparseMatrixd columns(SparseMatrixd a,const VectorXd& d) const {
        for(int k=0;k<a.outerSize();++k)for(SparseMatrixd::InnerIterator t(a,k);t;++t)t.valueRef()*=d[t.col()];
        return a;
    }
    ArclengthSystem system(const ArclengthSystem& old) const {
        const auto b=*this;ArclengthSystem s=old;
        s.residual=[b,old](const VectorXd& x,Real v){return old.residual(b.physical(x),v*b.voltage());};
        s.parameterDerivative=[b,old](const VectorXd& x,Real v)->VectorXd{return old.parameterDerivative(b.physical(x),v*b.voltage())*b.voltage();};
        const auto jac=[b](auto original){return [b,original](const VectorXd& x,Real v){return b.columns(original(b.physical(x),v*b.voltage()),b.scales_.device);};};
        s.jacobian=jac(old.jacobian);
        if(old.bootstrapJacobian)s.bootstrapJacobian=jac(old.bootstrapJacobian);
        if(old.protectedJacobian)s.protectedJacobian=jac(old.protectedJacobian);
        const auto direct=[b](auto original){return [b,original](const VectorXd& x,Real v,const VectorXd& d){return original(b.physical(x),v*b.voltage(),b.physical(d));};};
        if(old.jacobianInColumns)s.jacobianInColumns=direct(old.jacobianInColumns);
        if(old.bootstrapJacobianInColumns)s.bootstrapJacobianInColumns=direct(old.bootstrapJacobianInColumns);
        if(old.protectedJacobianInColumns)s.protectedJacobianInColumns=direct(old.protectedJacobianInColumns);
        s.solveJacobian=[b,old](const VectorXd& x,Real v,const VectorXd& rhs,VectorXd& answer){VectorXd p;const bool ok=old.solveJacobian(b.physical(x),v*b.voltage(),rhs,p);if(ok)answer=b.scaled(p);return ok;};
        if(old.solveJacobianPair)s.solveJacobianPair=[b,old](const VectorXd& x,Real v,const VectorXd& a,const VectorXd& z,VectorXd& aa,VectorXd& zz){VectorXd p,q;const bool ok=old.solveJacobianPair(b.physical(x),v*b.voltage(),a,z,p,q);if(ok){aa=b.scaled(p);zz=b.scaled(q);}return ok;};
        if(old.limitUpdate)s.limitUpdate=[b,old](const VectorXd& x,VectorXd& dx,Real& dv){VectorXd p=b.physical(dx);Real v=dv*b.voltage();old.limitUpdate(b.physical(x),p,v);dx=b.scaled(p);dv=v/b.voltage();};
        if(old.makeTrial)s.makeTrial=[b,old](const ArclengthState& x,const VectorXd& dx,Real dv,Real alpha){auto trial=old.makeTrial(b.physical(x),b.physical(dx),dv*b.voltage(),alpha);trial.state=b.scaled(trial.state);return trial;};
        if(old.predictState)s.predictState=[b,old](const VectorXd& a,const VectorXd& z,Real r){return b.scaled(old.predictState(b.physical(a),b.physical(z),r));};
        if(old.prescribedRows)s.prescribedRows=[b,old](const VectorXd& x,Real v){return old.prescribedRows(b.physical(x),v*b.voltage());};
        if(old.convergenceSatisfied)s.convergenceSatisfied=[b,old](const VectorXd& x,Real v){return old.convergenceSatisfied(b.physical(x),v*b.voltage());};
        if(old.columnScales)s.columnScales=[b,old](const ArclengthState& x)->VectorXd{VectorXd d=old.columnScales(b.physical(x));d.head(b.scales_.device.size()).array()/=b.scales_.device.array();d[b.scales_.device.size()]/=b.voltage();return d;};
        return s;
    }
    ArclengthScalarFunctional port(const ArclengthScalarFunctional& old) const {
        const auto b=*this;ArclengthScalarFunctional p;
        p.value=[b,old](const VectorXd& x,Real v){return old.value(b.physical(x),v*b.voltage())/b.current();};
        p.linearize=[b,old](const VectorXd& x,Real v){auto a=old.linearize(b.physical(x),v*b.voltage());a.value/=b.current();a.stateDerivative.array()*=b.scales_.device.array()/b.current();a.parameterDerivative*=b.voltage()/b.current();return a;};
        if(old.linearizeInColumns)p.linearizeInColumns=[b,old](const VectorXd& x,Real v,const VectorXd& d){auto a=old.linearizeInColumns(b.physical(x),v*b.voltage(),b.physical(d));a.value/=b.current();a.stateDerivative/=b.current();a.parameterDerivative*=b.voltage()/b.current();return a;};
        return p;
    }
    DensityUpdateErrorEvaluator metric(const DensityUpdateErrorEvaluator& old) const {
        if(!old)return {};const auto b=*this;
        return [b,old](const VectorXd& before,const VectorXd& raw,const VectorXd& after,Real v,Real dv){return old(b.physical(before),b.physical(raw),b.physical(after),v*b.voltage(),dv*b.voltage());};
    }
    CoupledLoadLineResult physical(CoupledLoadLineResult result) const {
        result.state=physical(result.state);result.directedCurrent_A_per_um*=current();result.loadLineResidual_V*=voltage();
        for(auto& t:result.trace) {
            t.innerVoltage_V*=voltage();t.loadLineResidual_V*=voltage();t.proposedVoltageUpdate_V*=voltage();t.rawVoltageUpdate_V*=voltage();
            t.bestTrialLoadLineResidual_V*=voltage();t.candidateLoadResidual_V*=voltage();t.circuitLinearizationResidual*=voltage();
            t.analyticCurrentDirectionalDerivative*=current();t.finiteDifferenceCurrentDirectionalDerivative*=current();
            if(t.rawStateUpdate.size()==scales_.device.size()) {t.rawStateUpdate=physical(t.rawStateUpdate);t.rawStateUpdateNorm=t.rawStateUpdate.cwiseAbs().maxCoeff(&t.rawStateUpdateIndex);}
            for(auto& a:t.trials) {if(a.appliedStateUpdate.size()==scales_.device.size())a.appliedStateUpdate=physical(a.appliedStateUpdate);a.appliedVoltageUpdate_V*=voltage();a.predictedLoadResidual_V*=voltage();a.actualLoadResidual_V*=voltage();}
        }
        return result;
    }
private:
    PrimaryVariableScales scales_;
};

} // namespace vela
