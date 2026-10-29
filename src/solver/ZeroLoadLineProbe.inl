NewtonLoadLineIterationEvaluation NewtonSolver::evaluateZeroLoadLineIteration(
    const DDSolution& state,const std::string& contactName,
    Real resistance,Real currentScale,const std::string& backend) const
{
    for(const auto& [name,bias]:contactBiases_)
        if(bias!=0)throw std::invalid_argument("zero load-line probe requires all contacts at 0 V");
    return evaluateZeroSupplyLoadLineIteration(state,contactName,resistance,currentScale,backend);
}

NewtonLoadLineIterationEvaluation NewtonSolver::evaluateZeroSupplyLoadLineIteration(
    const DDSolution& state,const std::string& contactName,
    Real resistance,Real currentScale,const std::string& backend) const
{
    if(!contactBiases_.contains(contactName) ||
       !std::isfinite(contactBiases_.at(contactName)) ||
       std::abs(contactBiases_.at(contactName))>0.01)
        throw std::invalid_argument("zero-supply diagnostic internal voltage must stay within 10 mV");
    return evaluateBoundedSupplyLoadLineIteration(
        state,contactName,resistance,currentScale,0.,0.01,backend);
}

NewtonLoadLineIterationEvaluation NewtonSolver::evaluateBoundedSupplyLoadLineIteration(
    const DDSolution& state,const std::string& contactName,
    Real resistance,Real currentScale,Real supply,Real maxAbsVoltage,
    const std::string& backend,const std::string& trialBoundaryPolicy,
    bool traceWeakRows,bool pairedTrialAudit,bool constantRelativeCoordinates) const
{
    const bool projectBoundary=trialBoundaryPolicy=="rounded_voltage_projection";
    const bool splitVoltage=trialBoundaryPolicy=="split_contact_voltage";
    if(trialBoundaryPolicy!="linear_direction" && !projectBoundary && !splitVoltage)
        throw std::invalid_argument("invalid trial boundary policy");
    if(projectBoundary && !cfg_.splitDDState)
        throw std::invalid_argument("rounded voltage projection requires split state");
    if(pairedTrialAudit && !cfg_.splitDDState)
        throw std::invalid_argument("paired trial audit requires split state");
    if(constantRelativeCoordinates && (!pairedTrialAudit || !cfg_.splitDDState || cfg_.carrierRegularizationScale!=0))
        throw std::invalid_argument("constant-relative coordinates require a split paired audit with zero regularization");
    if(splitVoltage && (!cfg_.splitDDState || cfg_.contactBoundaryMinorityElectronRelaxation))
        throw std::invalid_argument("split contact voltage needs unrelaxed split Dirichlet targets");
    if(!splitVoltage && !cfg_.splitContactBiasLow_V.empty())
        throw std::invalid_argument("contact low parts require split_contact_voltage policy in circuit probe");
    if(!contactBiases_.contains(contactName) || !std::isfinite(currentScale))
        throw std::invalid_argument("invalid zero load-line terminal");
    if(!std::isfinite(supply) || std::abs(supply)>0.2 ||
       !std::isfinite(maxAbsVoltage) || maxAbsVoltage<=0 || maxAbsVoltage>0.2 ||
       !std::isfinite(resistance) || resistance<=0)
        throw std::invalid_argument("bounded supply diagnostic requires positive resistance and a voltage box <= 0.2 V");
    const Real voltage=contactBiases_.at(contactName);
    const Real voltageLow=cfg_.splitContactBiasLow_V.contains(contactName)?cfg_.splitContactBiasLow_V.at(contactName):0.;
    const numerics::SplitContactVoltage initialVoltage{{voltage,voltageLow}};
    const auto circuit=[&](Real hi,Real lo,Real id) {
        return splitVoltage?numerics::SplitContactVoltage{{hi,lo}}.loadLine(resistance,id,supply):hi+resistance*id-supply;
    };
    if(!std::isfinite(voltage) || std::abs(voltage)>maxAbsVoltage)
        throw std::invalid_argument("initial contact voltage exceeds bounded diagnostic voltage box");
    for(const auto& [name,bias]:contactBiases_)
        if(name!=contactName && bias!=0)
            throw std::invalid_argument("zero-supply diagnostic other terminals must be zero");
    auto assembler=makeArclengthAssembler(&state);
    auto& a=*assembler;
    const VectorXd x=packReferencedSolution(a,state,{});
    const int N=mesh_.numNodes();
    VectorXd weights=VectorXd::Zero(3*N);
    for(const auto& ct:mesh_.contacts())if(ct.name==contactName)
        for(auto node:ct.node_ids) {
            weights[N+node]-=constants::q;
            weights[2*N+node]+=constants::q;
        }
    const Real scale=currentScale*a.continuityResidualScale();
    const auto biases=[&](Real v) {auto b=contactBiases_;b[contactName]=v;return b;};
    const auto boundaryAt=[&](Real v,Real lo) {
        auto lows=cfg_.splitContactBiasLow_V;lows[contactName]=lo;
        return splitVoltage?buildBoundaryConditions(a,biases(v),&lows):buildBoundaryConditions(a,biases(v));
    };
    const auto bcs=buildBoundaryConditions(a,biases(voltage));
    const VectorXd f=a.residual(x,bcs);
    NewtonLoadLineIterationEvaluation out;
    out.deviceJacobian=a.assembleJacobian(x,bcs);
    const VectorXd unbounded=a.residual(x,{});
    const auto unboundedJ=a.assembleJacobian(x,{});
    const Real current=scale*weights.dot(unbounded);
    const VectorXd gradient=scale*(unboundedJ.transpose()*weights);
    constexpr Real h=1e-4;
    VectorXd fv=(a.residual(x,buildBoundaryConditions(a,biases(voltage+h)))-
        a.residual(x,buildBoundaryConditions(a,biases(voltage-h))))/(2*h);
    if(cfg_.splitDDState && !cfg_.contactBoundaryMinorityElectronRelaxation) {
        // Bias changes only the prescribed Ohmic/metal target, not its row Jacobian.
        // Use its exact derivative even when Vbias is below the built-in ULP.
        fv.setZero();const Real invScale=1/(a.usesScaledState()?a.potentialScale():1.);
        for(const auto& ct:mesh_.contacts())if(ct.name==contactName)
            for(auto i:ct.node_ids) {
                if(bcs.psi.contains(i))fv[i]=-invScale;
                if(bcs.phin.contains(i))fv[N+i]=-invScale;
                if(bcs.phip.contains(i))fv[2*N+i]=-invScale;
            }
    }
    const auto physical=assembleCarrierReducedLoadLine(out.deviceJacobian,f,fv,gradient,
        0,current,voltage,supply,resistance,N,splitVoltage?std::optional<Real>(circuit(voltage,voltageLow,current)):std::nullopt);
    out.physicalBordered=physical.full;out.physicalRhs=physical.rhs;
    SparseMatrixd regularized=out.deviceJacobian;
    addCarrierRowRegularization(regularized,N,cfg_.carrierRegularizationScale);
    SparseMatrixd solveJ=regularized;
    VectorXd solveRhs=-f;
    if(cfg_.globalContinuityClosure.linearSystemStabilization=="replace_free_row_sum") {
        std::vector<Index> ec,hc;
        for(const auto& [node,value]:bcs.phin)ec.push_back(node);
        for(const auto& [node,value]:bcs.phip)hc.push_back(node);
        // Apply exactly the same row transformation to the voltage column.
        // Transforming just the device block would change the coupled equation.
        stabilizeGlobalContinuityLinearSystem(solveJ,solveRhs,N,ec,hc);
        SparseMatrixd duplicate=regularized;
        stabilizeGlobalContinuityLinearSystem(duplicate,fv,N,ec,hc);
    } else if(cfg_.globalContinuityClosure.linearSystemStabilization!="off")
        throw std::invalid_argument("unsupported coupled diagnostic stabilization");
    out.linearSystem=assembleCarrierReducedLoadLine(solveJ,-solveRhs,fv,gradient,
        0,current,voltage,supply,resistance,N,splitVoltage?std::optional<Real>(circuit(voltage,voltageLow,current)):std::nullopt);
    auto runtime=a.splitRuntime();
    if(!constantRelativeCoordinates)out.rawDirection=out.linearSystem.solve(backend);
    else {
        const auto eligible=runtime->constantTransportDirection(x,VectorXd::Zero(3*N),bcs);
        std::vector<std::vector<int>> groups;std::vector<VectorXd> columns;
        std::vector<bool> ec(N,false),hc(N,false);
        for(const auto& [i,v]:bcs.phin)ec[i]=true;
        for(const auto& [i,v]:bcs.phip)hc[i]=true;
        const auto replacement=selectGlobalContinuityReplacementRows(out.deviceJacobian,N,ec,hc);
        for(int block=1;block<=2;++block) {
            VectorXd w=VectorXd::Zero(3*N);std::vector<int> group;
            for(const auto& row:eligible) {
                const int i=row.at("row");
                if(i/N!=block || out.linearSystem.map[i]<0)continue;
                // Circuit-current derivatives from constant material contacts
                // are not covered by this independent diagnostic.
                if(out.linearSystem.full.coeff(3*N,i)!=0)
                    throw std::runtime_error("Constant relative block reaches circuit reaction row");
                w[i]=1;group.push_back(out.linearSystem.map[i]);
            }
            if(group.size()<2)continue;
            VectorXd action=out.deviceJacobian*w;
            for(const auto& row:runtime->constantTransportDirection(x,w,bcs))
                action[row.at("row").get<int>()]=row.at("action").get<Real>();
            if(cfg_.globalContinuityClosure.linearSystemStabilization=="replace_free_row_sum") {
                long double en=0,hp=0;
                for(int i=0;i<N;++i){if(!ec[i])en+=action[N+i];if(!hc[i])hp+=action[2*N+i];}
                if(replacement[0]>=0)action[replacement[0]]=static_cast<Real>(en);
                if(replacement[1]>=0)action[replacement[1]]=static_cast<Real>(hp);
            }
            VectorXd reduced=VectorXd::Zero(out.linearSystem.reduced.rows());
            for(int i=0;i<3*N;++i)if(out.linearSystem.map[i]>=0)reduced[out.linearSystem.map[i]]=action[i];
            groups.push_back(group);columns.push_back(reduced);
        }
        if(groups.empty())throw std::runtime_error("No eligible constant-material coordinate group");
        const auto solved=solveCommonRelative(out.linearSystem.reduced,out.linearSystem.reducedRhs,groups,columns,backend);
        out.coordinateMatrix=solved.matrix;out.coordinateBackend=solved.backendMatrix;
        out.coordinateDirection=solved.direction;out.coordinateRowScale=solved.rowScale;
        out.coordinateColumnScale=solved.columnScale;out.coordinateScaledRhs=solved.scaledRhs;
        out.directionRelative=out.linearSystem.lift(solved.relative);
        out.directionShared=VectorXd::Zero(3*N+1);
        for(int i=0;i<3*N+1;++i)if(out.linearSystem.map[i]>=0)out.directionShared[i]=solved.shared[out.linearSystem.map[i]];
        out.rawDirection=out.directionRelative+out.directionShared;
    }
    VectorXd step=out.rawDirection.head(3*N);
    Real dv=out.rawDirection[3*N];
    const Real potentialScale=a.usesScaledState()?a.potentialScale():1;
    const Real cap=configuredQuasiFermiUniformStepScale(step,cfg_,N,potentialScale,doping_);
    step*=cap;dv*=cap;
    const VectorXd baseLow=runtime?runtime->low():VectorXd{};
    const auto inf=[](const VectorXd& q){return q.lpNorm<Eigen::Infinity>();};
    const auto merit=[&](const VectorXd& q,Real g){return std::max(inf(q)/cfg_.abstol,std::abs(g)/1e-6);};
    const Real initialMerit=merit(f,circuit(voltage,voltageLow,current));
    const auto initialTerms=traceWeakRows ? a.carrierContinuityEquationTermDiagnostics(x,bcs)
        : std::vector<CoupledDDCarrierTermDiagnostic>{};
    const auto initialRows=evaluateCarrierRowConvergence(initialTerms,cfg_.carrierRowConvergence);
    const VectorXd physicalPredictedChange=out.physicalBordered*(cap*out.rawDirection);
    auto weakSamples=[&](const VectorXd& xx,const CoupledDDBoundaryConditions& bc,Real alpha) {
        auto samples=nlohmann::json::array();
        if(!traceWeakRows || initialRows.violations.empty())return samples;
        const auto terms=a.carrierContinuityEquationTermDiagnostics(xx,bc);
        const VectorXd prediction=f+alpha*physicalPredictedChange.head(3*N);
        for(const auto& bad:initialRows.violations) {
            const auto it=std::find_if(terms.begin(),terms.end(),[&](const auto& t){return t.nodeId==bad.nodeId;});
            if(it==terms.end())throw std::runtime_error("missing tracked weak row");
            const bool electron=bad.carrier=="electron";
            const Real residual=electron?it->electronResidual:it->holeResidual;
            const Real flux=electron?it->electronFlux:it->holeFlux;
            const Real fluxAbs=electron?it->electronFluxAbsSum:it->holeFluxAbsSum;
            const Real source=electron?it->electronRecombination:it->holeRecombination;
            const Real impact=electron?it->electronImpact:it->holeImpact;
            const Real denominator=std::max({std::abs(fluxAbs),std::abs(source),std::abs(impact),cfg_.carrierRowConvergence.scaleFloor});
            const int row=(electron?N:2*N)+bad.nodeId;
            samples.push_back({{"node",bad.nodeId},{"carrier",bad.carrier},{"initial_residual",bad.residual},
                {"initial_scale",bad.scale},{"linear_prediction_residual",prediction[row]},
                {"residual",residual},{"flux",flux},{"flux_abs_sum",fluxAbs},
                {"recombination",source},{"impact",impact},{"ratio",denominator>0?std::abs(residual)/denominator:0.},
                {"density_internal_label_m3",electron?it->electronDensity_m3:it->holeDensity_m3}});
        }
        return samples;
    };
    auto trials=nlohmann::json::array();
    bool accepted=false;
    VectorXd finalX=x;
    Real finalV=voltage,finalVLow=voltageLow;
    const std::vector<std::string> policies=pairedTrialAudit && !splitVoltage
        ? std::vector<std::string>{"linear_direction","rounded_voltage_projection"}
        : std::vector<std::string>{trialBoundaryPolicy};
    for(int k=0;k<=12;++k) {
        const Real alpha=std::ldexp(1.,-k);
        for(const auto& policy:policies) {
        try {
            VectorXd trial;
            if(constantRelativeCoordinates) {
                trial=runtime->candidate(x,VectorXd(cap*out.directionShared.head(3*N)),alpha,baseLow);
                trial=runtime->candidate(trial,VectorXd(cap*out.directionRelative.head(3*N)),alpha,runtime->low());
            } else trial=runtime?runtime->candidate(x,step,alpha,baseLow):VectorXd(x+alpha*step);
            const auto splitV=initialVoltage.shifted(alpha,dv);
            const Real v=splitVoltage?splitV.value.hi:voltage+alpha*dv;
            const Real vLow=splitVoltage?splitV.value.lo:0.;
            if(std::abs(v)>maxAbsVoltage)throw std::runtime_error("trial exceeds bounded diagnostic voltage box");
            const auto trialBcs=boundaryAt(v,vLow);
            if(policy=="rounded_voltage_projection")runtime->applyBoundary(trial,trialBcs);
            const VectorXd highDelta=trial-x;
            const VectorXd fullDelta=runtime?runtime->coordinateDelta(trial,x,baseLow):highDelta;
            int highZeroFullNonzero=0;
            for(int i=0;i<3*N;++i)if(highDelta[i]==0 && fullDelta[i]!=0)++highZeroFullNonzero;
            const auto r=a.residual(trial,trialBcs);
            const auto id=scale*weights.dot(a.residual(trial,{}));
            const auto g=circuit(v,vLow,id);
            const auto m=merit(r,g);
            const bool ok=trial.allFinite() && std::isfinite(m) && m<initialMerit;
            Real boundaryInf=0;
            for(const auto& [i,z]:trialBcs.psi)boundaryInf=std::max(boundaryInf,std::abs(r[i]));
            for(const auto& [i,z]:trialBcs.phin)boundaryInf=std::max(boundaryInf,std::abs(r[N+i]));
            for(const auto& [i,z]:trialBcs.phip)boundaryInf=std::max(boundaryInf,std::abs(r[2*N+i]));
            trials.push_back({{"alpha",alpha},{"rhs_inf",inf(r)},{"load_line_V",g},{"Id_A_per_um",id},{"merit",m},{"descent_accepted",ok},
                {"requested_voltage_delta_V",alpha*dv},{"applied_voltage_delta_V",v-voltage},
                {"voltage_low_V",vLow},{"applied_voltage_low_delta_V",vLow-voltageLow},
                {"high_zero_full_nonzero_components",highZeroFullNonzero},
                {"actual_split_delta_inf_internal",inf(fullDelta)},
                {"actual_high_delta_inf_internal",inf(highDelta)},
                {"boundary_rhs_inf",boundaryInf},{"policy",policy},
                {"weak_rows",weakSamples(trial,trialBcs,alpha)}});
            if(pairedTrialAudit) {
                const auto allRows=evaluateCarrierRowConvergence(
                    a.carrierContinuityEquationTermDiagnostics(trial,trialBcs),cfg_.carrierRowConvergence);
                auto ledger=trials.back();ledger["index"]=out.auditCandidates.size();
                ledger["voltage_V"]=v;ledger["all_row_violations"]=allRows.violations.size();
                ledger["all_row_max_ratio"]=allRows.maxRatio;
                ledger["all_rows_satisfied"]=allRows.satisfied;
                ledger["bad_rows"]=nlohmann::json::array();
                for(const auto& bad:allRows.violations)
                    ledger["bad_rows"].push_back({{"node",bad.nodeId},{"carrier",bad.carrier},
                        {"residual",bad.residual},{"scale",bad.scale}});
                const auto raw=a.residual(trial,{});
                ledger["raw_reaction_ports_A_per_um"]=nlohmann::json::object();
                for(const auto& ct:mesh_.contacts()) {
                    long double en=0,hp=0;
                    for(auto node:ct.node_ids){en-=raw[N+node];hp+=raw[2*N+node];}
                    ledger["raw_reaction_ports_A_per_um"][ct.name]={
                        {"electron",static_cast<Real>(en*constants::q*scale)},
                        {"hole",static_cast<Real>(hp*constants::q*scale)},
                        {"total",static_cast<Real>((en+hp)*constants::q*scale)}};
                }
                ledger["qualified"]=false;
                out.auditCandidates.push_back(makeSolution(a,trial,0));
                out.auditCandidateReports.push_back(std::move(ledger));
            } else if(ok){finalX=trial;finalV=v;finalVLow=vLow;accepted=true;break;}
        } catch(const std::exception& e){trials.push_back({{"alpha",alpha},{"policy",policy},{"error",e.what()}});}
        }
        if(accepted)break;
    }
    if(!accepted && runtime)runtime->setLow(baseLow);
    out.trialSolution=makeSolution(a,finalX,accepted?1:0);
    const auto finalBcs=boundaryAt(finalV,finalVLow);
    const VectorXd finalResidual=a.residual(finalX,finalBcs);
    const VectorXd finalRaw=a.residual(finalX,{});
    const auto finalTerms=a.carrierContinuityEquationTermDiagnostics(finalX,finalBcs);
    const auto rows=evaluateCarrierRowConvergence(finalTerms,cfg_.carrierRowConvergence);
    nlohmann::json ports=nlohmann::json::object();
    for(const auto& ct:mesh_.contacts()) {
        long double sum=0;
        for(auto node:ct.node_ids)sum+=static_cast<long double>(finalRaw[2*N+node])-finalRaw[N+node];
        ports[ct.name]=static_cast<Real>(sum*constants::q*scale);
    }
    // One generic direction-difference check audits the contact row independently
    // of the matrix product used to construct it. Restore the same split tail.
    VectorXd auditDirection=VectorXd::Zero(3*N);
    for(int i=0;i<3*N;++i)auditDirection[i]=std::sin(.37*(i+1));
    auto fd=nlohmann::json::array();
    const VectorXd finalLow=runtime?runtime->low():VectorXd{};
    for(Real eps:{1e-6,5e-7}) {
        const VectorXd xp=runtime?runtime->candidate(x,auditDirection,eps,baseLow):VectorXd(x+eps*auditDirection);
        const Real ip=scale*weights.dot(a.residual(xp,{}));
        const VectorXd xm=runtime?runtime->candidate(x,auditDirection,-eps,baseLow):VectorXd(x-eps*auditDirection);
        const Real im=scale*weights.dot(a.residual(xm,{}));
        fd.push_back({{"scaled_step",eps},{"analytic",gradient.dot(auditDirection)},{"central_difference",(ip-im)/(2*eps)}});
    }
    if(runtime)runtime->setLow(finalLow);
    out.report={{"schema",supply==0. && maxAbsVoltage==0.01 ? "zero_supply_bordered_iteration_v2" : "bounded_supply_bordered_iteration_v1"},{"qualified",false},{"new_sweep_points",0},
        {"accepted_updates",accepted?1:0},{"iteration_attempts",1},{"backend",backend},{"split_state",cfg_.splitDDState},
        {"device_dimension",3*N},{"backend_dimension",out.linearSystem.reduced.rows()},
        {"eliminated_carrier_constraints",3*N+1-out.linearSystem.reduced.rows()},
        {"potential_scale_V",potentialScale},{"current_scale_A_per_um",currentScale},
        {"carrier_regularization_scale",cfg_.carrierRegularizationScale},
        {"linear_system_stabilization",cfg_.globalContinuityClosure.linearSystemStabilization},
        {"supply_voltage_V",supply},{"max_abs_device_voltage_V",maxAbsVoltage},{"resistance_ohm_um",resistance},
        {"trial_boundary_policy",trialBoundaryPolicy},{"trace_weak_rows",traceWeakRows},
        {"paired_trial_audit",pairedTrialAudit},{"linear_solve_count",1},
        {"constant_relative_coordinates",constantRelativeCoordinates},
        {"archived_candidate_count",out.auditCandidates.size()},
        {"tracked_initial_bad_rows",initialRows.violations.size()},
        {"initial_voltage_V",voltage},{"initial_voltage_low_V",voltageLow},{"initial_rhs_inf",inf(f)},{"initial_Id_A_per_um",current},{"initial_load_line_V",circuit(voltage,voltageLow,current)},
        {"final_load_line_V",circuit(finalV,finalVLow,scale*weights.dot(finalRaw))},
        {"initial_merit",initialMerit},{"uniform_cap",cap},{"raw_voltage_direction_V",out.rawDirection[3*N]},
        {"limited_voltage_direction_V",dv},{"final_voltage_V",finalV},{"final_voltage_low_V",finalVLow},{"final_rhs_inf",inf(finalResidual)},
        {"row_violations",rows.violations.size()},{"row_max_ratio",rows.maxRatio},{"rows_satisfied",rows.satisfied},
        {"ports_A_per_um",ports},{"contact_gradient_fd",fd},{"trials",trials},
        {"stop_reason",accepted?"one_update_diagnostic_limit":"no_descending_trial"}};
    if(pairedTrialAudit && runtime && !constantRelativeCoordinates) {
        runtime->setLow(baseLow);
        out.report["constant_transport_direction"]=runtime->constantTransportDirection(x,step,bcs);
        for(auto& row:out.report["constant_transport_direction"]) {
            const int i=row.at("row");row["initial_residual"]=f[i];row["matrix_action"]=physicalPredictedChange[i];
        }
        runtime->setLow(finalLow);
    }
    return out;
}
