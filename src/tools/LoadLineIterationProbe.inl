// Bounded diagnostic only: one coupled linearization and at most
// one accepted update. No sweep, no qualification or gate override.
nlohmann::json runLoadLineIterationProbe(const std::string& configFile,
                                        const nlohmann::json& cfg)
{
    const auto cfgDir=configDirectory(configFile);
    NewtonProblem problem=loadNewtonProblem(configFile,cfg);
    DDStateCliScope scope(problem,cfg,cfgDir);
    const auto state=readExternalState(cfgDir,cfg,problem.mesh.numNodes());
    const auto solver=makeNewtonSolver(problem);
    const std::string contact=canonicalContactName(problem.mesh,cfg.value("contact",std::string("drain")));
    const vela::Real voltage=problem.biases.at(contact);
    const vela::Real supply=cfg.value("supply_voltage_V",0.0);
    const bool bounded=cfg.at("simulation_type")=="newton_bounded_supply_load_line_iteration_probe";
    const bool followup=cfg.at("simulation_type")=="newton_zero_supply_load_line_iteration_probe";
    if (!bounded && ((!followup && voltage!=0) || supply!=0))
        throw std::invalid_argument("load-line iteration probe is restricted to zero bias");
    const vela::Real resistance=cfg.value("resistance_ohm_um",1e12);
    const auto& units=problem.newton.inputScaling.unitSystem();
    const vela::Real scale=units.internalCurrentPerDeviceDepthToAPerUm(
        units.currentDensityAM2PerInternal()*units.lengthMPerInternal());
    const std::string backend=cfg.value("bordered_linear_backend",std::string("pardiso_mkl"));
    const std::filesystem::path directory=resolvePath(cfgDir,cfg.at("linear_system_output_dir").get<std::string>());
    if (!std::filesystem::create_directory(directory)) throw std::runtime_error("probe output directory must be new");
    const auto evaluation=bounded
        ? solver.evaluateBoundedSupplyLoadLineIteration(state,contact,resistance,scale,supply,
            cfg.value("max_abs_device_voltage_V",0.2),backend,
            cfg.value("trial_boundary_policy",std::string("linear_direction")),
            cfg.value("trace_weak_rows",false),cfg.value("paired_trial_audit",false),
            cfg.value("constant_relative_coordinates",false))
        : followup
        ? solver.evaluateZeroSupplyLoadLineIteration(state,contact,resistance,scale,backend)
        : solver.evaluateZeroLoadLineIteration(state,contact,resistance,scale,backend);
    const auto& bordered=evaluation.linearSystem;
    const auto& j=evaluation.deviceJacobian;
    const auto matrix=[&](const char* name,const vela::SparseMatrixd& a) {
        std::ofstream out(directory/name); out.exceptions(std::ios::badbit|std::ios::failbit);
        out<<std::setprecision(17)<<"%%MatrixMarket matrix coordinate real general\n"
           <<a.rows()<<' '<<a.cols()<<' '<<a.nonZeros()<<'\n';
        for(int c=0;c<a.outerSize();++c)
            for(vela::SparseMatrixd::InnerIterator e(a,c);e;++e)
                out<<e.row()+1<<' '<<e.col()+1<<' '<<e.value()<<'\n';
    };
    matrix("device.mtx",j); matrix("physical_bordered.mtx",evaluation.physicalBordered);
    matrix("bordered.mtx",bordered.full);
    matrix("reduced.mtx",bordered.reduced); matrix("backend.mtx",bordered.equilibrated);
    if(evaluation.coordinateMatrix.rows()) {
        matrix("coordinate.mtx",evaluation.coordinateMatrix);matrix("coordinate_backend.mtx",evaluation.coordinateBackend);
        std::ofstream coords(directory/"coordinate_vectors.csv");
        coords<<std::setprecision(17)<<"index,rhs,direction,row_scale,column_scale,scaled_rhs\n";
        for(int i=0;i<evaluation.coordinateDirection.size();++i)
            coords<<i<<','<<bordered.reducedRhs[i]<<','<<evaluation.coordinateDirection[i]<<','<<evaluation.coordinateRowScale[i]<<','<<evaluation.coordinateColumnScale[i]<<','<<evaluation.coordinateScaledRhs[i]<<'\n';
        std::ofstream parts(directory/"direction_parts.csv");parts<<std::setprecision(17)<<"index,shared,relative\n";
        for(int i=0;i<evaluation.directionShared.size();++i)parts<<i<<','<<evaluation.directionShared[i]<<','<<evaluation.directionRelative[i]<<'\n';
    }
    const auto& direction=evaluation.rawDirection;
    const vela::VectorXd defect=bordered.full*direction-bordered.rhs;
    std::ofstream vectors(directory/"vectors.csv");
    vectors<<std::setprecision(17)<<"index,reduced_index,rhs,direction,linear_defect,physical_rhs\n";
    for(int i=0;i<direction.size();++i)
        vectors<<i<<','<<bordered.map[i]<<','<<bordered.rhs[i]<<','<<direction[i]<<','<<defect[i]<<','<<evaluation.physicalRhs[i]<<'\n';
    std::ofstream scales(directory/"scales.csv");
    scales<<std::setprecision(17)<<"index,rhs,row_scale,column_scale,scaled_rhs\n";
    for(int i=0;i<bordered.reducedRhs.size();++i)
        scales<<i<<','<<bordered.reducedRhs[i]<<','<<bordered.rowScale[i]<<','<<bordered.columnScale[i]<<','<<bordered.scaledRhs[i]<<'\n';
    const auto& finalState=evaluation.trialSolution;
    {
        auto metadata=*vela::activeDDStateArchiveMetadata();
        metadata["state_role"]=bounded ? "bounded_supply_bordered_iteration_diagnostic" : "zero_bias_bordered_iteration_diagnostic";
        metadata["qualified"]=false;
        metadata["bias_contact"]=contact;
        metadata["bias_V"]=evaluation.report.at("final_voltage_V");
        metadata["contact_biases_V"][contact]=evaluation.report.at("final_voltage_V");
        metadata["contact_biases_low_V"][contact]=evaluation.report.value("final_voltage_low_V",0.);
        vela::DDStateArchiveScope trialScope(std::move(metadata));
        vela::writeDDSolutionState(directory/"trial_state.h5",finalState,problem.newton.inputScaling);
    }
    auto finalProblem=problem;
    finalProblem.biases[contact]=evaluation.report.at("final_voltage_V").get<vela::Real>();
    if(cfg.value("trial_boundary_policy",std::string("linear_direction"))=="split_contact_voltage")
        finalProblem.newton.splitContactBiasLow_V[contact]=evaluation.report.value("final_voltage_low_V",0.);
    const auto finalSolver=makeNewtonSolver(finalProblem);
    auto summary=evaluation.report;
    if(evaluation.coordinateMatrix.rows()) {
        summary["coordinate_linear_defect_inf"]=(evaluation.coordinateMatrix*evaluation.coordinateDirection-bordered.reducedRhs).lpNorm<Eigen::Infinity>();
        const vela::VectorXd scaledY=evaluation.coordinateDirection.array()/evaluation.coordinateColumnScale.array();
        summary["coordinate_scaled_linear_defect_inf"]=(evaluation.coordinateBackend*scaledY-evaluation.coordinateScaledRhs).lpNorm<Eigen::Infinity>();
        summary["actual_backend_matrix_file"]="coordinate_backend.mtx";
        summary["linear_defect_scope"]="legacy matrix and rounded combined direction; coordinate defect separately reported";
    }
    summary["audit_candidates"]=nlohmann::json::array();
    for(std::size_t index=0;index<evaluation.auditCandidates.size();++index) {
        auto record=evaluation.auditCandidateReports.at(index);
        const std::string name="candidate_"+std::to_string(index)+".h5";
        auto metadata=*vela::activeDDStateArchiveMetadata();
        metadata["state_role"]="same_direction_rejected_or_unaccepted_candidate";
        metadata["qualified"]=false;
        metadata["bias_contact"]=contact;metadata["bias_V"]=record.at("voltage_V");
        metadata["contact_biases_V"][contact]=record.at("voltage_V");
        metadata["contact_biases_low_V"][contact]=record.value("voltage_low_V",0.);
        vela::DDStateArchiveScope candidateScope(std::move(metadata));
        vela::writeDDSolutionState(directory/name,evaluation.auditCandidates[index],problem.newton.inputScaling);
        record["state_file"]=name;summary["audit_candidates"].push_back(std::move(record));
    }
    summary["linear_defect_inf"]=defect.lpNorm<Eigen::Infinity>();
    summary["final_block_residuals"]=blockResidualsJson(finalSolver.evaluateResidual(finalState).blockNorms);
    summary["closure"]=globalContinuityClosureJson(finalSolver.evaluateFrozenGlobalContinuityClosure(finalState));
    std::ofstream out(directory/"summary.json"); out<<summary.dump(2)<<'\n';
    return summary;
}
