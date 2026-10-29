#include "vela/solver/SlotPopulationHold.h"
// Opt-in bounded SLOT population loop. Never changes production Newton defaults.
// Included inside the runner's private namespace to share its checked IO helpers.
nlohmann::json runSlotPopulationMemory(const std::string& configFile,const nlohmann::json& cfg) {
 using namespace vela;
 using Json=nlohmann::json;
 const auto started=std::chrono::steady_clock::now();
 auto seconds=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();};
 const auto dir=configDirectory(configFile);
 const auto opt=cfg.at("population_memory");
 const bool updateDiagnostics=opt.value("update_diagnostics",false);
 const bool activityDiagnostic=opt.value("diagnostic_population_active_set",false);
 const bool mobilityDiagnostic=opt.value("diagnostic_freeze_transport_qf_mobility",false);
 if(activityDiagnostic && mobilityDiagnostic)throw std::invalid_argument("separate diagnostic arms required");
 const auto holdOpt=opt.value("diagnostic_density_hold",Json::object());
 const bool densityHold=holdOpt.value("enabled",false);
 const double holdThreshold=holdOpt.value("threshold_cm3",1e-6);
 const int holdRepeated=holdOpt.value("consecutive_nonpositive",2);
 const bool continueAfterPrimary=opt.value("diagnostic_continue_after_primary",false);
 const auto captureSteps=opt.value("capture_steps",std::vector<int>{});
 if(densityHold && (activityDiagnostic || !mobilityDiagnostic || !(holdThreshold>0) || !std::isfinite(holdThreshold) || holdRepeated<2))
  throw std::invalid_argument("density hold requires mobility diagnostic, positive threshold and at least two proposals");
 if(continueAfterPrimary && !(activityDiagnostic||mobilityDiagnostic||densityHold))
  throw std::invalid_argument("continue after primary is diagnostic only");
 const bool diagnosticOnly=activityDiagnostic||mobilityDiagnostic||densityHold;
 const bool diagnosticStopOnPrimary=opt.value("diagnostic_stop_on_primary",false);
 if(diagnosticStopOnPrimary && !mobilityDiagnostic)
  throw std::invalid_argument("diagnostic primary stop requires transport mobility diagnostic");
 const bool projectionDiagnostics=opt.value("projection_diagnostics",false)||diagnosticOnly||updateDiagnostics;
 const int refinementSweeps=opt.value("linear_row_refinement_sweeps",0);
 const int factorSteps=opt.value("linear_factor_refinement_steps",0);
 const bool shared=opt.value("shared_diagnostics",false);
 const bool stableFiniteUpdate=opt.value("stable_finite_update",false);
 const bool diluteLogInverse=opt.value("dilute_log_inverse",false);
 const bool consistentPsiWriteback=opt.value("consistent_psi_writeback",false);
 if(consistentPsiWriteback && !stableFiniteUpdate)
  throw std::invalid_argument("consistent psi writeback requires stable finite update");
 const int sourceMode=opt.value("diagnostic_source_linearization",0);
 if(sourceMode<0 || sourceMode>2 || (sourceMode && (!shared || !mobilityDiagnostic)))
  throw std::invalid_argument("source diagnostic requires shared mobility pair, mode 0..2");
 const bool compact=opt.value("compact_diagnostics",false);
 if(refinementSweeps<0 || refinementSweeps>4)throw std::invalid_argument("linear_row_refinement_sweeps must be 0..4");
 if(factorSteps<0 || factorSteps>3 || (factorSteps && refinementSweeps))throw std::invalid_argument("factor refinement requires 0..3 and separate row arm");
 if(shared && !mobilityDiagnostic)throw std::invalid_argument("shared diagnostics requires mobility pair");
 const auto acceptanceName=opt.value("acceptance_policy",std::string("primary_only"));
 if(acceptanceName!="primary_only" && acceptanceName!="primary_with_guards")
  throw std::invalid_argument("unknown SLOT population acceptance policy");
 const bool reportOnly=acceptanceName=="primary_only";
 const auto acceptancePolicy=reportOnly ? slot_population::AcceptancePolicy::PrimaryOnly
                                       : slot_population::AcceptancePolicy::PrimaryWithGuards;
 const int maxUpdates=opt.value("max_updates",40);
 const double limit=opt.value("soft_seconds",1080.);
 if(maxUpdates<1 || maxUpdates>40 || !std::isfinite(limit) || limit<=0 || limit>1080)
  throw std::invalid_argument("population memory bounds: 1..40 updates, (0,1080] seconds");
 const std::filesystem::path out=resolvePath(dir,opt.at("output_directory").get<std::string>());
 if(!std::filesystem::create_directory(out))throw std::runtime_error("population memory output must be new");
 auto put=[&](const auto& p,const Json& v){std::ofstream f(p);f.exceptions(std::ios::badbit|std::ios::failbit);f<<std::setw(2)<<v<<'\n';};
 put(out/"config.json",cfg);
 NewtonProblem problem=loadNewtonProblem(configFile,cfg);
 DDStateCliScope scope(problem,cfg,dir);
 auto state=readExternalState(dir,cfg,problem.mesh.numNodes());const auto seed=state;
 const auto solver=makeNewtonSolver(problem);
 const auto& units=problem.newton.inputScaling.unitSystem();
 const int N=problem.mesh.numNodes();
 if(!problem.newton.inputScaling.isUnitScaling() || problem.newton.temperature_K!=300. ||
    problem.newton.abstol!=1e-10 || problem.newton.lineSearch ||
    problem.newton.quasiFermiUpdateLimit_V!=0 || problem.newton.quasiFermiUpdateLimitMinority_V!=0)
  throw std::invalid_argument("population memory requires validated SLOT 300 K unit scaling and uncapped update policy");
 if(!state.hasReferencedElectronQuasiFermi() || !state.hasReferencedHoleQuasiFermi() ||
    state.packedLow.size()!=0 || problem.newton.carrierStatistics.model!="fermi_dirac" ||
    problem.newton.carrierRegularizationScale!=0 ||
    problem.newton.globalContinuityClosure.linearSystemStabilization!="off")
  throw std::invalid_argument("population memory requires referenced Fermi state and disabled stabilization");
 // Preserve the established physical update-error population support explicitly.
 // Binding a support file prevents silently changing the RMS denominator.
 Json support;{std::ifstream f(resolvePath(dir,opt.at("support_file").get<std::string>()));f>>support;}
 put(out/"support.json",support);
 const auto ids=support.get<std::vector<std::vector<int>>>();
 if(ids.size()!=3 || ids[0].size()!=static_cast<size_t>(N))throw std::invalid_argument("invalid population support");
 for(const auto& group:ids){std::vector<bool> seen(N);for(int i:group){if(i<0||i>=N||seen[i])throw std::invalid_argument("invalid support node");seen[i]=true;}}
 std::vector<bool> contact(N),configured(N);
 for(const auto& ct:problem.mesh.contacts())for(auto i:ct.node_ids){contact[i]=true;if(problem.biases.count(ct.name))configured[i]=true;}
 const auto poisson=solver.evaluatePoissonTerms(state);
 double initialBoundaryResidual=0.;
 for(const auto& row:poisson.rows)if(row.constrained){
  if(!std::isfinite(row.productionResidual))throw std::invalid_argument("nonfinite initial boundary residual");
  initialBoundaryResidual=std::max(initialBoundaryResidual,std::abs(row.productionResidual));
 }
 if(!reportOnly && initialBoundaryResidual>=1e-10)
  throw std::invalid_argument("initial seed boundary residual guard");
 std::vector<double> dos[2];dos[0].resize(N);dos[1].resize(N);
 for(const auto& row:poisson.rows){dos[0][row.nodeId]=row.electronDensityOfStates;dos[1][row.nodeId]=row.holeDensityOfStates;}
 for(int b=0;b<2;++b)for(int i:ids[b+1])if(!(dos[b][i]>0))throw std::invalid_argument("population support includes zero DOS");
 auto terms=solver.evaluateCarrierTermDiagnostics(state,true);
 auto populations=[&](const auto& t){std::vector<VectorXd> v{VectorXd(N),VectorXd(N)};for(const auto& r:t.rows){v[0][r.nodeId]=r.electronDensity_m3;v[1][r.nodeId]=r.holeDensity_m3;}return v;};
 auto density=populations(terms);
 const auto initialDensity=density;
 std::vector<std::vector<bool>> latched(2,std::vector<bool>(N)),previousProjected=latched;
 std::vector<std::vector<int>> projectionStreak(2,std::vector<int>(N));
 std::vector<std::vector<SlotPopulationHoldState>> holdState(2,std::vector<SlotPopulationHoldState>(N));
 LinearSolver constrainedSolver;
 Json result={{"converged",false},{"qualified",false},{"history",Json::array()},{"role","SLOT bounded opt-in population solve"}};
 result["acceptance_policy"]=acceptanceName;
 result["stable_finite_update"]=stableFiniteUpdate;
 result["dilute_log_inverse"]=diluteLogInverse;
 result["consistent_psi_writeback"]=consistentPsiWriteback;
 result["diagnostic_density_hold"]=holdOpt;
 result["diagnostic_continue_after_primary"]=continueAfterPrimary;
 result["diagnostic_population_active_set"]=activityDiagnostic;
 result["diagnostic_freeze_transport_qf_mobility"]=mobilityDiagnostic;
 result["diagnostic_stop_on_primary"]=diagnosticStopOnPrimary;
 if(activityDiagnostic)result["role"]="diagnostic population constraint; never a production parent";
 if(mobilityDiagnostic)result["role"]="diagnostic transport QF mobility linearization; never a production parent";
 result["closure_boundary_report_only"]=reportOnly;
 result["initial_boundary_residual_max"]=initialBoundaryResidual;
 NewtonLinearSystemEvaluation linear;
 SparseMatrixd actualMatrix;
 VectorXd actualRhs,actualDirection,unconstrainedActualDirection;
 std::vector<int> capturedHeldRows;
 DDSolution matrixParentState=seed;
 auto snapshot=[&](const std::string& label,const DDSolution& s){
  writeDDSolutionState(out/(label+".h5"),s,problem.newton.inputScaling);
  if(linear.direction.size()==0)return;
  writeDDSolutionState(out/(label+"_matrix_parent.h5"),matrixParentState,problem.newton.inputScaling);
  auto d=out/(label+"_matrix");std::filesystem::create_directory(d);
  auto matrix=[&](const char* name,const auto& a){std::ofstream f(d/name);f.exceptions(std::ios::badbit|std::ios::failbit);f<<std::setprecision(17)<<"%%MatrixMarket matrix coordinate real general\n"<<a.rows()<<' '<<a.cols()<<' '<<a.nonZeros()<<'\n';
   for(int j=0;j<a.outerSize();++j)for(SparseMatrixd::InnerIterator it(a,j);it;++it)f<<it.row()+1<<' '<<it.col()+1<<' '<<it.value()<<'\n';};
  if(compact){
   matrix("stabilized.mtx",linear.solveMatrix);
   // This entry already rejects both forms of stabilization. Preserve explicit
   // aliases instead of writing three identical MatrixMarket files.
   put(d/"matrix_aliases.json",{{"physical.mtx","stabilized.mtx"},{"regularized.mtx","stabilized.mtx"}});
  }else{matrix("physical.mtx",linear.physicalJacobian);matrix("regularized.mtx",linear.regularizedJacobian);matrix("stabilized.mtx",linear.solveMatrix);}
  if(actualDirection.size()==linear.direction.size()){
   matrix("actual.mtx",actualMatrix);
   std::ofstream af(d/"actual_vectors.csv");af<<std::setprecision(17)<<"index,solve_rhs,direction_internal\n";
   std::ofstream uf(d/"unconstrained_actual_direction.csv");uf<<std::setprecision(17)<<"index,direction_internal\n";
   for(int i=0;i<actualDirection.size();++i){af<<i<<','<<actualRhs[i]<<','<<actualDirection[i]<<'\n';uf<<i<<','<<unconstrainedActualDirection[i]<<'\n';}
   put(d/"held_rows.json",capturedHeldRows);
  }
  std::ofstream f(d/"vectors.csv");f.exceptions(std::ios::badbit|std::ios::failbit);f<<std::setprecision(17)<<"index,residual,solve_rhs,direction_internal\n";
  for(int i=0;i<linear.direction.size();++i)f<<i<<','<<linear.residual[i]<<','<<linear.solveRhs[i]<<','<<linear.direction[i]<<'\n';
 };
 result["initialization_s"]=seconds();
 result["diagnostic_source_linearization"]=sourceMode;
 try {
 for(int step=1;step<=maxUpdates;++step){
  if(seconds()>=limit){result["stop_reason"]="soft_time_limit";break;}
  const auto begin=seconds();matrixParentState=state;Json timing=Json::object();actualDirection.resize(0);capturedHeldRows.clear();
  NewtonLinearSystemEvaluation paired;
  PerformanceProfiler linearProfile({true,"unused.json"});
  {ActivePerformanceProfilerScope profileScope(&linearProfile);linear=solver.evaluateSourceDiagnosticPair(state,shared?&paired:nullptr,factorSteps,sourceMode);}
  timing[shared?"paired_systems_s":"original_system_s"]=seconds()-begin;
  for(const auto& [name,value]:linear.timings)timing[name]=value;
  const double vt=linear.potentialScale;
  auto qt=seconds();
  const auto refinement=slot_population::refineMixedRows(linear.solveMatrix,linear.solveRhs,linear.direction,vt,N,refinementSweeps);
  auto quality=slot_population::linearQuality(linear.solveMatrix,linear.solveRhs,linear.direction,vt,N);
  timing["original_quality_refinement_s"]=seconds()-qt;
  Json q={{"accepted",quality.accepted},{"mixed_ratio",quality.mixedRatio},{"component_backward",quality.componentBackward},{"constraint_error_V",quality.constraintV},{"block_normwise",quality.blockNormwise}};
  q["row_refinement"]={{"before",refinement.before.mixedRatio},{"after",refinement.after.mixedRatio},{"sweeps",refinement.sweeps},{"corrected_rows",refinement.correctedRows}};
  const auto observations=linearProfile.toJson().at("observations");
  q["profile"]=linearProfile.toJson();
  for(const char* key:{"linear.eliminated_dofs","linear.reduced_dimension","linear.retained_prefix","linear.mkl_refinement_steps","linear.mkl_perturbed_pivots","linear.mkl_threads"})
   if(observations.contains(key))q[key]=observations.at(key).at("last");
  if(shared)q["backend_observation_scope"]="last solve of shared pair";
  if(!quality.accepted){result["linear_quality"]=q;throw std::runtime_error("mixed_linear_guard");}
  // Always keep the original direction for U and for matrix evidence. A held
  // equation is replaced only in this separate diagnostic direction solve.
  VectorXd usedDirection=linear.direction;actualMatrix=linear.solveMatrix;actualRhs=linear.solveRhs;
  if(mobilityDiagnostic){
   auto dt=seconds();
   auto diagnosticProblem=problem;
   diagnosticProblem.newton.diagnosticFreezeTransportQfMobility=true;
   auto altered=shared?std::move(paired):makeNewtonSolver(diagnosticProblem).evaluateLinearSystem(state,nullptr,factorSteps);
   if(!shared)timing["diagnostic_system_s"]=seconds()-dt;
   dt=seconds();
   const bool sameResidual=(altered.residual.array()==linear.residual.array()).all() &&
                           (altered.solveRhs.array()==linear.solveRhs.array()).all();
   const auto mr=slot_population::refineMixedRows(altered.solveMatrix,altered.solveRhs,altered.direction,vt,N,refinementSweeps);
   const auto mq=slot_population::linearQuality(altered.solveMatrix,altered.solveRhs,altered.direction,vt,N);
   q["mobility_diagnostic"]={{"rhs_bitwise_equal",sameResidual},{"accepted",mq.accepted},{"mixed_ratio",mq.mixedRatio},{"block_normwise",mq.blockNormwise}};
   q["mobility_diagnostic"]["row_refinement"]={{"before",mr.before.mixedRatio},{"after",mr.after.mixedRatio},{"sweeps",mr.sweeps},{"corrected_rows",mr.correctedRows}};
   q["mobility_diagnostic"]["system_metrics"]=altered.timings;
   timing["diagnostic_quality_refinement_s"]=seconds()-dt;
   if(!sameResidual || !mq.accepted){result["linear_quality"]=q;throw std::runtime_error("diagnostic_mobility_linear_guard");}
   usedDirection=altered.direction;actualMatrix=std::move(altered.solveMatrix);actualRhs=altered.solveRhs;
  }
  unconstrainedActualDirection=usedDirection;
  std::vector<int> heldRows;int heldCount[2]={0,0},released[2]={0,0};
  std::vector<VectorXd> rawRelative{VectorXd::Zero(N),VectorXd::Zero(N)};
  auto unconstrainedRelative=rawRelative;
  for(int b=0;b<2;++b)for(int i:ids[b+1]){
   const double eta=slot_population::populationEta(density[b][i],dos[b][i]);
   const double g=slot_population::populationLogSlope(eta);
   rawRelative[b][i]=(b==0?-1.:1.)*g*(linear.direction[(b+1)*N+i]-linear.direction[i]);
   unconstrainedRelative[b][i]=(b==0?-1.:1.)*g*(unconstrainedActualDirection[(b+1)*N+i]-unconstrainedActualDirection[i]);
   if(densityHold && !contact[i]){
    const double proposal=unconstrainedRelative[b][i];
    if(!std::isfinite(proposal))throw std::runtime_error("nonfinite density hold proposal");
    auto& hs=holdState[b][i];const int oldReleases=hs.releases;
    hs.observe(units.internalConcentrationToM3(density[b][i])*1e-6,holdThreshold,proposal,holdRepeated);
    released[b]+=hs.releases-oldReleases;latched[b][i]=hs.held;
    if(hs.held){heldRows.push_back((b+1)*N+i);++heldCount[b];}
   }
   if(activityDiagnostic && !contact[i]){
    const bool hold=slot_population::retainPopulationHold(latched[b][i],rawRelative[b][i]);
    if(latched[b][i]&&!hold)++released[b];latched[b][i]=hold;
    if(hold){heldRows.push_back((b+1)*N+i);++heldCount[b];}
   }
  }
  std::sort(heldRows.begin(),heldRows.end());
  if(!heldRows.empty()){
   const auto a=slot_population::populationConstraintMatrix(actualMatrix,heldRows,N);
   auto b=actualRhs;for(int i:heldRows)b[i]=0.;
   PerformanceProfiler heldProfile({true,"unused.json"});
   {ActivePerformanceProfilerScope profileScope(&heldProfile);usedDirection=constrainedSolver.solve(a,b);}
   const auto cr=slot_population::refineWithFactors(a,b,usedDirection,vt,N,factorSteps,[&](const VectorXd& rhs){return constrainedSolver.solve(a,rhs);});
   const auto cq=slot_population::linearQuality(a,b,usedDirection,vt,N);
   actualMatrix=a;actualRhs=b;
   q["constrained_factor_refinement_steps"]=cr.sweeps;
   double constraint=0.;for(int i:heldRows)constraint=std::max(constraint,std::abs(vt*(usedDirection[i]-usedDirection[i%N])));
   q["constrained"]={{"accepted",cq.accepted},{"mixed_ratio",cq.mixedRatio},{"block_normwise",cq.blockNormwise},{"population_constraint_V",constraint},{"observations",heldProfile.toJson().at("observations")}};
   actualDirection=usedDirection;capturedHeldRows=heldRows;
   if(!cq.accepted || constraint>1e-14){result["linear_quality"]=q;throw std::runtime_error("diagnostic_constrained_linear_guard");}
  }
  population_step::FractionToBoundary populationBound;
  if(problem.newton.populationStepFraction>0.) {
   VectorXd rn=VectorXd::Zero(N),rp=VectorXd::Zero(N);std::vector<bool> active(N);
   for(int i=0;i<N;++i)if(!contact[i])for(int b=0;b<2;++b) {
    const bool enabled=b==0?terms.rows[i].electronContinuityActive:terms.rows[i].holeContinuityActive;
    if(!enabled||!(dos[b][i]>0.))continue;
    const double eta=slot_population::populationEta(density[b][i],dos[b][i]);
    const double rel=(b==0?-1.:1.)*slot_population::populationLogSlope(eta)*
     (usedDirection[(b+1)*N+i]-usedDirection[i]);
    (b==0?rn:rp)[i]=rel;active[i]=true;
   }
   populationBound=population_step::relativeFractionToBoundary(rn,rp,active,problem.newton.populationStepFraction);
   usedDirection*=populationBound.alpha;
  }
  actualDirection=usedDirection;capturedHeldRows=heldRows;
  const double matrixTime=seconds()-begin;
  state.packedState.resize(0);state.packedLow.resize(0);
  std::vector<VectorXd> delta{vt*linear.direction.head(N),VectorXd::Zero(N),VectorXd::Zero(N)};
  VectorXd applied=vt*usedDirection.head(N);int damped=0,projected[2]={0,0},overlap[2]={0,0};double excludedDirectionMax=0.;
  if(problem.newton.populationStepFraction==0.)for(int i=0;i<N;++i)if(std::abs(applied[i])>.0259){applied[i]=std::copysign(.0259*std::log1p(std::abs(applied[i])/.0259),applied[i]);++damped;}
  const VectorXd psiBefore=state.psi;
  std::vector<VectorXd> finiteQf,incrementWritten,qfHighWritten;
  if(updateDiagnostics){
   finiteQf.assign(2,VectorXd::Zero(N));incrementWritten=finiteQf;qfHighWritten=finiteQf;
  }
  state.psi+=applied;
  const VectorXd writtenPsi=state.psi-psiBefore;
  const double psiWritebackError=(writtenPsi-applied).cwiseAbs().maxCoeff();
  int stableFiniteCount=0;double finiteUpdateChange=0.;
  int diluteLogSelected=0,diluteLogApplied=0;
  auto targets=density;
  std::vector<VectorXd> usedRelative{VectorXd::Zero(N),VectorXd::Zero(N)};
  std::vector<std::vector<bool>> projectedNow(2,std::vector<bool>(N));
  std::vector<std::vector<bool>> updated(2,std::vector<bool>(N));
  for(int b=0;b<2;++b){const double sign=b==0?-1.:1.;
   for(int i:ids[b+1]){
    double eta=slot_population::populationEta(density[b][i],dos[b][i]);
    double g=slot_population::populationLogSlope(eta);
    delta[b+1][i]=sign*density[b][i]*g*(linear.direction[(b+1)*N+i]-linear.direction[i]);
   }
   for(int i=0;i<N;++i){const auto& row=terms.rows[i];
    if(row.nodeId!=i)throw std::runtime_error("carrier diagnostic node order");
    bool active=b==0?row.electronContinuityActive:row.holeContinuityActive;
    if(!active || dos[b][i]<=0 || contact[i]){
     excludedDirectionMax=std::max(excludedDirectionMax,std::abs(vt*linear.direction[(b+1)*N+i]));
     if(!reportOnly && excludedDirectionMax>1e-12)throw std::runtime_error("excluded QF direction guard");continue;}
    double rel=delta[b+1][i]/density[b][i];
    if(diagnosticOnly||problem.newton.populationStepFraction>0.){
     const double eta=slot_population::populationEta(density[b][i],dos[b][i]);
     const double g=slot_population::populationLogSlope(eta);
     // Match the baseline's arithmetic order as well as its formula. The
     // first unconstrained step must be bitwise identical in paired arms.
     const double populationDelta=sign*density[b][i]*g*(usedDirection[(b+1)*N+i]-usedDirection[i]);
     rel=latched[b][i]?0.:(problem.newton.populationStepFraction>0.
      ?sign*g*(usedDirection[(b+1)*N+i]-usedDirection[i]):populationDelta/density[b][i]);
    }
    usedRelative[b][i]=rel;
    if(problem.newton.populationStepFraction>0.&&!(1+rel>0.))throw std::runtime_error("global population step lost positivity");
    if(1+rel<=0){rel=-.999;++projected[b];projectedNow[b][i]=true;
     if(previousProjected[b][i])++overlap[b];
     if(activityDiagnostic)latched[b][i]=true;
    }
    const bool logSelected=diluteLogInverse && slot_population::usesDiluteLogInversePopulation(density[b][i],dos[b][i],rel);
    if(logSelected){++diluteLogSelected;++diluteLogApplied;}
    const double de=slot_population::relativeInversePopulation(density[b][i],dos[b][i],rel,diluteLogInverse);
    double dq=applied[i]+sign*vt*de;
    if(stableFiniteUpdate){
     const double psiStep=consistentPsiWriteback?writtenPsi[i]:applied[i];
     const double xp=usedDirection[i],xq=usedDirection[(b+1)*N+i];
     const double t=sign*(xq-xp);double repaired;
     if(!latched[b][i] && !projectedNow[b][i] && std::abs(t)<=1e-4){
      const double eta=slot_population::populationEta(density[b][i],dos[b][i]);
      repaired=slot_population::tangentQuasiFermiStep(eta,sign,vt,xp,xq,psiStep);
      ++stableFiniteCount;
      if(logSelected)--diluteLogApplied; // Existing small-tangent arithmetic wins.
     }else repaired=(double)((long double)psiStep+(long double)sign*vt*de);
     finiteUpdateChange=std::max(finiteUpdateChange,std::abs(repaired-dq));dq=repaired;
    }
    auto& inc=b==0?state.phinIncrement:state.phipIncrement;auto& qf=b==0?state.phin:state.phip;
    const double incBefore=updateDiagnostics?inc[i]:0.,qfBefore=updateDiagnostics?qf[i]:0.;
    inc[i]+=dq;qf[i]=(b==0?state.electronQuasiFermiReferenceAt(i):state.holeQuasiFermiReferenceAt(i))+inc[i];
    if(updateDiagnostics){
     finiteQf[b][i]=dq;
     incrementWritten[b][i]=static_cast<double>(update_diagnostics::coordinateChange(0.,incBefore,0.,inc[i]));
     qfHighWritten[b][i]=static_cast<double>((long double)qf[i]-qfBefore);
    }
    targets[b][i]*=1+rel;(b==0?state.n:state.p)[i]=targets[b][i];updated[b][i]=true;
   }
  }
  if(!state.psi.allFinite()||!state.phin.allFinite()||!state.phip.allFinite()||!state.n.allFinite()||!state.p.allFinite())throw std::runtime_error("nonfinite state");
  const double updateTime=seconds()-begin-matrixTime;
  PerformanceProfiler evaluationProfile({true,"unused.json"});
  ActivePerformanceProfilerScope evaluationScope(&evaluationProfile);
  auto metricStart=seconds();
  auto next=solver.evaluateCarrierTermDiagnostics(state,true);auto nd=populations(next);
  timing["carrier_evaluation_s"]=seconds()-metricStart;metricStart=seconds();
  const auto closure=solver.evaluateFrozenGlobalContinuityClosure(state);
  timing["closure_evaluation_s"]=seconds()-metricStart;
  if(!closure.enabled)throw std::runtime_error("closure diagnostics required");
  double mapping=0;Json mappingWorst=Json::object();bool boundary=true;long double errorSum=0;size_t count=0;
  for(int b=0;b<3;++b)for(int i:ids[b]){
   const double value=b==0?state.psi[i]:nd[b-1][i];
   // Same physical concentration error reference as the Python driver (cm^-3).
   const double ref=b==0?.025851995266484913:1e16/units.concentrationM3PerInternal();
   const double v=delta[b][i]/(std::abs(value)+ref)*1e5;errorSum+=(long double)v*v;++count;
   if(b>0 && !(value>0))throw std::runtime_error("population positivity guard");
  }
  for(int i=0;i<N;++i){
   if(configured[i] && std::abs(state.psi[i]-seed.psi[i])>1e-10)boundary=false;
   if(contact[i] && (state.phinIncrement[i]!=seed.phinIncrement[i] || state.phipIncrement[i]!=seed.phipIncrement[i]))boundary=false;
   for(int b=0;b<2;++b)if(updated[b][i]){
    const double err=std::abs(nd[b][i]/targets[b][i]-1);
    if(!std::isfinite(err) || err>mapping){mapping=err;mappingWorst={{"node",i},{"carrier",b==0?"electron":"hole"},{"relative_error",err},{"density_before_internal",density[b][i]},{"target_internal",targets[b][i]},{"rebuilt_internal",nd[b][i]},{"dos_internal",dos[b][i]},{"density_before_cm3",units.internalConcentrationToM3(density[b][i])*1e-6},{"target_cm3",units.internalConcentrationToM3(targets[b][i])*1e-6},{"rebuilt_cm3",units.internalConcentrationToM3(nd[b][i])*1e-6},{"used_relative",usedRelative[b][i]},{"held",latched[b][i]}};}
   }
  }
  if(!(mapping<1e-7)){
   result["inverse_guard"]={{"step",step},{"worst",mappingWorst},{"linear_quality",q},{"held",{heldCount[0],heldCount[1]}},{"projected",{projected[0],projected[1]}}};
   put(out/"inverse_failure.json",result["inverse_guard"]);
   writeNewtonCarrierTermProbeCsv(out/"inverse_failure_rows.csv",problem.mesh,problem.doping,next,1.,1.);
   throw std::runtime_error("population inverse guard");
  }
  const double error=std::sqrt((double)(errorSum/count)),rho=next.residual.blockNorms.combined;
  bool qualified=slot_population::accepts(rho,error,closure.satisfied,boundary,acceptancePolicy);
  const bool originalPolicySatisfied=qualified;
  if(diagnosticOnly)qualified=false; // A diagnostic direction is not a root certificate.
  Json rowCounts=Json::array();
  for(int b=0;b<2;++b){int bad=0;double worst=0;for(const auto& r:next.rows){if(!(b==0?r.electronContinuityActive:r.holeContinuityActive))continue;
   const double scale=std::max({b==0?r.electronFluxAbsSum:r.holeFluxAbsSum,std::abs(b==0?r.electronRecombination:r.holeRecombination),std::abs(b==0?r.electronImpact:r.holeImpact),1e-30});
   const double ratio=std::abs(b==0?r.electronResidual:r.holeResidual)/scale;worst=std::max(worst,ratio);if(ratio>.001)++bad;}
   rowCounts.push_back({{"failed",bad},{"max_ratio",worst}});}
  Json h={{"step",step},{"rho",rho},{"raw_error",error},{"qualified",qualified},{"boundary_preserved",boundary},{"mapping_error",mapping},{"rows",rowCounts},{"linear_quality",q},{"closure",globalContinuityClosureJson(closure)},{"projected",{projected[0],projected[1]}},{"psi_damped",damped},{"matrix_s",matrixTime},{"update_s",updateTime},{"wall_s",seconds()-begin}};
  Json totals=Json::array(),holdMetrics=Json::array();
  for(int b=0;b<2;++b){
   long double flux=0.,recombination=0.,impact=0.,residual=0.;double heldWorst=0.;int heldBad=0,activations=0,releaseTotal=0;
   for(const auto& r:next.rows){
    const bool active=b==0?r.electronContinuityActive:r.holeContinuityActive;if(!active)continue;
    const double f=b==0?r.electronFlux:r.holeFlux;
    const double sr=b==0?r.electronRecombination:r.holeRecombination;
    const double im=b==0?r.electronImpact:r.holeImpact;
    const double rr=b==0?r.electronResidual:r.holeResidual;
    flux+=f;recombination+=sr;impact+=im;residual+=rr;
    if(std::binary_search(heldRows.begin(),heldRows.end(),(b+1)*N+r.nodeId)){
     const double scale=std::max({b==0?r.electronFluxAbsSum:r.holeFluxAbsSum,std::abs(sr),std::abs(im),1e-30});
     const double ratio=std::abs(rr)/scale;heldWorst=std::max(heldWorst,ratio);if(ratio>.001)++heldBad;
    }
   }
   for(int i:ids[b+1]){activations+=holdState[b][i].activations;releaseTotal+=holdState[b][i].releases;}
   totals.push_back({{"active_row_flux_sum",(double)flux},{"active_row_recombination_sum",(double)recombination},{"active_row_impact_sum",(double)impact},{"active_row_residual_sum",(double)residual}});
   holdMetrics.push_back({{"held_failed_rows",heldBad},{"held_max_true_ratio",heldWorst},{"cumulative_activations",activations},{"cumulative_releases",releaseTotal}});
  }
  h["physical_active_row_sums_internal"]=totals;h["density_hold_metrics"]=holdMetrics;
  if(problem.newton.populationStepFraction>0.)h["global_population_step"]={
   {"fraction",problem.newton.populationStepFraction},{"alpha",populationBound.alpha},
   {"limiting_node",populationBound.limitingNode},{"limiting_carrier",populationBound.limitingCarrier},
   {"raw_U_policy_unchanged",true},{"local_projection_used",false}};
  h["mapping_worst"]=mappingWorst;
  if(updateDiagnostics){
   update_diagnostics::Rms total;Json blocks=Json::array();
   std::ofstream f(out/("update_stages_"+std::to_string(step)+".csv"));
   f.exceptions(std::ios::badbit|std::ios::failbit);
   f<<std::setprecision(17)<<"block,node_id,raw_normalized,actual_normalized,denominator,raw_delta,actual_delta,used_qf_direction_V,finite_qf_requested_V,qf_increment_written_V,qf_physical_high_written_V,qf_writeback_error_V,updated,projected,held,target_density_internal,rebuilt_density_internal\n";
   for(int b=0;b<3;++b){update_diagnostics::Rms block;
    for(int i:ids[b]){
     const double value=b==0?state.psi[i]:nd[b-1][i];
     const double ref=b==0?.025851995266484913:1e16/units.concentrationM3PerInternal();
     const double den=std::abs(value)+ref;
     const long double actual=b==0?(long double)state.psi[i]-psiBefore[i]:(long double)nd[b-1][i]-density[b-1][i];
     block.add(delta[b][i],actual,den,i);
     const double requested=b==0?applied[i]:finiteQf[b-1][i];
     const double written=b==0?(double)actual:incrementWritten[b-1][i];
     f<<b<<','<<i<<','<<delta[b][i]/den*1e5<<','<<(double)(actual/den*1e5)<<','<<den<<','<<delta[b][i]<<','<<(double)actual<<','
      <<vt*usedDirection[b*N+i]<<','<<requested<<','<<written<<','<<(b==0?(double)actual:qfHighWritten[b-1][i])<<','
      <<(double)((long double)written-requested)<<','<<(b==0||updated[b-1][i])<<','<<(b>0&&projectedNow[b-1][i])<<','<<(b>0&&latched[b-1][i])<<','
      <<(b==0?0.:targets[b-1][i])<<','<<(b==0?0.:nd[b-1][i])<<'\n';
    }
    total.merge(block);blocks.push_back(block.json());
   }
   h["update_diagnostic"]=total.json();h["update_diagnostic"]["blocks"]=blocks;
   h["update_diagnostic"]["raw_rms_recompute_difference"]=total.json()["raw_rms"].get<double>()-error;
   h["update_diagnostic"]["same_denominator_and_support"]=true;
   h["update_diagnostic"]["units"]="psi and QF: V; population and denominator: internal concentration";
   h["update_diagnostic"]["concentration_to_cm3"]=units.concentrationM3PerInternal()*1e-6;
   h["update_diagnostic"]["actual_density_note"]="regenerated binary64 populations; unresolved changes are not proof of convergence";
   h["convergence_report"]={{"schema","vela.convergence_report.v1"},
     {"numerical_stop",update_diagnostics::stopLayers(rho,error,originalPolicySatisfied)},
     {"device_checks",{{"closure",globalContinuityClosureJson(closure)},{"boundary_preserved",boundary},
       {"all_active_rows",rowCounts},{"row_report_is_production_qualification",false},
       {"load_line",{{"evaluated",false},{"reason","owned by external circuit controller"}}},
       {"file_readback",{{"evaluated",false},{"reason","owned by restart controller"}}}}},
     {"numerical_guards",{{"mapping_error",mapping},{"linear_quality",q}}}};
   h["convergence_report"]["device_checks"]["configured_row_gate"]=
    writeWeakRowLedger(out/("weak_rows_"+std::to_string(step)+".csv"),problem,next);
  }
  h["finite_update_arithmetic"]={{"stable_nodes",stableFiniteCount},{"max_qf_change_from_legacy_V",finiteUpdateChange},{"max_psi_writeback_error_V",psiWritebackError},{"consistent_psi_writeback",consistentPsiWriteback}};
  h["finite_update_arithmetic"]["dilute_log_selected"]=diluteLogSelected;
  h["finite_update_arithmetic"]["dilute_log_applied"]=diluteLogApplied;
  h["original_policy_satisfied"]=originalPolicySatisfied;h["held"]={heldCount[0],heldCount[1]};h["released"]={released[0],released[1]};h["projection_overlap"]={overlap[0],overlap[1]};
  const bool detail=!compact || originalPolicySatisfied || step==maxUpdates || seconds()>=limit || std::find(captureSteps.begin(),captureSteps.end(),step)!=captureSteps.end();
  const auto diagnosticsStart=seconds();
  if(projectionDiagnostics){
   std::ofstream f;if(detail){f.open(out/("population_"+std::to_string(step)+".csv"));f.exceptions(std::ios::badbit|std::ios::failbit);
   f<<std::setprecision(17)<<"carrier,node_id,raw_relative,used_relative,projected,held,projection_streak,density_before_cm3,target_cm3,density_after_cm3,qf_after_V,log10_change_from_seed,unconstrained_actual_relative,nonpositive_proposal_streak,hold_activations,hold_releases\n";}
   Json extrema=Json::array();
   for(int b=0;b<2;++b){double minDensity=std::numeric_limits<double>::infinity(),minQf=minDensity,maxQf=-minDensity,minLog=0.;int outside=0;
    for(int i:ids[b+1]){
     projectionStreak[b][i]=projectedNow[b][i]?projectionStreak[b][i]+1:0;
     const double pop=units.internalConcentrationToM3(nd[b][i])*1e-6;
     const double qf=b==0?state.phin[i]:state.phip[i];
     const double logDrop=std::log10(nd[b][i])-std::log10(initialDensity[b][i]);
     minDensity=std::min(minDensity,pop);minQf=std::min(minQf,qf);maxQf=std::max(maxQf,qf);minLog=std::min(minLog,logDrop);if(qf<-.5 || qf>5.5)++outside;
     const bool held=(activityDiagnostic||densityHold) && std::binary_search(heldRows.begin(),heldRows.end(),(b+1)*N+i);
     if(detail)f<<(b==0?"electron":"hole")<<','<<i<<','<<rawRelative[b][i]<<','<<usedRelative[b][i]<<','<<projectedNow[b][i]<<','<<held<<','<<projectionStreak[b][i]<<','<<units.internalConcentrationToM3(density[b][i])*1e-6<<','<<units.internalConcentrationToM3(targets[b][i])*1e-6<<','<<pop<<','<<qf<<','<<logDrop<<','<<unconstrainedRelative[b][i]<<','<<holdState[b][i].consecutiveNonpositive<<','<<holdState[b][i].activations<<','<<holdState[b][i].releases<<'\n';
    }
    extrema.push_back({{"min_density_cm3",minDensity},{"min_qf_V",minQf},{"max_qf_V",maxQf},{"min_log10_change_from_seed",minLog},{"outside_shadow_box_minus0p5_to5p5",outside}});
   }
   h["population_extrema"]=extrema;
   if(detail)writeNewtonCarrierTermProbeCsv(out/("rows_"+std::to_string(step)+".csv"),problem.mesh,problem.doping,next,1.,1.);
   timing["diagnostic_files_extrema_s"]=seconds()-diagnosticsStart;auto portStart=seconds();
   const auto packed=solver.packArclengthState(state);
   const double scale=units.internalCurrentPerDeviceDepthToAPerUm(units.currentDensityAM2PerInternal()*units.lengthMPerInternal());
   if(shared){const auto all=solver.evaluatePackedContactCurrents(packed,scale);for(const auto& [name,bias]:problem.biases)h["contact_currents_A_per_um"][name]=all.at(name);}
   else for(const auto& [name,bias]:problem.biases){auto current=solver.makeArclengthContactCurrentFunctional(name,scale);h["contact_currents_A_per_um"][name]=current.value(packed,bias);}
   timing["contact_currents_s"]=seconds()-portStart;
  }
  if(updateDiagnostics && h.contains("contact_currents_A_per_um"))
   h["convergence_report"]["device_checks"]["contact_currents_A_per_um"]=h["contact_currents_A_per_um"];
  h["evaluation_profile"]=evaluationProfile.toJson();
  h["timing"]=timing;h["wall_s"]=seconds()-begin;
  previousProjected=projectedNow;
  result["history"].push_back(h);result["final_residual"]=rho;result["matrix_parent_iteration"]=step-1;
  result["history"].back()["excluded_qf_direction_max_V"]=excludedDirectionMax;
  const auto outputStart=seconds();put(out/"progress.json",result);
  if((opt.value("validation_snapshots",false) && step<=opt.value("validation_snapshot_max_updates",40)) || std::find(captureSteps.begin(),captureSteps.end(),step)!=captureSteps.end())
   snapshot("iteration_"+std::to_string(step),state);
  result["history"].back()["timing"]["progress_snapshot_s"]=seconds()-outputStart;
  result["history"].back()["wall_s"]=seconds()-begin;
  if(qualified){
   const auto finishStart=seconds();
   const auto bp=solver.evaluatePoissonTerms(state);double boundaryR=0;for(const auto& r:bp.rows)if(r.constrained){
    if(!std::isfinite(r.productionResidual))throw std::runtime_error("nonfinite boundary residual");
    boundaryR=std::max(boundaryR,std::abs(r.productionResidual));
   }
   if(!std::isfinite(boundaryR))throw std::runtime_error("nonfinite boundary residual");
   if(!reportOnly && !(boundaryR<1e-10))throw std::runtime_error("boundary residual guard");
   const auto saveStart=seconds();snapshot("final",state);
   result["history"].back()["timing"]["accepted_snapshot_s"]=seconds()-saveStart;
   const auto readbackStart=seconds();
   Json reloadCfg=cfg;reloadCfg["state_file"]=(out/"final.h5").string();
   const auto reloaded=readExternalState(dir,reloadCfg,N);const auto readback=solver.evaluateResidual(reloaded);
   const double rr=readback.blockNorms.combined;bool consistent=std::isfinite(rr)&&std::abs(rr-rho)<=std::max(1e-15,1e-6*std::abs(rho));
   result["readback_audit"]={{"rhs_pass",rr<=1e-10},{"residual",rr},{"consistent",consistent}};
   result["history"].back()["timing"]["readback_s"]=seconds()-readbackStart;
   if(!consistent)throw std::runtime_error("saved state residual mismatch");
   Json currents=Json::object(), components=Json::object();const auto packed=solver.packArclengthState(state);
    vela::DDScalingSpec ddScaling;
    ddScaling.regionResolvedInterfaceAssembly =
        problem.newton.regionResolvedInterfaceAssembly;
    if (problem.newton.inputScaling.isUnitScaling()) {
        const vela::UnitScalingSystem sc = vela::UnitScalingSystem::fromInputs(
            problem.newton.temperature_K,
            (11.7 * vela::constants::eps0),
            vela::UnitScalingSystem::autoInputsFrom(
                problem.mesh, problem.doping, problem.matdb, 1e10),
            vela::UnitScalingReferenceConfig{},
            problem.newton.inputScaling.unitSystem());
        ddScaling.enabled = true;
        ddScaling.V0 = sc.V0();
        ddScaling.C0 = sc.C0();
        ddScaling.mu0 = sc.mu0();
        ddScaling.D0 = sc.D0();
        ddScaling.L0 = sc.L0();
        ddScaling.permittivityReference_F_per_m =
            11.7 * vela::constants::eps0;
        ddScaling.unitSystem = problem.newton.inputScaling.unitSystem();
        ddScaling.chargeAreaFactor = ddScaling.unitSystem.chargeAreaFactor();
        ddScaling.chargeLineFactor = ddScaling.unitSystem.chargeLineFactor();
        ddScaling.fieldFromCoordinateDeltaFactor =
            ddScaling.unitSystem.fieldFromCoordinateDeltaFactor();
        ddScaling.currentDensityLineIntegralFactor =
            ddScaling.unitSystem.currentDensityAM2PerInternal() *
            ddScaling.unitSystem.lengthMPerInternal();
    }
    const vela::ContactCurrent contactCurrent(
        problem.mesh, problem.matdb, problem.doping,
        problem.newton.mobility, problem.newton.temperature_K,
        ddScaling, problem.newton.bandgapNarrowing,
        problem.newton.carrierStatistics,
        problem.newton.electronQuantumPotential);

   double scale=units.internalCurrentPerDeviceDepthToAPerUm(units.currentDensityAM2PerInternal()*units.lengthMPerInternal());
   for(const auto& [name,bias]:problem.biases){auto f=solver.makeArclengthContactCurrentFunctional(name,scale);double value=f.value(packed,bias);if(!std::isfinite(value))throw std::runtime_error("nonfinite terminal current");currents[name]=value;
    const auto detail=contactCurrent.computeDetailed(state,name);
    components[name]={{"electron_flux_A_per_um",detail.totals.electronCurrent*1e-6},{"hole_flux_A_per_um",detail.totals.holeCurrent*1e-6},{"operator_minus_flux_A_per_um",value-detail.totals.totalCurrent*1e-6}};
   }
   result["contact_currents_A_per_um"]=currents;result["contact_flux_components"]=components;
   writeNewtonCarrierTermProbeCsv(out/"final_rows.csv",problem.mesh,problem.doping,next,1.,1.);result["boundary_residual_max"]=boundaryR;
   result["strict_qualification_satisfied"]=closure.satisfied&&boundary&&boundaryR<1e-10&&initialBoundaryResidual<1e-10&&excludedDirectionMax<=1e-12;
   result["rhs_satisfied"]=rho<=1e-10;result["raw_update_satisfied"]=error<1;
   result["qualified"]=true;result["converged"]=true;
   result["stop_reason"]=reportOnly?"SLOT_rhs_or_raw_update":"SLOT_primary_and_guards";
   result["history"].back()["timing"]["accepted_finalization_s"]=seconds()-finishStart;
   result["history"].back()["wall_s"]=seconds()-begin;break;
  }
  // Stop the explicitly authorized diagnostic continuation at its original
  // R/U criterion, without granting production qualification. The controller
  // must still audit saved-state consistency and every numerical guard.
  if(diagnosticStopOnPrimary && !continueAfterPrimary && originalPolicySatisfied){
   result["stop_reason"]="diagnostic_primary_satisfied";break;
  }
  terms=std::move(next);density=std::move(nd);
  if(rho>10){result["stop_reason"]="divergence";break;}
 }
 if(!result.contains("stop_reason"))result["stop_reason"]="max_updates";
 } catch(const std::exception& e){result["stop_reason"]="guard_or_runtime_failure";result["error"]=e.what();}
 double completedStepSeconds=0.;for(const auto& h:result["history"])completedStepSeconds+=h.at("wall_s").get<double>();
 result["unassigned_or_incomplete_step_s"]=std::max(0.,seconds()-result["initialization_s"].get<double>()-completedStepSeconds);
 const auto finalOutputStart=seconds();
 if(!result["qualified"].get<bool>())snapshot("failed",state);
 result["final_snapshot_s"]=seconds()-finalOutputStart;
 result["elapsed_s"]=seconds();put(out/"result.json",result);return result;
}
