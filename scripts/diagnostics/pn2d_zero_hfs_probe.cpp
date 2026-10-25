#define main embeddedRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
int main(int argc,char**argv){try{
 using namespace vela;
 if(argc!=5)throw std::runtime_error("config output inner outer required");
 nlohmann::json cfg;{std::ifstream f(argv[1]);f>>cfg;}cfg=canonicalizeDeck(cfg);
 const double v=std::stod(argv[3]),u=std::stod(argv[4]);
 for(auto& c:cfg["contacts"])if(c["name"]=="Anode")c["bias"]=v;
 auto problem=loadNewtonProblem(argv[1],cfg);DDStateCliScope scope(problem,cfg,configDirectory(argv[1]));
 auto state=readExternalState(configDirectory(argv[1]),cfg,problem.mesh.numNodes());
 auto owner=makeNewtonSolver(problem);
 if(cfg.contains("native_population_file")){
 nlohmann::json pop;{std::ifstream f(cfg["native_population_file"].get<std::string>());f>>pop;}
 auto units=problem.newton.inputScaling.unitSystem();
 state.packedState.resize(0);state.packedLow.resize(0);state.electronQfReference.resize(0);state.holeQfReference.resize(0);state.phinIncrement.resize(0);state.phipIncrement.resize(0);
 for(const auto& a:pop){int i=a[0];state.psi[i]=a[1];state.n[i]=units.m3ToInternalConcentration(a[2].get<double>()*1e6);state.p[i]=units.m3ToInternalConcentration(a[3].get<double>()*1e6);}
 for(const auto& t:owner.evaluatePoissonTerms(state).rows){state.phin[t.nodeId]=t.electronRequiredQuasiFermiPotential;state.phip[t.nodeId]=t.holeRequiredQuasiFermiPotential;}
 }
 auto sys=owner.makeArclengthSystem("Anode");auto x=owner.packArclengthState(state);
 auto J=sys.jacobian(x,v);auto F=sys.residual(x,v);auto fv=sys.parameterDerivative(x,v);
 if(cfg.value("diagnostic_jacobian",std::string{})=="frozen_transport") {
   auto tangentCfg=cfg;tangentCfg["solver"]["impact_ionization"]["model"]="none";
   tangentCfg["solver"]["mobility"]["jacobian_field_derivatives"]=false;
   auto tangentProblem=loadNewtonProblem(argv[1],tangentCfg);auto tangentOwner=makeNewtonSolver(tangentProblem);
   auto tangentSystem=tangentOwner.makeArclengthSystem("Anode");auto tangentX=tangentOwner.packArclengthState(state);
   J=tangentSystem.jacobian(tangentX,v);
   // F, fv, current, state and acceptance still use the original full model.
 }

 if(cfg.contains("diagnostic_jv")) {
   nlohmann::json checks=nlohmann::json::array();int N=problem.mesh.numNodes();
   for(int block=0;block<3;++block)for(double step:{1e-4,1e-5,1e-6,1e-7}) {
     VectorXd direction=VectorXd::Zero(x.size());for(int i=0;i<N;++i)direction[block*N+i]=std::sin(0.37*(i+1));
     VectorXd fd=(sys.residual(x+step*direction,v)-sys.residual(x-step*direction,v))/(2*step);VectorXd product=J*direction;
     nlohmann::json entry={{"input_block",block},{"step_scaled",step}};
     for(int output=0;output<3;++output){auto a=product.segment(output*N,N);auto b=fd.segment(output*N,N);entry[std::to_string(output)]=(a-b).norm()/std::max(a.norm(),b.norm());}
     checks.push_back(entry);
   }
   std::filesystem::create_directories(argv[2]);std::ofstream(std::filesystem::path(argv[2])/"jv.json")<<checks.dump(2);
 }
 auto unit=problem.newton.inputScaling.unitSystem();const double cs=unit.currentDensityAM2PerInternal()*unit.lengthMPerInternal()*1e-6;
 auto port=owner.makeArclengthContactCurrentFunctional("Anode",cs).linearize(x,v);
 const int n=problem.mesh.numNodes(),d=3*n;std::filesystem::path out=argv[2];std::filesystem::create_directories(out);
 std::ofstream m(out/"augmented.csv"),r(out/"rhs.csv");m<<std::setprecision(17)<<"row,col,value\n";r<<std::setprecision(17)<<"row,rhs,x\n";
 for(int k=0;k<J.outerSize();++k)for(SparseMatrixd::InnerIterator it(J,k);it;++it)m<<it.row()<<','<<it.col()<<','<<it.value()<<'\n';
 for(int i=0;i<d;++i){if(fv[i]!=0)m<<i<<','<<d<<','<<fv[i]<<'\n';if(port.stateDerivative[i]!=0)m<<d<<','<<i<<','<<1e7*port.stateDerivative[i]<<'\n';r<<i<<','<<-F[i]<<','<<x[i]<<'\n';}
 m<<d<<','<<d<<','<<1+1e7*port.parameterDerivative<<'\n';r<<d<<','<<-(v+1e7*port.value-u)<<','<<v<<'\n';
 nlohmann::json contact=nlohmann::json::array();for(const auto& c:problem.mesh.contacts())for(int id:c.node_ids)contact.push_back(id);
 auto gate=owner.makeArclengthConvergenceEvaluator("Anode")(x,v);
 nlohmann::json j={{"nodes",n},{"device_dimension",d},{"augmented_dimension",d+1},{"contact_nodes",contact},{"inner_V",v},{"outer_V",u},{"current_A_per_um",port.value},{"potential_scale_V",owner.evaluateResidual(state).potentialScale},{"residual_inf",F.lpNorm<Eigen::Infinity>()}};
 auto actual=owner.makeArclengthStateUnpacker()(x);
 std::ofstream states(out/"physical_state.csv");states<<std::setprecision(17)<<"node,psi_V,fn_V,fp_V,n_cm3,p_cm3,ni_cm3,ni_supplied_cm3,psi_boundary_target_V,fn_boundary_target_V,fp_boundary_target_V\n";
 for(const auto& t:owner.evaluatePoissonTerms(state).rows){int i=t.nodeId;states<<i<<','<<actual.psi[i]<<','<<actual.phin[i]<<','<<actual.phip[i]<<','<<unit.internalConcentrationToM3(actual.n[i])/1e6<<','<<unit.internalConcentrationToM3(actual.p[i])/1e6<<','<<unit.internalConcentrationToM3(t.intrinsicDensity)/1e6<<','<<std::sqrt(unit.internalConcentrationToM3(state.n[i])*unit.internalConcentrationToM3(state.p[i]))/1e6<<','<<actual.psi[i]-F[i]*j["potential_scale_V"].get<double>()<<','<<actual.phin[i]-F[n+i]*j["potential_scale_V"].get<double>()<<','<<actual.phip[i]-F[2*n+i]*j["potential_scale_V"].get<double>()<<'\n';}
 std::ofstream edges(out/"edges.csv");edges<<std::setprecision(17)<<"edge,node0,node1,mu_n,mu_p,doping_internal,field_n_internal,fn0,fn1,fp0,fp1\n";
 for(const auto& e:owner.evaluateSgEdgeFluxDiagnostics(state))edges<<e.edgeId<<','<<e.node0<<','<<e.node1<<','<<e.electronMobility_m2_V_s<<','<<e.holeMobility_m2_V_s<<','<<e.netDopingAvg_m3<<','<<e.electronMobilityField_V_m<<','<<e.phin0_V<<','<<e.phin1_V<<','<<e.phip0_V<<','<<e.phip1_V<<'\n';
 std::ofstream(out/"identity.json")<<j.dump(2);std::cout<<j.dump()<<'\n';
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
