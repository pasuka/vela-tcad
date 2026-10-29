// Uninstalled finite-coordinate audit. Never changes production stopping or
// updates; the optional 1% rule is an observed native control, not a new default.
#define main embeddedExampleRunnerMain
#include "../../src/tools/vela_example_runner.cpp"
#undef main
#include "vela/simulation/CoupledLoadLine.h"
#include "vela/physics/experimental/JoyceDixonStatistics.h"

int main(int argc,char**argv){try{
 using namespace vela;using J=nlohmann::json;namespace fs=std::filesystem;
 if(argc!=3)throw std::invalid_argument("manifest output-directory");
 J m,cfg;{std::ifstream f(argv[1]);f>>m;}{std::ifstream f(m.at("config").get<std::string>());f>>cfg;}
 cfg=canonicalizeDeck(cfg);const fs::path out=argv[2];fs::create_directories(out);
 const double inner=m.at("inner_V"),outer=m.at("outer_V");
 for(auto&c:cfg["contacts"])if(c["name"]=="drain")c["bias"]=inner;
 const auto configPath=m.at("config").get<std::string>();
 auto p=loadNewtonProblem(configPath,cfg);DDStateCliScope scope(p,cfg,configDirectory(configPath));
 auto solver=makeNewtonSolver(p);auto unpack=solver.makeArclengthStateUnpacker();auto sys=solver.makeArclengthSystem("drain");
 DDSolution seed;
 if(m.contains("fields")){
  const auto&f=m.at("fields");const int n=p.mesh.numNodes();seed.psi.resize(n);seed.phin.resize(n);seed.phip.resize(n);seed.n=seed.p=VectorXd::Zero(n);
  for(int i=0;i<n;++i)seed.psi[i]=f.at("psi_V")[i];
  for(const auto&v:f.at("electron"))seed.n[v[0].get<int>()]=p.newton.inputScaling.unitSystem().m3ToInternalConcentration(v[1].get<double>()*1e6);
  for(const auto&v:f.at("hole"))seed.p[v[0].get<int>()]=p.newton.inputScaling.unitSystem().m3ToInternalConcentration(v[1].get<double>()*1e6);
  for(const auto&t:solver.evaluatePoissonTerms(seed).rows){seed.phin[t.nodeId]=t.electronRequiredQuasiFermiPotential;seed.phip[t.nodeId]=t.holeRequiredQuasiFermiPotential;}
 }else seed=readDDSolutionState(m.at("state_file").get<std::string>(),p.mesh.numNodes(),p.newton.inputScaling);
 const VectorXd x=solver.packArclengthState(seed);seed=unpack(x);const int n=p.mesh.numNodes();
 const auto&u=p.newton.inputScaling.unitSystem();const double line=u.currentDensityAM2PerInternal()*u.lengthMPerInternal()*1e-6;
 auto port=solver.makeArclengthContactCurrentFunctional("drain",line);
 CoupledLoadLineConfig c;c.maxIterations=1;c.equationTolerance=1e-30;c.loadLineTolerance_V=1e-30;
 c.linearSolver="direct_bordered";c.lineSearchMode="full_step";c.linearizationAuditEnabled=true;c.applyDeviceUpdateLimit=false;c.maxVoltageUpdate_V=1e100;
 const auto result=CoupledLoadLineNewton(sys,port.value,c,port).solve({x,inner},outer,1e7);
 if(result.trace.empty()||result.trace[0].rawStateUpdate.size()!=3*n)throw std::runtime_error("raw direction missing");
 const auto&t=result.trace[0];const auto delta=t.rawStateUpdate;const double dv=t.rawVoltageUpdate_V;
 const double v0=solver.evaluateResidual(seed).potentialScale,vt=constants::kb*300/constants::q;
 auto cm=[&](double d){return u.internalConcentrationToM3(d)/1e6;};auto internal=[&](double d){return u.m3ToInternalConcentration(d*1e6);};
 VectorXd dn=VectorXd::Zero(n),dp=dn,nc=dn,nv=dn;std::vector<bool> active(n,false),nf(n,false),pf(n,false);
 for(const auto&v:solver.evaluateCarrierTermDiagnostics(seed,true).rows){nf[v.nodeId]=v.electronContinuityActive;pf[v.nodeId]=v.holeContinuityActive;}
 for(const auto&v:solver.evaluatePoissonTerms(seed).rows){int i=v.nodeId;nc[i]=cm(v.electronDensityOfStates);nv[i]=cm(v.holeDensityOfStates);active[i]=v.intrinsicDensity>0.;
  if(active[i]){dn[i]=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(cm(seed.n[i]),nc[i])*(delta[i]-delta[n+i])*v0/vt;
   dp[i]=experimental::JoyceDixonStatistics::densitySlopeFromPopulation(cm(seed.p[i]),nv[i])*(delta[2*n+i]-delta[i])*v0/vt;}}
 std::ofstream f(out/"direction.csv");f<<std::setprecision(17)<<"node,psi_V,fn_V,fp_V,n_cm3,p_cm3,Nc_cm3,Nv_cm3,e_free,h_free,dpsi_V,dfn_V,dfp_V,dn_cm3,dp_cm3\n";
 for(int i=0;i<n;++i)f<<i<<','<<seed.psi[i]<<','<<seed.phin[i]<<','<<seed.phip[i]<<','<<cm(seed.n[i])<<','<<cm(seed.p[i])<<','<<nc[i]<<','<<nv[i]<<','<<nf[i]<<','<<pf[i]<<','<<delta[i]*v0<<','<<delta[n+i]*v0<<','<<delta[2*n+i]*v0<<','<<dn[i]<<','<<dp[i]<<'\n';
 J report={{"v0_V",v0},{"thermal_voltage_V",vt},{"inner_V",inner},{"outer_V",outer},{"raw_inner_V",dv},{"linear_backward_error",t.componentwiseBackwardError},{"candidates",J::array()}};
 auto inspect=[&](const std::string&name,const VectorXd& z,double bias){J j={{"name",name},{"inner_V",bias}};
  try{const auto s=unpack(z);const auto r=sys.residual(z,bias);j["residual_finite"]=r.allFinite();j["equation_inf"]=r.lpNorm<Eigen::Infinity>();j["load_V"]=bias+1e7*port.value(z,bias)-outer;
   j["focus"]=J::array();for(int i:{145,152,308,321,415,418,1120})j["focus"].push_back({{"node",i},{"n_cm3",cm(s.n[i])},{"p_cm3",cm(s.p[i])},{"fn_V",s.phin[i]},{"fp_V",s.phip[i]}});
   writeDDSolutionState(out/(name+".h5"),s,p.newton.inputScaling);
  }catch(const std::exception&e){j["error"]=e.what();}
  report["candidates"].push_back(j);};
 inspect("seed",x,inner);inspect("qf_full",x+delta,inner+dv);
 DDSolution density=seed;density.packedState.resize(0);density.packedLow.resize(0);density.electronQfReference.resize(0);density.holeQfReference.resize(0);density.phinIncrement.resize(0);density.phipIncrement.resize(0);
 int negativeN=0,negativeP=0;for(int i=0;i<n;++i){density.psi[i]+=delta[i]*v0;
  if(active[i]){double a=cm(seed.n[i])+dn[i],b=cm(seed.p[i])+dp[i];negativeN+=a<=0.;negativeP+=b<=0.;
   density.n[i]=internal(a>0?a:.01*cm(seed.n[i]));density.p[i]=internal(b>0?b:.01*cm(seed.p[i]));}}
 report["linear_nonpositive_electron_count"]=negativeN;report["linear_nonpositive_hole_count"]=negativeP;
 for(const auto&v:solver.evaluatePoissonTerms(density).rows){density.phin[v.nodeId]=v.electronRequiredQuasiFermiPotential;density.phip[v.nodeId]=v.holeRequiredQuasiFermiPotential;}
 inspect("population_observed_1pct_control",solver.packArclengthState(density),inner+dv);
 std::ofstream(out/"summary.json")<<report.dump(2)<<'\n';std::cout<<report.dump(2)<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
