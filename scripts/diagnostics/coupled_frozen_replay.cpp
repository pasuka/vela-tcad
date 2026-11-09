#include "vela/simulation/CoupledLoadLine.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <cstdlib>
using namespace vela;static double num(const std::string&s){return std::strtod(s.c_str(),nullptr);}
int main(int argc,char**argv){try{
 if(argc!=3)throw std::invalid_argument("Usage: vela_coupled_frozen_replay <snapshot-directory> <report.json>");
 std::filesystem::path d=argv[1];std::ifstream v(d/"last_matrix_state.csv");std::string line;if(!v)throw std::runtime_error("Cannot read last_matrix_state.csv");std::getline(v,line);std::vector<double> bv,sv;
 while(std::getline(v,line)){std::stringstream in(line);std::vector<std::string> c;std::string t;while(std::getline(in,t,','))c.push_back(t);if(c.size()<4||std::stoi(c[0])!=static_cast<int>(bv.size()))throw std::runtime_error("Invalid ordered snapshot state");bv.push_back(num(c[2]));sv.push_back(num(c[3]));}
 if(bv.size()<2)throw std::runtime_error("Empty snapshot state");
 const int n=bv.size()-1;VectorXd rhs=Eigen::Map<VectorXd>(bv.data(),n+1),scales=Eigen::Map<VectorXd>(sv.data(),n+1);
 std::ifstream m(d/"last_matrix_matrix.csv");if(!m)throw std::runtime_error("Cannot read last_matrix_matrix.csv");std::getline(m,line);std::vector<Eigen::Triplet<double>> entries,dev;VectorXd fv=VectorXd::Zero(n),port=VectorXd::Zero(n);double pv=0.;
 while(std::getline(m,line)){std::stringstream in(line);std::string a,b,c;std::getline(in,a,',');std::getline(in,b,',');std::getline(in,c);int i=std::stoi(a),j=std::stoi(b);double value=num(c);if(i<0||j<0||i>n||j>n||!std::isfinite(value))throw std::runtime_error("Invalid finite snapshot matrix");entries.emplace_back(i,j,value);if(i<n&&j<n)dev.emplace_back(i,j,value);else if(i<n)fv[i]=value;else if(j<n)port[j]=value;else pv=value;}
 if(!rhs.allFinite()||!scales.allFinite()||(scales.array()<=0.).any())throw std::runtime_error("Invalid snapshot RHS or column scales");
 SparseMatrixd a(n+1,n+1),j(n,n);a.setFromTriplets(entries.begin(),entries.end());j.setFromTriplets(dev.begin(),dev.end());
 ArclengthSystem s;s.residual=[&](const VectorXd&x,double lambda){return VectorXd(j*x+fv*lambda-rhs.head(n));};s.jacobian=[&](const VectorXd&,double){return j;};s.parameterDerivative=[&](const VectorXd&,double){return fv;};s.solveJacobian=[](const VectorXd&,double,const VectorXd&,VectorXd&){return false;};s.columnScales=[&](const ArclengthState&){return scales;};
 ArclengthScalarFunctional p;p.value=[&](const VectorXd& x,double v){return port.dot(x)+(pv-1.)*v-rhs[n];};p.linearize=[&](const VectorXd&x,double v){ArclengthScalarLinearization l;l.value=p.value(x,v);l.stateDerivative=port;l.parameterDerivative=pv-1.;return l;};
 CoupledLoadLineConfig c;c.maxIterations=1;c.equationTolerance=1e-30;c.loadLineTolerance_V=1e-30;c.linearSolver="direct_bordered";c.lineSearchMode="full_step";c.maxVoltageUpdate_V=1e100;c.applyDeviceUpdateLimit=false;c.directRefinementSteps=1;
 const auto r=CoupledLoadLineNewton(s,p.value,c,p).solve({VectorXd::Zero(n),0.},0.,1.);VectorXd dx(n+1);dx.head(n)=r.state.x;dx[n]=r.state.lambda;
 auto defect=linear_refinement::highPrecisionDefect(a,-rhs,dx);
 nlohmann::json result={{"finite",dx.allFinite()},{"backward_error",defect.backwardMax},{"failure",r.failureReason},{"accepted_alpha",r.trace.empty()?0.:r.trace.front().acceptedAlpha}};
 std::ofstream(argv[2])<<result.dump(2);std::cout<<result.dump()<<'\n';return dx.allFinite()&&defect.backwardMax<1e-10&&!r.trace.empty()&&r.trace.front().acceptedAlpha>0.?0:2;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
