"""Prepare an isolated Canali PMI diagnostic; no deployment or simulations.

This is not a production HFS implementation or a claim of native equivalence.
Native baseline/observer identity must qualify before interpreting samples.
"""
import argparse,json,subprocess,os
from pathlib import Path
from decimal import Decimal,localcontext
import csv,io
import calibrate_simplemos_enormal_20260912 as e
a,d=e.a,e.d
L=e.R/'build-release/hfs_observer_preflight_20260912'
O=e.R/'reference_tcad/simplemos_sentaurus2022/hfs_observer_preflight_20260912'

HEADER=r'''#pragma once
#include <cmath>
#include <stdexcept>
struct Canali { double mu,dm,dF,dT; };
inline Canali canali(double m,double F,double T,bool electron) {
    if(m<0 || F<0 || T<=0 || !std::isfinite(m+F+T))
        throw std::invalid_argument("Canali diagnostic input domain");
    if(m==0)return {0,1,0,0};
    if(F==0)return {m,1,0,0};
    const double b0=electron?1.109:1.213,bt=electron?.66:.17;
    const double v0=electron?1.07e7:8.37e6,vt=electron?.87:.52;
    const double b=b0*std::pow(T/300.,bt),v=v0*std::pow(T/300.,-vt);
    const double x=m*F/v,y=std::pow(x,b),logden=std::log1p(y);
    const double mu=m*std::exp(-logden/b),dm=std::exp(-(1+1/b)*logden);
    const double df=-m*m/v*std::pow(x,b-1)*dm;
    const double dt=mu/T*(bt/b*logden-y/(1+y)*(bt*std::log(x)+vt));
    return {mu,dm,df,dt};
}
'''

PMI=r'''// Diagnostic replacement, NOT zero-additive. Must pass native identity.
#include "PMIModels.h"
#include "canali_diagnostic.h"
#include <array>
#include <map>
#include <tuple>
#include <mutex>
#include <fstream>
#include <iomanip>
#include <cctype>
#include <string>
class HfsObserver : public PMI_HighFieldMobility {
    bool electron_;
    std::string filename_;
    std::map<std::tuple<double,double,double>,std::array<double,14> > samples_;
    std::mutex lock_;
    unsigned long long sequence_=0;
public:
    HfsObserver(const PMI_Environment& e,PMI_HighFieldDrivingForce f,
        PMI_AnisotropyType a,bool electron):PMI_HighFieldMobility(e,f,a),electron_(electron) {
        std::string region=ReadRegionName();
        for(std::size_t i=0;i<region.size();++i)
            if(!std::isalnum(static_cast<unsigned char>(region[i])) && region[i]!='_')region[i]='_';
        filename_=std::string("hfs_observer_")+(electron?"e_":"h_")+region+".csv";
    }
    ~HfsObserver() {
        std::ofstream out(filename_.c_str());
        out<<"x_um,y_um,z_um,potential_V,n_cm3,p_cm3,temperature_K,carrier_temperature_K,mulow_cm2_V_s,F_V_cm,mu_cm2_V_s,driving_force_enum,last_call,interface_distance_um\n"<<std::setprecision(17);
        for(const auto& row:samples_) {
            for(std::size_t k=0;k<row.second.size();++k)out<<(k?",":"")<<row.second[k];
            out<<'\n';
        }
    }
    void Compute_mu(double pot,double n,double p,double t,double ct,double m,double F,double& mu) override {
        mu=canali(m,F,t,electron_).mu;
        double x,y,z;ReadCoordinate(x,y,z);
        const double distance=ReadDistanceFromSemiconductorInsulatorInterface();
        std::lock_guard<std::mutex> guard(lock_);
        samples_[std::make_tuple(x,y,z)]={{x,y,z,pot,n,p,t,ct,m,F,mu,
            static_cast<double>(HighFieldDrivingForce()),static_cast<double>(++sequence_),distance}};
    }
#define ZERO(name) void name(double,double,double,double,double,double,double,double& r) override {r=0;}
    ZERO(Compute_dmudpot) ZERO(Compute_dmudn) ZERO(Compute_dmudp)
    ZERO(Compute_dmudct) ZERO(Compute_dmudNa) ZERO(Compute_dmudNd)
#undef ZERO
#define PART(name,field) void name(double,double,double,double t,double,double m,double F,double& r) override {r=canali(m,F,t,electron_).field;}
    PART(Compute_dmudt,dT) PART(Compute_dmudmulow,dm) PART(Compute_dmudF,dF)
#undef PART
};
extern "C" PMI_HighFieldMobility* new_PMI_HighField_e_Mobility(const PMI_Environment& e,const PMI_HighFieldDrivingForce f,const PMI_AnisotropyType a) {return new HfsObserver(e,f,a,true);}
extern "C" PMI_HighFieldMobility* new_PMI_HighField_h_Mobility(const PMI_Environment& e,const PMI_HighFieldDrivingForce f,const PMI_AnisotropyType a) {return new HfsObserver(e,f,a,false);}
'''

TEST=r'''#include "canali_diagnostic.h"
#include <iostream>
#include <iomanip>
int main() {
    std::cout<<"e,m,F,T,mu,dm,dF,dT\n"<<std::setprecision(17);
    for(bool e:{false,true})for(double m:{0.,10.,500.,1500.})
      for(double F:{0.,1e-16,1.,1e3,1e5,1e8})for(double T:{250.,300.,350.}) {
        const auto r=canali(m,F,T,e);
        std::cout<<e<<','<<m<<','<<F<<','<<T<<','<<r.mu<<','<<r.dm<<','<<r.dF<<','<<r.dT<<'\n';
      }
}
'''

def main():
    assert not (O/'evidence.json').exists()
    L.mkdir(parents=True,exist_ok=True)
    for name,text in [('canali_diagnostic.h',HEADER),('pmi_vela_hfs_observer.C',PMI),('canali_diagnostic_test.cpp',TEST)]:
        (L/name).write_text(text,encoding='utf-8',newline='\n')
    env=os.environ.copy();env['PATH']='D:/msys64/ucrt64/bin'+os.pathsep+env['PATH']
    r=subprocess.run(['D:/msys64/ucrt64/bin/g++.exe','-std=c++20','-O2',str(L/'canali_diagnostic_test.cpp'),'-o',str(L/'canali_diagnostic_test.exe')],env=env,capture_output=True,text=True)
    (L/'compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    r=subprocess.run([str(L/'canali_diagnostic_test.exe')],env=env,capture_output=True,text=True)
    (L/'values.csv').write_text(r.stdout);assert r.returncode==0
    rows=[]
    with localcontext() as ctx:
        ctx.prec=90
        D=lambda x:Decimal.from_float(float(x))
        def mu(m,F,T,e):
            if not m or not F:return m
            b=D(1.109 if e else 1.213)*(T/300)**D(.66 if e else .17)
            v=D(1.07e7 if e else 8.37e6)*(T/300)**(-D(.87 if e else .52))
            return m/(1+(m*F/v)**b)**(1/b)
        for item in csv.DictReader(io.StringIO(r.stdout)):
            m,F,T=map(D,(item['m'],item['F'],item['T']));electron=item['e']=='1'
            ref=mu(m,F,T,electron)
            refs={'mu':ref};steps=[]
            if m and F:
                for variable,name in enumerate(('dm','dF','dT')):
                    values=[m,F,T];slopes=[]
                    for relative in ('1e-12','5e-13'):
                        h=values[variable]*Decimal(relative);plus=values.copy();minus=values.copy()
                        plus[variable]+=h;minus[variable]-=h
                        slopes.append((mu(*plus,electron)-mu(*minus,electron))/(2*h))
                    refs[name]=slopes[-1];steps.append(float(abs(slopes[0]/slopes[1]-1)))
            else:refs.update(dm=Decimal(1),dF=Decimal(0),dT=Decimal(0))
            errors={key:float(abs(D(item[key])/val-1)) if val else abs(float(item[key])) for key,val in refs.items()}
            rows.append(dict(electron=electron,m=float(m),F=float(F),T=float(T),**errors,step_error=max(steps,default=0.),qualified=max(errors.values())<=1e-10 and max(steps,default=0.)<=1e-10))
    a.write_csv(O/'checks.csv',rows)
    summary=dict(checks=len(rows),qualified=sum(r['qualified'] for r in rows),max_error=max(max(r[k] for k in ('mu','dm','dF','dT')) for r in rows),scope='Standalone alpha=0 Canali formula and partials, original default Si parameters. No native ABI/build/identity qualification, no production changes, no HFS simulation.',sample_scope='Observer keeps last cell-context call per vertex. Final-state identity and cell matching are required; missing or ambiguous samples are not filled.',all_qualified=all(r['qualified'] for r in rows))
    a.write(O/'summary.json',summary)
    d.matrix.freeze(O/'evidence.json',[Path(__file__).resolve(),e.L/'vendor_reference/enormal_native_models.par',O/'checks.csv',O/'summary.json']+[p for p in L.iterdir() if p.is_file()])
    print(summary,flush=True);assert summary['all_qualified']

if __name__=='__main__':main()
