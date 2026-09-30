"""Compare the diagnostic C++ Poisson kernel with independent Decimal arithmetic."""
from decimal import Decimal as D,localcontext
from pathlib import Path
import subprocess
import validate_simplemos_poisson_precision_20260910 as c
import analyze_simplemos_e1089_precision_20260909 as a
s=c.s;LOCAL=c.LOCAL/'kernel_check';OUT=c.OUT/'kernel_check'
CPP=r'''
#include "simplemos_poisson_precision.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
    if(argc!=3)return 2;
    nlohmann::json jobs;std::ifstream(argv[1])>>jobs;
    std::ofstream out(argv[2]);out<<std::setprecision(40)<<"case,mode,node,value\n";
    for(const auto& job:jobs) {
        std::vector<simplemos_poisson_precision::Node> nodes;
        std::vector<simplemos_poisson_precision::Edge> edges;
        for(const auto& n:job["nodes"])nodes.push_back({n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7],n[8],n[9],n[10],n[11],n[12],n[13],n[14],n[15],n[16],n[17],n[18]});
        for(const auto& e:job["edges"])edges.push_back({e[0],e[1],e[2]});
        for(bool packed:{false,true}) {
            const auto r=simplemos_poisson_precision::evaluate(nodes,edges,packed);
            for(std::size_t i=0;i<r.size();++i)out<<job["key"].get<std::string>()<<','<<(packed?"packed":"kernel")<<','<<i<<','<<r[i]<<'\n';
        }
    }
}
'''

def run():
    s.a.verify(c.v.OUT/'replay_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    (LOCAL/'check.cpp').write_text(CPP)
    cmd=['D:/msys64/ucrt64/bin/c++.exe','-std=c++20','-O2','-I'+str(c.HEADER.parent),str(LOCAL/'check.cpp'),'-o',str(LOCAL/'check.exe')]
    s.a.write(LOCAL/'compile.json',cmd);r=subprocess.run(cmd,env=s.V.environment(),capture_output=True,text=True);(LOCAL/'compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    jobs=[];expected={};scales={}
    with localcontext() as ctx:
        ctx.prec=100
        for j in s.a.read(c.v.OUT/'contract.json')['jobs']:
            root=Path(j['dest'])
            for label in ('base','trial_0','trial_5','trial_12'):
                key=j['key']+'_'+j['label']+'_'+label;pn=a.records(root/(label+'_poisson_nodes.csv'));nd=a.records(root/(label+'_nodes.csv'));pe=a.records(root/(label+'_poisson_edges.csv'));nodes=[]
                for n,d in zip(pn,nd):
                    nodes.append([float(x) for x in (n['psi'],n['n'],n['p'],n['xpsi'],n['xn'],n['xp'],n['potential_scale'],d['eref'],d['href'],d['ni'],d['Vt'],n['net_doping'],n['voln'],n['volp'],n['vold'],n['q'],n['area_factor'],n['poisson_scale'],n['interface_rhs'])])
                jobs.append(dict(key=key,nodes=nodes,edges=[[int(e['node0']),int(e['node1']),float(e['G'])] for e in pe]))
                for mode in ('kernel','packed'):
                    psi=[n['psi'] if mode=='kernel' else n['xpsi']*n['potential_scale'] for n in pn];rr=[D(0) for n in pn];ss=[D(0) for n in pn]
                    for e in pe:
                        i,jj=int(e['node0']),int(e['node1']);f=e['G']*(psi[i]-psi[jj]);rr[i]+=f;rr[jj]-=f;ss[i]+=abs(f);ss[jj]+=abs(f)
                    for i,(n,d) in enumerate(zip(pn,nd)):
                        ne,ho=n['n'],n['p']
                        if mode=='packed':
                            ne=d['ni']*a.hp.clamp((psi[i]-d['eref']-n['xn']*n['potential_scale'])/d['Vt']).exp() if d['ni']>0 else D(0)
                            ho=d['ni']*a.hp.clamp((d['href']+n['xp']*n['potential_scale']-psi[i])/d['Vt']).exp() if d['ni']>0 else D(0)
                        qn=n['q']*ne*n['voln']*n['area_factor'];qp=n['q']*ho*n['volp']*n['area_factor'];qd=n['q']*n['net_doping']*n['vold']*n['area_factor']
                        expected[key,mode,i]=(rr[i]+qn-qp-qd-n['interface_rhs'])/n['poisson_scale']
                        scales[key,mode,i]=max((ss[i]+abs(qn)+abs(qp)+abs(qd)+abs(n['interface_rhs']))/n['poisson_scale'],D('1e-300'))
        s.a.write(LOCAL/'inputs.json',jobs)
        r=subprocess.run([str(LOCAL/'check.exe'),str(LOCAL/'inputs.json'),str(LOCAL/'results.csv')],env=s.V.environment(),capture_output=True,text=True);assert r.returncode==0,r.stderr
        results=s.a.rows(LOCAL/'results.csv');assert len(results)==len(expected)
        errors={}
        for row in results:
            k=(row['case'],row['mode'],int(row['node']));err=abs(D(row['value'])-expected[k])/scales[k];assert err<D('1e-25'),(k,err)
            errors[row['mode']]=max(errors.get(row['mode'],D(0)),err)
    summary=dict(checks=len(results),passed=True,max_scaled_error={k:float(e) for k,e in errors.items()},gate=1e-25,reference='Independent 100-digit Decimal reconstruction at 16 frozen actual states, all Poisson nodes including unconstrained boundary rows, two modes.')
    s.a.write(OUT/'summary.json',summary);s.d.matrix.freeze(OUT/'evidence.json',[Path(__file__).resolve(),c.HEADER,c.v.OUT/'replay_evidence.json',OUT/'summary.json']+[p for p in LOCAL.iterdir() if p.is_file()]);print(summary,flush=True)

if __name__=='__main__':run()
