"""Replay native PMI inputs through its exact C++ header and independent Decimal FD.

These checks qualify local constitutive partials, not the native assembled
Jacobian or state-to-driving-field chain. All zero-field samples are checked
separately by their 300 K one-sided limit.
"""
import argparse
import csv
import io
import os
import subprocess
from decimal import Decimal as D, localcontext
from pathlib import Path
import check_simplemos_hfs_native_identity_20260912 as n

a,d=n.a,n.d
L=n.L/'derivatives_20260914';O=n.O/'derivatives_20260914'

CPP=r'''#include "canali_diagnostic.h"
#include <fstream>
#include <iostream>
#include <iomanip>
int main(int argc,char**argv) {
 if(argc!=2)return 2;
 std::ifstream in(argv[1]); int id,e;double m,F,T;
 std::cout<<"id,mu,dm,dF,dT\n"<<std::setprecision(17);
 while(in>>id>>e>>m>>F>>T){auto r=canali(m,F,T,e!=0);
  std::cout<<id<<','<<r.mu<<','<<r.dm<<','<<r.dF<<','<<r.dT<<'\n';}
}
'''


def mu(m,F,T,e):
    if m==0 or F==0:return m
    b=D.from_float(1.109 if e else 1.213)*(T/300)**D.from_float(.66 if e else .17)
    v=D.from_float(1.07e7 if e else 8.37e6)*(T/300)**(-D.from_float(.87 if e else .52))
    return m/(1+(m*F/v)**b)**(1/b)


def main(stage):
    assert not (O/f'{stage}_evidence.json').exists()
    a.verify(n.O/f'{stage}_observer_evidence.json')
    assert a.read(n.O/f'{stage}_observer_summary.json')['all_observed_final_states_qualified']
    a.write(O/f'{stage}_contract.json',dict(decimal_precision=90,
        relative_steps=['1e-12','5e-13'],relative_gate=1e-10,
        scope='Every actual last-call PMI sample. Exact compiled Canali header; local m,F,T partials only. No device Jacobian certification.'))
    L.mkdir(parents=True,exist_ok=True)
    header=n.L/'pmi/canali_diagnostic.h'
    src=L/f'{stage}_probe.cpp';src.write_text(CPP,newline='\n')
    exe=L/f'{stage}_probe.exe'
    env=os.environ.copy();env['PATH']='D:/msys64/ucrt64/bin'+os.pathsep+env['PATH']
    r=subprocess.run(['D:/msys64/ucrt64/bin/g++.exe','-std=c++20','-O2','-I'+str(header.parent),str(src),'-o',str(exe)],env=env,capture_output=True,text=True)
    (L/f'{stage}_compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    samples=a.rows(n.O/f'{stage}_observer_samples.csv')
    inputs=L/f'{stage}_input.txt'
    inputs.write_text(''.join(f"{i} {int(s['carrier']=='e')} {s['mulow_cm2_V_s']} {s['F_V_cm']} 300\n" for i,s in enumerate(samples)),newline='\n')
    r=subprocess.run([str(exe),str(inputs)],env=env,capture_output=True,text=True)
    (L/f'{stage}_cpp.csv').write_text(r.stdout);assert r.returncode==0,r.stderr
    values=list(csv.DictReader(io.StringIO(r.stdout)));assert len(values)==len(samples)
    results=[]
    with localcontext() as ctx:
        ctx.prec=90
        for i,(s,v) in enumerate(zip(samples,values)):
            assert int(v['id'])==i
            m,F=D.from_float(float(s['mulow_cm2_V_s'])),D.from_float(float(s['F_V_cm']))
            e=s['carrier']=='e';xx=[m,F,D(300)]
            refs={'mu':mu(*xx,e)};step_error=D(0)
            if m and F:
                for k,name in enumerate(('dm','dF','dT')):
                    slopes=[]
                    for hrel in ('1e-12','5e-13'):
                        h=xx[k]*D(hrel);xp=xx.copy();xm=xx.copy();xp[k]+=h;xm[k]-=h
                        slopes.append((mu(*xp,e)-mu(*xm,e))/(2*h))
                    refs[name]=slopes[-1]
                    step_error=max(step_error,abs(slopes[0]/slopes[1]-1) if slopes[1] else abs(slopes[0]))
            else:refs.update(dm=D(1),dF=D(0),dT=D(0))
            errors={k:float(abs(D.from_float(float(v[k]))/r-1)) if r else abs(float(v[k])) for k,r in refs.items()}
            observer_error=abs(float(v['mu'])/float(s['mu_cm2_V_s'])-1)
            results.append(dict(key=s['key'],carrier=s['carrier'],node=s['node'],F=float(F),m=float(m),
                zero_field=not F,**{k+'_relative':v for k,v in errors.items()},
                step_relative=float(step_error),observer_mu_relative=observer_error,
                qualified=max(*errors.values(),float(step_error),observer_error)<=1e-10))
            if i%1000==0:print(stage,i,'/',len(samples),flush=True)
    a.write_csv(O/f'{stage}_checks.csv',results)
    summary=dict(samples=len(results),qualified=sum(r['qualified'] for r in results),
        zero_field_samples=sum(r['zero_field'] for r in results),
        max_error=max(max(r[k] for k in ('mu_relative','dm_relative','dF_relative','dT_relative','step_relative','observer_mu_relative')) for r in results),
        all_local_partials_qualified=all(r['qualified'] for r in results),device_chain_qualified=False,production_changed=False)
    a.write(O/f'{stage}_summary.json',summary)
    d.matrix.freeze(O/f'{stage}_evidence.json',[Path(__file__).resolve(),header,n.O/f'{stage}_observer_evidence.json',
        O/f'{stage}_contract.json',O/f'{stage}_checks.csv',O/f'{stage}_summary.json']+list(L.glob(f'{stage}_*')))
    print(summary,flush=True);assert summary['all_local_partials_qualified']


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=('pilot','rest'),required=True)
    main(p.parse_args().stage)
