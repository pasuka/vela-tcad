"""Three-block Poisson directional derivative audit at the two frozen bases."""
from pathlib import Path
import copy
import subprocess
import numpy as np
from decimal import Decimal as D,localcontext
import calibrate_simplemos_packed_poisson_response_20260910 as r
import analyze_simplemos_e1089_precision_20260909 as a
s=r.s;LOCAL=r.LOCAL/'jvp';OUT=r.OUT/'jvp'

def run():
    s.a.verify(r.OUT/'preflight_evidence.json');s.a.verify(r.c.OUT/'kernel_check/evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    fixtures=s.a.read(r.c.LOCAL/'kernel_check/inputs.json');jobs=[];tests=[];refs={};matrices={};masks={}
    for case in s.a.read(r.OUT/'contract.json')['cases']:
        label='minus_full' if case['device']=='n19' else 'plus_full';key=case['key']+'_'+label+'_base'
        fixture=next(x for x in fixtures if x['key']==key);nodes=fixture['nodes'];N=len(nodes)
        pn=s.a.rows(r.c.v.LOCAL/case['key']/label/'base_poisson_nodes.csv');mask=np.array([int(n['bcpsi'])==0 for n in pn]);masks[case['key']]=mask
        J=s.matrix(r.LOCAL/'fixed'/case['key']/'zero/jacobian.csv',N);old=s.matrix(r.c.v.LOCAL/case['key']/label/'jacobian.csv',N)
        diff=J-old;error=np.linalg.norm(diff.data)/np.linalg.norm(old.data);assert error<1e-12,error
        matrices[case['key']]=J;refs[case['key']]=fixture
        for block,field in enumerate(('psi','n','p')):
            direction=np.sin(np.arange(N)*.713+1.17);direction[np.array([bool(int(n['bc'+field])) for n in pn])]=0
            for step in (1e-4,1e-5,1e-6):
                pair=[]
                for sign in (-1,1):
                    job=copy.deepcopy(fixture);job['key']=case['key']+'_'+field+'_'+str(step)+'_'+str(sign)
                    for i,n in enumerate(job['nodes']):n[3+block]+=sign*step*direction[i]
                    jobs.append(job);pair.append(job)
                tests.append(dict(key=case['key'],field=field,block=block,step=step,pair=pair,J_identity_relative=error))
    s.a.write(LOCAL/'inputs.json',jobs)
    cmd=[str(r.c.LOCAL/'kernel_check/check.exe'),str(LOCAL/'inputs.json'),str(LOCAL/'results.csv')];s.a.write(LOCAL/'command.json',cmd)
    result=subprocess.run(cmd,env=s.V.environment(),capture_output=True,text=True);assert result.returncode==0,result.stderr
    values={}
    for row in s.a.rows(LOCAL/'results.csv'):
        if row['mode']=='packed':values.setdefault(row['case'],[]).append(D(row['value']))
    rows=[]
    with localcontext() as ctx:
        ctx.prec=80
        for test in tests:
            minus,plus=test['pair'];N=len(plus['nodes']);h=D.from_float(test['step']);block=test['block'];key=test['key'];mask=masks[key]
            actual=np.zeros(3*N)
            for i in range(N):actual[block*N+i]=float((D.from_float(plus['nodes'][i][3+block])-D.from_float(minus['nodes'][i][3+block]))/(2*h))
            predicted=np.asarray(matrices[key]@actual)[:N][mask]
            difference=np.array([float((y-x)/(2*h)) for x,y in zip(values[minus['key']],values[plus['key']])])[mask]
            error=np.linalg.norm(difference-predicted)/np.linalg.norm(predicted);assert error<=1e-6,(key,test['field'],test['step'],error)
            rows.append(dict(key=key,field=test['field'],step=test['step'],relative=float(error),J_identity_relative=test['J_identity_relative'],qualified=True))
    s.a.write_csv(OUT/'checks.csv',rows);summary=dict(checks=len(rows),qualified=len(rows),max_relative=max(x['relative'] for x in rows),gate=1e-6,scope='Poisson block only, three state-direction blocks and three step sizes at two frozen low-Vd bases; actual representable central displacement used. Does not claim arbitrary carrier cross-block audit or production integration.')
    s.a.write(OUT/'summary.json',summary);s.d.matrix.freeze(OUT/'evidence.json',[Path(__file__).resolve(),r.OUT/'preflight_evidence.json',r.c.OUT/'kernel_check/evidence.json',OUT/'checks.csv',OUT/'summary.json']+[p for p in LOCAL.iterdir() if p.is_file()]);print(summary,flush=True)

if __name__=='__main__':run()
