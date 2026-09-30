"""Fresh full-J source preflight and two-amplitude response with packed Poisson."""
import argparse
import subprocess
from pathlib import Path
import numpy as np
import validate_simplemos_poisson_half_source_20260910 as h
c=h.c;s=c.s;w=c.m.w
LOCAL=c.LOCAL/'response';OUT=c.OUT/'response';RUNNER=LOCAL/'jacobian_runner.exe'

def prepare():
    s.a.verify(c.OUT/'build_evidence.json');s.a.verify(c.OUT/'kernel_check/evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    link=s.a.read(c.LOCAL/'link_command.json');old=next(x for x in link if 'vela_example_runner.cpp.obj' in x);link[link.index(old)]=str(s.LOCAL/'jacobian_adapter/runner.o');link[link.index('-o')+1]=str(RUNNER)
    s.a.write(LOCAL/'link_command.json',link);r=subprocess.run(link,env=s.V.environment(),capture_output=True,text=True);(LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    cases=s.a.read(w.OUT/'settled_source/contract.json')['cases'];files=[Path(__file__).resolve(),c.OUT/'build_evidence.json',c.OUT/'kernel_check/evidence.json',s.OUT/'jacobian_build_evidence.json',LOCAL/'link_command.json',LOCAL/'link.log',RUNNER]
    for case in cases:
        base=Path(case['baseline']);cfg=s.a.read(base/'config.json')
        for label in ('zero','unit'):
            dest=LOCAL/'fixed'/case['key']/label;files+=s.V.probes(cfg,dest,base/'state.csv')
            deck=s.a.read(dest/'functional.json');deck.update(simulation_type='parameter_jacobian',output_csv=str(dest/'jacobian.csv'));s.a.write(dest/'jacobian.json',deck);files.append(dest/'jacobian.json')
        files += [Path(case['source_file']),base/'state.csv',base/'all_row.csv']
    s.a.write(OUT/'contract.json',dict(cases=cases,gates=s.a.read(w.OUT/'settled_source/contract.json')['gates'],scope='Fresh full Jacobian and fixed source insertion on the identical original settled baselines, using the same packed-Poisson candidate binary. Original full/half DC outputs retained in their own directories. No fitted direction or relaxed thresholds.'))
    s.a.write_csv(OUT/'scope.csv',s.a.rows(w.OUT/'settled_source/scope.csv'));s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json',OUT/'scope.csv'])
    def execute(path,case,alpha):
        env=s.env(case,alpha);env.update(VELA_STABLE_MERIT='1',VELA_POISSON_PRECISION='packed')
        return s.V.execute(path,RUNNER if path.stem=='jacobian' else c.RUNNER,env)
    s.LOCAL=LOCAL;s.OUT=OUT;s.execute=execute;s.preflight()

def analyze():
    for p in (OUT/'preflight_evidence.json',c.OUT/'dc_evidence.json',h.OUT/'dc_evidence.json'):s.a.verify(p)
    assert all(r['qualified']=='True' for r in s.a.rows(OUT/'preflight.csv'))
    dc=[r for r in s.a.rows(c.OUT/'dc.csv') if r['mode']=='packed']+s.a.rows(h.OUT/'dc.csv')
    tangent=s.a.rows(OUT/'tangent.csv');native=s.a.rows(s.OUT/'native/comparison/terminal_derivatives.csv');results=[];ports=[]
    for case in s.a.read(OUT/'contract.json')['cases']:
        geo,mask=s.fields.support(case);ids=np.flatnonzero(mask)
        def state(path):
            rows=s.d.ordered(path,geo.count);return [[s.fields.physical(rows[i],field) for i in ids] for field in ('psi','phin','phip')]
        def delta(x,y):return np.array([[float(a-b) for a,b in zip(xx,yy)] for xx,yy in zip(x,y)])
        base=state(Path(case['baseline'])/'state.csv');zero=state(c.LOCAL/'packed'/case['key']/'zero/state.csv');noise=np.linalg.norm(delta(zero,base));ds={};odd={}
        for amp,root in (('full',c.LOCAL/'packed'),('half',h.LOCAL)):
            for sign in ('plus','minus'):ds[sign,amp]=delta(state(root/case['key']/(sign+'_'+amp)/'state.csv'),zero)
            odd[amp]=(ds['plus',amp]-ds['minus',amp])/2
        tr=[r for r in tangent if r['key']==case['key']];tv=np.array([[float(r[f+'_V']) for r in tr] for f in ('psi','phin','phip')]);linearity=np.linalg.norm(odd['full']-2*odd['half'])/np.linalg.norm(odd['full'])
        dc_ok=all(r['qualified']=='True' for r in dc if r['key']==case['key'])
        assert sum(r['key']==case['key'] for r in dc)==5
        for amp,alpha,root in (('full',.001,c.LOCAL/'packed'),('half',.0005,h.LOCAL)):
            pred=np.linalg.norm(odd[amp]-alpha*tv)/np.linalg.norm(alpha*tv);even=np.linalg.norm((ds['plus',amp]+ds['minus',amp])/2)/np.linalg.norm(odd[amp]);snr=np.linalg.norm(odd[amp])/max(noise,1e-300)
            good=dc_ok and pred<=1e-3 and linearity<=1e-3 and even<=.01 and snr>=100
            results.append(dict(key=case['key'],amplitude=amp,DC_qualified=dc_ok,prediction_relative=float(pred),two_amplitude_relative=float(linearity),even_over_odd=float(even),signal_over_drift=float(snr),qualified=bool(good)))
            plus=s.a.read(root/case['key']/('plus_'+amp)/'config.status.json')['contact_currents_A_per_um'];minus=s.a.read(root/case['key']/('minus_'+amp)/'config.status.json')['contact_currents_A_per_um']
            for port in ('drain','substrate'):
                derivative=(plus[port]-minus[port])/(2*alpha);nr=next(r for r in native if r['key']==case['key'] and r['amplitude']==amp and r['port']==port);target=float(nr['derivative_A_per_um']);error=abs(derivative/target-1)
                ports.append(dict(key=case['key'],amplitude=amp,port=port,vela_derivative_A_per_um=derivative,native_derivative_A_per_um=target,relative=error,qualified=bool(dc_ok and error<=1e-3)))
    for name,data in (('response',results),('ports',ports)):s.a.write_csv(OUT/(name+'.csv'),data)
    summary=dict(DC=len(dc),qualified_DC=sum(r['qualified']=='True' for r in dc),field_responses=len(results),qualified_fields=sum(r['qualified'] for r in results),port_responses=len(ports),qualified_ports=sum(r['qualified'] for r in ports),native_new_runs=0,production_changed=False)
    s.a.write(OUT/'summary.json',summary);s.d.matrix.freeze(OUT/'evidence.json',[OUT/'preflight_evidence.json',c.OUT/'dc_evidence.json',h.OUT/'dc_evidence.json',s.OUT/'native/comparison/evidence.json',OUT/'response.csv',OUT/'ports.csv',OUT/'summary.json']);print(summary,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','analyze'));globals()[p.parse_args().action]()
