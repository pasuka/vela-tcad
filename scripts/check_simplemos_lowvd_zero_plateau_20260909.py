"""Two predeclared zero-source restarts on each low-Vd control; no rebasing claim."""
from pathlib import Path
import math
import numpy as np
import validate_simplemos_restart_source_followup_20260909 as w
s=w.s;LOCAL=w.LOCAL/'zero_plateau';OUT=w.OUT/'zero_plateau'

def drift(first,second,c):
    geo,mask=s.fields.support(c);a=s.d.ordered(first,geo.count);b=s.d.ordered(second,geo.count)
    return math.sqrt(math.fsum(float(s.fields.physical(x,k)-s.fields.physical(y,k))**2 for x,y,keep in zip(a,b,mask) if keep for k in ('psi','phin','phip')))

def run():
    s.a.verify(w.OUT/'evidence.json');cases=[c for c in s.a.read(w.OUT/'contract.json')['cases'] if c['vd']==.05]
    s.a.write(OUT/'contract.json',dict(cases=cases,restarts=2,scope='Follow zero-source trajectories without changing source, amplitude or gates. Retain original 6/8 calibration; this is a stationarity diagnostic before any new baseline experiment.',stationarity_gate='Each of the two successive free-Si field drifts <= original half-amplitude signal / 100.',original_gates=s.a.read(w.OUT/'contract.json')['gates']))
    s.d.matrix.freeze(OUT/'freeze.json',[Path(__file__).resolve(),OUT/'contract.json',w.OUT/'evidence.json',w.RUNNER])
    rows=[]
    for c in cases:
        zero=w.LOCAL/'dc'/c['key']/'zero';response=next(r for r in s.a.rows(w.OUT/'response.csv') if r['key']==c['key'] and r['amplitude']=='half')
        original_drift=drift(zero/'state.csv',Path(c['baseline'])/'state.csv',c);limit=original_drift*float(response['signal_over_drift'])/100
        base=zero
        for iteration in range(1,3):
            dest=LOCAL/c['key']/str(iteration);cfg=s.a.read(base/'config.json');cfg.update(state_file=str(base/'state.csv'),output_state_file=str(dest/'state.csv'))
            s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest);s.d.matrix.freeze(dest/'freeze.json',[base/'state.csv',dest/'config.json',OUT/'freeze.json'])
            status=s.V.execute(dest/'config.json',w.RUNNER,s.env(c,0.))
            for name in ('all_row','acceptance_edges'):s.V.execute(dest/(name+'.json'),w.RUNNER,s.env(c,0.))
            q=s.prior.q.run.w.old.prior.old.qualify(c,dest);delta=drift(dest/'state.csv',base/'state.csv',c)
            row=dict(key=c['key'],restart=iteration,drift_V=delta,limit_V=limit,stationary=delta<=limit,**q);rows.append(row);print(row,flush=True);base=dest
    s.a.write_csv(OUT/'results.csv',rows);s.d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'results.csv']+[p for p in LOCAL.rglob('*') if p.is_file()])

if __name__=='__main__':run()
