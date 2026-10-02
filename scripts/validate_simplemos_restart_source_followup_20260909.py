"""Isolate cancellation-free restarts in the already calibrated source experiment."""
import argparse
import copy
import subprocess
from pathlib import Path
import calibrate_simplemos_phumob_local_source_20260909 as s

OLD=s.LOCAL/'scaled_source';OLDOUT=s.OUT/'scaled_source'
LOCAL=s.REPO/'build-release/phumob_restart_followup_20260909'
OUT=s.REPO/'reference_tcad/simplemos_sentaurus2022/phumob_restart_followup_20260909'
RUNNER=LOCAL/'source_runner.exe'

def prepare():
    for p in (OLDOUT/'evidence.json',s.OUT/'restart_coordinates/build_evidence.json',s.OUT/'native/comparison/evidence.json'):s.a.verify(p)
    LOCAL.mkdir(exist_ok=False)
    cmd=s.a.read(OLD/'link_command.json');cmd.insert(1,str(s.LOCAL/'restart_coordinates/newton.o'));cmd[cmd.index('-o')+1]=str(RUNNER)
    s.a.write(LOCAL/'link_command.json',cmd)
    r=subprocess.run(cmd,env=s.V.environment(),capture_output=True,text=True);(LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    files=[Path(__file__).resolve(),s.OUT/'restart_coordinates/build_evidence.json',OLDOUT/'evidence.json',s.OUT/'native/comparison/evidence.json']+[Path(x) for x in cmd[1:] if Path(x).is_file() and Path(x).is_relative_to(s.REPO)]
    cases=s.a.read(OLDOUT/'contract.json')['cases']
    for c in cases:
        for job in c['jobs']:
            original=Path(job['dest']);dest=LOCAL/'dc'/c['key']/job['label'];cfg=s.a.read(original/'config.json')
            cfg['output_state_file']=str(dest/'state.csv');s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest);job['dest']=str(dest)
            files += [original/'config.json',Path(cfg['state_file'])]+list(dest.glob('*.json'))
        dest=LOCAL/'identity'/c['key'];cfg=s.a.read(Path(c['jobs'][0]['dest'])/'config.json');cfg['output_state_file']=str(dest/'state.csv');cfg['solver']['max_iter']=0;cfg['solver']['carrier_row_convergence']['min_newton_max_iter']=0
        s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest);files+=list(dest.glob('*.json'))
    s.a.write(OUT/'contract.json',dict(cases=cases,change='Link the previously tested cancellation-free reference subtraction object with the identical calibrated source assembler. No other solver or source changes.',gates=s.a.read(OLDOUT/'contract.json')['gates'],reuse_tangent=str(OLDOUT/'tangent.csv'),production_modified=False,finite_volume_replacement=False))
    files += [OUT/'contract.json',OLDOUT/'tangent.csv',OLDOUT/'port_calibration/ports.csv',LOCAL/'link_command.json',LOCAL/'link.log']
    s.d.matrix.freeze(OUT/'freeze.json',files)

def identity():
    s.a.verify(OUT/'freeze.json');rows=[]
    for c in s.a.read(OUT/'contract.json')['cases']:
        dest=LOCAL/'identity'/c['key'];geo,mask=s.fields.support(c)
        status=s.V.execute(dest/'config.json',RUNNER,s.env(c,0.));assert status['iterations']==0
        s.V.execute(dest/'all_row.json',RUNNER,s.env(c,0.))
        before=s.d.ordered(Path(c['baseline'])/'all_row.csv',geo.count);after=s.d.ordered(dest/'all_row.csv',geo.count)
        residual=all(x[k]==y[k] for x,y in zip(before,after) for k in ('electron_residual','hole_residual','electron_recombination','hole_recombination'))
        before=s.d.ordered(Path(c['baseline'])/'state.csv',geo.count);after=s.d.ordered(dest/'state.csv',geo.count)
        delta=s.fields.delta_states(after,before,mask)
        identity=all(x[k]==y[k] for x,y in zip(before,after) for k in ('psi','electron_qf_increment_V','hole_qf_increment_V'))
        rows.append(dict(key=c['key'],residual_identity=residual,state_identity=identity,**delta,qualified=residual and identity));print(rows[-1],flush=True)
    s.a.write_csv(OUT/'identity.csv',rows);s.d.matrix.freeze(OUT/'identity_evidence.json',[OUT/'freeze.json',OUT/'identity.csv']+[p for p in (LOCAL/'identity').rglob('*') if p.is_file()])

def bind():
    s.LOCAL=LOCAL;s.OUT=OUT;s.RUNNER=RUNNER

def run():
    s.a.verify(OUT/'identity_evidence.json');assert all(r['qualified']=='True' for r in s.a.rows(OUT/'identity.csv'))
    # The unchanged source Jacobian and four bitwise state/residual identities
    # justify reusing the frozen baseline tangent, not fitting a new direction.
    s.a.write_csv(OUT/'preflight.csv',s.a.rows(OUT/'identity.csv'))
    s.a.write_csv(OUT/'tangent.csv',s.a.rows(OLDOUT/'tangent.csv'))
    s.d.matrix.freeze(OUT/'preflight_evidence.json',[OUT/'identity_evidence.json',OUT/'preflight.csv',OUT/'tangent.csv',OLDOUT/'preflight_evidence.json'])
    bind();s.run();s.analyze()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','identity','run'));globals()[p.parse_args().action]()
