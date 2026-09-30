"""Step-convergence audit for the retained finite-amplitude HFS Jv failures.

Does not change C++ source, matrix tolerances, DC gates or old failure files.
Only the directions that failed at 1e-6 V are repeated at 1e-10/1e-12 V.
The original 1e-20 V results must also pass, with every weak entry retained.
"""
import argparse,math,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_hfs_candidate_20260914 as h
a,d=h.a,h.d
L=h.L/'step_convergence';O=h.O/'step_convergence'

def entries(path):
    out=[]
    for r in a.rows(path):
        an=float(r['analytic_derivative']);fd=float(r['finite_difference_derivative']);scale=max(abs(an),abs(fd));rel=abs(an-fd)/scale if scale else 0.
        out.append(dict(**r,true_relative=rel,qualified=math.isfinite(rel) and rel<=1e-4))
    return out

def prepare():
    assert not (O/'contract.json').exists()
    a.write(O/'contract.json',dict(steps_V=[1e-10,1e-12],row_relative_gate=1e-4,original_microstep_V=1e-20,weak_entries_excluded=0,
        scope='Resolve finite-amplitude challenges by same-matrix branch-local step convergence. Preserve the initial 1e-6 V failure and do not certify finite-amplitude response from a tangent check.',
        required='All original 1e-20 entries pass. Every failed direction is repeated with all rows at both smaller amplitudes. Both new amplitudes pass the original relative tolerance. The unchanged source/binary and all DC gates remain frozen.',
        production_change=False,acceptance_tolerance_change=False))
    d.matrix.freeze(O/'contract_evidence.json',[O/'contract.json',h.O/'candidate_freeze.json',Path(__file__).resolve()])

def run_case(cc):
    source=h.L/'jvp'/cc['key'];assert (source/'config.status.json').exists()
    old=entries(source/'rows.csv');bad={r['direction'] for r in old if not r['qualified']}
    assert all(r['qualified'] for r in old if float(r['amplitude_V'])<1e-10),'Microstep matrix failure requires implementation investigation'
    cfg=a.read(source/'config.json');directions=[r for r in cfg['directions'] if r['name'] in bad]
    if not directions:return dict(key=cc['key'],repeated_directions=0,qualified=True)
    dest=L/cc['key'];out=O/cc['key']
    assert not (out/'evidence.json').exists()
    cfg.update(output_csv=str(dest/'jvp.csv'),row_output_csv=str(dest/'rows.csv'),directions=[dict(r,name=r['name']+'_refine_'+str(step),amplitude_V=step) for r in directions for step in (1e-10,1e-12)])
    a.write(dest/'config.json',cfg);d.matrix.freeze(out/'input_evidence.json',[dest/'config.json',source/'config.json',source/'rows.csv',O/'contract_evidence.json',h.O/'candidate_freeze.json'])
    s=h.V.execute(dest/'config.json',h.RUNNER,h.V.environment());assert s['exit_code']==0,s
    rows=entries(dest/'rows.csv');bad=[r for r in rows if not r['qualified']]
    a.write(out/'failures.json',bad)
    summary=dict(key=cc['key'],repeated_directions=len(directions),checks=len(rows),nonzero=sum(float(r['analytic_derivative'])!=0 or float(r['finite_difference_derivative'])!=0 for r in rows),max_relative=max(r['true_relative'] for r in rows),failures=len(bad),qualified=not bad)
    a.write(out/'summary.json',summary);d.matrix.freeze(out/'evidence.json',[out/'input_evidence.json',out/'failures.json',out/'summary.json']+[f for f in dest.iterdir() if f.is_file()]);print(summary,flush=True);return summary

def run(keys):
    a.verify(O/'contract_evidence.json');a.verify(h.O/'candidate_freeze.json')
    cc=[c for c in a.read(h.O/'inputs.json')['cases'] if c['key'] in keys]
    assert len(cc)==len(keys)
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(run_case,cc))
    assert all(r['qualified'] for r in rows),rows

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','run'));p.add_argument('--keys',nargs='*');args=p.parse_args();prepare() if args.action=='prepare' else run(args.keys)
