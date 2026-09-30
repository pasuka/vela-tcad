"""Pipeline a qualified case's response while other DC cases finish.

Uses exactly the later full-cohort paths, templates, seeds and solver protocol.
The cohort driver rechecks/reuses these cached attempts, never replaces them.
"""
import argparse,math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_enormal_same_response_20260912 as r
q,a,d=r.q,r.a,r.d

def main(device,vd):
    a.verify(q.O/'candidate_freeze.json');a.verify(q.O/'real_column_evidence.json')
    assert a.read(q.O/'real_column_summary.json')['qualified']
    cc=next(c for c in a.read(q.O/'inputs.json')['cases'] if c['device']==device and c['vd']==vd and c['index']==40)
    contract=a.read(r.n.O/'native_contract.json');stage='pilot' if cc['case']==contract['pilot_case'] else 'rest'
    a.verify(r.n.O/f'{stage}_evidence.json');assert a.read(r.n.O/f'{stage}_summary.json')['all_qualified']
    selected=[];files=[]
    for arm in ('phumob_seed','native_seed'):
        paths=sorted((q.L/'vela'/cc['case']/arm/'vg_040').glob('attempt_*/result.json'))
        eligible=[(p,a.read(p)) for p in paths if a.read(p)['qualified']];assert eligible,(device,vd,arm)
        p,v=eligible[0];selected.append(v);files += [p,Path(v['dest'])/'state.csv']
    geo,mask=q.V.m.previous.prior.support(cc)
    delta=q.V.m.previous.prior.delta_states(*[d.ordered(Path(v['dest'])/'state.csv',geo.count) for v in selected],mask)
    error=abs(selected[0]['current_A_per_um']/selected[1]['current_A_per_um']-1)
    assert q.v.dual_qualified(True,True,delta,error)
    jobs=[];seed=Path(selected[0]['dest'])/'state.csv'
    for label,value in [('zero',0.),('minus_large',-.001),('minus_small',-.0005),('plus_small',.0005),('plus_large',.001)]:
        cfg=q.config_case(cc);surface=cfg['solver']['mobility']['surface']
        surface['acoustic_factor']=surface['roughness_factor']=1+value
        path=r.L/'inputs'/cc['key']/label/'template.json';q.v.write_same(path,cfg)
        jobs.append(dict(**{k:v for k,v in cc.items() if k!='template'},template=str(path),label=label,seed=str(seed),delta=value))
        files.append(path)
    freeze=r.O/'case_prerequisites'/f'{cc["key"]}.json'
    files += [Path(__file__).resolve(),q.O/'candidate_freeze.json',q.O/'real_column_evidence.json',r.n.O/f'{stage}_evidence.json']
    if not freeze.exists():d.matrix.freeze(freeze,files)
    a.verify(freeze);q.configure();q.v.LOCAL=r.L;q.v.OUT=r.O
    def run(j):return q.v.solve_target(j,j['index'],j['label'],Path(j['seed']))
    attempts=run(jobs[0]);zero=next((v for v in attempts if v['qualified']),attempts[-1])
    drift=abs(zero.get('current_A_per_um',math.nan)/selected[0]['current_A_per_um']-1)
    assert zero['qualified'] and drift<=q.v.GATES['initialization_Id_relative'],zero
    with ThreadPoolExecutor(max_workers=2) as pool:attempts += [v for group in pool.map(run,jobs[1:]) for v in group]
    out=r.O/'case_attempts'/f'{cc["key"]}.csv';q.v.csv_union(out,attempts)
    d.matrix.freeze(r.O/'case_evidence'/f'{cc["key"]}.json',[freeze,out]+[f for f in (r.L/'vela'/cc['case']).rglob('*') if f.is_file()])
    print(dict(case=cc['key'],attempts=len(attempts),failed=sum(not v['qualified'] for v in attempts),release='Case response attempts only. Full-cohort comparison remains gated.'),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--device',required=True);ap.add_argument('--vd',type=float,required=True);args=ap.parse_args();main(args.device,args.vd)
