"""Native-initialized 0 -> 1/2 -> 1 constant continuation, original gates.

The two direct native-seed failures stay frozen. This records a separate
initialization path and never changes a failed run's acceptance status.
"""
import argparse,math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_constants_finite_20260908 as c
a=c.a;d=c.d;v=c.v;LOCAL=c.LOCAL/'continuation';OUT=c.OUT/'continuation'

def env(alpha):
    e=v.environment()
    for key in c.NATIVE:e['VELA_CANDIDATE_'+key.upper()+'_RELATIVE']=format(alpha*(c.NATIVE[key]/c.BASE[key]-1),'.17g')
    return e

def prepare():
    a.verify(c.OUT/'validation_evidence.json');cs=c.cases();files=[Path(__file__).resolve(),c.OUT/'validation_evidence.json',c.RUNNER]
    for case in cs:
        case['continuation_jobs']=[];seed=v.m.previous.LOCAL/case['key']/'native_referenced_joint/state.csv';files.append(seed)
        for label,alpha in [('zero',0.),('half',.5),('full',1.)]:
            dest=LOCAL/case['key']/label;cfg=c.old.cfg(case);cfg.update(state_file=str(seed),output_state_file=str(dest/'state.csv'))
            a.write(dest/'config.json',cfg);v.post_config(cfg,dest)
            files += [dest/(n+'.json') for n in ('config','all_row','acceptance_edges')]+v.probes(cfg,dest/'post',dest/'state.csv')
            case['continuation_jobs'].append(dict(label=label,alpha=alpha,dest=str(dest)));seed=dest/'state.csv'
    a.write(OUT/'contract.json',dict(cases=cs,DC=12,route=[0,.5,1],seed='Original native-referenced state, then each previously qualified continuation state.',
        retained_failure='Direct native initialization has two carrier-row line-search rejections; no failed status or original gate changed.',
        gates='Identical DC gates; final dual phi 1e-6 V, density relative 1e-4, Id relative 1e-6 against the direct qualified Vela-seeded endpoint.'))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])

def run():
    a.verify(OUT/'freeze.json')
    def one(case):
        rows=[]
        for job in case['continuation_jobs']:
            dest=Path(job['dest']);e=env(job['alpha']);status=v.execute(dest/'config.json',c.RUNNER,e)
            r=dict(key=case['key'],label=job['label'],alpha=job['alpha'],qualified=False,failure=status.get('failure_reason',''))
            if (dest/'state.csv').exists():
                for n in ('all_row','acceptance_edges'):v.execute(dest/(n+'.json'),c.RUNNER,e)
                r.update(c.s.prior.old.qualify(case,dest))
                if job['label']=='full':
                    for n in ('functional','edges','terms'):v.execute(dest/'post'/(n+'.json'),c.RUNNER,e)
            rows.append(r);print(case['key'],job['label'],r['qualified'],r.get('max_row_ratio'),flush=True)
            if not r['qualified']:break
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/'dc.csv',[{k:r.get(k,'') for k in keys} for r in rows])

def analyze():
    dc=a.rows(OUT/'dc.csv');duals=[]
    for case in c.cases():
        geo,mask=v.m.previous.prior.support(case);root=LOCAL/case['key']/'full';base=c.LOCAL/case['key']/'manual_constants'
        if not (root/'state.csv').exists():continue
        delta=v.m.previous.prior.delta_states(d.ordered(root/'state.csv',geo.count),d.ordered(base/'state.csv',geo.count),mask)
        row=next(r for r in dc if r['key']==case['key'] and r['label']=='full');Id=float(row['current_A_per_um']);original=a.read(base/'config.status.json')['contact_currents_A_per_um']['drain'];err=abs(Id/original-1)
        duals.append(dict(key=case['key'],**delta,Id_relative=err,native_Id_relative=Id/case['native_Id_A_per_um']-1,qualified=row['qualified']=='True' and max(delta[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and delta['density_max_relative']<=1e-4 and err<=1e-6))
    a.write_csv(OUT/'dual.csv',duals)
    summary=dict(DC=len(dc),qualified_DC=sum(r['qualified']=='True' for r in dc),qualified_dual=sum(r['qualified'] for r in duals),qualified=len(dc)==12 and len(duals)==4 and all(r['qualified']=='True' for r in dc) and all(r['qualified'] for r in duals),direct_native_failures_retained=2)
    a.write(OUT/'summary.json',summary);d.matrix.freeze(OUT/'validation_evidence.json',[OUT/'freeze.json',OUT/'summary.json']+list(OUT.glob('*.csv'))+[x for x in LOCAL.rglob('*') if x.is_file()]);print(summary,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run','analyze'));globals()[parser.parse_args().action]()
