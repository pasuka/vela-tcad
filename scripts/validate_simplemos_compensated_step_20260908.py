"""Psi carry A/B within each comparator; no residual/J or SRH changes."""
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import build_simplemos_compensated_step_20260908 as b
import validate_simplemos_stable_merit_20260908 as first

a=b.a;d=b.d;t=b.t;c=b.c;v=b.v;p=b.p;OUT=b.OUT;LOCAL=b.LOCAL/'validation';RUNNER=b.RUNNER


def execute(path,job):
    e=t.env(job['case'],job['axis'],1.)
    if job['mode']=='stable':e.update(VELA_STABLE_MERIT='1',VELA_STABLE_MERIT_TRACE=str(Path(job['dest'])/'merit_trace.csv'))
    if job['carry']:e['VELA_COMPENSATED_PSI_STEP']='1'
    return v.execute(path,RUNNER,e)


def prepare():
    a.verify(OUT/'build_evidence.json');a.verify(first.OUT/'freeze.json');jobs=[];files=[Path(__file__).resolve(),OUT/'build_evidence.json',first.OUT/'freeze.json',RUNNER]
    for previous in a.read(first.OUT/'contract.json')['jobs']:
        if previous['name'] not in ('failure_original_restart','control_finite'):continue
        carry_options=[True]
        if previous['case']['device']=='n23' and ((previous['name']=='failure_original_restart' and previous['mode']=='stable') or (previous['name']=='control_finite' and previous['case']['vd']==.05 and previous['mode']=='legacy')):carry_options.append(False)
        for carry in carry_options:
            source=Path(previous['dest']);dest=LOCAL/previous['case']['key']/previous['name']/previous['mode']/('carry' if carry else 'identity')
            cfg=a.read(source/'config.json');cfg['output_state_file']=str(dest/'state.csv');cfg['solver']['local_update_diagnostics']['csv_file']=str(dest/'updates.csv')
            a.write(dest/'config.json',cfg);v.post_config(cfg,dest);files += [dest/(name+'.json') for name in ('config','all_row','acceptance_edges')]+v.probes(cfg,dest/'post',dest/'state.csv')+[Path(cfg['state_file']),source/'config.json']
            jobs.append(dict(case=previous['case'],axis=previous['axis'],name=previous['name'],mode=previous['mode'],carry=carry,dest=str(dest),baseline=str(source)))
    assert len(jobs)==14
    a.write(OUT/'contract.json',dict(jobs=jobs,DC=14,candidate='TwoSum product/carry + TwoSum potential update, commit psi remainder only for accepted candidate; reset before each coupled Newton solve. All residual/J/state exports use the actual represented candidate.',
        isolation='Two original failed histories + four qualified joint-control states; carry on under legacy and stable comparator independently. Two carry-off identity runs. No source volumes or physical fields/formulas changed.',
        gates='Same original 1814-row 1e-6, KCL/Id 1e-8, independent global closure, 200 iterations, state 1e-6 V, density 1e-4, Id 1e-6.',
        boundary='This accumulates lost updates within a solve. It does not implement persistent electrostatic reference coordinates or high precision residual evaluation. Remainders are not exported as new physical state. No automatic restart or default promotion.'))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json']);print('Frozen 14 psi-carry DC',flush=True)


def run():
    a.verify(OUT/'freeze.json')
    def one(job):
        dest=Path(job['dest']);status=execute(dest/'config.json',job);row=dict(key=job['case']['key'],name=job['name'],mode=job['mode'],carry=job['carry'],qualified=False,failure=status.get('failure_reason',''))
        if (dest/'state.csv').exists():
            for n in ('all_row','acceptance_edges'):execute(dest/(n+'.json'),job)
            row.update(c.s.prior.old.qualify(job['case'],dest))
            for n in ('functional','edges','terms'):execute(dest/'post'/(n+'.json'),job)
        print(row['key'],row['name'],row['mode'],row['carry'],row['qualified'],row.get('max_row_ratio'),status.get('iterations'),flush=True);return row
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,a.read(OUT/'contract.json')['jobs']))
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/'dc.csv',[{k:r.get(k,'') for k in keys} for r in rows])


def analyze():
    a.verify(first.OUT/'evidence.json');a.verify(OUT/'freeze.json');rows=a.rows(OUT/'dc.csv');old=a.rows(first.OUT/'dc.csv');checks=[];ports=[];identities=[]
    for job in a.read(OUT/'contract.json')['jobs']:
        dest=Path(job['dest']);base=Path(job['baseline']);case=job['case'];geo,mask=v.m.previous.prior.support(case)
        row=next(r for r in rows if r['key']==case['key'] and r['name']==job['name'] and r['mode']==job['mode'] and (r['carry']=='True')==job['carry']);prior=next(r for r in old if r['key']==case['key'] and r['name']==job['name'] and r['mode']==job['mode'])
        if not job['carry']:
            q=a.read(dest/'config.status.json');o=a.read(base/'config.status.json');identity=a.sha(dest/'state.csv')==a.sha(base/'state.csv') and all(q[k]==o[k] for k in ('iterations','converged','failure_reason'))
            identities.append(dict(key=case['key'],name=job['name'],mode=job['mode'],identity=identity))
        delta=v.m.previous.prior.delta_states(d.ordered(base/'state.csv',geo.count),d.ordered(dest/'state.csv',geo.count),mask);ie=abs(float(row['current_A_per_um'])/float(prior['current_A_per_um'])-1)
        close=max(delta[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and delta['density_max_relative']<=1e-4 and ie<=1e-6
        checks.append(dict(key=case['key'],name=job['name'],mode=job['mode'],carry=job['carry'],**delta,Id_relative=ie,baseline_qualified=prior['qualified']=='True',candidate_qualified=row['qualified']=='True',state_close=close))
        f=a.read(dest/'post/functional.status.json');err=max(abs(float(row['current_A_per_um'])/f['current_A_per_um']-1),abs(f['current_A_per_um']/f['contact_current_extractor_A_per_um']-1));ports.append(dict(key=case['key'],name=job['name'],mode=job['mode'],carry=job['carry'],relative=err,qualified=err<=1e-8))
    for name,data in [('comparison',checks),('ports',ports),('identity',identities)]:a.write_csv(OUT/(name+'.csv'),data)
    summary=dict(DC=len(rows),carry_qualified=sum(r['qualified']=='True' and r['carry']=='True' for r in rows),carry_total=12,identity_passed=sum(r['identity'] for r in identities),
        carry_failures_qualified=sum(r['qualified']=='True' and r['carry']=='True' and r['name'].startswith('failure') for r in rows),failure_cases=4,
        carry_controls_qualified=sum(r['qualified']=='True' and r['carry']=='True' and r['name'].startswith('control') for r in rows),control_cases=8,ports_passed=sum(r['qualified'] for r in ports),production_changes=False)
    a.write(OUT/'summary.json',summary);d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',first.OUT/'evidence.json',OUT/'summary.json']+list(OUT.glob('*.csv'))+[f for f in LOCAL.rglob('*') if f.is_file()]);print(summary,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run','analyze'));globals()[parser.parse_args().action]()
