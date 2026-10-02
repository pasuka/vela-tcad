"""Frozen 24-DC A/B: only stable norm-increment comparison changes."""
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import build_simplemos_stable_merit_20260908 as b

t=b.t;c=b.c;a=b.a;d=b.d;v=b.v;p=b.p;LOCAL=b.LOCAL/'validation';OUT=b.OUT;RUNNER=b.RUNNER


def execute(path,job):
    e=t.env(job['case'],job['axis'],1.)
    if job['mode']=='stable':e.update(VELA_STABLE_MERIT='1',VELA_STABLE_MERIT_TRACE=str(Path(job['dest'])/'merit_trace.csv'))
    return v.execute(path,RUNNER,e)


def prepare():
    a.verify(OUT/'build_evidence.json');a.verify(b.audit.OUT/'final_evidence.json');jobs=[]
    files=[Path(__file__).resolve(),OUT/'build_evidence.json',b.audit.OUT/'final_evidence.json',RUNNER]
    def add(case,axis,name,seed,source,reference):
        for mode in ('legacy','stable'):
            dest=LOCAL/case['key']/name/mode;cfg=a.read(source/'config.json');cfg.update(state_file=str(seed),output_state_file=str(dest/'state.csv'))
            cfg['solver'].pop('local_update_diagnostics',None)
            cfg['solver']['local_update_diagnostics']=dict(enabled=True,csv_file=str(dest/'updates.csv'),nodes=[r['node'] for r in case['mapped_nodes'].values()],first_iterations=200,every_iterations=1)
            assert cfg['solver']['max_iter']==200 and cfg['solver']['carrier_row_convergence']['eps_row']==1e-6 and cfg['solver']['global_continuity_closure']['mode']=='off'
            a.write(dest/'config.json',cfg);v.post_config(cfg,dest)
            files.extend([seed,source/'config.json',dest/'config.json',dest/'all_row.json',dest/'acceptance_edges.json']+v.probes(cfg,dest/'post',dest/'state.csv'))
            jobs.append(dict(case=case,axis=axis,name=name,mode=mode,dest=str(dest),reference=str(reference)))
    for j in a.read(b.audit.OUT/'contract.json')['jobs']:
        case=j['case'];source=Path(j['source'])
        add(case,j['axis'],'failure_original_restart',Path(a.read(source/'config.json')['state_file']),source,source)
        add(case,j['axis'],'failure_saved_state',source/'state.csv',source,source)
    cases={c['key']:c for c in a.read(t.OUT/'contract.json')['cases']}
    for row in a.rows(b.audit.prior.OUT/'selected_states.csv'):
        if row['axis']!='joint':continue
        assert row['qualified']=='True';source=p.REPO/row['path'];case=cases[row['key']]
        add(case,'joint','control_'+row['label'],source/'state.csv',source,source)
    assert len(jobs)==24
    a.write(OUT/'contract.json',dict(jobs=jobs,DC=24,
        candidate='Guarded stable sum((Rtrial-Rbase)*(Rtrial+Rbase)) with compensated long double and exact binary64 integer fallback. Same mathematical strict L2 descent rule and zero-zero exception.',
        isolation='Only decrease comparison changes; no residual/Jacobian/state/volume/limits/acceptance changes. Two failures: original restart input and exported failed state, each legacy/stable. Eight qualified joint endpoints: each legacy/stable.',
        gates=dict(all_free_carrier_rows=1e-6,kcl_over_Id=1e-8,state_potential_V=1e-6,state_density_relative=1e-4,Id_relative=1e-6,port_agreement=1e-8),
        bounded='One 200-iteration solve per job; retain failures. No automatic repeated restart. No production promotion from comparison tests alone.'))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json']);print('Frozen',len(jobs),'merit-only DC',flush=True)


def run():
    a.verify(OUT/'freeze.json')
    def one(job):
        dest=Path(job['dest']);s=execute(dest/'config.json',job);row=dict(key=job['case']['key'],name=job['name'],mode=job['mode'],qualified=False,failure=s.get('failure_reason',''))
        if (dest/'state.csv').exists():
            for name in ('all_row','acceptance_edges'):execute(dest/(name+'.json'),job)
            row.update(c.s.prior.old.qualify(job['case'],dest))
            for name in ('functional','edges','terms'):execute(dest/'post'/(name+'.json'),job)
        print(row['key'],row['name'],row['mode'],row['qualified'],row.get('max_row_ratio'),s.get('iterations'),flush=True);return row
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,a.read(OUT/'contract.json')['jobs']))
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/'dc.csv',[{k:r.get(k,'') for k in keys} for r in rows])


def analyze():
    a.verify(OUT/'freeze.json');jobs=a.read(OUT/'contract.json')['jobs'];dc=a.rows(OUT/'dc.csv');comparisons=[];ports=[];identities=[];trace=[]
    for j in jobs:
        dest=Path(j['dest']);row=next(r for r in dc if r['key']==j['case']['key'] and r['name']==j['name'] and r['mode']==j['mode']);f=a.read(dest/'post/functional.status.json');I=float(row['current_A_per_um'])
        rel=max(abs(I/f['current_A_per_um']-1),abs(f['current_A_per_um']/f['contact_current_extractor_A_per_um']-1));ports.append(dict(key=row['key'],name=j['name'],mode=j['mode'],relative=rel,qualified=rel<=1e-8))
        if j['mode']=='legacy' and j['name']=='failure_original_restart':
            old=a.read(Path(j['reference'])/'config.status.json');new=a.read(dest/'config.status.json')
            identities.append(dict(key=row['key'],identity=a.sha(dest/'state.csv')==a.sha(Path(j['reference'])/'state.csv') and all(old[k]==new[k] for k in ('iterations','converged','failure_reason'))))
        if j['mode']=='stable':
            hist=a.rows(dest/'merit_trace.csv') if (dest/'merit_trace.csv').exists() else []
            trace.append(dict(key=row['key'],name=j['name'],attempts=len(hist),accepted=sum(r['accepted']=='1' for r in hist),exact_fallback=sum(r['exact_fallback']=='1' for r in hist)))
    for j in (x for x in jobs if x['mode']=='stable'):
        case=j['case'];geo,mask=v.m.previous.prior.support(case);stable=Path(j['dest']);legacy=stable.parent/'legacy'
        delta=v.m.previous.prior.delta_states(d.ordered(legacy/'state.csv',geo.count),d.ordered(stable/'state.csv',geo.count),mask)
        rows={mode:next(r for r in dc if r['key']==case['key'] and r['name']==j['name'] and r['mode']==mode) for mode in ('legacy','stable')}
        ie=abs(float(rows['stable']['current_A_per_um'])/float(rows['legacy']['current_A_per_um'])-1)
        close=max(delta[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and delta['density_max_relative']<=1e-4 and ie<=1e-6
        comparisons.append(dict(key=case['key'],name=j['name'],**delta,Id_relative=ie,legacy_qualified=rows['legacy']['qualified']=='True',stable_qualified=rows['stable']['qualified']=='True',state_close=close,
            strict_pair_qualified=close and all(r['qualified']=='True' for r in rows.values())))
    for name,rows in [('comparison',comparisons),('ports',ports),('legacy_identity',identities),('trace_summary',trace)]:a.write_csv(OUT/(name+'.csv'),rows)
    stable=[r for r in dc if r['mode']=='stable'];summary=dict(DC=len(dc),stable_qualified=sum(r['qualified']=='True' for r in stable),stable_total=len(stable),
        stable_failed_path_qualified=sum(r['qualified']=='True' for r in stable if r['name'].startswith('failure')),stable_failed_path_total=4,
        stable_controls_qualified=sum(r['qualified']=='True' for r in stable if r['name'].startswith('control')),stable_controls_total=8,
        legacy_identity=all(r['identity'] for r in identities),ports_qualified=sum(r['qualified'] for r in ports),all_stable_qualified=all(r['qualified']=='True' for r in stable),production_changes=False)
    a.write(OUT/'summary.json',summary);d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'summary.json']+list(OUT.glob('*.csv'))+[x for x in LOCAL.rglob('*') if x.is_file()]);print(summary,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run','analyze'));globals()[parser.parse_args().action]()
