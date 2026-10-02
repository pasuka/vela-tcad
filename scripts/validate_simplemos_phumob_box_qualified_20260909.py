"""Repeat failed fixed-state preflight using qualified, material-consistent states."""
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import validate_simplemos_phumob_box_20260909 as v
import analyze_simplemos_phumob_box_20260909 as report

a,d,q=v.a,v.d,v.q
PARENT_OUT,PARENT_LOCAL=v.OUT,v.LOCAL
OUT,LOCAL=PARENT_OUT/'qualified_state',PARENT_LOCAL/'qualified_state'


def prepare():
    for name in ('fixed_analysis_evidence.json','probe_evidence.json'):
        a.verify(PARENT_OUT/name)
    original=a.read(PARENT_OUT/'contract.json');contract=copy.deepcopy(original);files=[]
    failures=[]
    for job in contract['fixed']:
        old=Path(job['dir']);dest=LOCAL/'fixed'/job['key']
        source=q.LOCAL/'dc/old_slotboom'/job['case']/'vela'/f"vg_{job['index']:03d}"/'attempt_0/state.csv'
        assert source.exists(),source
        oldcfg=a.read(old/'edges.json');stored={int(r['node_id']):r for r in a.rows(Path(oldcfg['state_file']))}
        ratios=[]
        for row in a.rows(old/'edges.csv'):
            for k in (0,1):
                node=int(row[f'node{k}'])
                for carrier,column in [('electron','electrons_m3'),('hole','holes_m3')]:
                    value=float(row[f'{carrier}_density{k}_m3']);raw=float(stored[node][column])
                    if value>0 and raw>0:ratios.append(value/raw)
        failures.append(dict(key=job['key'],old_seed=str(oldcfg['state_file']),minimum_density_ratio=min(ratios),maximum_density_ratio=max(ratios),
            diagnosis='Initial seed retained older intrinsic density. Newton/SG recomputed populations under current matched-ni material; raw cell/edge-mobility probes used stored populations.',replacement=str(source)))
        for name in ('cells','edges','mobility','jvp'):
            cfg=a.read(old/(name+'.json'));cfg['state_file']=str(source);cfg['output_csv']=str(dest/(name+'.csv'))
            if name=='jvp':cfg['row_output_csv']=str(dest/'rows.csv')
            a.write(dest/(name+'.json'),cfg);files += [old/(name+'.json'),dest/(name+'.json')]
        job['dir']=str(dest);job['seed']=str(source);files.append(source)
    contract['amendment']='Only fixed-state audit inputs change to prior qualified final OldSlotboom states. DC initialization jobs and all physical/numerical gates are unchanged. Initial failed outputs remain in parent directory.'
    contract['density_preflight_relative_gate']=1e-10
    a.write(OUT/'contract.json',contract);a.write_csv(OUT/'initial_state_failure_ledger.csv',failures)
    d.matrix.freeze(OUT/'input_freeze.json',files+[OUT/'contract.json',OUT/'initial_state_failure_ledger.csv',Path(__file__).resolve(),
        PARENT_OUT/'input_freeze.json',PARENT_OUT/'fixed_analysis_evidence.json',PARENT_OUT/'probe_evidence.json'])
    print(failures,flush=True)


def fixed():
    a.verify(OUT/'input_freeze.json')
    d.matrix.freeze(OUT/'fixed_freeze.json',[OUT/'input_freeze.json',q.run.RUNNER,Path(__file__).resolve(),
        v.REPO/'src/equation/CoupledDDAssembler.cpp',v.REPO/'include/vela/equation/AssemblerUtils.h',v.REPO/'src/physics/MobilityModel.cpp'])
    def one(job):
        dest=Path(job['dir']);rows=[]
        for name in ('cells','edges','mobility'):
            status=q.run.V.execute(dest/(name+'.json'),q.run.RUNNER,q.run.V.environment())
            rows.append(dict(key=job['key'],probe=name,**status))
            if status['exit_code']!=0:return rows
        cfg=a.read(dest/'edges.json');stored={int(r['node_id']):r for r in a.rows(Path(cfg['state_file']))}
        differences=[]
        for row in a.rows(dest/'edges.csv'):
            for k in (0,1):
                node=int(row[f'node{k}'])
                for carrier,column in [('electron','electrons_m3'),('hole','holes_m3')]:
                    value=float(row[f'{carrier}_density{k}_m3']);raw=float(stored[node][column])
                    if value>0 and raw>0:differences.append(abs(value/raw-1))
                    else:assert value==raw==0
        preflight=dict(max_density_relative=max(differences),gate=1e-10,qualified=max(differences)<=1e-10)
        a.write(dest/'density_preflight.json',preflight)
        if not preflight['qualified']:return rows
        status=q.run.V.execute(dest/'jvp.json',q.run.RUNNER,q.run.V.environment())
        rows.append(dict(key=job['key'],probe='jvp',**status));return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['fixed']) for r in group]
    a.write(OUT/'fixed_execution.json',dict(results=rows))
    d.matrix.freeze(OUT/'fixed_raw_evidence.json',[OUT/'fixed_freeze.json',OUT/'fixed_execution.json']+[p for p in (LOCAL/'fixed').rglob('*') if p.is_file()])
    assert len(rows)==32 and all(r['exit_code']==0 for r in rows)
    print('Qualified input preflight and all 32 probes completed',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','fixed','analyze_fixed','probes','legacy','candidate','compare','post_jvp'))
    action=parser.parse_args().action
    v.OUT,v.LOCAL=OUT,LOCAL;report.OUT,report.LOCAL=OUT,LOCAL
    if action in ('prepare','fixed'):globals()[action]()
    elif action in ('probes','compare','post_jvp'):getattr(report,action)()
    elif action in ('legacy','candidate'):
        a.verify(OUT/'probe_evidence.json');assert a.read(OUT/'probe_summary.json')['failed']==0
        v.solve(action)
    else:v.analyze_fixed()
