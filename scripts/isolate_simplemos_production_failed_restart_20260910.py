"""Same frozen failed state, one-axis numerical controls; never cohort recovery."""
import copy
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import validate_simplemos_production_weighted_20260910 as v

a,d,run=v.a,v.d,v.run
LOCAL=v.LOCAL/'failed_restart_isolation';OUT=v.OUT/'failed_restart_isolation'


def main():
    v.configure()
    base=v.LOCAL/'dc/phumob/m65_n23_vd_1p000000_endpoint/native/vg_000/attempt_1'
    cfg=a.read(base/'config.json');jobs=[];files=[Path(__file__).resolve(),v.OUT/'validation_freeze.json',base/'state.csv',base/'config.status.json']
    for name,changes in (('combined',{}),('double_poisson',{'poisson_residual_precision':'double'}),
                         ('ordinary_merit',{'stable_merit_comparison':False}),('unprojected_contacts',{'exact_dirichlet_updates':False})):
        dest=LOCAL/name;deck=copy.deepcopy(cfg);deck['solver'].update(changes)
        deck['solver']['local_update_diagnostics']=dict(enabled=True,nodes=[800,982,786],csv_file=str(dest/'updates.csv'),first_iterations=200,every_iterations=1)
        deck['output_state_file']=str(dest/'state.csv');p=dest/'original.json';a.write(p,deck)
        jobs.append(dict(case='m65_n23_vd_1p000000_endpoint',device='n23',vd=1.,vg=0.,index=0,model=name,arm='diagnostic',original_config=str(p),seed=str(base/'state.csv')))
        files.append(p)
    a.write(OUT/'contract.json',dict(jobs=jobs,scope='Additional same-state reload diagnostics only. Not a second allowed cohort reload; no 16-point qualification credit.',gate=1e-6))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])
    run.LOCAL=LOCAL;run.OUT=v.OUT
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(lambda j:run.attempt(j,Path(j['seed']),0),jobs))
    run.v.csv_union(OUT/'results.csv',rows)
    d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'results.csv']+[p for p in LOCAL.rglob('*') if p.is_file()])
    print(rows,flush=True)


if __name__=='__main__':main()
