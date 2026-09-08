"""Gate perturbation-invariance runs on both unshifted states passing all-row checks."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import validate_simplemos_minority_invariance_20260906 as m

a=m.a;OUT=m.OUT


def prepare():
    a.verify(OUT/'freeze.json')
    trigger=m.LOCAL/'m65_n19_vd_0p050000_endpoint/vg_007/vela/result.json'
    assert not a.read(trigger)['all_row_qualified']
    groups={}
    for job in a.read(OUT/'contract.json')['jobs']:
        groups.setdefault((job['case'],job['vg']),[]).append(job)
    cases=[dict(case=key[0],vg=key[1],jobs=jobs) for key,jobs in groups.items()]
    a.write(OUT/'screen_contract.json',dict(status='frozen_before_remaining_screen',cases=cases,
        trigger='First unshifted Vela state failed the all-row gate; its hole +/- initializations have not completed after prolonged iterations.',
        order='Run existing Vela and coherent native initializations at every target bias. Only if both pass the frozen all-row/legacy/global/KCL gates, run +/- shifted initializations.',
        interrupted='Already running shifted initializations in the first failed target are explicitly interrupted, not called nonconverged or used as final-state evidence.',
        acceptance='All tolerances, max_iter, physical equations and original 32 initial configurations unchanged.',
        interpretation='Failure of the unshifted state prerequisite blocks qualification of initialization invariance at that bias; omitted shifted solves cannot establish nonuniqueness.'))
    m.d.matrix.freeze(OUT/'screen_freeze.json',[Path(__file__).resolve(),OUT/'screen_contract.json',OUT/'freeze.json',trigger])


def run_one(job):
    result=Path(job['config']).parent/'result.json'
    return a.read(result) if result.exists() else m.one(job)


def one(case):
    jobs={j['initialization']:j for j in case['jobs']};done=[];deferred=[]
    for label in ('vela','native'):done.append(run_one(jobs[label]))
    if all(r['all_row_qualified'] for r in done):
        for label in ('hole_plus','hole_minus'):done.append(run_one(jobs[label]))
    else:
        for label in ('hole_plus','hole_minus'):
            deferred.append(dict(case=case['case'],vg=case['vg'],initialization=label,status='not_qualified_prerequisite',
                final_state_used=False,reason='At least one unshifted initialization failed the frozen all-row gate.'))
    print(case['case'],case['vg'],'screen complete',sum(r['all_row_qualified'] for r in done),'/',len(done),flush=True)
    return done,deferred


def run():
    a.verify(OUT/'screen_freeze.json')
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(one,a.read(OUT/'screen_contract.json')['cases']))
    a.write_csv(OUT/'screen_runs.csv',[r for g,_ in results for r in g])
    deferred=[r for _,g in results for r in g]
    if deferred:a.write_csv(OUT/'deferred.csv',deferred)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','run'));globals()[p.parse_args().action]()
