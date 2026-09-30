"""Qualify the unchanged production matrix after the retained step study.

Original finite-amplitude failures stay failed. The branch-local derivative
gate requires two additional amplitudes for every failed direction, plus the
original sub-ULP check. No numerical tolerance, source, or DC gate is changed.
"""
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import refine_simplemos_hfs_jvp_20260914 as r
h=r.h;a=h.a;d=h.d

def gate():
    a.verify(h.O/'candidate_freeze.json');a.verify(h.O/'jvp_evidence.json');a.verify(h.O/'branch_steps_evidence.json');a.verify(r.O/'contract_evidence.json')
    assert not (h.O/'tangent_evidence.json').exists()
    files=[h.O/'jvp_evidence.json',h.O/'branch_steps_evidence.json',r.O/'contract_evidence.json',Path(__file__).resolve()]
    cases=[];original_failures=0
    for cc in a.read(h.O/'inputs.json')['cases']:
        old=r.entries(h.L/'jvp'/cc['key']/'rows.csv');bad={x['direction'] for x in old if not x['qualified']}
        original_failures+=sum(not x['qualified'] for x in old)
        micro=[x for x in old if float(x['amplitude_V'])<1e-10];assert micro and all(x['qualified'] for x in micro)
        new=[]
        if bad:
            out=r.O/cc['key'];a.verify(out/'evidence.json');s=a.read(out/'summary.json')
            assert s['qualified'] and s['repeated_directions']==len(bad)
            new=r.entries(r.L/cc['key']/'rows.csv');assert all(x['qualified'] for x in new)
            for direction in bad:
                for step in ('1e-10','1e-12'):
                    count=sum(x['direction']==direction+'_refine_'+step for x in new)
                    assert count==len(micro)//6,(cc['key'],direction,step,count,len(micro))
            files.append(out/'evidence.json')
        cases.append(dict(key=cc['key'],initial_failed_directions=len(bad),micro_checks=len(micro),additional_checks=len(new),max_micro_relative=max(x['true_relative'] for x in micro),max_refined_relative=max((x['true_relative'] for x in new),default=0.),qualified=True))
    branches=a.read(h.O/'branch_steps_summary.json');assert all(not x['hfs_cross'] and not x['enormal_sign_cross'] for x in branches if x['step_V']<1e-6)
    summary=dict(cases=cases,qualified=True,original_failed_entries=original_failures,original_challenge_passed=False,
        row_relative_gate=1e-4,weak_entries_excluded=0,acceptance_tolerance_changed=False,source_changed_since_freeze=False,
        claim='Branch-local assembled Jv validated by all original microsteps and two smaller same-branch steps for every failed direction. The original 1e-6 finite-amplitude challenge remains failed; no claim of global smoothness or finite-candidate qualification.')
    a.write(h.O/'tangent_qualification.json',summary);d.matrix.freeze(h.O/'tangent_evidence.json',files+[h.O/'tangent_qualification.json']);print(summary,flush=True)

def dc():
    a.verify(h.O/'tangent_evidence.json');assert a.read(h.O/'tangent_qualification.json')['qualified'];a.verify(h.O/'candidate_freeze.json')
    h.v.LOCAL=h.L;h.v.OUT=h.O;h.v.RUNNER=h.RUNNER;h.V.post_config=h.c.qualification.previous.post_config
    jobs=[(cc,arm,cc[arm]) for cc in a.read(h.O/'inputs.json')['cases'] for arm in ('enormal_seed','native_seed')]
    def run(job):cc,arm,seed=job;return h.v.solve_target(cc,cc['index'],arm,Path(seed))
    with ThreadPoolExecutor(max_workers=2) as pool:attempts=[x for group in pool.map(run,jobs) for x in group]
    h.v.csv_union(h.O/'attempts.csv',attempts);d.matrix.freeze(h.O/'dc_evidence.json',[h.O/'candidate_freeze.json',h.O/'tangent_evidence.json',h.O/'attempts.csv',Path(__file__).resolve()]+[f for f in (h.L/'vela').rglob('*') if f.is_file()])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('gate','dc','summarize'));s=p.parse_args().action
    h.summarize() if s=='summarize' else globals()[s]()
