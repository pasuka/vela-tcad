"""Summarize the frozen same-matrix control without interpreting trials as DC passes."""
from pathlib import Path
import check_simplemos_minority_linear_20260906 as l

a=l.a


def run():
    a.verify(l.OUT/'freeze.json');a.verify(l.v.OUT/'trace_result.json')
    previous=a.rows(l.v.OUT/'last_updates.csv');precision=a.rows(l.v.OUT/'precision_summary.csv')
    summaries=[];identities=[];files=[Path(__file__).resolve(),l.OUT/'freeze.json',l.v.OUT/'trace_result.json',
        l.v.OUT/'last_updates.csv',l.v.OUT/'precision_summary.csv']
    for job in a.read(l.OUT/'contract.json')['jobs']:
        path=Path(job['config']).parent;rows=a.rows(path/'linear.csv')
        old={(r['carrier'],int(r['node_id'])):r for r in previous if r['tag']==job['tag']}
        initial=[r for r in rows if r['profile']=='legacy' and r['refinement']=='0']
        p=next(r for r in precision if r['tag']==job['tag']);node=int(p['worst_hole_node'])
        assert len(initial)==1814
        defect=max(abs(float(r['original_residual'])-float(old[r['carrier'],int(r['node'])]['residual']))/float(r['row_scale']) for r in initial)
        step=max(abs(float(r['step_V'])-float(old[r['carrier'],int(r['node'])]['raw_linear_step_V'])) for r in initial)
        weight=max(abs(float(r['row_weight'])-float(old[r['carrier'],int(r['node'])]['row_weight'])) for r in initial)
        assert defect==step==weight==0,(job,defect,step,weight)
        identities.append(dict(tag=job['tag'],rows=len(initial),max_residual_difference_scaled=defect,max_step_difference_V=step,max_row_weight_difference=weight))
        for profile in ('legacy','equilibrated'):
            for k in range(5):
                group=[r for r in rows if r['profile']==profile and int(r['refinement'])==k]
                assert len(group)==1814
                target=next(r for r in group if r['carrier']=='hole' and int(r['node'])==node)
                worst=max(group,key=lambda r:float(r['linear_ratio']));trial=max(group,key=lambda r:float(r['trial_ratio']))
                summary=dict(tag=job['tag'],profile=profile,refinement=k,rows=len(group),
                    baseline_violations=int(p['baseline_violations']),baseline_qualified=int(p['baseline_violations'])==0,
                    max_linear_ratio=float(worst['linear_ratio']),worst_linear_carrier=worst['carrier'],worst_linear_node=worst['node'],
                    linear_rows_above_1e6=sum(float(r['linear_ratio'])>1e-6 for r in group),
                    max_trial_ratio=float(trial['trial_ratio']),worst_trial_carrier=trial['carrier'],worst_trial_node=trial['node'],
                    trial_rows_above_1e6=sum(float(r['trial_ratio'])>1e-6 for r in group),
                    worst_hole_node=node,worst_hole_step_V=float(target['step_V']),worst_hole_linear_ratio=float(target['linear_ratio']),
                    worst_hole_trial_ratio=float(target['trial_ratio']),worst_hole_row_scale=float(target['row_scale']))
                summaries.append(summary)
                if k in (0,1,4):print(job['tag'],profile,k,summary['max_linear_ratio'],summary['max_trial_ratio'],summary['trial_rows_above_1e6'],flush=True)
        files += [Path(job['config']),path/'linear.csv',path/'config.status.json',path/'config.stdout.txt',path/'config.stderr.txt']
    a.write_csv(l.OUT/'summary.csv',summaries);a.write_csv(l.OUT/'legacy_identity.csv',identities)
    a.write(l.OUT/'result.json',dict(input_hashes={a.rel(p):a.sha(p) for p in sorted(set(files))},
        states=len(identities),all_legacy_rows_and_raw_steps_identical=True,groups=len(summaries),
        rows=sum(r['rows'] for r in summaries),nonlinear_solves=0,gate_changed=False,
        interpretation='Linear defect uses 100-digit dot accumulation on double J,F,dx. Trial uses uncapped x+dx and fixed ORIGINAL row scales; trial counts are not DC qualification, positivity/KCL/Poisson or initialization invariance certification. Column/row equilibration alone is not uniformly better.',
        observations='One refinement recovers n19 node 41 from an exactly zero raw step and removes its fixed-scale trial defect. High-Vd node 973 remains a finite-amplitude nonlinear counterexample after accurate linear solution.'))


if __name__=='__main__':run()
