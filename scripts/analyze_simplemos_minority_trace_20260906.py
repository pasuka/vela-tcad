"""Separate local linear accuracy, update constraints and scalar merit rejection."""
from pathlib import Path
import csv
import math
import numpy as np
import audit_simplemos_minority_residual_20260906 as v

a=v.a


def run():
    a.verify(v.OUT/'freeze.json');summaries=[];trials=[];updates=[];files=[Path(__file__).resolve(),v.OUT/'freeze.json',v.OUT/'traces.csv',v.OUT/'precision_summary.csv']
    worst={r['tag']:int(r['worst_hole_node']) for r in a.rows(v.OUT/'precision_summary.csv')}
    for job in a.read(v.OUT/'contract.json')['traces']:
        path=Path(job['config']).parent;s=a.read(path/'config.status.json');history=s['newton_trace'];raw=a.rows(path/'updates.csv')
        it=max(int(r['iteration']) for r in raw);current=[r for r in raw if int(r['iteration'])==it]
        trace=next(t for t in history if t['iteration']==it);prev=max((t for t in history if t['iteration']<it),key=lambda t:t['iteration'])
        scales=np.array([max(history[0]['blocks'][key],1.) for key in ('psi','phin','phip')]);base=np.array([prev['blocks'][key] for key in ('psi','phin','phip')])
        base_norm=float(np.linalg.norm(base/scales));node=worst[job['tag']]
        target=next(r for r in current if int(r['node_id'])==node and r['carrier']=='hole')
        bad=a.read(Path(job['source'])/'all_row.status.json')['carrier_row_convergence']['violations'];badset={(int(r['node_id']),r['carrier']) for r in bad}
        zero_bad=sum(float(r['raw_linear_step_V'])==0 and (int(r['node_id']),r['carrier']) in badset for r in current)
        summary=dict(tag=job['tag'],mode=job['mode'],iteration=it,converged=s['converged'],failure=s['failure_reason'],
            raw_linear_global_relative=float(current[0]['raw_linear_residual_l2'])/max(prev['blocks']['combined'],1e-300),
            capped_rows=sum(float(r['raw_linear_step_V'])!=float(r['capped_step_V']) for r in current),zero_update_previous_bad_rows=zero_bad,
            baseline_merit=base_norm,poisson_merit_fraction=float((base[0]/scales[0])**2)/max(base_norm**2,1e-300),
            worst_hole_node=node,worst_raw_step_V=float(target['raw_linear_step_V']),worst_capped_step_V=float(target['capped_step_V']),
            worst_row_weight=float(target['row_weight']),worst_linear_defect_over_rhs=abs(float(target['raw_linear_residual']))/max(abs(float(target['residual'])),1e-300),
            attempts=len(trace['trials']),accepted=trace['accepted'])
        summary.update(worst_baseline_residual=float(target['residual']),worst_linear_residual=float(target['raw_linear_residual']))
        groups={}
        with (path/'trials.csv').open(newline='') as f:
            for r in csv.DictReader(f):
                if int(r['iteration'])!=it:continue
                groups.setdefault(int(r['attempt']),[]).append(r)
        for trial in trace['trials']:
            nodes=groups[trial['attempt']];assert len(nodes)==s['nodes']
            matrix=np.array([[float(r[k]) for k in ('rpsi','rn','rp')] for r in nodes]);blocks=np.linalg.norm(matrix,axis=0);norm=float(np.linalg.norm(blocks/scales))
            assert math.isclose(norm,trial['residual'],rel_tol=1e-12,abs_tol=1e-290),(job,norm,trial)
            assert math.isclose(base_norm*(1-1e-4*trial['damping']),trial['target'],rel_tol=1e-12)
            value=next(float(r['rp']) for r in nodes if int(r['node'])==node)
            trials.append(dict(tag=job['tag'],iteration=it,**trial,poisson_l2=float(blocks[0]),electron_l2=float(blocks[1]),hole_l2=float(blocks[2]),
                electron_block_ratio=float(blocks[1])/max(base[1],1e-300),hole_block_ratio=float(blocks[2])/max(base[2],1e-300),
                worst_hole_node=node,worst_hole_residual=value,worst_hole_residual_ratio=abs(value)/max(abs(float(target['residual'])),1e-300)))
        summaries.append(summary);updates.extend(dict(tag=job['tag'],**r) for r in current)
        files += [path/n for n in ('config.status.json','updates.csv','trials.csv','state.csv')]
        print(job['tag'],summary['capped_rows'],summary['worst_raw_step_V'],summary['worst_linear_defect_over_rhs'],flush=True)
    a.write_csv(v.OUT/'update_summary.csv',summaries);a.write_csv(v.OUT/'last_updates.csv',updates);a.write_csv(v.OUT/'last_trials.csv',trials)
    a.write(v.OUT/'trace_result.json',dict(input_hashes={a.rel(p):a.sha(p) for p in sorted(set(files))},identity_replays=6,single_step_diagnostics=1,
        trials=len(trials),all_trial_norms_reconstructed=True,boundary='Last failed trial is explained only after replay identity; the saved 200-iteration failure has a new single-step diagnostic, not a historical-step reconstruction.'))


if __name__=='__main__':run()
