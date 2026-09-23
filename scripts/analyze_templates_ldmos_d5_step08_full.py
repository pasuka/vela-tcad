"""Read-only qualification and cost summary of a completed D5 full-curve pair."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
from analyze_templates_ldmos_d5_newton_cost import analyze_curve


def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',type=Path,required=True)
    p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();sys.path.insert(0,str(a.package/'scripts'))
    import run_templates_ldmos_linked_d5 as linked
    batch=read(a.batch/'summary.json')
    if batch['status']!='completed':raise ValueError('Batch is not complete')
    if not all(x['pass_all'] for x in batch['joint'].values()):raise ValueError('Joint gate failed')
    report=dict(batch_sha256=sha(a.batch/'summary.json'),cases=[],paired_states={},timing={})
    ledgers={}
    for case in batch['cases']:
        root=Path(case['directory']);ledger=read(root/'fixed/ledger.json');score=read(root/'score/summary.json')
        if ledger['status']!='completed' or len(ledger['exact_points'])!=31:raise ValueError('Incomplete curve')
        counters=Counter();profile=Counter();recoveries=[];guards=Counter();regions={k:Counter() for k in ('low','high')}
        for run in ledger['runs']:
            dest=root/run['case'];status=read(dest/'status.json')
            row=status['curve'][0];pred=run.get('outer_predictor')
            counts=dict(updates=run['Newton_updates'],services=1,
                failed_attempts=sum(x['status']!='accepted' for x in status['attempts']),
                internal_recovery_points=int(row.get('carrier_row_recovery_attempted',0)),
                internal_recovery_cycles=int(row.get('carrier_row_recovery_cycles',0)),
                density_passes=int(row.get('carrier_row_recovery_density_passes',0)))
            counters.update(counts);regions['low' if run['target_V']<=28/3+1e-10 else 'high'].update(counts)
            profile.update(read(dest/'performance_profile.json')['counters'])
            guards[run.get('predictor_guard_reason') or 'used']+=1
            if counts['internal_recovery_points']:recoveries.append(dict(case=run['case'],target_V=run['target_V'],**counts))
            if pred:
                for key in ('previous','current','predicted'):
                    if sha(Path(pred[key+'_state']))!=pred[key+'_sha256']:raise ValueError('Prediction lineage changed')
                if not 0<pred['ratio']<=2+1e-10:raise ValueError('Unprotected prediction')
        if counters['updates']!=ledger['total_Newton_updates']:raise ValueError('Update accounting mismatch')
        report['cases'].append(dict(**case,counters=dict(counters),regions=regions,
            profile_counters=dict(profile),recoveries=recoveries,guards=dict(guards),
            score=score,ledger_sha256=sha(root/'fixed/ledger.json'),
            trace_diagnostics=analyze_curve(root),
            max_actual_step=max(t['accepted_step_V'] for t in ledger['transfers'])))
        ledgers[(case['mode'],case['gate'])]=ledger
    for gate in (4,8):
        comparisons=[]
        for x,y in zip(ledgers[('candidate',gate)]['exact_points'],ledgers[('baseline',gate)]['exact_points']):
            if x['bias_V']!=y['bias_V']:raise ValueError('Exact target mismatch')
            diff=linked.state_difference(Path(x['state']),Path(y['state']))
            error=max(diff[k]['max_absolute'] for k in ('psi','phin','phip'))
            comparisons.append(dict(bias_V=x['bias_V'],max_potential_V=error,state_difference=diff))
            if error>1e-8:raise ValueError('Fresh paired state comparison failed')
        report['paired_states'][str(gate)]=comparisons
        get=lambda mode:next(c for c in report['cases'] if c['mode']==mode and c['gate']==gate)
        candidate,baseline=get('candidate'),get('baseline')
        report['timing'][str(gate)]=dict(candidate_seconds=candidate['wall_seconds'],
            baseline_seconds=baseline['wall_seconds'],ratio=candidate['wall_seconds']/baseline['wall_seconds'],
            candidate_updates=candidate['updates'],baseline_updates=baseline['updates'])
    report['status']='pass';report['joint']=batch['joint']
    a.output.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report['timing']),flush=True)


if __name__=='__main__':main()
