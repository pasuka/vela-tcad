"""Classify every recorded Newton update; missing diagnostics remain explicit."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re


def analyze(result):
    counts=Counter(); records=[];floor_blockers=Counter();damped_limiters=Counter()
    linear_samples=[]
    for h in result['history']:
        before=h['scaled_l2_before'];after=h['scaled_l2_after']
        trace=h.get('iteration_trace',{})
        initial=h.get('initial_alpha',trace.get('initial_alpha'))
        ratio=after/before if before else None
        if not h['accepted']:phase='rejected'
        elif before<1e-9:phase='floor'
        elif initial is None:phase='missing_initial_alpha'
        elif initial<1e-2:phase='stall'
        elif h['alpha']<1:phase='damped'
        elif ratio is not None and ratio<=.1:phase='quadratic'
        else:phase='linear'
        counts[phase]+=1
        if phase=='floor':
            gates=trace.get('before_gates')
            if gates is None:floor_blockers['missing_trace']+=1
            else:
                for name,block in zip(('psi','electron','hole'),gates['blocks']):
                    if not block['satisfied']:floor_blockers[name]+=1
                if not gates['row']['satisfied']:floor_blockers['row']+=1
        if phase=='damped':
            component=h.get('limiter_component',trace.get('raw_limiter_component'))
            key='missing' if component is None else ('none' if initial>=1 or component<0 else ('psi','fn','fp','T')[component])
            damped_limiters[key]+=1
        if phase=='linear':
            for sample in trace.get('carrier_samples',[]):
                if 'relative_change_by_component' not in sample:continue
                parts=sample['relative_change_by_component'];total=sum(abs(v) for v in parts)
                linear_samples.append(dict(sample,iteration=h['iteration'],
                    qf_fraction_of_absolute_linear_change=abs(parts[sample['component']])/total if total else None))
        compact_trace={k:v for k,v in trace.items() if k not in ('before_gates','after_gates')}
        for stage in ('before_gates','after_gates'):
            if stage in trace:
                gates=trace[stage]
                compact_trace[stage]=dict(blocks=gates['blocks'],row={k:v for k,v in gates['row'].items() if k!='violations'},
                                         violation_count=len(gates['row']['violations']))
        records.append(dict(iteration=h['iteration'],phase=phase,merit_ratio=ratio,
            initial_alpha=initial,alpha=h['alpha'],trace=compact_trace or None,
            line_search_trials=h.get('line_search_trials',trace.get('line_search_trials'))))
    return dict(newton_updates=result['newton_updates'],stop=result['diagnostic_stop'],
                phases=dict(counts),floor_blockers=dict(floor_blockers),
                damped_limiters=dict(damped_limiters),linear_samples=linear_samples,records=records)


def summarize(runs):
    """Named frozen vg4/vg8 replays; counts overlap where multiple gates fail."""
    groups={}
    for gate in (4,8):
        selected=[r for r in runs if re.search(r'[/\\]vg%d_\d+[/\\]'%gate,r['path'])]
        if not selected:continue
        phases=Counter();blockers=Counter();limiters=Counter();dominant=positive=matched=sampled=0
        for run in selected:
            phases.update(run['phases']);blockers.update(run['floor_blockers']);limiters.update(run['damped_limiters'])
            by_iteration={}
            for q in run['linear_samples']:by_iteration.setdefault(q['iteration'],[]).append(q)
            for candidates in by_iteration.values():
                q=max(candidates,key=lambda v:abs(v['qf_direction_over_Vt']));sampled+=1
                if q['qf_fraction_of_absolute_linear_change']<.9:continue
                dominant+=1
                if not -1<q['relative_linear_change']<0:continue
                positive+=1
                if abs(q['actual_density_ratio']/math.exp(q['relative_linear_change'])-1)<.01:matched+=1
        groups[str(gate)]=dict(attempts=len(selected),phases=dict(phases),floor_blockers=dict(blockers),
            damped_raw_limiters=dict(limiters),linear_sampled=sampled,qf_fraction_ge_90_percent=dominant,
            positive_decreasing_linear_target=positive,actual_density_ratio_matches_exp_within_1_percent=matched)
    return groups


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--from-analysis',action='store_true',help='Summarize a saved analysis JSON instead of scanning raw point outputs')
    args=p.parse_args();totals=Counter();reports=[]
    if args.from_analysis:
        reports=json.loads(args.root.read_text(encoding='utf-8'))['runs']
        args.output.write_text(json.dumps(dict(source_sha256=hashlib.sha256(args.root.read_bytes()).hexdigest(),
            scope='For each linear update choose the larger recorded QF/Vt maximum; component dominance threshold 90%; density-ratio check 1%',
            groups=summarize(reports)),indent=2),encoding='utf-8')
        return
    for path in sorted(args.root.rglob('output.json')):
        result=json.loads(path.read_text(encoding='utf-8'))
        if 'history' not in result:continue
        report=analyze(result);totals.update(report['phases'])
        report.update(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        reports.append(report)
    args.output.write_text(json.dumps(dict(schema='vela.newton_phases.v1',
        scope='Recorded updates including failed attempts; phase names are diagnostic bins, not convergence-order proofs',
        phases=dict(totals),groups=summarize(reports),runs=reports),indent=2),encoding='utf-8')


if __name__=='__main__':main()
