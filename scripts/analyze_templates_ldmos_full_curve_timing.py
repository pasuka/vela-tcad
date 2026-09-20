"""Read-only aggregation of full-curve timing, including partial batch status."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import statistics

from run_templates_ldmos_linked_d5 import read, write, digest
from run_templates_ldmos_analysis_reuse import aggregate


def repeat_statistics(values):
    if not values or any(v <= 0 for v in values): raise ValueError('Positive timings required')
    return dict(count=len(values),values=values,median=statistics.median(values),
        minimum=min(values),maximum=max(values),
        relative_range_percent=100*(max(values)-min(values))/statistics.mean(values))


def profile_breakdown(directory):
    ledger=read(directory/'fixed/ledger.json')
    groups={key:dict(requests=0,updates=0,counters=Counter(),seconds=Counter(),
                    main_only_profiles=0,main_only_analyses=0,mixed_profiles=0)
            for key in ('prefix_0_to_9p333333','remaining_high_voltage')}
    for run in ledger['runs']:
        group=groups['prefix_0_to_9p333333' if run['target_V']<=9.33333333333333+1e-12 else 'remaining_high_voltage']
        profile=read(directory/run['case']/'performance_profile.json')
        group['requests']+=1;group['updates']+=run['Newton_updates']
        group['counters'].update(profile['counters'])
        for stage in profile['stages']: group['seconds'][stage['name']]+=stage['total_ns']/1e9
        main=next((stage['calls'] for stage in profile['stages'] if stage['name']=='newton.linear_l2_row_column'),0)
        calls=profile['counters'].get('linear.solve_calls',0)
        if calls and main==calls:
            group['main_only_profiles']+=1
            group['main_only_analyses']+=profile['counters'].get('linear.analyze_calls',0)
        elif calls: group['mixed_profiles']+=1
    return groups


def summarize(root):
    report=read(root/'summary.json');cases=[];groups=defaultdict(list)
    for case in report['cases']:
        if case['status']!='pass': continue
        directory=root/case['name'];ledger=read(directory/'fixed/ledger.json')
        totals=aggregate(directory)
        if totals!=case['totals']: raise ValueError('Recorded aggregates changed: '+case['name'])
        samples=case['system_cpu_samples'];duration=sum(x['interval_seconds'] for x in samples)
        if duration<=0: raise ValueError('Missing load samples')
        cpu=case['cpu_seconds']['total']+case['audit']['child_cpu_seconds']
        row={k:case[k] for k in ('key','name','stage','profile','points','round','gate','config','wall_seconds')}
        row.update(process_tree_cpu_seconds=cpu,updates=ledger['total_Newton_updates'],
            requests=len(ledger['runs']),rollbacks=len(ledger['rollbacks']),totals=totals,
            segments=profile_breakdown(directory),
            mean_system_busy_percent=sum(x['busy_percent']*x['interval_seconds'] for x in samples)/duration,
            background_cpu_seconds_estimate=max(0.,sum(x['busy_cpu_seconds'] for x in samples)-cpu))
        cases.append(row);groups[(row['profile'],row['points'],row['config'],row['gate'])].append(row)
    repeats=[]
    for (profile,points,config,gate),rows in sorted(groups.items()):
        rows.sort(key=lambda x:x['round'])
        repeats.append(dict(profile=profile,points=points,config=config,gate=gate,
            rounds=[c['round'] for c in rows],
            wall_seconds=repeat_statistics([c['wall_seconds'] for c in rows]),
            cpu_seconds=repeat_statistics([c['process_tree_cpu_seconds'] for c in rows]),
            updates=[c['updates'] for c in rows],rollbacks=[c['rollbacks'] for c in rows]))
    paired=[]
    lookup={(c['profile'],c['points'],c['round'],c['gate'],c['config']):c for c in cases}
    for c in cases:
        for control in ('U0','U1'):
            if c['config']==control: continue
            other=lookup.get((c['profile'],c['points'],c['round'],c['gate'],control))
            if other:
                paired.append(dict(profile=c['profile'],points=c['points'],round=c['round'],
                    gate=c['gate'],control=control,candidate=c['config'],
                    wall_saving_percent=100*(1-c['wall_seconds']/other['wall_seconds']),
                    cpu_saving_percent=100*(1-c['process_tree_cpu_seconds']/other['process_tree_cpu_seconds'])))
    return dict(status=report['status'],completed_cases=len(cases),planned_cases=len(report['schedule']),
        completed_exact_points=sum(c['points'] for c in cases),cases=cases,repeats=repeats,pairs=paired,
        joint_qualifications=report['joint_qualifications'],comparisons=report['comparisons'],
        repeat_comparisons=report['repeat_comparisons'],
        unsuccessful_attempts=[c for c in report['cases'] if c['status'] not in ('pass','running')],
        source_summary_sha256=digest(root/'summary.json'),analyzer_sha256=digest(Path(__file__)),
        timing_scope='Controller wall including in-run scoring; excludes freeze and independent post-audit/joint score.',
        cpu_scope='Controller plus point-service CPU; worker shutdown tail excluded.',
        segment_scope='Requests grouped by target voltage; nested stage seconds are not additive wall time.',
        memory_scope='Solver-process lifetime working-set maximum, not simultaneous process-tree usage.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    args=parser.parse_args();result=summarize(args.root.resolve());write(args.root/'analysis.json',result)
    print(result['status'],result['completed_cases'],result['completed_exact_points'])
