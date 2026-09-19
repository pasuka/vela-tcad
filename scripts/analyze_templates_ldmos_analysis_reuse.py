"""Summarize paired cache controls without confusing auxiliary and Newton solves."""
import argparse
from pathlib import Path
from run_templates_ldmos_linked_d5 import read, write
from run_templates_ldmos_analysis_reuse import aggregate, compare


def summarize(root):
    report = read(root/'summary.json')
    if report['status'] != 'pass':
        raise ValueError('Require a completed study')
    cases = []
    for case in report['cases']:
        if case['status'] != 'pass':
            raise ValueError('Failed case in completed study')
        directory = root/case['name']; ledger = read(directory/'fixed/ledger.json')
        main_only, mixed, inactive = [], [], []
        for run in ledger['runs']:
            profile = read(directory/run['case']/'performance_profile.json')
            main = next((x['calls'] for x in profile['stages']
                         if x['name']=='newton.linear_l2_row_column'), 0)
            total = profile['counters'].get('linear.solve_calls', 0)
            item = dict(bias_V=run['target_V'], main_calls=main, all_calls=total,
                        analyses=profile['counters'].get('linear.analyze_calls', 0))
            if total == 0: inactive.append(item)
            elif main == total: main_only.append(item)
            else: mixed.append(item)
        totals = aggregate(directory)
        if totals != case['totals']:
            raise ValueError('Aggregates changed: '+case['name'])
        samples = case['system_cpu_samples']
        duration = sum(x['interval_seconds'] for x in samples)
        cpu = case['cpu_seconds']['total'] + case['audit']['child_cpu_seconds']
        cases.append(dict(name=case['name'], mode=case['mode'], gate=case['gate'], round=case['round'],
            wall_seconds=case['wall_seconds'], process_tree_cpu_seconds=cpu,
            updates=ledger['total_Newton_updates'], requests=len(ledger['runs']),
            rollbacks=len(ledger['rollbacks']), totals=totals,
            main_only_profiles=main_only, mixed_profiles=mixed, inactive_profiles=inactive,
            mean_system_busy_percent=sum(x['busy_percent']*x['interval_seconds'] for x in samples)/duration,
            background_cpu_seconds_estimate=max(0.,sum(x['busy_cpu_seconds'] for x in samples)-cpu)))
    pairs = []
    for candidate in cases:
        if candidate['mode'] != 'reuse': continue
        group = {c['mode']:c for c in cases if (c['gate'],c['round'])==(candidate['gate'],candidate['round'])}
        for mode in ('subprocess','worker'):
            if mode not in group: continue
            control = group[mode]
            state = compare(root/control['name'], root/candidate['name'])
            if not state['equivalent']: raise ValueError('State mismatch')
            pairs.append(dict(gate=candidate['gate'], round=candidate['round'], control=mode,
                wall_saving_percent=100*(1-candidate['wall_seconds']/control['wall_seconds']),
                cpu_saving_percent=100*(1-candidate['process_tree_cpu_seconds']/control['process_tree_cpu_seconds']),
                max_state_difference_V=state['max_state_difference_V'],
                same_target_sequence=state['same_target_sequence'], same_update_sequence=state['same_update_sequence']))
    return dict(status='pass', backend=report.get('backend','strumpack'), cases=cases, pairs=pairs,
        analysis_scope='Counters include main Newton, Poisson recorrection and recovery systems; mixed profiles cannot assign each analysis to an operator.',
        timing_scope='Curve controller wall; excludes fixture freezing and independent post-run study audit.',
        memory_scope='Peak of solver-process lifetime high-water marks; not simultaneous whole-process-tree memory.')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    args=parser.parse_args();write(args.root/'analysis.json',summarize(args.root))
