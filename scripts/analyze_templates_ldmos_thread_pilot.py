"""Verify exported exact states and compare each thread setting with its own control."""
import argparse
from pathlib import Path
import statistics

from run_templates_ldmos_linked_d5 import read, write, digest, state_difference
from analyze_templates_ldmos_host_pilot import verified_point


def analyze(root):
    report = read(root/'summary.json')
    if report['status'] not in ('pass', 'completed_with_rejections'):
        raise ValueError('Batch is not complete')
    cases = {c['key']: c for c in report['cases']}
    expected = {c['key'] for c in report['configurations']}
    if set(cases) != expected or len(cases) != len(report['cases']):
        raise ValueError('Missing or duplicate configurations')
    rows = []
    for config in report['configurations']:
        case = cases[config['key']]
        screens = [s for s in report['screen'] if s['key'] == config['key']]
        row = dict(**config, status=case['status'], reason=case.get('reason'),
                   matrix_rounds=screens)
        if case['status'] != 'pass':
            rows.append(row)
            continue
        if (len(screens) != 3 or {s['round'] for s in screens} != {0, 1, 2}
                or any(s['status'] != 'pass' for s in screens)):
            raise ValueError('Curve lacks three qualified matrix rounds')
        if not case['thread_audit']['verified']:
            raise ValueError('Thread audit missing')
        control = cases[config['backend']+'_t1']
        if control['status'] != 'pass':
            raise ValueError('Control is not qualified')
        directory = root/config['key']
        baseline = root/control['key']
        left, right = (read(d/'fixed/ledger.json') for d in (baseline, directory))
        if any(l['status'] != 'completed' or len(l['exact_points']) != 8 for l in (left, right)):
            raise ValueError('Incomplete exact points')
        differences = []
        for a, b in zip(left['exact_points'], right['exact_points']):
            if a['bias_V'] != b['bias_V']:
                raise ValueError('Voltage mismatch')
            differences.append(state_difference(verified_point(baseline, a), verified_point(directory, b)))
        maximum = max(d[k]['max_absolute'] for d in differences for k in ('psi', 'phin', 'phip'))
        if maximum > 1e-8:
            raise ValueError('State equivalence failed')
        wall = case['wall_seconds']
        samples = case['system_cpu_samples']
        duration = sum(s['interval_seconds'] for s in samples)
        row.update(wall_seconds=wall, speed_ratio=control['wall_seconds']/wall,
                   reduction_percent=100*(1-wall/control['wall_seconds']),
                   cpu_seconds=case['process_tree_cpu_seconds'],
                   average_cpu_cores=case['process_tree_cpu_seconds']/wall,
                   mean_system_busy_percent=sum(s['busy_percent']*s['interval_seconds'] for s in samples)/duration,
                   matrix_median_seconds=statistics.median(s['wall_seconds'] for s in screens),
                   updates=right['total_Newton_updates'],
                   same_targets=[r['target_V'] for r in left['runs']]==[r['target_V'] for r in right['runs']],
                   same_updates=[r['Newton_updates'] for r in left['runs']]==[r['Newton_updates'] for r in right['runs']],
                   max_state_difference_V=maximum, audit=case['audit'], totals=case['totals'])
        rows.append(row)
    return dict(status='verified', batch_status=report['status'],
                scope='One Vg8 first-eight-point D5 curve per configuration; matrix screening has three rounds',
                source_summary_sha256=digest(root/'summary.json'), runner_sha256=report['runner_sha256'],
                cases=rows)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args=parser.parse_args()
    result=analyze(args.root.resolve())
    write(args.output, result)
    for row in result['cases']:
        print(row['key'], row['status'], round(row.get('wall_seconds', 0), 3),
              'speed', round(row.get('speed_ratio', 0), 3), 'updates', row.get('updates'))
