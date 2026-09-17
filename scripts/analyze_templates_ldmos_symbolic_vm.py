"""Audit completed symbolic-HFS controls, explicitly retaining partial matrices.

Native field identity and D0 physical-field gates are separate audits. This
report cannot promote a candidate even when numerical and timing checks pass.
"""
import argparse
import json
from pathlib import Path

from analyze_templates_ldmos_neutral_cache import signature
from evidence_paths import candidate_path
from run_templates_ldmos_contact_sweep import read, sha
from run_templates_ldmos_electrothermal_curve import state_gate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix', type=Path, required=True)
    p.add_argument('--path-map', nargs=2)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--allow-partial', action='store_true')
    args = p.parse_args()
    summary = read(args.matrix/'summary.json')
    if not args.allow_partial:
        assert summary['status'] == 'complete'
    assert sha(args.matrix/'vela_example_runner') == summary['runner_sha256']
    resolve = lambda path: candidate_path(path, args.path_map) if args.path_map else Path(path)
    report = dict(status='partial', matrix_status=summary['status'],
                  summary_sha256=sha(args.matrix/'summary.json'), runs=[], pairs=[], repeats=[],
                  physical_and_native_field_audits_required=True, promoted=False)
    report['prior_interrupted_runs'] = summary.get('prior_interrupted_runs', [])
    report['continuation_summary_sha256'] = summary.get('continuation_summary_sha256')
    rows = {}
    trajectories = {}
    for row in summary['runs']:
        key = (row['repeat'], row['gate_V'], row['variant'])
        assert key not in rows
        result = {k: row[k] for k in ('repeat', 'gate_V', 'variant', 'exit_code',
                                     'external_wall_seconds', 'child_cpu_seconds')}
        if row.get('interrupted_at_deadline') or row['exit_code'] != 0:
            result['complete'] = False
            report['runs'].append(result)
            continue
        directory = resolve(row['directory'])
        result['complete'] = True
        if row['variant'] != 'native':
            assert row['qualification']['pass_gate']
            for name in ('input', 'deck'):
                assert sha(directory/(name+'.json')) == row[name+'_sha256']
            assert sha(directory/'results/ledger.json') == row['qualification']['ledger_sha256']
            ledger = read(directory/'results/ledger.json')
            assert ledger['status'] == 'complete' and len(ledger['exact_points']) == 31
            signatures, raw_hashes = [], []
            totals = dict(drain_updates=0, initialization_updates=0, trials=0, failed_attempts=0)
            for initial in (True, False):
                for attempt in ledger['initialization_runs' if initial else 'runs']:
                    path = resolve(attempt['result']) if initial else resolve(attempt['directory'])/'output.json'
                    value = read(path)
                    assert value['performance']['linear_solver'] == 'umfpack'
                    if not initial:
                        passed = state_gate(value, attempt['bias_V'])['pass_gate']
                        assert passed == attempt['gate']['pass_gate']
                        totals['failed_attempts'] += not passed
                    totals['initialization_updates' if initial else 'drain_updates'] += value['newton_updates']
                    totals['trials'] += sum(h['line_search_trials'] for h in value['history'])
                    signatures.append((initial, attempt.get('bias_V'), attempt.get('parent_bias_V'),
                                       signature(attempt.get('prediction')), signature(value)))
                    raw_hashes.append(sha(path))
            assert raw_hashes == row['qualification']['output_sha256']
            assert all(v == row['qualification']['costs'][k] for k, v in totals.items())
            for point in ledger['exact_points']:
                assert state_gate(read(resolve(point['result'])), point['bias_V'])['pass_gate']
            trajectories[key] = signatures
            result.update(costs=totals, max_state_delta=row['qualification']['max_state_delta'])
        report['runs'].append(result)
        rows[key] = result
    expected = {(r, g, v) for r in range(summary['repeats']) for g in (4, 8)
                for v in ('native', 'baseline', 'high')}
    report['complete_runs'] = len(rows)
    report['expected_runs'] = len(expected)
    for r in range(summary['repeats']):
        for gate in (4, 8):
            group = {v: rows.get((r, gate, v)) for v in ('baseline', 'high', 'native')}
            if group['baseline'] and group['high']:
                baseline, high = group['baseline'], group['high']
                pair = dict(repeat=r, gate_V=gate,
                            wall_reduction_fraction=1-high['external_wall_seconds']/baseline['external_wall_seconds'],
                            cpu_reduction_fraction=1-high['child_cpu_seconds']/baseline['child_cpu_seconds'])
                if group['native']:
                    pair['high_wall_ratio_to_native'] = high['external_wall_seconds']/group['native']['external_wall_seconds']
                    pair['high_cpu_ratio_to_native'] = high['child_cpu_seconds']/group['native']['child_cpu_seconds']
                    pair['wall_gate'] = pair['high_wall_ratio_to_native'] <= 1.5
                report['pairs'].append(pair)
    for gate in (4, 8):
        for variant in ('baseline', 'high'):
            keys = [(r, gate, variant) for r in range(summary['repeats'])]
            if all(k in trajectories for k in keys):
                exact = all(trajectories[k] == trajectories[keys[0]] for k in keys[1:])
                report['repeats'].append(dict(gate_V=gate, variant=variant, exact_non_timing_trajectory=exact))
    if set(rows) == expected and summary['status'] == 'complete':
        report['status'] = 'complete_numerical_and_timing_audit'
    report['numerical_repeat_pass'] = len(report['repeats']) == 4 and all(r['exact_non_timing_trajectory'] for r in report['repeats'])
    report['wall_gate_pass'] = len(report['pairs']) == 2*summary['repeats'] and all(r.get('wall_gate', False) for r in report['pairs'])
    report['faster_every_pair'] = len(report['pairs']) == 2*summary['repeats'] and all(r['wall_reduction_fraction'] > 0 for r in report['pairs'])
    args.output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
