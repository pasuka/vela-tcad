"""Audit full cache controls, including initialization and rejected trajectories.

Run against original VM paths. Exact numerical signatures avoid retaining every
large field simultaneously. Native fields are audited separately after copying.
"""
import argparse
import hashlib
import json
from pathlib import Path

from analyze_templates_ldmos_contact_sweep import numerical
from run_templates_ldmos_contact_sweep import read, sha
from run_templates_ldmos_electrothermal_curve import state_gate


def signature(value):
    return hashlib.sha256(json.dumps(numerical(value), sort_keys=True,
                                    separators=(',', ':')).encode()).hexdigest()


def inspect(directory):
    ledger = read(directory / 'results/ledger.json')
    assert ledger['status'] == 'complete' and len(ledger['exact_points']) == 31
    costs = dict(updates=0, initialization_updates=0, floor_updates=0,
                 trials=0, assemblies=0, assembly_seconds=0.,
                 root_solves=0, root_hits=0, factorizations=0,
                 factorization_seconds=0., linear_solve_seconds=0.)
    signatures = []
    for stage, rows in (('initialization', ledger['initialization_runs']),
                        ('drain', ledger['runs'])):
        for row in rows:
            path = Path(row['result']) if stage == 'initialization' else Path(row['directory'])/'output.json'
            data = read(path)
            if stage == 'drain':
                assert state_gate(data, row['bias_V'])['pass_gate'] == row['gate']['pass_gate']
                costs['floor_updates'] += sum(h['scaled_l2_before'] < 1e-9 for h in data['history'])
            costs['updates' if stage == 'drain' else 'initialization_updates'] += data['newton_updates']
            costs['trials'] += sum(h['line_search_trials'] for h in data['history'])
            perf = data['performance']
            for key, output in (('assembly_calls', 'assemblies'), ('assembly_seconds', 'assembly_seconds'),
                                ('factorizations', 'factorizations'), ('factorization_seconds', 'factorization_seconds'),
                                ('linear_solve_seconds', 'linear_solve_seconds')):
                costs[output] += perf[key]
            roots = perf.get('neutral_root_counts', [0, 0])
            costs['root_solves'] += roots[0]; costs['root_hits'] += roots[1]
            signatures.append((stage, row.get('bias_V'), row.get('parent_bias_V'),
                               signature(row.get('prediction')), signature(data)))
    for point in ledger['exact_points']:
        assert state_gate(read(point['result']), point['bias_V'])['pass_gate']
    costs['drain_attempts'] = len(ledger['runs'])
    costs['failed_attempts'] = sum(not r['gate']['pass_gate'] for r in ledger['runs'])
    return signatures, costs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix', type=Path, required=True)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); summary = read(a.matrix/'summary.json')
    assert summary['status'] == 'complete'
    assert sha(a.matrix/'vela_example_runner') == summary['runner_sha256']
    rows = summary['runs']
    assert len(rows) == 6 and {(r['gate_V'], r['variant']) for r in rows} == {
        (g, v) for g in (4, 8) for v in ('native', 'baseline', 'roots')}
    report = dict(status='pass', scope='One full timing pair per gate; not repeated timing stability.',
                  runner_sha256=summary['runner_sha256'], gates={})
    for gate in (4, 8):
        by = {r['variant']:r for r in rows if r['gate_V'] == gate}
        prior, _ = inspect(a.reference/('r0_vg%d_contact' % gate))
        result = {}
        for variant in ('baseline', 'roots'):
            row = by[variant]; directory = Path(row['directory'])
            assert sha(directory/'input.json') == row['input_sha256']
            assert sha(directory/'results/ledger.json') == row['ledger_sha256']
            trajectory, costs = inspect(directory)
            assert trajectory == prior, 'Changed numerical trajectory: '+str(directory)
            ratio = row['external_wall_seconds']/by['native']['external_wall_seconds']
            result[variant] = dict(costs=costs, exact_to_qualified_R9=True,
                wall_seconds=row['external_wall_seconds'], cpu_seconds=row['child_cpu_seconds'],
                wall_ratio_to_native=ratio, cpu_ratio_to_native=row['child_cpu_seconds']/by['native']['child_cpu_seconds'])
            if ratio > 1.5: report['status'] = 'performance_gate_failed'
        result['native'] = {k:by['native'][k] for k in ('external_wall_seconds','child_cpu_seconds')}
        report['gates'][str(gate)] = result
    a.output.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
