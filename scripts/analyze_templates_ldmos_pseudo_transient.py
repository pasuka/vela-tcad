"""Summarize frozen PTC point failures and independent first-state directions."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_matrix(root):
    summary = read(root / 'summary.json')
    if summary['status'] != 'completed_diagnostic_matrix':
        raise ValueError('Incomplete diagnostic matrix: ' + str(root))
    if sha(root / 'electrothermal_probe.exe') != summary['binary_sha256']:
        raise ValueError('Frozen executable changed')
    for name, digest in summary['runtime_sha256'].items():
        if sha(root / name) != digest:
            raise ValueError('Frozen DLL changed: ' + name)
    outputs = {}
    for path in root.glob('*/output.json'):
        digest = sha(path)
        if digest in outputs:
            raise ValueError('Ambiguous duplicate output digest')
        outputs[digest] = path
    rows = []
    for run in summary['runs']:
        path = outputs[run['output_sha256']]
        if sha(path.parent / 'input.json') != run['input_sha256']:
            raise ValueError('Frozen input changed')
        rows.append((run, read(path)))
    return summary, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', type=Path, required=True)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    matrix, runs = checked_matrix(args.matrix)
    audit, audits = checked_matrix(args.audit)
    if matrix.get('direction_audit_only', False) or not audit['direction_audit_only']:
        raise ValueError('Matrix/audit scopes confused')
    if len(runs) != 12 or len(audits) != 2:
        raise ValueError('Expected twelve point controls and two first-state audits')
    result = dict(schema='vela.pseudo_transient_findings.v1',
                  matrix_sha256=sha(args.matrix / 'summary.json'),
                  audit_sha256=sha(args.audit / 'summary.json'), points=[], directions=[])
    for run, data in runs:
        history = data['history']
        result['points'].append(dict(
            variant=run['variant'], vg=run['vg'], qualified=run['qualified'],
            updates=data['newton_updates'], attempts=len(history), stop=data['diagnostic_stop'],
            trials=run['trials'], fallbacks=run['fallbacks'], wall_seconds=run['wall_seconds'],
            assemblies=run['performance']['assembly_calls'], factorizations=run['performance']['factorizations'],
            blocks=data['electrical_block_gates'], row_max_ratio=data['carrier_row_gate']['max_ratio'],
            mass_seconds=data.get('pseudo_transient', {}).get('mass_seconds', 0.),
            initial_tau_s=data.get('pseudo_transient', {}).get('initial_tau_s'),
            # QF-only controls in the first frozen executable did not measure
            # density targets; their placeholder zero must not imply positivity.
            first_nonpositive_targets=(history[0].get('pseudo_transient', {}).get('nonpositive_density_targets')
                                       if run['variant'].startswith('ptc_v1_') else None)))
    for run, data in audits:
        directions = data['pseudo_direction_audit']['directions']
        if len(directions) != 8:
            raise ValueError('Incomplete coefficient isolation')
        for direction in directions:
            if 'error' in direction or not direction['steady_residual_exact']:
                raise ValueError('Failed or residual-changing direction audit')
            slope = direction['true_steady_merit_directional_slope']
            if not math.isfinite(slope):
                raise ValueError('Nonfinite merit slope')
            result['directions'].append(dict(vg=run['vg'], **direction))
    result['max_relative_linear_residual'] = max(
        d['relative_linear_residual'] for d in result['directions'])
    result['scope'] = 'Diagnostic counts only; failed solves are not speedups or physical qualification.'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(points=len(runs), directions=len(result['directions']),
                          max_relative_linear_residual=result['max_relative_linear_residual'])))


if __name__ == '__main__':
    main()
