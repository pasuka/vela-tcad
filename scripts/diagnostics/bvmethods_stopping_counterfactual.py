"""Read archived iteration-entry rows under coupled OR/AND scalar stopping rules.

The first divergence is a counterfactual only. Per-row KCL and other independent
physical gates are generally unavailable; this tool never labels such a row a
verified accepted state or manufactures a new IV trajectory. Candidate diagnostics
are not substituted for the report from an already-applied update.
"""
import argparse
import copy
import csv
import hashlib
import json
import math
from pathlib import Path


def configure_mode(config, mode='enforce'):
    if mode not in ('enforce', 'enforce_joint'):
        raise ValueError('expected enforce (OR) or enforce_joint (AND)')
    result = copy.deepcopy(config)
    circuit = result['sweep']['external_circuit']
    if circuit.get('solver') != 'coupled_newton':
        raise ValueError('coupled stopping requires coupled_newton')
    circuit.setdefault('coupled_update_convergence', {})['mode'] = mode
    return result


def finite(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def field(row, *names):
    for name in names:
        if name in row:
            return row[name]
    return None


def evaluate(row, equation_tolerance=1e-8, load_tolerance=1e-6):
    if not all(math.isfinite(x) and x > 0 for x in [equation_tolerance, load_tolerance]):
        raise ValueError('stopping tolerances must be positive and finite')
    f = finite(field(row, 'device_residual_inf', 'device_residual'))
    load = finite(field(row, 'load_residual_V', 'load_V'))
    rms = finite(field(row, 'update_rms', 'rms'))
    valid = str(field(row, 'update_valid')).lower() in ('1', 'true')
    residual = f is not None and load is not None and f <= equation_tolerance and abs(load) <= load_tolerance
    update = valid and rms is not None and rms < 1.
    load_gate = load is not None and abs(load) <= load_tolerance
    numerical_finite = f is not None and load is not None
    candidate = finite(field(row, 'candidate_update_rms', 'candidate_rms'))
    return {'device_residual': f, 'load_V': load, 'rms': rms, 'update_valid': valid,
            'residual_branch': residual, 'update_branch': update,
            'or_scalar_eligible': numerical_finite and load_gate and (residual or update),
            'and_scalar_eligible': numerical_finite and load_gate and residual and update,
            'candidate_rms_lt_one': str(field(row, 'candidate_update_valid', 'candidate_valid')).lower() in ('1', 'true') and candidate is not None and candidate < 1.,
            'outgoing_candidate_rejected': str(field(row, 'recovery_rejected', 'rejected')).lower() in ('1', 'true')}


def analyze(path, equation_tolerance=1e-8, load_tolerance=1e-6):
    path = Path(path)
    with path.open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.DictReader(stream))
    attempts = []
    previous_key = None
    previous_iteration = -1
    for index, row in enumerate(rows):
        key = tuple(field(row, *names) for names in [('outer_V',), ('attempted_step_V', 'step_V'), ('retry',)])
        iteration = int(row['iteration'])
        if key != previous_key or iteration < previous_iteration:
            attempts.append({'outer_V': finite(key[0]), 'step_V': finite(key[1]),
                             'retry': finite(key[2]), 'first_OR': None, 'first_AND': None,
                             'candidate_only_rows': 0, 'rows': 0})
        attempt = attempts[-1]
        e = evaluate(row, equation_tolerance, load_tolerance)
        record = {'csv_line': index+2, 'iteration': iteration,
                  'mode': field(row, 'derivative_mode', 'mode'), **e}
        attempt['rows'] += 1
        if e['candidate_rms_lt_one'] and not e['or_scalar_eligible']:
            attempt['candidate_only_rows'] += 1
        # Zero-bias initialization has a separate anchor contract, not a valid update sample.
        zero_anchor = finite(key[0]) == 0. and iteration == 0 and e['residual_branch']
        if e['or_scalar_eligible'] and attempt['first_OR'] is None:
            attempt['first_OR'] = {**record, 'zero_anchor': zero_anchor}
        if (e['and_scalar_eligible'] or zero_anchor) and attempt['first_AND'] is None:
            attempt['first_AND'] = {**record, 'zero_anchor': zero_anchor}
        previous_key, previous_iteration = key, iteration
    earlier = [a for a in attempts if a['first_OR'] and not a['first_OR']['zero_anchor']
               and (not a['first_AND'] or a['first_OR']['csv_line'] < a['first_AND']['csv_line'])]
    return {'input_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'equation_tolerance': equation_tolerance, 'load_tolerance_V': load_tolerance,
            'attempts': len(attempts), 'rows': len(rows),
            'OR_scalar_eligible_attempts': sum(a['first_OR'] is not None for a in attempts),
            'AND_scalar_eligible_attempts': sum(a['first_AND'] is not None for a in attempts),
            'OR_stops_earlier_attempts': len(earlier),
            'OR_eligible_before_outgoing_rejection': sum(a['first_OR']['outgoing_candidate_rejected'] for a in earlier),
            'candidate_only_rows': sum(a['candidate_only_rows'] for a in attempts),
            'earlier_examples': earlier,
            'qualification': 'scalar counterfactual; independent physical gates need replay; later archived steps are not an OR trajectory'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    cfg = commands.add_parser('configure')
    cfg.add_argument('input', type=Path)
    cfg.add_argument('output', type=Path)
    cfg.add_argument('--mode', choices=['enforce', 'enforce_joint'], default='enforce')
    data = commands.add_parser('analyze')
    data.add_argument('input', type=Path)
    data.add_argument('output', type=Path)
    data.add_argument('--equation-tolerance', type=float, default=1e-8)
    data.add_argument('--load-tolerance', type=float, default=1e-6)
    args = parser.parse_args()
    if args.command == 'configure':
        result = configure_mode(json.loads(args.input.read_text(encoding='utf-8-sig')), args.mode)
    else:
        result = analyze(args.input, args.equation_tolerance, args.load_tolerance)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
