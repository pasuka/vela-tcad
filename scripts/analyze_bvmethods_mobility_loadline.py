#!/usr/bin/env python3
"""Audit diagnostic signal and accepted load-line points without relaxing physical gates."""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

CONTACTS = ('drain', 'source', 'gate', 'substrate')

def jvp_block_metrics(row: dict, block: str) -> dict:
    """Undo the legacy max(1, FD-norm) display floor before comparing small blocks."""
    analytic = float(row['analytic_' + block + '_norm'])
    finite_difference = float(row['finite_difference_' + block + '_norm'])
    error = float(row[block + '_relative_error']) * max(1.0, finite_difference)
    if not all(math.isfinite(x) and x >= 0 for x in (analytic, finite_difference, error)):
        raise ValueError('Derivative norms must be finite and nonnegative')
    signal = max(analytic, finite_difference)
    return {'analytic_change_norm': analytic, 'finite_difference_change_norm': finite_difference,
            'absolute_error_norm': error, 'relative_error': error / signal if signal > 0 else None,
            'finite_difference_signal_resolved': finite_difference > 0,
            'qualification': 'diagnostic only; unresolved FD is not a derivative pass'}

def endpoint_probe_eligible(config: dict) -> bool:
    """The existing endpoint decomposition explicitly evaluates Boltzmann SG flux."""
    statistics = config['solver'].get('carrier_statistics', 'boltzmann')
    if isinstance(statistics, dict):
        statistics = statistics.get('model', 'boltzmann')
    return statistics == 'boltzmann'

def audit_accepted_point(row: dict) -> dict:
    if row.get('converged') != '1':
        raise ValueError('Rejected sweep CSV currents can be placeholders; use a saved-state functional')
    inner, outer = float(row['inner_voltage_V']), float(row['outer_voltage_V'])
    drain = float(row['current_total_A_per_um'])
    load = outer - inner - 1e7 * drain
    currents = {c: float(row['contact_' + c + '_total_A_per_um']) for c in CONTACTS}
    kcl = math.fsum(currents.values())
    kcl_limit = max(1e-12, 1e-5 * abs(drain))
    finite = all(math.isfinite(x) for x in (inner, outer, load, drain, kcl, *currents.values()))
    initial_zero = inner == 0 and outer == 0
    reused_seed = row.get('solver_method') == 'boundary_seed_state'
    coupled = row.get('solver_method') == 'coupled_load_line_newton'
    try:
        proof = json.loads(row.get('validation_diagnostics') or '{}')
    except (ValueError, TypeError):
        proof = {}
    checked_coupled = coupled and isinstance(proof, dict) and proof.get('coupled_acceptance_checked') is True
    # Older coupled records leave carrier/closure diagnostics at defaults.
    measured = not initial_zero and not reused_seed and (not coupled or checked_coupled)
    electron = float(row['global_electron_continuity_closure_ratio']) if measured else None
    hole = float(row['global_hole_continuity_closure_ratio']) if measured else None
    carrier = float(row['carrier_row_max_ratio']) if measured else None
    checks = {'finite': finite, 'load_line': abs(load) <= 1e-6,
              'terminal_kcl_csv': abs(kcl) <= kcl_limit,
              'carrier_rows': 0 <= carrier <= .001 if carrier is not None else None,
              'global_continuity': (0 <= electron <= .01 and 0 <= hole <= .01) if measured else None}
    if coupled:
        checks['coupled_enforced_physical_gates'] = checked_coupled and all(proof.get(k) is True for k in (
            'physical_residuals_finite', 'actual_load_line_satisfied', 'carrier_row_enforced',
            'carrier_row_satisfied', 'global_continuity_enforced', 'global_continuity_satisfied', 'terminal_kcl_satisfied'))
    return {'inner_V': inner, 'outer_V': outer, 'current_A_per_um': drain,
            'load_line_residual_recomputed_V': load, 'carrier_row_max_ratio': carrier,
            'electron_closure_ratio': electron, 'hole_closure_ratio': hole,
            'kcl_csv_A_per_um': kcl, 'kcl_limit_A_per_um': kcl_limit,
            'checks': checks, 'pass': all(checks.values()),
            'seed_gates_require_parent_evidence': reused_seed or initial_zero,
            'coupled_carrier_and_continuity_metrics_unmeasured': coupled and not checked_coupled,
            'independent_saved_state_kcl_required': True}

def audit_saved_state(terminal: dict, outer: float) -> dict:
    """Recompute gates from actual four-contact functionals, not sweep placeholders."""
    currents = {c: float(terminal['totals_A_per_um'][c]) for c in CONTACTS}
    kcl = math.fsum(currents.values())
    inner = float(terminal['inner_bias_V'])
    load = outer - inner - 1e7 * currents['drain']
    components = {}
    for contact in CONTACTS:
        q = terminal['components'][contact]
        electron = float(q['contact_current_extractor_electron_A_per_um'])
        hole = float(q['contact_current_extractor_hole_A_per_um'])
        # The extractor reports raw carrier fluxes: conventional total = electron - hole.
        error = currents[contact] - (electron - hole)
        components[contact] = math.isfinite(error) and abs(error) <= 1e-12 * max(abs(electron) + abs(hole), 1e-30)
    checks = {'finite': all(math.isfinite(x) for x in (outer, inner, load, kcl, *currents.values())),
              'load_line': abs(load) <= 1e-6,
              'terminal_kcl': abs(kcl) <= max(1e-12, 1e-5 * abs(currents['drain'])),
              'component_identity': all(components.values())}
    return {'state_sha256': terminal['state_sha256'], 'outer_V': outer, 'inner_V': inner,
            'load_line_residual_recomputed_V': load, 'kcl_recomputed_A_per_um': kcl,
            'checks': checks, 'pass': all(checks.values())}


def summarize(audit_root: Path, followup_root: Path, restore_root: Path | None = None) -> dict:
    first = json.loads((audit_root / 'status.json').read_text())
    follow = json.loads((followup_root / 'status.json').read_text())
    endpoints = {}
    restore = json.loads((restore_root / 'status.json').read_text()) if restore_root else None
    statuses = [first, follow] + ([restore] if restore else [])
    for status in statuses:
        for step in status['steps']:
            if 'points' not in step or not step['points'] or 'outer_voltage_V' not in step['points'][0]:
                continue
            endpoints[step['case']] = {'process_exit': step['exit'],
                'points': [audit_accepted_point(row) for row in step['points'] if row.get('converged') == '1']}
    probes = {}
    paths = list((audit_root / 'fixed_direction').glob('*/jvp.csv'))
    paths += list(followup_root.glob('transport_frozen_mobility_retry/jvp.csv'))
    for path in paths:
        with path.open() as stream:
            probes[path.parent.name] = [{
                'direction': row['direction'], 'amplitude_V': float(row['amplitude_V']),
                'blocks': {b: jvp_block_metrics(row, b) for b in ('phin', 'phip')}}
                for row in csv.DictReader(stream)]
    edge_configs = list((audit_root / 'edge_transport').glob('*/config.json'))
    independent = {}
    for status in statuses[1:]:
        for name, terminal in status.get('terminal_audits', {}).items():
            case = name.removesuffix('_last_accepted')
            if case in endpoints and endpoints[case]['points']:
                independent[case] = audit_saved_state(terminal, endpoints[case]['points'][-1]['outer_V'])
    return {'schema': 'vela.bvmethods.mobility_loadline.audit.v1',
            'original_full_qs_qualified': False,
            'states': first['states'], 'runner_sha256': first['runner_sha256'],
            'newton_direction_audits': first['audit'], 'fixed_direction_probes': probes,
            'endpoint_probe_eligible': bool(edge_configs) and all(endpoint_probe_eligible(json.loads(p.read_text())) for p in edge_configs),
            'endpoint_probe_exclusion': 'Fermi source case: legacy endpoint ledger uses Boltzmann formulas; excluded from production derivative accuracy claims',
            'external_cases': endpoints, 'independent_terminal_audits': follow.get('terminal_audits', {}),
            'independent_load_line_and_kcl': independent,
            'restore_controls': restore,
            'initial_invalid_frozen_probe_retained': True,
            'completed': all(s['phase'] == 'completed' for s in statuses)}

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-root', type=Path, required=True)
    parser.add_argument('--followup-root', type=Path, required=True)
    parser.add_argument('--restore-root', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.audit_root, args.followup_root, args.restore_root)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'completed': result['completed'], 'endpoint_probe_eligible': result['endpoint_probe_eligible'],
                     'new_nonzero_points': {k: sum(p['outer_V'] > 0 and not p['seed_gates_require_parent_evidence'] for p in v['points']) for k, v in result['external_cases'].items()}}))
    return 0 if result['completed'] else 1

if __name__ == '__main__':
    raise SystemExit(main())
