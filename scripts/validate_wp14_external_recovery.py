#!/usr/bin/env python3
"""Independently check native fixed-resistor recovery, beyond a Good Bye log.

This qualifies a transient acquisition and its subsequent DC branch only. It
does not qualify the original zero-bias quasistationary acquisition trajectory.
Inputs are a byte-verified C7 archive; raw native files remain outside Git.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
from reference_artifacts import log_exit_status, lf_sha256
from vela_sdevice.plt.document import read_plt_document

RESISTANCE = 1e7  # ohm*um, 2-D
STOP = 1.443e-3  # A/um
CONTACTS = ('drain', 'source', 'gate', 'substrate')


def check_state(row: dict, *, require_stop: bool = True) -> dict:
    drain = row['drain']
    current = drain['total_current_A_per_um']
    outer, inner = drain['outer_voltage_V'], drain['inner_voltage_V']
    displacement = drain['DisplacementCurrent']
    conduction = drain['e_current_A_per_um'] + drain['h_current_A_per_um']
    scale = max(abs(current), 1e-20)
    load_error = abs(outer - inner - RESISTANCE * current)
    kcl = abs(sum(row[name]['total_current_A_per_um'] for name in CONTACTS)) / scale
    checks = {
        'declared_stop': not require_stop or abs(current) >= STOP,
        'fixed_resistor_load_line': load_error <= 1e-6,
        'terminal_kcl': kcl <= 1e-5,
        'negligible_displacement': abs(displacement) / scale <= 1e-6,
        'current_components': abs(current - conduction - displacement) / scale <= 1e-8,
        'finite_positive_voltage': outer > 0 and inner > 0,
    }
    return {'status': 'pass' if all(checks.values()) else 'fail', 'checks': checks,
            'outer_voltage_V': outer, 'inner_voltage_V': inner,
            'total_current_A_per_um': current, 'conduction_current_A_per_um': conduction,
            'displacement_current_A_per_um': displacement,
            'load_line_residual_V': load_error, 'kcl_relative': kcl}


def validate(archive: Path) -> dict:
    results, rows_by_job = {}, {}
    jobs = ('external_transient_saved', 'external_dc_saved', 'external_recovery_saved')
    fingerprints = {}
    baseline = Path(__file__).resolve().parents[1] / 'reference_tcad/bvmethods_sentaurus2018/source/external_resistor_sdevice.cmd'
    physics_pattern = r'\bPhysics\s*\{([^}]*)\}'
    expected_physics = re.sub(r'\s+', '', re.search(physics_pattern, baseline.read_text()).group(1))
    parameter = archive / 'source/models.par'
    if lf_sha256(parameter.read_bytes()) != lf_sha256((baseline.parent / 'models.par').read_bytes()):
        raise ValueError('original repository material parameters are required')
    fingerprints['source/models.par'] = hashlib.sha256(parameter.read_bytes()).hexdigest()
    for job in jobs:
        cmd = archive / 'source' / (job + '_sdevice.cmd')
        plt = archive / 'plt' / (job + '_des.plt')
        log = archive / 'log' / (job + '_des.log')
        source = cmd.read_text()
        physics = re.sub(r'\s+', '', re.search(physics_pattern, source).group(1))
        if physics != expected_physics or 'Resistor=1e7' not in source:
            raise ValueError(job + ': original physical models and fixed resistance are required')
        if job != 'external_dc_saved' and 'Current(Contact="drain" AbsVal=1.443e-3)' not in source:
            raise ValueError(job + ': original declared stop is required')
        if job == 'external_recovery_saved' and 'Goal { Name="drain" Voltage=100000.0 }' not in source:
            raise ValueError(job + ': original final voltage goal is required')
        names, rows = read_plt_document(plt.read_text())
        if set(names) != set(CONTACTS) or not rows:
            raise ValueError(job + ': complete terminal data required')
        results[job] = check_state(rows[-1])
        results[job]['native_exit'] = log_exit_status(log.read_text())
        if results[job]['native_exit'] != 'Good Bye': results[job]['status'] = 'fail'
        rows_by_job[job] = rows
        for path in (cmd, plt, log):
            fingerprints[path.relative_to(archive).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    transient, dc = (results[job] for job in jobs[:2])
    steady_comparison = {
        'relative_current_difference': (abs(dc['total_current_A_per_um'] / transient['total_current_A_per_um'] - 1)
                                        if transient['total_current_A_per_um'] != 0 else None),
        'inner_voltage_difference_V': abs(dc['inner_voltage_V'] - transient['inner_voltage_V']),
        'outer_voltage_difference_V': abs(dc['outer_voltage_V'] - transient['outer_voltage_V']),
    }
    comparison_pass = (steady_comparison['relative_current_difference'] is not None
                       and steady_comparison['relative_current_difference'] <= 1e-4
                       and steady_comparison['inner_voltage_difference_V'] <= 1e-3
                       and steady_comparison['outer_voltage_difference_V'] <= 1e-6)
    recovery_rows = rows_by_job[jobs[2]]
    anchors = [row for row in recovery_rows if abs(row['drain']['outer_voltage_V'] - 14300) <= 1e-6]
    anchor_proof = bool(anchors) and all(abs(row['drain']['total_current_A_per_um']) < STOP
                                      and check_state(row, require_stop=False)['status'] == 'pass' for row in anchors)
    branch = [check_state(row, require_stop=False) for row in recovery_rows
              if row['drain']['outer_voltage_V'] >= 14300]
    branch_pass = bool(branch) and all(state['status'] == 'pass' for state in branch)
    for seed in ('external_transient_seed_des.sav', 'external_transient_seed_circuit_des.sav'):
        path = archive / 'tdr' / seed
        if not path.is_file(): raise ValueError('complete device and circuit Save seed required: ' + seed)
        fingerprints[path.relative_to(archive).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return {'schema': 'vela.wp14.external_recovery_validation.v1',
            'status': 'pass' if all(value['status'] == 'pass' for value in results.values())
                          and comparison_pass and anchor_proof and branch_pass else 'fail',
            'qualified_scope': 'same-model fixed-resistor transient acquisition, independent DC endpoint, and DC continuation from below the original stop; original zero-bias QS trajectory excluded',
            'resistance_ohm_um': RESISTANCE, 'stop_current_A_per_um': STOP,
            'states': results, 'steady_comparison': steady_comparison,
            'steady_comparison_pass': comparison_pass, 'below_stop_14300V_anchor_pass': anchor_proof,
            'dc_branch_states_checked': len(branch), 'dc_branch_pass': branch_pass,
            'native_fingerprints_sha256': fingerprints}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = validate(args.archive)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('fixed-resistor recovery: ' + result['status'])
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
