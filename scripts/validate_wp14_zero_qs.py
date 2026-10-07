#!/usr/bin/env python3
"""Check a native zero-bias, fixed-resistor voltage QS acquisition independently.

This checker distinguishes acquisition from a restart or a transient seed. It
does not establish Vela/native field agreement or stability of a branch.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

from reference_artifacts import lf_sha256, log_exit_status
from vela_sdevice.cmd_parser.lexer import strip_comments
from vela_sdevice.plt.document import read_plt_document
from validate_wp14_external_recovery import CONTACTS, RESISTANCE, STOP

BASELINE = Path(__file__).resolve().parents[1] / 'reference_tcad/bvmethods_sentaurus2018/source'


def block(text: str, name: str) -> str:
    """Extract a balanced native block, including its nested declarations."""
    match = re.search(r'\b' + name + r'\s*\{', text)
    if match is None:
        raise ValueError('missing ' + name + ' block')
    start, depth = match.end(), 1
    for index in range(start, len(text)):
        depth += (text[index] == '{') - (text[index] == '}')
        if depth == 0:
            return text[start:index]
    raise ValueError('unclosed ' + name + ' block')


def source_checks(source: str, baseline: str) -> dict:
    source, baseline = strip_comments(source), strip_comments(baseline)
    compact = lambda value: re.sub(r'\s+', '', value)
    solve = block(source, 'Solve')
    goals = re.findall(r'Goal\s*\{\s*Name\s*=\s*"drain"\s+Voltage\s*=\s*([\d.eE+-]+)\s*\}', solve)
    return {
        'same_physics': compact(block(source, 'Physics')) == compact(block(baseline, 'Physics')),
        'same_electrodes_and_fixed_resistor': compact(block(source, 'Electrode')) == compact(block(baseline, 'Electrode')),
        'original_stop': bool(re.search(r'Current\s*\(\s*Contact\s*=\s*"drain"\s+AbsVal\s*=\s*1\.443e-3\s*\)', source)),
        'voltage_qs_only': ('Quasistationary' in solve
                            and not re.search(r'\b(?:Transient|Continuation|Load|Set|Unset|Current)\b', solve)
                            and len(goals) == len(re.findall(r'\bGoal\s*\{', solve))),
        'original_final_voltage_goal': bool(goals) and float(goals[-1]) == 100000.0,
        'original_avalanche_current_approximation': not re.search(r'\bAvalDensGradQF\b', block(source, 'Math')),
    }


def trajectory_checks(rows: list[dict], log: str) -> dict:
    if not rows:
        raise ValueError('terminal trajectory is empty')
    drain = [row['drain'] for row in rows]
    finite = all(math.isfinite(value) for row in rows for name in CONTACTS for value in row[name].values())
    load = [abs(d['outer_voltage_V'] - d['inner_voltage_V'] - RESISTANCE * d['total_current_A_per_um']) for d in drain]
    # At zero current a relative KCL measure is undefined. The absolute floor
    # applies only to low-current states, preserving the recovery's 1e-5
    # relative budget at every high-current state.
    kcl_checks = []
    kcl_ratios = []
    displacement = []
    component = []
    for row, d in zip(rows, drain):
        current = abs(d['total_current_A_per_um'])
        residual = abs(sum(row[name]['total_current_A_per_um'] for name in CONTACTS))
        kcl_checks.append(residual <= max(1e-12, 1e-5 * current))
        if current >= 1e-7:
            kcl_ratios.append(residual / current)
        displacement.append(abs(d['DisplacementCurrent']) <= max(1e-20, 1e-6 * current))
        component.append(abs(d['total_current_A_per_um'] - d['e_current_A_per_um'] - d['h_current_A_per_um'] - d['DisplacementCurrent']) <= max(1e-20, 1e-8 * current))
    checks = {
        'finite': finite,
        'starts_at_zero_voltage': abs(drain[0]['outer_voltage_V']) <= 1e-12,
        'starts_below_original_stop': abs(drain[0]['total_current_A_per_um']) < STOP,
        'original_fixed_bias_contacts_all_states': all(
            abs(row[name]['outer_voltage_V']) <= 1e-12
            for row in rows for name in ('source', 'gate', 'substrate')),
        'monotone_outer_voltage': all(b['outer_voltage_V'] >= a['outer_voltage_V'] - 1e-8 for a, b in zip(drain, drain[1:])),
        'fixed_resistor_all_states': max(load) <= 1e-6,
        'terminal_kcl_all_states': all(kcl_checks),
        'negligible_displacement_all_states': all(displacement),
        'current_components_all_states': all(component),
        'original_stop_reached': abs(drain[-1]['total_current_A_per_um']) >= STOP,
        'native_current_stop_reported': 'Contact current exceeded limit:' in log,
        'native_success': log_exit_status(log) == 'Good Bye',
    }
    return {'checks': checks, 'status': 'pass' if all(checks.values()) else 'fail',
            'states_checked': len(rows), 'last': drain[-1],
            'max_load_line_residual_V': max(load),
            'max_kcl_relative_above_1e_minus7_A_per_um': max(kcl_ratios, default=None),
            'low_current_kcl_absolute_floor_A_per_um': 1e-12}


def validate(directory: Path) -> dict:
    cmd = directory / 'external_resistor_sdevice.cmd'
    par = directory / 'models.par'
    plt = directory / 'external_resistor_des.plt'
    log = directory / 'external_resistor_des.log'
    source = source_checks(cmd.read_text(encoding='utf-8'), (BASELINE / cmd.name).read_text(encoding='utf-8'))
    source['same_parameters'] = lf_sha256(par.read_bytes()) == lf_sha256((BASELINE / par.name).read_bytes())
    names, rows = read_plt_document(plt.read_text(encoding='utf-8'))
    if set(names) != set(CONTACTS):
        raise ValueError('all four original contacts are required')
    result = trajectory_checks(rows, log.read_text(encoding='utf-8'))
    result['source_checks'] = source
    result['status'] = 'pass' if result['status'] == 'pass' and all(source.values()) else 'fail'
    result['schema'] = 'vela.wp14.zero_qs_validation.v1'
    result['qualified_scope'] = 'native fixed-resistor voltage QS acquisition from zero; field/frontend agreement and branch stability excluded'
    result['sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (cmd, par, plt, log)}
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = validate(args.directory)
    except (OSError, ValueError, KeyError) as exc:
        result = {'schema': 'vela.wp14.zero_qs_validation.v1',
                  'status': 'fail', 'input_error': str(exc)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print('zero-bias pure QS: ' + result['status'])
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
