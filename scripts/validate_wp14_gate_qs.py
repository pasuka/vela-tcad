#!/usr/bin/env python3
"""Qualify the separately authorized 0 -> 0.2 -> 0 V auxiliary gate QS path.

Native acquisition only: this does not qualify the original fixed-gate path,
dynamic stability, or Vela/native field agreement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

from reference_artifacts import lf_sha256
from vela_sdevice.cmd_parser.lexer import strip_comments
from vela_sdevice.plt.document import read_plt_document
from validate_wp14_zero_qs import BASELINE, block, source_checks, trajectory_checks
from validate_wp14_external_recovery import CONTACTS, STOP

MESH_SHA256 = '7cd3968f5dd0e3cd63ef62ed0fcced465dd33bea465167f0be8139c16025902e'
FILES = ('external_resistor_des.plt', 'gate_seed_external_resistor_des.plt',
         'final_external_resistor_des.plt')
SCHEMA = 'vela.wp14.auxiliary_gate_qs_validation.v1'


def gate_source_checks(source: str, baseline: str) -> dict:
    checks = source_checks(source, baseline)
    del checks['voltage_qs_only']
    solve = block(strip_comments(source), 'Solve')
    goals = re.findall(r'Goal\s*\{\s*Name\s*=\s*"(\w+)"\s+Voltage\s*=\s*([\d.eE+-]+)\s*\}', solve)
    checks['authorized_goal_sequence'] = ([(name, float(value)) for name, value in goals]
        == [('drain', 13000.0), ('gate', 0.2), ('gate', 0.0), ('drain', 100000.0)]
        and len(goals) == len(re.findall(r'\bGoal\s*\{', solve)))
    checks['qs_only_no_loaded_seed'] = (len(re.findall(r'\bQuasistationary\s*\(', solve)) == 4
        and not re.search(r'\b(?:Transient|Continuation|Load|Set|Unset|Current)\b', solve))
    checks['declared_phase_outputs'] = bool(re.search(r'NewCurrentPrefix\s*=\s*"gate_seed_"', solve)
        and re.search(r'NewCurrentPrefix\s*=\s*"final_"', solve))
    return checks


def _same_state(left: dict, right: dict) -> bool:
    for name in CONTACTS:
        a, b = left[name], right[name]
        for key in ('outer_voltage_V', 'inner_voltage_V'):
            if abs(a[key] - b[key]) > 1e-6:
                return False
        scale = max(abs(a['total_current_A_per_um']), abs(b['total_current_A_per_um']))
        if abs(a['total_current_A_per_um'] - b['total_current_A_per_um']) > max(1e-12, 1e-8 * scale):
            return False
    return True


def gate_trajectory_checks(phases: list[list[dict]], log: str) -> dict:
    if len(phases) != 3 or any(not phase for phase in phases):
        raise ValueError('three nonempty acquisition/gate/final trajectories required')
    acquisition, gate, final = phases
    rows = [row for phase in phases for row in phase]
    result = trajectory_checks(rows, log)
    checks = result['checks']
    del checks['original_fixed_bias_contacts_all_states']
    gate_voltage = [row['gate']['outer_voltage_V'] for row in gate]
    peak = max(range(len(gate_voltage)), key=gate_voltage.__getitem__)
    near = lambda value, target: abs(value - target) <= 1e-6
    checks.update({
        'source_and_substrate_fixed_zero_all_states': all(
            abs(row[name]['outer_voltage_V']) <= 1e-12 for row in rows for name in ('source', 'substrate')),
        'gate_zero_during_acquisition_and_final': all(
            abs(row['gate']['outer_voltage_V']) <= 1e-12 for row in acquisition + final),
        'acquisition_reached_13000V_below_stop': near(acquisition[-1]['drain']['outer_voltage_V'], 13000)
            and abs(acquisition[-1]['drain']['total_current_A_per_um']) < STOP,
        'gate_phase_drain_held_13000V': all(near(row['drain']['outer_voltage_V'], 13000) for row in gate),
        'gate_sweep_0_to_0_2_to_0': (abs(gate_voltage[0]) <= 1e-12
            and abs(gate_voltage[-1]) <= 1e-12 and abs(gate_voltage[peak] - 0.2) <= 1e-12
            and all(-1e-12 <= value <= 0.2 + 1e-12 for value in gate_voltage)
            and all(b >= a - 1e-12 for a, b in zip(gate_voltage[:peak], gate_voltage[1:peak+1]))
            and all(b <= a + 1e-12 for a, b in zip(gate_voltage[peak:], gate_voltage[peak+1:]))),
        'phase_states_continuous': _same_state(acquisition[-1], gate[0]) and _same_state(gate[-1], final[0]),
        'final_starts_at_13000V_below_stop': near(final[0]['drain']['outer_voltage_V'], 13000)
            and abs(final[0]['drain']['total_current_A_per_um']) < STOP,
    })
    result.update(status='pass' if all(checks.values()) else 'fail',
                  phase_state_counts=[len(phase) for phase in phases],
                  gate_peak_V=gate_voltage[peak],
                  before_gate=acquisition[-1]['drain'], after_gate=gate[-1]['drain'])
    return result


def validate(directory: Path) -> dict:
    cmd, par = directory / 'external_resistor_sdevice.cmd', directory / 'models.par'
    log, mesh = directory / 'external_resistor_des.log', directory / 'bvmethods_nmos_msh.tdr'
    checks = gate_source_checks(cmd.read_text(encoding='utf-8'),
                                (BASELINE / cmd.name).read_text(encoding='utf-8'))
    checks['same_parameters'] = lf_sha256(par.read_bytes()) == lf_sha256((BASELINE / par.name).read_bytes())
    checks['same_mesh'] = hashlib.sha256(mesh.read_bytes()).hexdigest() == MESH_SHA256
    phases = []
    for name in FILES:
        contacts, rows = read_plt_document((directory / name).read_text(encoding='utf-8'))
        if set(contacts) != set(CONTACTS):
            raise ValueError(name + ': original four contacts required')
        phases.append(rows)
    result = gate_trajectory_checks(phases, log.read_text(encoding='utf-8'))
    result.update(schema=SCHEMA, source_checks=checks,
                  qualified_scope='native zero-origin auxiliary gate 0->0.2->0V QS acquisition and final voltage QS; original fixed-gate path, stability, and Vela/field agreement excluded',
                  sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (cmd, par, mesh, log, *(directory / name for name in FILES))})
    if not all(checks.values()):
        result['status'] = 'fail'
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = validate(args.directory)
    except (OSError, ValueError, KeyError) as exc:
        result = {'schema': SCHEMA, 'status': 'fail', 'input_error': str(exc)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print('auxiliary gate pure QS: ' + result['status'])
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
