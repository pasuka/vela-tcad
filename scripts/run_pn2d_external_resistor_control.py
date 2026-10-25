"""Bounded PN2D load-line control using the current atomic BV template.

Runs from each arm's own equilibrium, never from a mapped commercial state.
Raw inputs and results belong in an ignored, fresh output directory.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

from generate_pn2d_config import render_named_template


def validate_mobility_basis(config: dict) -> None:
    basis = config['solver']['mobility']['doping_concentration_basis']
    if basis not in ('total_impurity', 'cell_reconstructed_total_impurity'):
        raise ValueError('PN2D impurity scattering requires total Nd + Na, not net doping')


def qualify_rows(rows: list[dict[str, str]], arm: str, resistance: float) -> dict:
    """Recompute physical balances rather than trust a successful process exit."""
    checks = []
    for row in rows:
        current = float(row['current_total_A_per_um'])
        kcl = sum(float(row[key]) for key in row
                  if key.startswith('contact_') and key.endswith('_total_A_per_um'))
        check = {'converged': row['converged'] == '1',
                 'carrier_rows': int(row['carrier_row_violations']) == 0,
                 'global_closure': row['global_continuity_closure_satisfied'] == '1',
                 'kcl_A_per_um': kcl,
                 'kcl_passed': math.isfinite(kcl) and
                     abs(kcl) <= max(1e-12, 1e-5 * abs(current))}
        if arm != 'direct':
            load = (float(row['inner_voltage_V']) + resistance * current
                    - float(row['outer_voltage_V']))
            check.update(load_V=load, load_passed=math.isfinite(load) and abs(load) <= 1e-6)
        check['passed'] = all(check[key] for key in
            ['converged', 'carrier_rows', 'global_closure', 'kcl_passed']) and check.get('load_passed', True)
        checks.append(check)
    voltage_column = 'bias_V' if arm == 'direct' else 'outer_voltage_V'
    targets = [0.0, -0.002, -0.004, -0.006, -0.008, -0.01]
    targets_present = all(any(abs(float(row[voltage_column]) - target) <= 1e-12
                              for row in rows) for target in targets)
    return {'passed': targets_present and all(check['passed'] for check in checks),
            'targets_present': targets_present, 'points': checks}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runner', type=Path, required=True)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--doping', type=Path, required=True)
    parser.add_argument('--materials', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resistance', type=float, default=1e7)
    parser.add_argument('--timeout', type=float, default=180)
    parser.add_argument('--arms', nargs='+', choices=['direct', 'nested_scalar', 'coupled_newton'],
                        default=['direct', 'nested_scalar', 'coupled_newton'])
    args = parser.parse_args()
    if not args.resistance > 0:
        parser.error('resistance must be positive')
    runner = args.runner.resolve(strict=True)
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    inputs = root / 'inputs'
    inputs.mkdir()
    identities = {}
    for name, path in [('mesh.json', args.mesh), ('doping.csv', args.doping),
                       ('materials.json', args.materials)]:
        path = path.resolve(strict=True)
        shutil.copyfile(path, inputs / name)
        identities[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    base, manifest = render_named_template('pn2d_bv', {
        'mesh_file': '../inputs/mesh.json',
        'node_doping_file': '../inputs/doping.csv',
        'materials_file': '../inputs/materials.json',
    })
    validate_mobility_basis(base)
    report = {'schema': 'vela.pn2d_external_resistor_control.v1',
              'input_sha256': identities, 'template': manifest,
              'runner_sha256': hashlib.sha256(runner.read_bytes()).hexdigest(),
              'resistance_ohm_um': args.resistance, 'arms': {}}
    repo = Path(__file__).resolve().parents[1]
    report['source_head'] = subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    report['source_sha256'] = {name: hashlib.sha256((repo / name).read_bytes()).hexdigest()
        for name in ['src/simulation/DCSweep.cpp', 'include/vela/simulation/CoupledLoadLine.h',
                     'configs/templates/pn2d_bv.template.json',
                     'scripts/run_pn2d_external_resistor_control.py']}
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               VELA_LINEAR_SOLVER='umfpack')
    passed = True
    for arm in args.arms:
        folder = root / arm
        folder.mkdir()
        config = copy.deepcopy(base)
        config['output_csv'] = 'points.csv'
        sweep = config['sweep']
        sweep.update(mode='iv', start=0.0, stop=-0.01, step=-0.002,
                     bias_points=[0.0, -0.002, -0.004, -0.006, -0.008, -0.01],
                     initial_step=0.002, max_step=0.002,
                     csv_all_contacts=True, write_state_file='state.h5',
                     write_state_every_point_prefix='accepted/state',
                     stop_on_failure=True)
        sweep['diagnostics']['newton_history']['rejected_state_directory'] = 'rejected'
        if arm != 'direct':
            sweep['external_circuit'] = {
                'mode': 'series_resistor', 'solver': arm,
                'resistance_ohm_um': args.resistance, 'current_direction': 1.0,
                'residual_tolerance_V': 1e-6,
                'coupled_equation_tolerance': 1e-8,
                'coupled_linear_solver': 'direct_bordered',
            }
        (folder / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
        start = time.monotonic()
        with (folder / 'run.log').open('w') as log:
            try:
                code = subprocess.run([str(runner), '--config', 'config.json'], cwd=folder,
                                      env=env, stdout=log, stderr=subprocess.STDOUT,
                                      timeout=args.timeout).returncode
            except subprocess.TimeoutExpired:
                code = 124
        rows = []
        if (folder / 'points.csv').exists():
            with (folder / 'points.csv').open(newline='') as stream:
                rows = list(csv.DictReader(stream))
        qualification = qualify_rows(rows, arm, args.resistance)
        report['arms'][arm] = {'exit_code': code, 'qualification': qualification,
                               'seconds': time.monotonic() - start,
                               'rows': len(rows), 'last_row': rows[-1] if rows else None}
        passed = passed and code == 0 and qualification['passed']
        (root / 'execution.json').write_text(json.dumps(report, indent=2) + '\n')
        print(arm, code, len(rows), flush=True)
    # These bounded physical checks do not qualify a full BV trajectory or
    # establish parity with a native external-resistor reference.
    return 0 if passed else 2


if __name__ == '__main__':
    raise SystemExit(main())
