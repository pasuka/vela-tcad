#!/usr/bin/env python3
"""Stage TDX DF-ISE inputs and qualify the full BVMethods sd-sim run."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

from compare_reference_curves import interpolate_at
from sdevice_acceptance import GATE_FIXTURES, REPO, csv_rows, finite, read_json, sha
from reference_artifacts import decks_equivalent
from vela_sdevice.orchestrator.support import write_json
from vela_sdevice.plt import read_plt_document

FAMILY = REPO / 'reference_tcad/bvmethods_sentaurus2018/u2022_12_20261109'
DECK = 'pp5_des.cmd'
CURVE_COLUMNS = ('outer_voltage_V', 'inner_voltage_V', 'drain_total_current_A_per_um')
MESH_KEY = 'bvmethods_sentaurus2018'
SOURCE_GRID = 'n1_msh.tdr'


def reference_grid():
    return read_json(GATE_FIXTURES / 'manifest.json')['grids'][MESH_KEY]


def stage_inputs(output, tdx=None, tdr=None):
    grid = reference_grid()
    inputs = output / 'inputs'
    inputs.mkdir(parents=True, exist_ok=True)
    # Validate both files before publishing either one.
    with tempfile.TemporaryDirectory(prefix='tdx-', dir=output) as folder:
        temp = Path(folder)
        source = GATE_FIXTURES
        if tdx is not None:
            if sha(tdr) != grid['source_tdr']['sha256']:
                raise ValueError('BVMethods source TDR SHA-256 mismatch')
            shutil.copyfile(tdr, temp / SOURCE_GRID)
            with (output / 'tdx.log').open('w', encoding='utf-8') as log:
                subprocess.run([str(tdx), '-dd', '-M', '0', '-S', '0', SOURCE_GRID, 'n1_msh'], cwd=temp, check=True,
                               stdout=log, stderr=subprocess.STDOUT)
            source = temp
        for name, digest in grid['files'].items():
            if sha(source / name) != digest:
                raise ValueError(f'{name}: TDX DF-ISE SHA-256 mismatch')
        for name in grid['files']:
            shutil.copyfile(source / name, inputs / name)

    source = FAMILY / 'source'
    original = (source / DECK).read_text(encoding='utf-8')
    deck, count = re.subn(r'(\bGrid\s*=\s*")n1_msh\.tdr(")', r'\g<1>n1_msh.grd\2', original)
    if count != 1:
        raise ValueError('Expected exactly one original BVMethods grid declaration')
    (inputs / DECK).write_text(deck, encoding='utf-8')
    for name in ('pp5_des.par', 'pp5_des.vela.json'):
        shutil.copyfile(source / name, inputs / name)
    write_json(output / 'staged.json', {'mesh_files': grid['files'], 'tdr_sha256': grid['source_tdr']['sha256'],
                                      'export': 'tdx -dd' if tdx else 'pinned tdx -dd export'})


def check_inputs(path):
    grid = reference_grid()
    for name, digest in grid['files'].items():
        if sha(path / name) != digest:
            raise ValueError(f'{name}: staged mesh changed')
    for name in ('pp5_des.par', 'pp5_des.vela.json'):
        if sha(path / name) != sha(FAMILY / 'source' / name):
            raise ValueError(f'{name}: staged input changed')
    original = (FAMILY / 'source' / DECK).read_bytes()
    staged = (path / DECK).read_bytes()
    if not decks_equivalent(original, staged.replace(b'"n1_msh.grd"', b'"n1_msh.tdr"')):
        raise ValueError('Staged BVMethods CMD differs from original')


def curve(path):
    return [tuple(finite(row[k]) for k in CURVE_COLUMNS) for row in csv_rows(path)]


def qualify(run, published):
    status = read_json(run / 'run_manifest.json')
    stages = status['stages']
    plan = read_json(run / 'stage_plan.json')['stages']
    if (status['exit_code'] != 0 or len(stages) != len(plan) or len(stages) < 2
            or any(stage['runner_exit_code'] != 0 or stage['accepted_points'] < 1 for stage in stages)
            or not all(stage['goal_reached'] for stage in stages[:-1])
            or stages[-1].get('stop_criterion_reached') != 'current'):
        raise ValueError('Expected a complete sd-sim stage chain ending at the drain current stop')
    for stage in plan:
        solver = read_json(run / stage['deck'])['solver']
        if (solver.get('state_variables') != 'potential_density'
                or solver.get('primary_variable_coordinates') != 'unit'):
            raise ValueError(f"{stage['deck']}: expected Unit population coordinates")
    mesh = read_json(run / 'mesh.json')
    geometry = reference_grid()['equivalence']
    if (len(mesh['nodes']) != geometry['nodes'] or len(mesh['triangles']) != geometry['triangles']
            or {c['name']: len(c['edge_node_ids']) for c in mesh['contacts']} != geometry['contact_edges']):
        raise ValueError('Imported BVMethods mesh does not match the qualified topology')
    rows = csv_rows(run / plan[-1]['csv'])
    if len(rows) != stages[-1]['accepted_points']:
        raise ValueError('Sweep CSV and manifest accepted counts differ')
    if any(row['converged'] not in ('1', 'true') for row in rows):
        raise ValueError('Nonconverged row in accepted BV data')
    points = [(finite(r['outer_voltage_V']), finite(r['inner_voltage_V']),
               finite(r['current_total_A_per_um'])) for r in rows]
    if points[0][0] != 0 or any(a[0] >= b[0] for a, b in zip(points, points[1:])):
        raise ValueError('BVMethods curve must start at zero and increase in outer voltage')
    reference = read_json(FAMILY / 'comparison.json')['contract']
    criterion, = plan[-1]['break']['current']
    stop = finite(criterion['abs_val_A_per_um'])
    if criterion['contact'] != 'drain' or stop != 1.3e-3:
        raise ValueError('Expected the original 1.3e-3 A/um drain stop')
    if abs(points[-1][2]) < stop or any(abs(p[2]) >= stop for p in points[:-1]):
        raise ValueError('Drain current did not stop the sweep at the first qualified crossing')
    threshold = 1e-4
    bracket = next(((a, b) for a, b in zip(points, points[1:]) if abs(a[2]) < threshold <= abs(b[2])), None)
    if bracket is None:
        raise ValueError('No accepted points bracket the BV measurement current')
    a, b = bracket
    bv = a[1] + (b[1] - a[1]) * (threshold - abs(a[2])) / (abs(b[2]) - abs(a[2]))
    bv_error = abs(bv / reference['existing_BV_acceptance_anchor_V'] - 1)

    checks = [json.loads(row['validation_diagnostics']) for row in rows]
    load_failures = sum(not d['actual_load_line_satisfied'] for d in checks)
    kcl_failures = sum(not d['terminal_kcl_satisfied'] for d in checks)
    nonfinite = sum(not d['physical_residuals_finite'] for d in checks)
    stopping_failures = sum(not (finite(d['coupled_device_residual_inf']) <= 1e-8
                                and d['density_update_valid'] and 0 <= finite(d['density_update_rms']) < 1)
                            for d in checks)
    native = curve(FAMILY / 'reference_curves/native_u2022_12_adaptive_0_100.csv')
    if not native or native[0][0] != 0 or native[-1][0] != 100:
        raise ValueError('Expected the bounded 0-100 V native curve')
    candidate = [(p[0], p[2]) for p in points]
    errors = []
    for outer, _, current in native:
        value = interpolate_at(candidate, outer)
        if value is None:
            raise ValueError(f'BVMethods curve does not cover the native {outer:g} V reference')
        errors.append((abs(value - current) / max(abs(current), 1e-12), outer))
    worst_curve_error, worst_curve_voltage = max(errors)
    _, published_rows = read_plt_document(published.read_text(encoding='utf-8'), sign_factor=1.0)
    expected_count = sum(stage['accepted_points'] for stage in stages)
    published_curve = [tuple(row['drain'][column] for column in
                             ('outer_voltage_V', 'inner_voltage_V', 'total_current_A_per_um'))
                       for row in published_rows[-len(points):]]
    if len(published_rows) != expected_count or published_curve != points:
        raise ValueError('Published drain PLT differs from the accepted CSV rows')
    passed = (bv_error <= reference['BV_relative_tolerance'] and worst_curve_error <= 0.02
              and not any((load_failures, kcl_failures, nonfinite, stopping_failures)))
    return dict(passed=passed, accepted_points=len(rows),
                newton_iterations=sum(int(r['iterations']) for r in rows),
                reached_current_stop=True, final_outer_voltage_V=points[-1][0],
                final_inner_voltage_V=points[-1][1], final_drain_current_A_per_um=points[-1][2],
                BV_inner_V=bv, BV_relative_error=bv_error, native_curve_points=len(errors),
                native_curve_worst_relative_error_with_floor=worst_curve_error,
                native_curve_worst_outer_voltage_V=worst_curve_voltage, load_failures=load_failures,
                kcl_failures=kcl_failures, nonfinite_points=nonfinite, stopping_failures=stopping_failures,
                carrier_report_failures=sum(not d['carrier_row_satisfied'] for d in checks),
                continuity_report_failures=sum(not d['global_continuity_satisfied'] for d in checks))


NUMERICAL_CONTROLS = {
    'coefficient-box': {'solver': {'bvmethods_numerics': {'avalanche_weights': 'coefficient_box'}}},
    'integral-statistics': {'bvmethods_statistics': 'fermi_dirac'},
    'effective-sg': {'bvmethods_statistics': 'fermi_dirac',
                     'solver': {'carrier_statistics': {'sg_formulation': 'effective_potential'}}},
    'gss-laux': {'solver': {'bvmethods_numerics': {'avalanche_current': 'gss_laux'}}},
    'sharp-switch': {'solver': {'bvmethods_numerics': {'avalanche_switch': 'sharp'}}},
}


def run_case(output, launcher, numerical_control=None):
    library = os.environ.get('VELA_MKL_RT_LIBRARY') or os.environ.get('VELA_MKL_RT_DLL')
    if not library or not Path(library).is_absolute() or not Path(library).is_file():
        raise ValueError('Set VELA_MKL_RT_LIBRARY to the absolute path of the installed oneMKL runtime')
    inputs = output / 'inputs'
    check_inputs(inputs)
    with tempfile.TemporaryDirectory(prefix='full-bv-', dir=output) as folder:
        work = Path(folder)
        for source in inputs.iterdir():
            shutil.copyfile(source, work / source.name)
        if numerical_control:
            sidecar = work / 'pp5_des.vela.json'
            data = read_json(sidecar)
            overrides = NUMERICAL_CONTROLS[numerical_control]
            data.setdefault('solver', {}).update(overrides.get('solver', {}))
            if 'bvmethods_statistics' in overrides:
                data['bvmethods_statistics'] = overrides['bvmethods_statistics']
            write_json(sidecar, data)
        run = work / 'pp5_des_vela'
        env = dict(os.environ, VELA_PYTHON=sys.executable, VELA_BIN=str(launcher.parent),
                   VELA_LINEAR_SOLVER='pardiso_mkl', VELA_LINEAR_THREADS='1', MKL_NUM_THREADS='1',
                   MKL_THREADING_LAYER='SEQUENTIAL', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
        try:
            with (output / 'full-bv.log').open('w', encoding='utf-8') as log:
                subprocess.run([str(launcher), '--physical-constants-profile', 'sentaurus_u2022_12',
                                '--threads', '1', '--allow-lossy', '--exit-on-failure', DECK],
                               cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            report = qualify(run, work / 'n5_des.plt')
        except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
            write_json(output / 'full-bv.json', {'passed': False, 'error': str(error)})
            raise
        finally:
            for file in (run / 'fail_report.json', run / 'run_manifest.json', run / 'stage_plan.json',
                         run / 'mesh_qualification.json', work / 'n5_des.log'):
                if file.is_file():
                    shutil.copyfile(file, output / file.name)
        write_json(output / 'full-bv.json', report)
        if not report['passed']:
            raise ValueError(f'BVMethods acceptance failed. See {output / "full-bv.json"}')
    print(f"BVMethods: {report['accepted_points']} accepted points, BV={report['BV_inner_V']:.9g} V")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('stage', 'run'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--launcher', type=Path)
    parser.add_argument('--numerical-control', choices=NUMERICAL_CONTROLS)
    parser.add_argument('--tdx', help='TDX executable name or path')
    parser.add_argument('--tdr', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.action == 'stage':
        if bool(args.tdx) != bool(args.tdr):
            parser.error('--tdx and --tdr must be supplied together')
        stage_inputs(args.output.resolve(), Path(shutil.which(args.tdx) or args.tdx).resolve() if args.tdx else None, args.tdr)
    else:
        if args.launcher is None:
            parser.error('--launcher is required for run')
        run_case(args.output.resolve(), args.launcher.resolve(), args.numerical_control)


if __name__ == '__main__':
    main()
