#!/usr/bin/env python3
"""Prepare, optionally qualify and continue SLOT in one explicit electrical frame.

CLI voltages are absolute physical drain voltages, not solver-frame coordinates.
Each invocation creates a new bundle. Resuming uses the last qualified archive
and its base config as the next invocation's explicit inputs.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import run_slot_ldmos_noavalanche_continuation as direct
import slot_reference_frame as frame
import slot_poly_equilibrium as poly


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(args):
    document = json.loads(args.base_config.read_text(encoding='utf-8-sig'))
    direct.require_direct_device_bias(document)
    solver = document.get('solver', {})
    if (solver.get('method') != 'newton' or
        solver.get('carrier_row_convergence', {}).get('mode') != 'enforce' or
        solver.get('global_continuity_closure', {}).get('mode') != 'enforce'):
        raise ValueError('requires full Newton, enforced carrier rows and global closure')
    targets = [frame.finite(t) for t in args.targets]
    bias = frame.finite(args.initial_drain_voltage_V)
    if targets != sorted(set(targets)) or any(t <= bias for t in targets):
        raise ValueError('targets must be strictly increasing above the initial physical drain voltage')
    if bool(args.previous_state) != (args.previous_drain_voltage_V is not None):
        raise ValueError('previous state and physical voltage must be supplied together')
    if args.previous_source_origin_V is not None and not args.previous_state:
        raise ValueError('previous origin requires a previous state')
    if args.previous_state and frame.finite(args.previous_drain_voltage_V) >= bias:
        raise ValueError('previous physical drain voltage must be below the current voltage')
    resources = {}
    def resolve(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if key in {'mesh_file', 'materials_file', 'node_doping_file',
                           'geometry_file', 'avalanche_geometry_file',
                           'external_averagebox_couples_file'}:
                    path = (args.base_config.parent / item).resolve()
                    resources[str(path)] = sha(path)
                    value[key] = str(path)
                else:
                    resolve(item)
        elif isinstance(value, list):
            for item in value:
                resolve(item)
    resolve(document)
    restore_poly = getattr(args, 'restore_poly_equilibrium', False)
    corrector_mode = getattr(args, 'carrier_corrector_mode', 'carrier_only')
    if corrector_mode not in ('carrier_only', 'carrier_density', 'density_then_qf'):
        raise ValueError('unsupported carrier corrector mode')
    if corrector_mode == 'density_then_qf' and args.carrier_corrector_passes < 2:
        raise ValueError('density_then_qf requires at least two candidate passes')
    if restore_poly:
        if args.carrier_corrector_passes < 1:
            raise ValueError('Poly equilibrium candidates require a carrier corrector')
        poly._profile(document)
        poly._components(document, args.base_config.parent,
                         bias-document.get('potential_origin_V', 0.))
    translated = frame.translate_document(document, args.target_origin_V)
    args.output.mkdir(parents=True, exist_ok=False)
    records = []
    seeds = [('current', args.initial_state, bias, document)]
    if args.previous_state:
        previous_document = frame.translate_document(document,
            document.get('potential_origin_V', 0.) if args.previous_source_origin_V is None
            else args.previous_source_origin_V)
        seeds.insert(0, ('previous', args.previous_state,
                        args.previous_drain_voltage_V, previous_document))
    for name, source, voltage, source_document in seeds:
        candidate = args.output / (name+'.h5')
        info = frame.translate_candidate(source.resolve(), candidate, source_document,
                                         args.base_config.parent, translated, voltage)
        records.append(dict(name=name, source=str(source.resolve()), source_sha256=sha(source),
                            candidate=str(candidate), candidate_sha256=sha(candidate), **info))
    base_path = args.output/'base.json'
    base_path.write_text(json.dumps(translated, indent=2), encoding='utf-8')
    manifest = dict(schema='vela.slot_ldmos.reference_frame.v1', completed=False,
        prepared=True, resources=resources, source_config_sha256=sha(args.base_config),
        tool_hashes={str(Path(p).resolve()): sha(p) for p in
                     (__file__, frame.__file__, frame.states.__file__, frame.states.archive.__file__, poly.__file__)},
        restore_poly_equilibrium=restore_poly,
        carrier_corrector_mode=corrector_mode,
        base_sha256=sha(base_path), potential_origin_V=args.target_origin_V,
        targets_physical_drain_V=targets, seeds=records, phases=[])
    save(args.output, manifest)
    return translated, manifest


def save(output, manifest):
    # Small durable record per completed phase; no running-status claim.
    target = output/'manifest.json';temporary = output/'manifest.tmp'
    with temporary.open('w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, allow_nan=False)
        f.flush();os.fsync(f.fileno())
    os.replace(temporary, target)


def execute(args, base, manifest):
    for path, digest in manifest['resources'].items():
        if sha(path) != digest:
            raise ValueError('physical resource changed since preparation: '+path)
    manifest['runner_sha256'] = sha(args.runner)
    manifest['driver_sha256'] = sha(Path(direct.__file__))
    previous = None
    for seed in manifest['seeds']:
        if sha(seed['candidate']) != seed['candidate_sha256']:
            raise ValueError('translated candidate changed')
        folder = args.output/seed['name'];folder.mkdir(exist_ok=False)
        deck = copy.deepcopy(base);direct.redirect_diagnostics(deck, folder)
        config = folder/'base.json';config.write_text(json.dumps(deck, indent=2), encoding='utf-8')
        origin = manifest['potential_origin_V'];bias = seed['solver_drain_voltage_V']
        targets = [bias]
        if seed['name'] == 'current':
            targets += [v-origin for v in manifest['targets_physical_drain_V']]
        command = [sys.executable, str(Path(direct.__file__).resolve()),
            '--repo-root', str(args.repo_root), '--runner', str(args.runner),
            '--base-config', str(config), '--initial-state', seed['candidate'],
            '--initial-bias', str(bias), '--targets', *map(str, targets),
            '--max-newton-iterations', str(args.max_newton_iterations),
            '--carrier-corrector-passes', str(args.carrier_corrector_passes),
            '--carrier-corrector-mode', manifest.get('carrier_corrector_mode', 'carrier_only'),
            '--audit-local-continuity', '--save-rejected-states', '--secant-predictor',
            '--min-step-V', str(args.min_step_V), '--max-attempts', str(args.max_attempts),
            '--timeout-seconds', str(args.timeout_seconds)]
        if previous:
            command += ['--secant-previous-state', previous['state'],
                        '--secant-previous-bias', str(previous['bias'])]
        if manifest.get('restore_poly_equilibrium', False):
            command.append('--restore-poly-equilibrium')
        (folder/'invocation.json').write_text(json.dumps(command, indent=2), encoding='utf-8')
        with (folder/'driver.log').open('w', encoding='utf-8') as log:
            process = subprocess.run(command, cwd=args.repo_root, stdout=log, stderr=subprocess.STDOUT)
        (folder/'exit_code.txt').write_text(str(process.returncode), encoding='utf-8')
        path = folder/'reports'/('adaptive_'+direct.physics_branch(base)+'_summary.json')
        summary = json.loads(path.read_text()) if path.exists() else None
        phase = dict(name=seed['name'], exit_code=process.returncode,
                     completed=bool(summary and summary['completed']), points=[])
        if summary:
            source_bias = next(c.get('bias', 0.) for c in base['contacts'] if c['name']=='source')
            for attempt in summary['attempts']:
                if attempt['converged']:
                    voltage = attempt['target_bias_V']
                    phase['points'].append(dict(solver_drain_voltage_V=voltage,
                        physical_drain_voltage_V=voltage+origin, Vds_V=voltage-source_bias,
                        current_A_per_um=float(attempt['current_A_per_um']),
                        state=attempt['state_file'], state_sha256=sha(attempt['state_file'])))
        manifest['phases'].append(phase)
        save(args.output, manifest)
        if process.returncode != 0 or not phase['completed']:
            return 1
        previous = dict(state=summary['last_successful_state'], bias=summary['last_successful_bias_V'])
    manifest['completed'] = True
    save(args.output, manifest)
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('base-config', 'initial-state', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--initial-drain-voltage-V', type=float, required=True)
    p.add_argument('--target-origin-V', type=float, required=True)
    p.add_argument('--previous-state', type=Path)
    p.add_argument('--previous-drain-voltage-V', type=float)
    p.add_argument('--previous-source-origin-V', type=float)
    p.add_argument('--targets', nargs='*', type=float, default=[])
    p.add_argument('--run', action='store_true', help='qualify anchors, then advance; preparation alone grants no qualification')
    p.add_argument('--runner', type=Path)
    p.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--max-newton-iterations', type=int, default=24)
    p.add_argument('--carrier-corrector-passes', type=int, default=4)
    p.add_argument('--carrier-corrector-mode', choices=['carrier_only', 'carrier_density', 'density_then_qf'],
                   default='carrier_only', help='candidate correction coordinates; full coupled gates remain required')
    p.add_argument('--restore-poly-equilibrium', action='store_true',
                   help='guarded source-free Poly candidate after carrier correction; full gates still apply')
    p.add_argument('--max-attempts', type=int, default=48)
    p.add_argument('--min-step-V', type=float, default=.00625)
    p.add_argument('--timeout-seconds', type=float, default=300.)
    args = p.parse_args()
    if args.run and not args.runner:
        p.error('--run requires --runner')
    if args.max_newton_iterations < 24 or args.max_attempts < 1 or args.carrier_corrector_passes < 0 or not (frame.finite(args.min_step_V)>0 and frame.finite(args.timeout_seconds)>0):
        p.error('invalid bounded execution settings (SLOT requires at least 24 Newton iterations)')
    for key in ('base_config', 'initial_state', 'output', 'repo_root', 'runner', 'previous_state'):
        value = getattr(args, key)
        if value is not None:
            setattr(args, key, value.resolve())
    base, manifest = prepare(args)
    code = execute(args, base, manifest) if args.run else 0
    print(json.dumps(dict(prepared=True, completed=manifest['completed'], manifest=str(args.output/'manifest.json'))))
    return code


if __name__ == '__main__':
    raise SystemExit(main())
