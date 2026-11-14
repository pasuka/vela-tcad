#!/usr/bin/env python3
"""Run the PN2D GRD/DAT reference suite through sd-sim."""
from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

import compare_pn2d_variants as variants
import pn2d_variants as campaign
import sdevice_acceptance as acceptance
import sentaurus_sde_cases as sde

ARCHIVE = campaign.FIXTURE / 'ci/pn2d-ci-reference-v1.tar.gz'
ARCHIVE_SHA256 = 'b7d6709b8988e4c5f3dbf2cb6bcbd9d5cbbdcc38917a7220cebc133699e1dd8b'
CONTACTS = ('Anode', 'Cathode')
GRID = re.compile(r'(\bGrid\s*=\s*")pn2d_msh\.(tdr|grd)(")')


def stage_reference(destination):
    if campaign.sha(ARCHIVE) != ARCHIVE_SHA256:
        raise ValueError('PN2D CI reference archive SHA-256 mismatch')
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(ARCHIVE) as archive:
        archive.extractall(destination, filter='data')


def stage_inputs(source, output, deck, spec, experiment, branch):
    text, count = GRID.subn(r'\g<1>pn2d_msh.grd\3', (source / deck).read_text())
    if count != 1:
        raise ValueError(f'{deck}: expected exactly one PN2D Grid declaration')
    native = GRID.sub(r'\g<1>pn2d_msh.tdr\3', text)
    if experiment is None:
        original = (campaign.FIXTURE / 'source' / deck).read_text()
        acceptance.augmented_source(original, native)
    else:
        cfg = {**spec['baseline'], **experiment['overrides']}
        mesh = {**spec['meshes']['M3'], **spec['sentaurus_numerics_by_stage'].get(experiment['stage'], {})}
        expected = campaign.sdevice_text(spec, cfg, branch, mesh)
        if not acceptance.ra.decks_equivalent(expected.encode(), native.encode()):
            raise ValueError(f'{experiment["id"]}/{branch}: archived CMD differs from the M3 generator')

    par = campaign.FIXTURE / 'source/models.par'
    if not acceptance.ra.decks_equivalent(par.read_bytes(), (source / par.name).read_bytes()):
        raise ValueError('Archived PN2D models.par differs from repository')
    for name in ('pn2d_msh.grd', 'pn2d_msh.dat', par.name):
        shutil.copyfile(source / name, output / name)
    (output / deck).write_text(text, encoding='utf-8')
    sidecar = Path(deck).with_suffix('.vela.json')
    if (source / sidecar).is_file():
        if deck == 'pn2d_bv_sdevice.cmd' and (
            acceptance.read_json(source / sidecar) != acceptance.read_json(acceptance.COMPENSATED_BV_SIDECAR)
        ):
            raise ValueError('Archived BV sidecar differs from the qualified compensated profile')
        shutil.copyfile(source / sidecar, output / sidecar)


def run_sdevice(launcher, output, deck, log, compensated):
    args = [str(launcher), '--physical-constants-profile', 'codata2018',
            '--threads', '1', '--allow-lossy', '--exit-on-failure']
    if compensated:
        args.append('--pn2d-population')
    env = dict(os.environ, VELA_PYTHON=sys.executable, VELA_LINEAR_SOLVER='sparselu',
               VELA_BIN=str(launcher.parent))
    with log.open('w', encoding='utf-8') as stream:
        subprocess.run([*args, deck], cwd=output, env=env, stdout=stream,
                       stderr=subprocess.STDOUT, check=True)
    run = output / (Path(deck).stem + '_vela')
    status = acceptance.read_json(run / 'run_manifest.json')
    plan = acceptance.read_json(run / 'stage_plan.json')['stages']
    if status['exit_code'] != 0 or len(status['stages']) != len(plan) or not all(
        stage['goal_reached'] and stage['runner_exit_code'] == 0 for stage in status['stages']
    ):
        raise ValueError('SDevice stages did not reach all goals')
    if compensated and not all(
        acceptance.read_json(run / stage['deck'])['solver'].get('compensated_population_state')
        for stage in plan
    ):
        raise ValueError('SDevice stage did not apply the compensated population state')
    return run, status, plan


def field_result(source, run, mesh, stage, bias, limits, digest):
    candidates = list(run.glob(f'stage{stage:02d}_*_{bias:g}V.vtu'))
    if len(candidates) != 1:
        raise ValueError(f'{bias:g} V: expected one field snapshot, found {len(candidates)}')
    native_fields = source / 'fields' / f'{bias:g}V'
    return {'bias_V': bias, **acceptance.field_metrics(
        native_fields, candidates[0], mesh, limits, bias, digest)}


def compare_base(case, source, deck, run, status, plan):
    contract = acceptance.read_json(acceptance.CONTRACT)
    limits = contract['cases'][f'pn2d_{case}']
    candidate = acceptance.raw_rows([run / stage['csv'] for stage in plan], 'Anode')
    native = acceptance.native_rows(source / 'reference.plt', 'Anode')
    targets = [0.] if case == '0v' else acceptance.requested_output_biases(deck, 'Anode')
    accepted = {round(row['bias_V'], 9) for row in candidate}
    missing = [v for v in targets if round(v, 9) not in accepted]
    expected_count = {'0v': 1, 'iv': 201, 'bv': 401}[case]
    coverage = not missing and len(targets) == expected_count
    curve = acceptance.curve_metrics(native, candidate, limits, contract,
                                     output_biases=targets if case == 'iv' else None)
    mesh = acceptance.read_json(run / 'mesh.json')
    fields = [field_result(source, run, mesh, 0, 0., contract['spatial'], status['sha256'])]
    if case == 'iv':
        fields.append(field_result(source, run, mesh, 1, .2, contract['spatial'], status['sha256']))
    if case == 'bv':
        profile = contract['biased_spatial_qualification']['pn2d_bv']
        fields.extend(field_result(source, run, mesh, 1, bias, profile['spatial'], status['sha256'])
                      for bias in profile['required_biases_V'])
    knee = (acceptance.breakdown_knee_metrics(native, candidate, limits['breakdown_check'])
            if case == 'bv' else {'status': 'not_required'})
    passed = (coverage and curve['status'] == 'pass' and knee['status'] != 'fail'
              and all(field['status'] == 'pass' for field in fields))
    return dict(passed=passed, target_count=len(targets), missing_biases_V=missing,
                curve=curve, fields=fields, knee=knee)


def compare_variant(source, run, status, plan, spec, branch):
    contract = campaign.read_json(campaign.CONTRACT)
    required = campaign.biases(spec['branches'][branch])
    decks = [campaign.read_json(run / stage['deck']) for stage in plan]
    declared = [float(deck['sweep']['stop']) for deck in decks]
    if len(required) != len(declared) or any(abs(a - b) > 1e-9 for a, b in zip(required, declared)):
        raise ValueError('Emitted stages differ from required exact biases')
    source_floor = decks[1]['solver']['global_continuity_closure']['source_floor']
    native = variants.plt_rows(source / 'reference.plt')
    direct = variants.rows(source / 'direct-terminal.csv')
    metrics = {'sentaurus': [], 'vela_json': []}
    conservation = []
    for index, voltage in enumerate(required):
        stage = max(index, 1)  # The first coupled sweep publishes the terminal state at 0 V.
        terminals = [r for r in variants.rows(run / f'stage{stage:02d}_terminal_balance.csv')
                     if r['converged'] in ('1', 'true')]
        raw = [r for r in variants.rows(run / plan[stage]['csv']) if r['converged'] in ('1', 'true')]
        states = {contact: variants.exact_row(terminals, voltage, contact=contact) for contact in CONTACTS}
        state = variants.exact_row(raw, voltage)
        oracle = variants.exact_row(native, voltage, 'Anode OuterVoltage')
        for contact in CONTACTS:
            reference = variants.exact_row(direct, voltage, contact=contact)
            for component, native_column in (('total', 'TotalCurrent'), ('electron', 'eCurrent'), ('hole', 'hCurrent')):
                column = f'current_{component}_A_per_um'
                sign = -1 if component == 'hole' else 1
                value = sign * float(states[contact][column])
                for key, expected in (('sentaurus', float(oracle[f'{contact} {native_column}'])),
                                      ('vela_json', sign * float(reference[column]))):
                    metrics[key].append(dict(bias_V=voltage, contact=contact, component=component,
                                             **variants.compare_current(expected, value, contract['current'])))
        currents = [float(states[c]['current_total_A_per_um']) for c in CONTACTS]
        closure = variants.continuity_satisfied(state, contract, states, source_floor)
        imbalance = abs(sum(currents))
        limit = max(contract['conservation']['absolute_floor_A_per_um'],
                    contract['conservation']['kcl_relative_max'] * max(map(abs, currents)))
        conservation.append(dict(bias_V=voltage, imbalance_A_per_um=imbalance,
                                 passed=closure and imbalance <= limit, closure=closure))

    rmse = {}
    for key, rows in metrics.items():
        errors = [row['log_error_dex'] for row in rows if row['log_error_dex'] is not None]
        rmse[key] = math.sqrt(sum(e * e for e in errors) / len(errors)) if errors else None
    spatial = {**contract['spatial'], 'coordinate_tolerance_um': contract['input']['coordinate_tolerance_um'],
               'junction_x_um': spec['geometry_um']['junction_x'], 'junction_half_width_um': .1}
    mesh = acceptance.read_json(run / 'mesh.json')
    fields = [field_result(source, run, mesh, required.index(bias), bias, spatial, status['sha256'])
              for bias in spec['saved_biases_V'] if (bias < 0) == (branch == 'reverse')]
    passed = (all(all(row['pass'] for row in rows) for rows in metrics.values())
              and all(v is not None and v <= contract['current']['log_rmse_max_dex'] for v in rmse.values())
              and all(row['passed'] for row in conservation)
              and all(field['status'] == 'pass' for field in fields))
    return dict(passed=passed, target_count=len(required), fields=fields, currents=metrics,
                current_log_rmse_dex=rmse, conservation=conservation)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case')
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--launcher', type=Path)
    args = parser.parse_args()
    case = args.case
    if case == 'stage':
        stage_reference(args.reference)
        return
    args.output.mkdir(parents=True, exist_ok=True)
    report = args.output / f'{case}.json'
    if case == 'sde':
        manifest = sde.load_manifest(sde.DEFAULT_MANIFEST)
        case_spec = next(c for c in manifest['cases'] if c['id'] == 'pn2d')
        with tempfile.TemporaryDirectory(prefix='sde-', dir=args.output) as folder:
            root = Path(folder)
            result = sde.run_case(case_spec, root / 'case', root / 'runtime', repeat=2)
        campaign.write_json(report, result)
        if result['status'] != 'generated':
            raise SystemExit(result.get('error', 'PN2D model generation failed'))
        return

    if args.launcher is None:
        parser.error('--launcher is required for SDevice cases')
    spec = campaign.read_json(campaign.SPEC)
    is_base = case in ('0v', 'iv', 'bv')
    experiment = branch = None
    if is_base:
        source = args.reference / 'base' / case
        deck = f'pn2d_{case}_sdevice.cmd'
    else:
        variant, branch = case.split('-')
        experiment = next(e for e in spec['experiments'] if e['id'] == variant)
        source = args.reference / 'variants' / variant / branch
        deck = f'{branch}_sdevice.cmd'

    log = args.output / f'{case}.log'
    with tempfile.TemporaryDirectory(prefix=f'{case}-', dir=args.output) as folder:
        output = Path(folder).resolve()
        stage_inputs(source, output, deck, spec, experiment, branch)
        try:
            run, status, plan = run_sdevice(args.launcher.resolve(), output, deck, log, case != '0v')
            if is_base:
                result = compare_base(case, source, output / deck, run, status, plan)
            else:
                result = compare_variant(source, run, status, plan, spec, branch)
        finally:
            work = output / (Path(deck).stem + '_vela')
            for name in ('fail_report.json', 'run_manifest.json'):
                if (work / name).is_file():
                    shutil.copyfile(work / name, args.output / f'{case}-{name}')
    campaign.write_json(report, result)
    if not result['passed']:
        raise SystemExit(f'PN2D {case}: acceptance failed. See {report}')
    print(f"PN2D {case}: {result['target_count']} targets, {len(result['fields'])} spatial checks passed.")


if __name__ == '__main__':
    main()
