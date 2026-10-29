#!/usr/bin/env python3
"""Prepare and audit isolated native model/Math ablations; never qualify the original deck."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

from analyze_sentaurus_bvmethods import first_upward_crossing
from reference_artifacts import parse_version_banner
from validate_wp14_external_recovery import check_state
from validate_wp14_zero_qs import block, trajectory_checks
from vela_sdevice.plt.document import read_plt_document

REPO = Path(__file__).resolve().parents[1]
SOURCE = REPO / 'reference_tcad/bvmethods_sentaurus2018/source'
MESH_SHA256 = '7cd3968f5dd0e3cd63ef62ed0fcced465dd33bea465167f0be8139c16025902e'
VARIANTS = {
    'baseline': ('control', '', ''),
    'no_btbt': ('Physics', 'Band2Band(E2)', ''),
    'no_enormal': ('Physics', ' Enormal', ''),
    'no_hfs': ('Physics', ' HighFieldsaturation(GradQuasiFermi)', ''),
    'no_bgn': ('Physics', 'EffectiveIntrinsicDensity(OldSlotboom)',
               'EffectiveIntrinsicDensity(NoBandGapNarrowing)'),
    'boltzmann': ('Physics', '  Fermi\n', ''),
    'constant_srh': ('Physics', 'SRH(DopingDep)', 'SRH'),
    'no_avalanche': ('Physics', 'Avalanche(Eparallel)', ''),
    'no_extrapolate': ('Math', '  Extrapolate\n', ''),
    'no_aval_derivatives': ('Math', '  AvalDerivatives\n', '  -AvalDerivatives\n'),
    'no_mobility_derivatives': ('Math', 'Math {\n', 'Math {\n  -Derivatives\n'),
}
ANCHORS = (1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 1.3e-3, 1.443e-3)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_variant(source: str, name: str) -> str:
    section, before, after = VARIANTS[name]
    if name == 'baseline':
        return source
    if source.count(before) != 1:
        raise ValueError(f'{name}: expected exactly one source token {before!r}')
    result = source.replace(before, after, 1)
    for key in ('File', 'Electrode', 'Plot', 'Solve'):
        if block(result, key) != block(source, key):
            raise ValueError(f'{name}: unexpected change to {key}')
    other = 'Math' if section == 'Physics' else 'Physics'
    if block(result, other) != block(source, other):
        raise ValueError(f'{name}: unexpected change to {other}')
    return result



def classify_vela_control(config: dict) -> str:
    """Use the actual model/Jacobian selectors, never an experiment directory label."""
    solver = config['solver']
    mobility = solver['mobility']
    controls = []
    if solver['band_to_band']['model'] == 'none':
        controls.append('no_btbt')
    if mobility['model'] == 'masetti_lombardi':
        controls.append('no_hfs')
    if not mobility.get('jacobian_field_derivatives', True):
        controls.append('frozen_hfs_jacobian')
    if len(controls) > 1:
        raise ValueError('Vela screen must change only one selector')
    return controls[0] if controls else 'baseline'


def prepare(root: Path, mesh: Path) -> dict:
    if sha(mesh) != MESH_SHA256:
        raise ValueError('original native mesh identity mismatch')
    root.mkdir(parents=True, exist_ok=False)
    text = (SOURCE / 'external_resistor_sdevice.cmd').read_text(encoding='utf-8')
    plan = {'schema': 'vela.bvmethods.native_ablation.plan.v1',
            'scope': 'Single-factor native sensitivity and failure diagnosis; not original-model qualification or Vela/native agreement',
            'mesh_sha256': MESH_SHA256, 'parameter_sha256': sha(SOURCE / 'models.par'),
            'baseline_raw_cmd_sha256': sha(SOURCE / 'external_resistor_sdevice.cmd'),
            'baseline_lf_cmd_sha256': hashlib.sha256(text.encode()).hexdigest(),
            'temperature_K': 300, 'resistance_ohm_um': 1e7,
            'outer_goal_V': 100000, 'original_stop_A_per_um': 1.443e-3,
            'bv_extraction_A_per_um': 1e-4, 'variants': {}}
    for name, (section, before, after) in VARIANTS.items():
        case = root / name
        case.mkdir()
        cmd = case / 'external_resistor_sdevice.cmd'
        cmd.write_text(make_variant(text, name), encoding='utf-8', newline='\n')
        (case / 'models.par').write_bytes((SOURCE / 'models.par').read_bytes())
        (case / 'bvmethods_nmos_msh.tdr').write_bytes(mesh.read_bytes())
        plan['variants'][name] = {'section': section, 'removed': before, 'replacement': after,
                                  'cmd_sha256': sha(cmd)}
    (root / 'plan.json').write_text(json.dumps(plan, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    return plan


def summarize(root: Path) -> dict:
    plan = json.loads((root / 'plan.json').read_text())
    plan['base_plan_sha256'] = sha(root / 'plan.json')
    followup = root / 'followup_plan.json'
    if followup.exists():
        extra = json.loads(followup.read_text())
        for key in ('mesh_sha256', 'parameter_sha256', 'baseline_lf_cmd_sha256'):
            if extra[key] != plan[key]:
                raise ValueError('follow-up input mismatch: ' + key)
        if set(extra['variants']) & set(plan['variants']):
            raise ValueError('follow-up cannot replace an initial variant')
        plan['variants'].update(extra['variants'])
        plan['followup_plan_sha256'] = sha(followup)
    cases = {}
    for name, spec in plan['variants'].items():
        case = root / name
        paths = [case / n for n in ('external_resistor_sdevice.cmd', 'models.par',
                  'bvmethods_nmos_msh.tdr', 'external_resistor_des.plt', 'external_resistor_des.log',
                  'run.out', 'exit_code.txt', 'started_utc.txt', 'finished_utc.txt')]
        missing = [p.name for p in paths if not p.is_file()]
        item = {'missing_files': missing, 'original_model_qualified': False}
        cases[name] = item
        if missing:
            item['status'] = 'incomplete'
            continue
        if sha(paths[0]) != spec['cmd_sha256'] or sha(paths[1]) != plan['parameter_sha256'] or sha(paths[2]) != plan['mesh_sha256']:
            raise ValueError(f'{name}: input identity mismatch')
        item['sha256'] = {p.name: sha(p) for p in paths}
        log = paths[4].read_text(errors='replace')
        item['native_version'] = parse_version_banner(log)
        item['thread_settings'] = re.findall(r'^.*Number of (?:assembly|solver) threads:.*$', log, re.MULTILINE)
        physics_start = log.find('  Physical models:')
        physics_end = log.find('  Plot variables:', physics_start)
        item['physical_model_log_excerpt'] = log[physics_start:physics_end if physics_end > physics_start else physics_start + 5000]
        item['numerical_model_log_excerpt'] = log[max(0, log.find('  Numerical parameters:')):physics_start]
        item['process_exit_code'] = int((case / 'exit_code.txt').read_text())
        item['started_utc'] = (case / 'started_utc.txt').read_text().strip()
        item['finished_utc'] = (case / 'finished_utc.txt').read_text().strip()
        try:
            contacts, rows = read_plt_document(paths[3].read_text())
            if set(contacts) != {'drain', 'source', 'gate', 'substrate'}:
                raise ValueError('four complete contacts required')
            audit = trajectory_checks(rows, log)
            item['trajectory'] = audit
            physical = {k: v for k, v in audit['checks'].items()
                        if k not in ('original_stop_reached', 'native_current_stop_reported', 'native_success')}
            item['accepted_trajectory_checks_pass'] = all(physical.values())
            item['minstep_failure'] = bool(re.search(r'Step-size (?:less than MinStep|is too small)', log))
            item['timeout'] = item['process_exit_code'] == 124
            item['native_stop_reached'] = (item['process_exit_code'] == 0 and all(
                audit['checks'][key] for key in ('original_stop_reached', 'native_current_stop_reported', 'native_success')))
            item['status'] = ('stop_reached' if audit['status'] == 'pass' and item['process_exit_code'] == 0
                              else 'native_stop_reached_audit_failed' if item['native_stop_reached']
                              else 'timeout' if item['timeout']
                              else 'minstep_failure' if item['minstep_failure'] else 'failed_or_goal_only')
            currents = [abs(r['drain']['total_current_A_per_um']) for r in rows]
            item['current_anchors'] = []
            for current in ANCHORS:
                anchor = {'current_abs_A_per_um': current}
                crossing = next((i for i in range(1, len(currents))
                                 if currents[i - 1] < current <= currents[i]), None)
                anchor['accepted_sample_bracket'] = ([{
                    'current_abs_A_per_um': currents[i],
                    'inner_voltage_V': rows[i]['drain']['inner_voltage_V'],
                    'outer_voltage_V': rows[i]['drain']['outer_voltage_V']}
                    for i in (crossing - 1, crossing)] if crossing is not None else None)
                anchor['bracket_state_checks'] = ([check_state(rows[i], require_stop=False)
                                                   for i in (crossing - 1, crossing)] if crossing is not None else None)
                anchor['bracket_state_checks_pass'] = (all(state['status'] == 'pass' for state in anchor['bracket_state_checks'])
                                                      if crossing is not None else False)
                for contact in contacts:
                    anchor[contact] = {key: first_upward_crossing([r[contact][key] for r in rows], currents, current)
                                       for key in ('inner_voltage_V', 'outer_voltage_V', 'total_current_A_per_um')}
                item['current_anchors'].append(anchor)
            item['criterion_voltage_V'] = item['current_anchors'][2]['drain']['inner_voltage_V']
            item['avalanche_enabled'] = name != 'no_avalanche'
            with (case / 'trajectory.csv').open('w', newline='', encoding='utf-8') as stream:
                fields = ['time'] + [c + ':' + k for c in contacts for k in rows[0][c]]
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                writer.writerows({'time': r['time'], **{c + ':' + k: v for c in contacts for k, v in r[c].items()}} for r in rows)
        except (ValueError, KeyError) as error:
            item.update(status='invalid_trajectory', error=str(error))
    baseline = cases.get('baseline', {}).get('criterion_voltage_V')
    for item in cases.values():
        value = item.get('criterion_voltage_V')
        item['criterion_delta_vs_baseline_V'] = value - baseline if value is not None and baseline is not None else None
    result = {'schema': 'vela.bvmethods.native_ablation.results.v1', 'plan': plan,
              'cases': cases, 'interpretation': 'A native-only ablation measures sensitivity; it does not prove the cause of Vela/native differences. Failed scans retain usable audited accepted points; reaching the current stop alone cannot qualify the original model.'}
    (root / 'results.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return result


def plot_results(root: Path, result: dict) -> None:
    """Plot native accepted points; axes use V and A/um, never external voltage as BV."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for group in ('Physics', 'Math'):
        fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
        for name, spec in result['plan']['variants'].items():
            if name != 'baseline' and spec['section'] != group:
                continue
            if not result['cases'][name].get('trajectory'):
                continue
            _, rows = read_plt_document((root / name / 'external_resistor_des.plt').read_text())
            pairs = [(row['drain']['inner_voltage_V'], abs(row['drain']['total_current_A_per_um']))
                     for row in rows if 1e-8 <= abs(row['drain']['total_current_A_per_um']) <= 3e-3]
            if pairs:
                voltage, current = zip(*pairs)
                ax.semilogy(voltage, current, '.-', linewidth=1.5, markersize=3, label=name)
        ax.axhline(1e-4, color='black', linestyle='--', linewidth=1, label='100 uA/um extraction')
        ax.axhline(1.443e-3, color='gray', linestyle=':', linewidth=1, label='Original current stop')
        ax.set(xlabel='Drain InnerVoltage (V)', ylabel='Absolute drain current (A/um)',
               title=f'Native {group} controls: accepted points, not original-model qualification',
               ylim=(1e-8, 3e-3))
        ax.grid(True, which='both', alpha=.2)
        ax.legend(fontsize=8, loc='best')
        for suffix in ('png', 'svg'):
            fig.savefig(root / f'{group.lower()}_iv.{suffix}', dpi=180)
        if group == 'Physics':
            ax.set_xlim(5.8, 7.5)
            ax.set_title('Native Physics controls: near-threshold zoom')
            for suffix in ('png', 'svg'):
                fig.savefig(root / f'physics_iv_zoom.{suffix}', dpi=180)
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'summarize'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--mesh', type=Path)
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    if args.action == 'prepare':
        if args.mesh is None:
            parser.error('--mesh is required for prepare')
        result = prepare(args.root, args.mesh)
        print(json.dumps({'variants': list(result['variants']), 'root': str(args.root)}))
    else:
        result = summarize(args.root)
        if args.plot:
            plot_results(args.root, result)
        print(json.dumps({k: {f: v.get(f) for f in ('status', 'criterion_voltage_V', 'criterion_delta_vs_baseline_V', 'minstep_failure', 'accepted_trajectory_checks_pass')} for k, v in result['cases'].items()}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
