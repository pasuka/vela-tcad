"""Four-point extension of the sealed, isolated vector-mobility Jacobian experiment."""
import argparse
import copy
import json
import math
import subprocess
from decimal import Decimal
from pathlib import Path
import numpy as np
import validate_simplemos_vector_chain as chain
import compare_simplemos_qualified_charge_response as spatial
import seal_simplemos_vector_chain as sealed

p = chain.p
a = p.audit
REPO = p.REPO
LOCAL = REPO / 'build-release/simplemos_vector_chain_matrix_20260906'
OUT = p.ROOT / 'vector_chain_matrix'
CONTRACT = OUT / 'contract.json'
FREEZE = OUT / 'freeze.json'
DOC = REPO / 'docs/validation/simplemos_vector_chain_matrix_2026-09-06.md'
REUSED = 'm65_n23_vd_0p050000_endpoint'


def freeze(path, files):
    a.write(path, {'input_hashes': {a.rel(f): a.sha(f) for f in sorted(set(files))}})


def execute(path, enabled=False, charge=None, runner=chain.RUNNER):
    status_path = path.with_suffix('.status.json')
    if status_path.exists():
        return a.read(status_path)
    env = chain.env(enabled, charge)
    env['VELA_LINEAR_SOLVER'] = 'sparselu'
    proc = subprocess.run([str(runner), '--config', str(path), '--log', 'off'],
                          env=env, capture_output=True, text=True)
    path.with_suffix('.stdout.txt').write_text(proc.stdout)
    path.with_suffix('.stderr.txt').write_text(proc.stderr)
    try:
        status = json.loads(proc.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError) as exc:
        raise RuntimeError(f'No numerical status: {path}, exit {proc.returncode}') from exc
    status['exit_code'] = proc.returncode
    a.write(status_path, status)
    return status


def prepare():
    sealed.verify()
    spatial.verify()
    cases = a.rows(spatial.precision.OUT / 'case_ledger.csv')
    assert len(cases) == 4 and all(r['qualified'] == 'True' for r in cases)
    cfg = a.read(p.CONTRACT)
    files = [Path(__file__).resolve(), chain.RUNNER, p.RUNNER, p.CONTRACT,
             chain.OUT / 'evidence.json', spatial.CONTRACT,
             spatial.precision.OUT / 'case_ledger.csv']
    nv, _ = spatial.native_arrays()
    for c in cases:
        if c['case'] == REUSED:
            continue
        key = c['case']
        geo = spatial.m73.Geometry(c['device'])
        masks, _, _ = spatial.old.m78.supports(c['device'], geo, .05)
        free = np.zeros(geo.count, dtype=bool)
        free[geo.free] = True
        volume = geo.volumes['signed_si'].copy()
        if c['device'] == 'n23':
            volume = np.array([nv.get(i, 0) * 1e-12 for i in range(geo.count)])
        selected = masks['channel'] & free
        base = spatial.precision.LOCAL / key
        deck = a.read(base / 'config.json')
        deck['state_file'] = str(base / 'state.csv')
        files += [base / 'config.json', base / 'state.csv']
        files += [Path(deck[k]) for k in ('mesh_file', 'materials_file', 'node_doping_file')]
        c['nodes'] = geo.count
        c['selected_nodes'] = int(selected.sum())
        c['area_m2'] = float(np.sum(volume * selected))
        # Freeze scaling from the independent, already qualified volume-source probe.
        states = {int(r['node_id']): r for r in a.rows(base / 'state.csv')}
        n = np.array([float(states[i]['electrons_m3']) for i in range(geo.count)])
        prior = spatial.LOCAL / key
        known = np.array([float(r['source_poisson']) for r in a.rows(prior / 'response.csv')])
        phys = spatial.m73.Q * n * (geo.volumes['barycentric_si'] - geo.volumes['all_cell'])
        active = geo.free[np.abs(phys[geo.free]) > np.max(np.abs(phys[geo.free])) * 1e-6]
        c['factor'] = float(np.median(known[active] / phys[active]))
        files += [prior / 'response.csv', prior / 'adjoint.csv', prior / 'config.json']
        for name, amp in p.CASES:
            charge = LOCAL / key / 'sources' / (name + '.csv')
            a.write_csv(charge, [{'node_id': i, 'charge_C_per_m': float(spatial.m73.Q * amp * 1e6 * volume[i] * selected[i])} for i in range(geo.count)])
            charge.write_text(charge.read_text(), newline='\n')
            files.append(charge)
            for mode in ('off', 'on'):
                dest = LOCAL / key / mode / name
                d = copy.deepcopy(deck)
                d['output_state_file'] = str(dest / 'state.csv')
                a.write(dest / 'config.json', d)
                files.append(dest / 'config.json')
        d = a.read(prior / 'config.json')
        d['output_csv'] = str(LOCAL / key / 'adjoint/adjoint.csv')
        d['electron_volume_response']['output_csv'] = str(LOCAL / key / 'adjoint/response.csv')
        a.write(LOCAL / key / 'adjoint/config.json', d)
        files.append(LOCAL / key / 'adjoint/config.json')
    a.write(CONTRACT, {'status': 'frozen_before_execution', 'cases': cases,
        'reused_case': REUSED, 'backend': 'sparselu', 'gates': cfg['gates'],
        'current_identity_relative': 1e-8, 'adjoint_relative_residual': 1e-10,
        'jvp_minimum_defect_reduction': .99,
        'global_audit_profile': cfg['global_audit_profile'],
        'scope': '30 new local DC reclosures, three adjoints, actual-direction Jv; reuse ten n23 low-Vd states. Same physical residual, strict gates and isolated binary. n19 geometric signed Si area; n23 native signed Si area. No production modification or new Sentaurus run.',
        'm82_released': False, 'm83_released': False})
    files.append(CONTRACT)
    freeze(FREEZE, files)
    print('Frozen 30 new DC reclosures and three adjoints', flush=True)


def new_cases():
    return [c for c in a.read(CONTRACT)['cases'] if c['case'] != REUSED]


def preflight():
    a.verify(FREEZE)
    checks = []
    for c in new_cases():
        root = LOCAL / c['case']
        tables = {}
        for mode, name in [('off', n) for n, _ in p.CASES] + [('on', 'zero')]:
            dest = root / 'preflight' / (mode + '_' + name)
            deck = a.read(root / mode / name / 'config.json')
            deck.pop('output_state_file')
            deck.update(simulation_type='newton_residual_probe', output_csv=str(dest / 'residual.csv'))
            a.write(dest / 'config.json', deck)
            status = execute(dest / 'config.json', mode == 'on', root / 'sources' / (name + '.csv'))
            assert status['exit_code'] == 0
            tables[mode, name] = a.rows(dest / 'residual.csv')
        assert tables['off', 'zero'] == tables['on', 'zero']
        for name, _ in p.CASES:
            charge = a.rows(root / 'sources' / (name + '.csv'))
            delta, expected = [], []
            for z, r, q in zip(tables['off', 'zero'], tables['off', name], charge, strict=True):
                assert z['node_id'] == r['node_id'] == q['node_id']
                assert all(r[k] == z[k] for k in ('phin_residual', 'phip_residual', 'donors_m3', 'acceptors_m3', 'net_doping_m3', 'ni_eff_m3'))
                delta.append(float(r['psi_residual']) - float(z['psi_residual']))
                expected.append(-c['factor'] * float(q['charge_C_per_m']))
            delta, expected = np.array(delta), np.array(expected)
            err = float(np.linalg.norm(delta - expected) / max(np.linalg.norm(expected), 1e-300))
            assert err <= a.read(CONTRACT)['gates']['source_relative_L2']
            assert np.all(delta[expected == 0] == 0)
            checks.append({'case': c['case'], 'name': name, 'source_relative_L2': err, 'passed': True})
        print(c['case'], 'Poisson-only source and residual identity verified', flush=True)
    a.write_csv(OUT / 'source_preflight.csv', checks)


def run():
    a.verify(FREEZE)
    assert len(a.rows(OUT / 'source_preflight.csv')) == 15
    results = []
    for c in new_cases():
        root = LOCAL / c['case']
        for mode in ('off', 'on'):
            for name, amp in p.CASES:
                dest = root / mode / name
                charge = root / 'sources' / (name + '.csv')
                s = execute(dest / 'config.json', mode == 'on', charge)
                check = a.read(dest / 'config.json')
                check.pop('output_state_file')
                check['state_file'] = str(dest / 'state.csv')
                check['solver']['carrier_row_convergence']['mode'] = 'report'
                check['solver']['global_continuity_closure'] = a.read(CONTRACT)['global_audit_profile']
                check.update(simulation_type='newton_carrier_term_probe', output_csv=str(dest / 'terms.csv'), carrier_term_probe={'solved_equation_terms': True})
                if not (dest / 'acceptance.json').exists():
                    a.write(dest / 'acceptance.json', check)
                audit = execute(dest / 'acceptance.json', False, charge, p.RUNNER)
                cc = s['contact_currents_A_per_um']
                current = cc['drain']
                kcl = abs(math.fsum(cc.values())) / max(abs(current), 1e-300)
                local, glob = audit['carrier_row_convergence'], audit['global_continuity_closure']
                r = {'case': c['case'], 'mode': mode, 'name': name, 'amplitude_cm_3': amp,
                     'current_A_per_um': current, 'iterations': s['iterations'], 'reason': s['convergence_reason'],
                     'exit_code': s['exit_code'], 'local_violations': local['violation_count'],
                     'local_max_ratio': local['max_ratio'], 'kcl_over_Id': kcl,
                     'electron_global_qualified': glob['electron']['qualified'], 'hole_global_qualified': glob['hole']['qualified'],
                     'qualified': s['exit_code'] == 0 and audit['exit_code'] == 0 and s['converged'] and local['satisfied'] and glob['satisfied'] and kcl <= 1e-8}
                results.append(r)
                print(json.dumps(r), flush=True)
        s = execute(root / 'adjoint/config.json', True)
        assert s['exit_code'] == 0 and s['adjoint_relative_residual'] <= 1e-10
    a.write_csv(OUT / 'dc_ledger.csv', results)


def score(currents, pred, base, qualified, gates):
    """Independent odd/even-response and restart-drift qualification."""
    required = {name for name, _ in p.CASES}
    if set(currents) != required or any(not math.isfinite(v) or v <= 0 for v in currents.values()) or not math.isfinite(pred) or pred == 0:
        raise ValueError('Five finite positive currents and a nonzero prediction required')
    zero = currents['zero']
    drift = abs(zero - base)
    full = (currents['plus_1e13'] - currents['minus_1e13']) / 2
    half = (currents['plus_5e12'] - currents['minus_5e12']) / 2
    if full == 0 or half == 0:
        raise ValueError('No resolved odd response')
    linearity = abs(2 * half / full - 1)
    out = []
    for suffix, amp, fd in (('1e13', 1e13, full), ('5e12', 5e12, half)):
        target = pred * amp / 1e13
        plus, minus = currents['plus_' + suffix], currents['minus_' + suffix]
        even = abs((plus + minus) / 2 - zero) / abs(fd)
        ratio = abs(fd) / max(drift, 1e-300)
        err = abs(fd / target - 1)
        zero_dex = abs(math.log10(zero / base))
        passed = qualified and (plus-zero)*target > 0 and (minus-zero)*target < 0 and err <= gates['fd_vs_adjoint_relative'] and even <= gates['even_nonlinear_fraction'] and linearity <= gates['two_amplitude_derivative_relative'] and ratio >= gates['minimum_signal_to_zero_drift'] and zero_dex <= gates['zero_drift_dex']
        out.append({'amplitude_cm_3': amp, 'central_response_A_per_um': fd,
                    'prediction_A_per_um': target, 'relative_error': err,
                    'two_amplitude_relative': linearity, 'even_nonlinear_fraction': even,
                    'signal_to_zero_drift': ratio, 'zero_drift_dex': zero_dex, 'passed': bool(passed)})
    return out


def analyze():
    a.verify(FREEZE)
    cfg = a.read(CONTRACT)
    ledger = a.rows(OUT / 'dc_ledger.csv')
    cal, identity = [], []
    for c in cfg['cases']:
        key = c['case']
        values = {}
        for mode in ('off', 'on'):
            if key == REUSED:
                records = a.rows((p.OUT if mode == 'off' else chain.OUT) / ('case_ledger.csv' if mode == 'off' else 'dc_ledger.csv'))
                pred = (a.read(p.CONTRACT)['vela_prediction_at_1e13_A_per_um'] if mode == 'off' else a.read(chain.OUT / 'fixed_state_result.json')['prediction_A_per_um_at_1e13'])
            else:
                records = [r for r in ledger if r['case'] == key and r['mode'] == mode]
                folder = spatial.LOCAL / key if mode == 'off' else LOCAL / key / 'adjoint'
                adj = {int(r['node_id']): r for r in a.rows(folder / 'adjoint.csv')}
                charge = a.rows(LOCAL / key / 'sources/plus_1e13.csv')
                pred = math.fsum(float(adj[int(q['node_id'])]['lambda_poisson']) * c['factor'] * float(q['charge_C_per_m']) for q in charge)
            assert len(records) == 5
            currents = {r['name']: float(r['current_A_per_um']) for r in records}
            values[mode] = currents
            qualified = all(r['qualified'] == 'True' for r in records)
            for r in score(currents, pred, float(c['current_A_per_um']), qualified, cfg['gates']):
                cal.append(dict(case=key, device=c['device'], vd=float(c['vd']), mode=mode, strict_states=sum(r['qualified'] == 'True' for r in records), **r))
        for name, _ in p.CASES:
            err = abs(values['on'][name] / values['off'][name] - 1)
            identity.append({'case': key, 'name': name, 'relative_change': err, 'passed': err <= cfg['current_identity_relative']})
    a.write_csv(OUT / 'calibration.csv', cal)
    a.write_csv(OUT / 'current_identity.csv', identity)
    result = {'new_nonlinear_solves': len(ledger), 'new_strict_states': sum(r['qualified'] == 'True' for r in ledger),
        'on_passed_amplitudes': sum(r['passed'] for r in cal if r['mode'] == 'on'),
        'off_passed_amplitudes': sum(r['passed'] for r in cal if r['mode'] == 'off'),
        'max_on_fd_adjoint_relative': max(r['relative_error'] for r in cal if r['mode'] == 'on'),
        'max_current_relative_change': max(r['relative_change'] for r in identity),
        'passed': all(r['passed'] for r in cal if r['mode'] == 'on') and all(r['passed'] for r in identity),
        'new_sentaurus_runs': 0, 'm82_released': False, 'm83_released': False}
    a.write(OUT / 'result.json', result)
    print(json.dumps(result, indent=2), flush=True)


def jvp():
    a.verify(FREEZE)
    checks = []
    for c in new_cases():
        root = LOCAL / c['case']
        adj = {int(r['node_id']): r for r in a.rows(spatial.LOCAL / c['case'] / 'adjoint.csv')}
        weight = {'psi': 'lambda_poisson', 'phin': 'lambda_electron', 'phip': 'lambda_hole'}
        for suffix in ('1e13', '5e12'):
            dest = root / 'jvp' / suffix
            plus = root / 'off' / ('plus_' + suffix) / 'state.csv'
            minus = root / 'off' / ('minus_' + suffix) / 'state.csv'
            dirs = {k: [] for k in weight}
            for u, v in zip(a.rows(plus), a.rows(minus), strict=True):
                assert u['node_id'] == v['node_id']
                for k, prefix in (('psi', None), ('phin', 'electron'), ('phip', 'hole')):
                    if prefix:
                        ref, inc = prefix + '_qf_reference_V', prefix + '_qf_increment_V'
                        value = (Decimal(u[ref]) - Decimal(v[ref]) + Decimal(u[inc]) - Decimal(v[inc])) / 2
                    else:
                        value = (Decimal(u[k]) - Decimal(v[k])) / 2
                    dirs[k].append({'node_id': int(u['node_id']), 'component0': float(value)})
            paths = {}
            for k, rows in dirs.items():
                path = dest / (k + '.csv')
                a.write_csv(path, rows)
                paths[k] = str(path)
            inputs = [plus, minus] + [Path(v) for v in paths.values()]
            for mode in ('off', 'on'):
                deck = a.read(root / mode / 'zero/config.json')
                deck.pop('output_state_file')
                deck.update(simulation_type='actual_direction_jvp', plus_state_file=str(plus), minus_state_file=str(minus), direction_files=paths, output_csv=str(dest / mode / 'jvp.csv'))
                a.write(dest / mode / 'config.json', deck)
                inputs.append(dest / mode / 'config.json')
            freeze(dest / 'freeze.json', inputs)
            tables = {}
            for mode in ('off', 'on'):
                s = execute(dest / mode / 'config.json', mode == 'on')
                assert s['exit_code'] == 0
                tables[mode] = a.rows(dest / mode / 'jvp.csv')
                assert len(tables[mode]) == c['nodes'] * 36
            before, after = tables['off'], tables['on']
            same_f = all(u['fd'] == v['fd'] and u['actual_endpoint_fd'] == v['actual_endpoint_fd'] for u, v in zip(before, after, strict=True))
            assert same_f
            for h in (1, .5, .25):
                sums = {}
                for mode, rows in tables.items():
                    group = [r for r in rows if r['mode'] == 'all' and float(r['step']) == h]
                    assert len(group) == 3 * c['nodes']
                    sums[mode] = math.fsum(float(adj[int(r['node_id'])][weight[r['row_block']]]) * (float(r['analytic']) - float(r['fd'])) for r in group)
                reduction = 1 - abs(sums['on'] / sums['off']) if sums['off'] != 0 else None
                checks.append({'case': c['case'], 'amplitude_cm_3': float(suffix), 'step': h,
                    'old_weighted_defect_A_per_um': sums['off'], 'new_weighted_defect_A_per_um': sums['on'],
                    'residual_identical': same_f, 'defect_reduction': reduction,
                    'passed': same_f and reduction is not None and reduction >= .99})
            print(c['case'], suffix, 'Jv completed', flush=True)
    a.write_csv(OUT / 'jvp_checks.csv', checks)
    print(json.dumps({'passed': sum(r['passed'] for r in checks), 'checks': len(checks), 'minimum_reduction': min(r['defect_reduction'] for r in checks)}), flush=True)


def seal():
    a.verify(FREEZE)
    freeze(OUT / 'evidence.json', [DOC] + [f for root in (LOCAL, OUT) for f in root.rglob('*') if f.is_file()])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'preflight', 'run', 'analyze', 'jvp', 'seal', 'verify'))
    action = parser.parse_args().action
    if action == 'verify':
        a.verify(FREEZE)
        a.verify(OUT / 'evidence.json')
        print('Four-point evidence verified')
    else:
        globals()[action]()
