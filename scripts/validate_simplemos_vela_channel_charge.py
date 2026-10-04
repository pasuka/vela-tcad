"""Isolated, frozen Poisson-only nodal charge FD on a strict Vela baseline."""
import argparse
import copy
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np
import validate_simplemos_native_channel_charge as native
import seal_simplemos_native_channel_charge as native_seal
import run_simplemos_convergence_audit as audit
import run_simplemos_strict_flux_precision as precision

REPO, ROOT = native.REPO, native.ROOT
LOCAL = REPO / 'build-release/simplemos_vela_channel_charge_fd'
OUT = ROOT / 'vela_channel_charge_fd'
RUNNER = LOCAL / 'charge_runner.exe'
CONTRACT = OUT / 'contract.json'
FREEZE = OUT / 'freeze.json'
SOURCE_ENV = 'VELA_DIAGNOSTIC_NODAL_CHARGE_C_PER_M_FILE'
BASE = precision.LOCAL / 'm65_n23_vd_0p050000_endpoint'
CASES = native.CASES

HOOK = r'''
    // Diagnostic overlay only: integrated physical fixed charge in C/m.
    // Constant in state, so it has no Jacobian, carrier, mobility or doping term.
    if (const char* filename = std::getenv("VELA_DIAGNOSTIC_NODAL_CHARGE_C_PER_M_FILE")) {
        std::ifstream input(filename);
        std::string line;
        if (!input || !std::getline(input, line) || line != "node_id,charge_C_per_m")
            throw std::runtime_error("Invalid diagnostic charge file/header");
        std::vector<bool> seen(mesh_.numNodes(), false);
        while (std::getline(input, line)) {
            if (line.empty()) continue;
            const auto comma = line.find(',');
            if (comma == std::string::npos) throw std::runtime_error("Invalid charge row");
            const auto node = std::stoll(line.substr(0, comma));
            const double charge = std::stod(line.substr(comma + 1));
            if (node < 0 || static_cast<std::size_t>(node) >= seen.size() ||
                seen[node] || !std::isfinite(charge))
                throw std::runtime_error("Invalid/duplicate diagnostic charge node");
            seen[node] = true;
            fixedInterfaceChargeRhs_(node) += charge;
        }
        if (!std::all_of(seen.begin(), seen.end(), [](bool v) { return v; }))
            throw std::runtime_error("Missing diagnostic charge node");
    }
'''


def env(charge=None):
    e = os.environ.copy()
    e['PATH'] = 'D:/msys64/ucrt64/bin' + os.pathsep + e['PATH']
    e.pop(SOURCE_ENV, None)
    if charge is not None:
        e[SOURCE_ENV] = str(charge)
    return e


def build():
    if RUNNER.exists():
        raise FileExistsError(RUNNER)
    LOCAL.mkdir(parents=True, exist_ok=True)
    source = (REPO / 'src/equation/CoupledDDAssembler.cpp').read_text()
    marker = '    carrierStatisticsModel_ = carrierStatisticsModel(carrierStatistics_);'
    assert source.count(marker) == 1
    target = LOCAL / 'CoupledDDAssembler.cpp'
    target.write_text('#include <cstdlib>\n#include <fstream>\n' + source.replace(marker, HOOK + '\n' + marker), newline='\n')
    args = audit.read(REPO / 'build-release/simplemos_convergence_audit/build_command.json')
    where = next(i for i, v in enumerate(args) if v.endswith('libvela_core.a'))
    args.insert(where, str(target))
    args[-1] = str(RUNNER)
    audit.write(LOCAL / 'build_command.json', args)
    p = subprocess.run(args, capture_output=True, text=True, env=env())
    (LOCAL / 'build.log').write_text(p.stdout + p.stderr)
    if p.returncode:
        raise RuntimeError(p.stderr[-3000:])
    print('Built isolated charge assembler with original Newton algorithm', flush=True)


def prepare():
    native_seal.verify()
    audit.verify(precision.FREEZE)
    if CONTRACT.exists():
        raise FileExistsError(CONTRACT)
    old = native.old
    geo = native.upstream.m73.Geometry('n23')
    volumes, _ = native.upstream.native_arrays()
    masks = old.read_csv(native.OUT / 'node_mask.csv')
    selected = {int(r['node_id']) for r in masks if r['selected'] == 'True'}
    assert len(selected) == 186 and not (selected & set(geo.contact_nodes))
    c = audit.read(BASE / 'config.json')
    c['state_file'] = str(BASE / 'state.csv')
    c.pop('output_state_file')
    # Preserve historical branch provenance; this is an explicitly requested isolated validation.
    paths = [Path(__file__), RUNNER, LOCAL / 'CoupledDDAssembler.cpp', LOCAL / 'build_command.json',
             BASE / 'config.json', BASE / 'state.csv', precision.FREEZE,
             native.CONTRACT, native.OUT / 'evidence.json', native.OUT / 'node_mask.csv',
             native.OUT / 'calibration.csv', REPO / 'build-release/libvela_core.a']
    for v in audit.read(LOCAL / 'build_command.json'):
        if Path(v).is_file(): paths.append(Path(v))
    for key in ('mesh_file', 'node_doping_file', 'materials_file'): paths.append(Path(c[key]))
    for name, amp in CASES:
        dest = LOCAL / name
        deck = copy.deepcopy(c)
        deck['output_state_file'] = str(dest / 'state.csv')
        audit.write(dest / 'config.json', deck)
        rows = [{'node_id': i, 'charge_C_per_m': native.upstream.m73.Q * amp * 1e6 * volumes.get(i, 0) * 1e-12 * (i in selected)} for i in range(geo.count)]
        audit.write_csv(dest / 'charge.csv', rows)
        # C++ text reader uses LF on every platform.
        p = dest / 'charge.csv'; p.write_text(p.read_text(), newline='\n')
        paths.extend((dest / 'config.json', p))
    prior = native.upstream.LOCAL / 'm65_n23_vd_0p050000_endpoint'
    response = old.read_csv(prior / 'response.csv')
    known = np.array([float(r['source_poisson']) for r in response])
    state = old.read_csv(BASE / 'state.csv')
    density = np.array([float(r['electrons_m3']) for r in state])
    physical = native.upstream.m73.Q * density * (geo.volumes['barycentric_si'] - geo.volumes['all_cell'])
    active = geo.free[np.abs(physical[geo.free]) > np.max(np.abs(physical[geo.free])) * 1e-6]
    factor = float(np.median(known[active] / physical[active]))
    paths += [prior / 'response.csv', prior / 'adjoint.csv']
    cfg = audit.read(native.CONTRACT)
    audit.write(CONTRACT, {'status': 'frozen_before_execution', 'device': 'n23', 'vd': .05, 'vg': .9,
        'cases': [{'name': n, 'charge_cm_3': a} for n, a in CASES], 'selected_nodes': 186,
        'source': 'Physical nodal Q=q*amplitude_cm^-3*1e6*native_signed_Si_area_m2. Add Q to existing fixedInterfaceChargeRhs before scaling and Dirichlet replacement. Constant Poisson-only source; same nodal mask and signed Si areas as native calibration and prior Vela adjoint.',
        'source_environment_variable': SOURCE_ENV, 'frozen_poisson_residual_per_physical_charge': factor,
        'vela_prediction_at_1e13_A_per_um': cfg['vela_qualified_prediction_A_per_um_at_1e13'],
        'solver': c['solver'], 'global_audit_profile': audit.read(precision.CONTRACT)['global_profile'],
        'gates': {'all_states_strict_qualified': True, 'kcl_over_Id': 1e-8, 'zero_drift_dex': 1e-5,
            'fd_vs_adjoint_relative': .001, 'two_amplitude_derivative_relative': .001,
            'even_nonlinear_fraction': .01, 'minimum_signal_to_zero_drift': 100,
            'source_relative_L2': 1e-8, 'carrier_source_absolute': 0},
        'scope': 'Five local nonlinear reclosures only; original acceptance unchanged, no production algorithm change, no new Sentaurus run or curve sweep.',
        'm82_released': False, 'm83_released': False})
    paths.append(CONTRACT)
    audit.write(FREEZE, {'input_hashes': {(audit.rel(p) if p.is_relative_to(REPO) else str(p)): audit.sha(p) for p in sorted(set(paths))}})
    print('Frozen five local reclosures and original strict acceptance gates', flush=True)


def execute(deck, path, charge=None, runner=RUNNER):
    audit.write(path, deck)
    p = subprocess.run([str(runner), '--config', str(path), '--log', 'off'], env=env(charge), capture_output=True, text=True)
    path.with_suffix('.stdout.txt').write_text(p.stdout)
    path.with_suffix('.stderr.txt').write_text(p.stderr)
    status = json.loads(p.stdout.strip().splitlines()[-1])
    status['exit_code'] = p.returncode
    audit.write(path.with_suffix('.status.json'), status)
    return status


def preflight():
    audit.verify(FREEZE)
    cfg = audit.read(CONTRACT)
    rows = {}
    for name, _ in CASES:
        dest = LOCAL / name
        deck = audit.read(dest / 'config.json'); deck.pop('output_state_file')
        deck.update(simulation_type='newton_residual_probe', output_csv=str(dest / 'preflight.csv'))
        execute(deck, dest / 'preflight.json', dest / 'charge.csv')
        rows[name] = audit.rows(dest / 'preflight.csv')
    dest = LOCAL / 'unmodified_reference'
    deck = audit.read(LOCAL / 'zero/config.json'); deck.pop('output_state_file')
    deck.update(simulation_type='newton_residual_probe', output_csv=str(dest / 'residual.csv'))
    execute(deck, dest / 'config.json', runner=precision.previous.first.PRODUCTION)
    reference = audit.rows(dest / 'residual.csv')
    assert reference == rows['zero'], 'Zero overlay differs from original production residual'
    checks = []
    for name, _ in CASES:
        charge = audit.rows(LOCAL / name / 'charge.csv')
        delta, expected = [], []
        for z, r, q in zip(rows['zero'], rows[name], charge):
            assert z['node_id'] == r['node_id'] == q['node_id']
            for k in ('phin_residual', 'phip_residual', 'donors_m3', 'acceptors_m3', 'net_doping_m3', 'ni_eff_m3'):
                assert r[k] == z[k], (name, k)
            delta.append(float(r['psi_residual']) - float(z['psi_residual']))
            expected.append(-cfg['frozen_poisson_residual_per_physical_charge'] * float(q['charge_C_per_m']))
        delta, expected = np.array(delta), np.array(expected)
        error = float(np.linalg.norm(delta-expected) / max(np.linalg.norm(expected), 1e-300))
        assert error <= cfg['gates']['source_relative_L2'], (name, error)
        assert np.all(delta[expected == 0] == 0), 'Charge outside target nodes'
        checks.append({'name': name, 'source_relative_L2': error, 'nonzero_nodes': int(np.count_nonzero(delta)), 'carrier_rows_unchanged': True})
    audit.write_csv(OUT / 'source_preflight.csv', checks)
    print('Source verified: unchanged zero residual, Poisson-only target nodes, correct physical scale', flush=True)


def run():
    audit.verify(FREEZE)
    assert (OUT / 'source_preflight.csv').exists()
    cfg = audit.read(CONTRACT); results = []
    for name, amp in CASES:
        dest = LOCAL / name; deck = audit.read(dest / 'config.json')
        s = execute(deck, dest / 'solve.json', dest / 'charge.csv')
        check = copy.deepcopy(deck); check.pop('output_state_file')
        check['state_file'] = str(dest / 'state.csv')
        check['solver']['carrier_row_convergence']['mode'] = 'report'
        check['solver']['global_continuity_closure'] = cfg['global_audit_profile']
        check.update(simulation_type='newton_carrier_term_probe', output_csv=str(dest / 'terms.csv'), carrier_term_probe={'solved_equation_terms': True})
        a = execute(check, dest / 'acceptance.json', dest / 'charge.csv')
        currents = s['contact_currents_A_per_um']; current = currents['drain']
        kcl = abs(math.fsum(currents.values())) / abs(current)
        local, glob = a['carrier_row_convergence'], a['global_continuity_closure']
        r = {'name': name, 'amplitude_cm_3': amp, 'current_A_per_um': current, 'exit_code': s['exit_code'],
            'native_converged': s['converged'], 'iterations': s['iterations'], 'convergence_reason': s['convergence_reason'],
            'failure_reason': s['failure_reason'], 'final_residual': s['final_residual'],
            'local_satisfied': local['satisfied'], 'local_violations': local['violation_count'], 'local_max_ratio': local['max_ratio'],
            'global_satisfied': glob['satisfied'], 'electron_global_qualified': glob['electron']['qualified'],
            'hole_global_qualified': glob['hole']['qualified'], 'kcl_over_Id': kcl,
            'qualified': bool(s['exit_code'] == 0 and s['converged'] and local['satisfied'] and glob['satisfied'] and kcl <= cfg['gates']['kcl_over_Id'])}
        results.append(r); print(json.dumps(r), flush=True)
    audit.write_csv(OUT / 'case_ledger.csv', results)


def analyze():
    audit.verify(FREEZE); cfg = audit.read(CONTRACT); cases = audit.rows(OUT / 'case_ledger.csv')
    currents = {r['name']: float(r['current_A_per_um']) for r in cases}
    baseline = next(r for r in audit.rows(precision.OUT / 'case_ledger.csv') if r['device'] == 'n23' and float(r['vd']) == .05)
    native_currents = {r['name']: float(r['current_A_per_um']) for r in audit.rows(native.OUT / 'case_ledger.csv')}
    all_states = all(r['qualified'] == 'True' for r in cases)
    cal = []
    for suffix, amp in (('1e13', 1e13), ('5e12', 5e12)):
        r = native.compare(currents['plus_'+suffix], currents['minus_'+suffix], currents['zero'], cfg['vela_prediction_at_1e13_A_per_um'], amp, float(baseline['current_A_per_um']))
        r['zero_drift_dex'] = abs(math.log10(currents['zero']/float(baseline['current_A_per_um'])))
        s = (native_currents['plus_'+suffix] - native_currents['minus_'+suffix]) / 2
        r['native_central_delta_A_per_um'] = s
        r['vela_to_native_absolute_ratio'] = r['central_delta_A_per_um']/s
        r['vela_to_native_normalized_ratio'] = (r['central_delta_A_per_um']/currents['zero'])/(s/native_currents['zero'])
        cal.append(r)
    change = abs(2*cal[1]['central_delta_A_per_um']/cal[0]['central_delta_A_per_um']-1)
    g = cfg['gates']
    for r in cal:
        r['two_amplitude_derivative_relative'] = change
        r['pass'] = bool(all_states and r['relative_error'] <= g['fd_vs_adjoint_relative'] and r['same_sign'] and r['positive_perturbation_sign_pass'] and r['negative_perturbation_sign_pass'] and r['zero_drift_dex'] <= g['zero_drift_dex'] and r['even_nonlinear_fraction'] <= g['even_nonlinear_fraction'] and change <= g['two_amplitude_derivative_relative'] and r['signal_to_zero_control_drift'] >= g['minimum_signal_to_zero_drift'])
    audit.write_csv(OUT / 'calibration.csv', cal)
    result = {'status': 'passed' if all(r['pass'] for r in cal) else 'failed', 'qualified_states': sum(r['qualified']=='True' for r in cases),
        'passed_amplitudes': sum(r['pass'] for r in cal), 'max_fd_adjoint_relative': max(r['relative_error'] for r in cal),
        'two_amplitude_derivative_relative': change, 'vela_native_normalized_ratios': [r['vela_to_native_normalized_ratio'] for r in cal],
        'new_nonlinear_solves': 5, 'new_sentaurus_runs': 0, 'm82_released': False, 'm83_released': False}
    audit.write(OUT / 'result.json', result); print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('build','prepare','preflight','run','analyze','verify'))
    {'build':build,'prepare':prepare,'preflight':preflight,'run':run,'analyze':analyze,'verify':lambda:audit.verify(FREEZE)}[p.parse_args().action]()
