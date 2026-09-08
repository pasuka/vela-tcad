"""Review and seal native channel FD outputs without changing frozen inputs."""
import argparse
import json
import math
from pathlib import Path

import validate_simplemos_native_channel_charge as run

old = run.old
EVIDENCE = run.OUT / 'evidence.json'
DOC = run.REPO / 'docs/validation/simplemos_native_channel_charge_fd_execution_2026-09-05.md'


def review():
    run.verify()
    cfg = old.read_json(run.CONTRACT)
    result = old.read_json(run.OUT / 'result.json')
    ledger = old.read_csv(run.OUT / 'case_ledger.csv')
    cal = old.read_csv(run.OUT / 'calibration.csv')
    support = old.read_csv(run.OUT / 'charge_support_audit.csv')
    expected_cases = dict(run.CASES)
    assert len(ledger) == len(expected_cases) == 5
    assert {r['name'] for r in ledger} == set(expected_cases)
    assert len(cal) == 2
    base = run.fields('zero')
    nodes = set(base['SpaceCharge'])
    assert len(nodes) == 942
    assert len(support) == 5 * len(nodes)
    assert len({(r['name'], int(r['node_id'])) for r in support}) == len(support)
    for name, charge in run.CASES:
        fields = run.fields(name)
        assert all(set(values) == nodes for values in fields.values())
        assert all(math.isfinite(v) for values in fields.values() for v in values.values())
        rows = [r for r in support if r['name'] == name]
        assert {int(r['node_id']) for r in rows} == nodes
        assert sum(r['selected'] == 'True' for r in rows) == 186
        for r in rows:
            for k in ('inferred_charge_cm_3', 'expected_charge_cm_3', 'error_cm_3', 'tolerance_cm_3'):
                assert math.isfinite(float(r[k]))
            assert float(r['expected_charge_cm_3']) == charge * (r['selected'] == 'True')
    currents = {r['name']: float(r['current_A_per_um']) for r in ledger}
    assert all(math.isfinite(v) and v > 0 for v in currents.values())
    assert all(int(r['exit_code']) == 0 for r in ledger)
    qualified = old.read_csv(run.upstream.ROOT / 'strict_flux_precision/case_ledger.csv')
    vela_rows = [r for r in qualified if r['device'] == 'n23' and float(r['vd']) == .05]
    assert len(vela_rows) == 1
    vela_current = float(vela_rows[0]['current_A_per_um'])
    comparisons = []
    for suffix, amplitude in (('1e13', 1e13), ('5e12', 5e12)):
        r = next(r for r in cal if float(r['amplitude_cm_3']) == amplitude)
        central = (currents['plus_' + suffix] - currents['minus_' + suffix]) / 2
        assert math.isclose(central, float(r['central_delta_A_per_um']), rel_tol=1e-14, abs_tol=0)
        vela = cfg['vela_qualified_prediction_A_per_um_at_1e13'] * amplitude / 1e13
        comparisons.append({'amplitude_cm_3': amplitude, 'native_central_delta_A_per_um': central,
            'vela_predicted_delta_A_per_um': vela, 'vela_to_native_absolute_response_ratio': vela / central,
            'vela_to_native_normalized_response_ratio': (vela / vela_current) / (central / currents['zero']),
            'native_fd_dlog10_Id_linear': central / currents['zero'] / math.log(10),
            'native_direction_qualified': r['pass'] == 'True'})
    old.write_csv(run.OUT / 'vela_native_comparison.csv', comparisons)
    assert result['passed_amplitudes'] == sum(r['pass'] == 'True' for r in cal)
    assert result['passed'] == all(r['pass'] == 'True' for r in cal)
    assert result['charge_support_failed_rows'] == sum(r['pass'] != 'True' for r in support)
    assert not result['m82_released'] and not result['m83_released']
    print(json.dumps({'reviewed': True, 'support_rows': len(support), 'comparisons': comparisons}, indent=2))


def seal():
    review()
    if EVIDENCE.exists():
        raise FileExistsError(EVIDENCE)
    paths = [Path(__file__), run.SCRIPT, DOC, run.CONTRACT,
             run.LOCAL / 'input.tgz', run.LOCAL / 'results.tgz']
    for folder in (run.OUT, run.LOCAL / 'raw', run.LOCAL / 'exports'):
        paths.extend(p for p in folder.rglob('*') if p.is_file())
    cfg = old.read_json(run.CONTRACT)
    paths.extend(run.REPO / p for p in cfg['input_hashes'])
    old.write_json(EVIDENCE, {'status': old.read_json(run.OUT / 'result.json')['status'],
        'native_dc_runs': 5, 'new_ac_runs': 0, 'new_bias_sweeps': 0,
        'authorization': 'User explicitly approved this five-case payload, destination, execution and retrieval.',
        'remote_root': run.REMOTE, 'm82_released': False, 'm83_released': False,
        'hashes': {old.portable(p): old.sha256(p) for p in sorted(set(paths))}})
    print('Sealed native channel FD evidence')


def verify():
    run.verify()
    for p, digest in old.read_json(EVIDENCE)['hashes'].items():
        assert old.sha256(run.REPO / p) == digest, p
    print('Native channel FD evidence verified')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('review', 'seal', 'verify'))
    {'review': review, 'seal': seal, 'verify': verify}[parser.parse_args().action]()
