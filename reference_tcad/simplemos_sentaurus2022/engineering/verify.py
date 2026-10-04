"""Read-only audit of the shipped SimpleMOS engineering package; never solves."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_hashes(root, hashes):
    root = root.resolve()
    for name, expected in hashes.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root) or not path.is_file() or sha(path) != expected:
            raise ValueError('Missing or changed file: ' + name)


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def verify(root, inputs=None):
    verify_hashes(root, read(root / 'sha256.json'))
    contract = read(root / 'engineering_contract_v1.json')
    expected = {(f'n{n}_vd_{vd}', i) for n in range(17, 25)
                for vd in ('0p05', '1') for i in range(51)}
    for n in range(17, 25):
        rel = f'n{n}/template.json'
        require(sha(root / 'inputs' / rel) == contract['frozen_input_hashes'][rel], rel)
        cfg = read(root / 'inputs' / rel)
        require({k: cfg[k] for k in contract['profile']} == contract['profile'], 'Profile mismatch')
    joined = root / 'evidence/joined'
    verify_hashes(joined, read(joined / 'seal.json'))
    data = rows(joined / 'comparison.csv')
    keys = {(r['case'], int(r['index'])) for r in data}
    require(len(data) == 816 and keys == expected, 'Missing/duplicate comparison points')
    limits = dict(psi_max_V=contract['dual_gates']['potential_max_V'],
                  phin_max_V=contract['dual_gates']['potential_max_V'],
                  phip_max_V=contract['dual_gates']['potential_max_V'],
                  density_relative=contract['dual_gates']['density_relative'],
                  Id_relative=contract['dual_gates']['Id_relative'])
    current_limit = contract['engineering_reference']['current_error_percent']
    errors = []
    for row in data:
        require(row['qualified'] == 'True', 'Unqualified comparison point')
        require(abs(float(row['vg']) - int(row['index']) * .05) <= 1e-12, 'Wrong gate grid')
        for name, limit in limits.items():
            value = float(row[name])
            require(math.isfinite(value) and 0 <= value <= limit, 'Dual gate: ' + name)
        for name in ('cold_Id_error_percent', 'native_Id_error_percent'):
            value = float(row[name]); errors.append(abs(value))
            require(math.isfinite(value) and abs(value) <= current_limit, 'Current gate')
    fields = rows(joined / 'physical_fields.csv')
    names = {'psi', 'phin', 'phip', 'electrons_m3', 'holes_m3', 'SRH'}
    require(len(fields) == 4896 and {(r['case'], int(r['index']), r['field']) for r in fields}
            == {(c, i, f) for c, i in expected for f in names}, 'Field coverage mismatch')
    for row in fields:
        unit = 'weighted_relative_L1' if row['field'] == 'SRH' else 'dex' if row['field'].endswith('_m3') else 'V'
        require(row['unit'] == unit and math.isfinite(float(row['max_abs']))
                and float(row['max_abs']) >= 0, 'Invalid descriptive field statistic')
    sources = read(joined / 'sources.json')
    require(len(sources) == 16 and {r['case'] for r in sources} == {c for c, _ in expected}
            and all(r['points'] == 51 for r in sources), 'Source coverage mismatch')
    summary = read(root / 'summary.json'); old = read(joined / 'summary.json')
    require(summary['numerical_passed'] is True and old['numerical_passed'] is True
            and summary['points'] == old['points'] == 816
            and summary['curves'] == old['curves'] == 16
            and summary['field_rows'] == old['field_rows'] == 4896, 'Summary coverage mismatch')
    require(max(errors) == summary['max_abs_Id_error_percent'] == old['max_abs_Id_error_percent'],
            'Summary maximum differs')
    regression = root / 'evidence/regression'; status = read(regression / 'status.json')
    require(summary['full_regression_passed'] is True and summary['ctest_passed'] == 993
            and status['phase'] == 'passed' and status['tests'] == 993 and status['failures'] == 0
            and status['binary_before'] == status['binary_after'], 'Regression summary mismatch')
    for name in ('job.exit', 'ctest.exit'):
        require((regression / name).read_text().strip() == '0', 'Regression exit mismatch')
    if inputs is not None:
        verify_hashes(inputs, contract['frozen_input_hashes'])
    return dict(passed=True, points=816, curves=16, field_rows=4896,
                recorded_ctest_passed=993, external_inputs_checked=inputs is not None,
                scope='Shipped evidence audit only; no new numerical qualification or solver run')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, help='Optional restored frozen data/inputs directory')
    args = parser.parse_args()
    print(json.dumps(verify(Path(__file__).resolve().parent, args.inputs), indent=2))
