"""Collect HFS pilot/rest evidence; gate a replacement observer before use.

No SDevice process is launched here. Numerical identity is distinct from
qualifying the observer's last-call samples or their cell interpretation.
"""
import argparse
import math
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import prepare_simplemos_hfs_native_20260912 as p
import calibrate_simplemos_enormal_20260912 as calibration

a, d, L, O = p.a, p.d, p.L, p.O


def jobs(stage):
    contract = a.read(O / 'contract.json')
    return [j for j in contract['jobs']
            if (j['name'] in contract['pilot']) == (stage == 'pilot')]


def collect(stage, expected):
    a.verify(O / 'input_evidence.json')
    if stage == 'rest':
        a.verify(O / 'pilot_identity_evidence.json')
        assert a.read(O / 'pilot_identity_summary.json')['all_identity_qualified']
    archive = L / f'{stage}_results.tgz'
    assert expected and a.sha(archive) == expected
    raw = L / f'{stage}_raw'
    assert not raw.exists()
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            assert (raw / member.name).resolve().is_relative_to(raw.resolve())
            assert member.isfile() or member.isdir()
        tar.extractall(raw, filter='data')
    for source in (L / 'pmi').iterdir():
        if source.is_file():
            assert a.sha(source) == a.sha(raw / 'pmi' / source.name)
    points, exports = [], []
    for job in jobs(stage):
        root = raw / 'bundle' / job['name']
        for source in (L / 'bundle' / job['name']).iterdir():
            assert source.is_file()
            assert a.sha(source) == a.sha(root / source.name)
        code = int((root / 'exit_code.txt').read_text())
        log = (root / 'console.log').read_text(errors='replace')
        values = p.e.c.n.exporter.pltrows(root / 'native_des.plt')
        assert len(values) == 1
        value = values[0]
        currents = [value[k + ' TotalCurrent']
                    for k in ('drain', 'source', 'gate', 'substrate')]
        finite = all(math.isfinite(x) for x in currents)
        current = currents[0]
        kcl = abs(math.fsum(currents)) / max(abs(current), 1e-300)
        bias = max(abs(value['gate OuterVoltage'] - job['vg']),
                   abs(value['drain OuterVoltage'] - job['vd']))
        good = (code == 0 and finite and current != 0 and
                'T-2022.03-SP2' in log and 'Good Bye' in log and
                bias <= 1e-10 and kcl <= 1e-8)
        points.append(dict(**job, exit_code=code, Id_A_per_um=current,
                           kcl_over_Id=kcl, bias_error_V=bias, qualified=good))
        exports.append(dict(case=job['case'], index=job['index'],
                            tdr=str(root / 'final_des.tdr'),
                            export=str(L / f'{stage}_exports' / job['name'])))
    a.write_csv(O / f'{stage}_points.csv', points)
    a.write(O / f'{stage}_export_contract.json', dict(jobs=exports))
    summary = dict(points=len(points), qualified=sum(r['qualified'] for r in points),
                   all_native_qualified=all(r['qualified'] for r in points))
    a.write(O / f'{stage}_native_summary.json', summary)
    d.matrix.freeze(O / f'{stage}_native_evidence.json',
                    [Path(__file__).resolve(), O / 'input_evidence.json', archive] +
                    [O / f'{stage}_{n}' for n in
                     ('points.csv', 'export_contract.json', 'native_summary.json')] +
                    [f for f in raw.rglob('*') if f.is_file()])
    print(summary, flush=True)
    assert summary['all_native_qualified']


def export(stage):
    a.verify(O / f'{stage}_native_evidence.json')
    assert a.read(O / f'{stage}_native_summary.json')['all_native_qualified']
    items = a.read(O / f'{stage}_export_contract.json')['jobs']
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(p.e.c.n.exporter.export_one, items))
    d.matrix.freeze(O / f'{stage}_export_evidence.json',
                    [Path(__file__).resolve(), O / f'{stage}_native_evidence.json'] +
                    [f for f in (L / f'{stage}_exports').rglob('*') if f.is_file()])


def identity(stage):
    a.verify(O / f'{stage}_export_evidence.json')
    points = {r['name']: r for r in a.rows(O / f'{stage}_points.csv')}
    gates = a.read(O / 'contract.json')['gates']
    fields, checks = [], []
    for job in jobs(stage):
        if job['arm'] != 'observed':
            continue
        base_name = job['key'] + '_baseline'
        base, observed = [L / f'{stage}_exports' / name
                          for name in (base_name, job['name'])]
        # A missing export must fail rather than silently shrink the comparison.
        assert {f.name for f in (base / 'fields').glob('*.csv')} == {
            f.name for f in (observed / 'fields').glob('*.csv')}
        rows = calibration.identity(base, observed)
        assert rows
        fields.extend(dict(key=job['key'], **row) for row in rows)
        potential = max(r['max_absolute'] for r in rows if 'Potential' in r['field'])
        density = max(r['max_relative'] for r in rows
                      if r['field'].startswith(('eDensity_', 'hDensity_')))
        mobility = max(r['max_relative'] for r in rows
                       if r['field'].startswith(('eMobility_', 'hMobility_')))
        current = abs(float(points[job['name']]['Id_A_per_um']) /
                      float(points[base_name]['Id_A_per_um']) - 1)
        passed = (current <= gates['observer_Id_relative'] and
                  potential <= gates['observer_psi_V'] and
                  density <= gates['observer_density_relative'] and
                  mobility <= gates['observer_mobility_relative'])
        checks.append(dict(key=job['key'], fields=len(rows),
                           current_relative=current, potential_max_V=potential,
                           density_relative=density, mobility_relative=mobility,
                           qualified=passed))
    assert checks
    a.write_csv(O / f'{stage}_identity_fields.csv', fields)
    a.write_csv(O / f'{stage}_identity.csv', checks)
    summary = dict(pairs=len(checks), qualified=sum(r['qualified'] for r in checks),
                   all_identity_qualified=all(r['qualified'] for r in checks),
                   samples_final_state_qualified=False,
                   cell_formation_qualified=False,
                   scope='Native replacement identity only; no production or field-formation promotion.')
    a.write(O / f'{stage}_identity_summary.json', summary)
    d.matrix.freeze(O / f'{stage}_identity_evidence.json',
                    [Path(__file__).resolve(), O / 'contract.json', O / f'{stage}_export_evidence.json'] +
                    [O / f'{stage}_{n}' for n in
                     ('identity_fields.csv', 'identity.csv', 'identity_summary.json')])
    print(summary, checks, flush=True)
    assert summary['all_identity_qualified']


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('collect', 'export', 'identity'))
    parser.add_argument('--stage', choices=('pilot', 'rest'), required=True)
    parser.add_argument('--sha')
    args = parser.parse_args()
    if args.action == 'collect':
        collect(args.stage, args.sha)
    elif args.action == 'export':
        export(args.stage)
    else:
        identity(args.stage)
