"""Join sealed cloud, T470p low-Vd, and retry evidence without changing gates.

This is a numerical evidence audit, not a reconstruction of missing process
exit codes and not a declaration that the repository regression suite passed.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

CASES = {f'n{n}_vd_{vd}' for n in range(17, 25) for vd in ('0p05', '1')}
FIELDS = {'psi', 'phin', 'phip', 'electrons_m3', 'holes_m3', 'SRH'}
LIMITS = dict(psi_max_V=1e-6, phin_max_V=1e-6, phip_max_V=1e-6,
              density_relative=1e-4, Id_relative=1e-6)


def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify_manifest(root):
    manifest = read(root / 'download_manifest.json')
    for name, expected in manifest.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root.resolve()) or sha(path) != expected:
            raise ValueError(f'Evidence identity mismatch: {name}')
    return len(manifest)


def bounded(value, limit):
    return type(value) in (int, float) and math.isfinite(value) and abs(value) <= limit


def validate_curve(case, data):
    rows = data.get('results', [])
    if case not in CASES or data.get('passed') is not True or data.get('points') != 51 or len(rows) != 51:
        raise ValueError(f'Incomplete curve: {case}')
    if any(type(r.get('index')) is not int for r in rows) or {r['index'] for r in rows} != set(range(51)):
        raise ValueError(f'Duplicate or missing point: {case}')
    for r in rows:
        if r.get('case') != case or not bounded(r.get('vg', math.inf) - r['index'] * .05, 1e-12):
            raise ValueError(f'Wrong join key or bias: {case}')
        if r.get('passed') is not True or r.get('reconstruction_passed') is not True or r['dual'].get('qualified') is not True:
            raise ValueError(f'Failed qualification: {case}/{r["index"]}')
        if any(not bounded(r['dual'].get(k), limit) for k, limit in LIMITS.items()):
            raise ValueError(f'Dual gate exceeded: {case}/{r["index"]}')
        if set(r['Id_error_percent']) != {'cold', 'native'} or any(not bounded(v, 2.) for v in r['Id_error_percent'].values()):
            raise ValueError(f'Current gate exceeded: {case}/{r["index"]}')
    return rows


def validate_fields(rows, keys):
    triples = {(r['case'], int(r['index']), r['field']) for r in rows}
    expected = {(case, index, field) for case, index in keys for field in FIELDS}
    if len(rows) != len(expected) or triples != expected:
        raise ValueError('Duplicate or incomplete physical-field coverage')
    for r in rows:
        expected_unit = 'weighted_relative_L1' if r['field'] == 'SRH' else 'dex' if r['field'].endswith('_m3') else 'V'
        if r['unit'] != expected_unit or not math.isfinite(float(r['max_abs'])) or float(r['max_abs']) < 0:
            raise ValueError('Invalid field units or nonfinite field difference')


def csvout(path, rows):
    names = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=names)
        writer.writeheader(); writer.writerows(rows)


def join(cloud, low, retry, core_manifest, output):
    roots = {'cloud': cloud, 'low': low, 'retry': retry}
    verified = {k: verify_manifest(v) for k, v in roots.items()}
    # Linux/Windows binaries differ; source identity is checked separately.
    expected_core = {k: v['raw'] for k, v in read(core_manifest).items()}
    if read(low / 'core_sha256.json') != expected_core:
        raise ValueError('Cross-platform core source mismatch')
    if (cloud / 'job.exit').read_text().strip() != '0':
        raise ValueError('Cloud did not exit successfully')
    if (low / 'retry_replay/replay_cached.exit').read_text().strip() != '0':
        raise ValueError('Independent frozen-driver replay did not exit successfully')
    # PowerShell redirects stdout to UTF-16 on this host.
    raw = (low / 'retry_replay/replay_cached.stdout.log').read_bytes()
    log = raw.decode('utf-16') if raw.startswith(b'\xff\xfe') else raw.decode('utf-8-sig')
    replay = json.loads(log.strip().splitlines()[-1])
    if replay.get('passed') is not True or replay.get('unchanged_summary_sha256') != sha(retry/'matrix/n17_vd_1/summary.json'):
        raise ValueError('Replay is not linked to retry summary')
    if read(cloud/'source_verified.json')['runner_sha256'] != read(cloud/'supervisor/seal.json')['runner']:
        raise ValueError('Cloud runner identity mismatch')
    rows = []; fields = []; sources = []; seen = set()
    for label, root in roots.items():
        for path in sorted((root/'matrix').glob('*/summary.json')):
            case = path.parent.name
            if case in seen: raise ValueError('Overlapping case ownership: ' + case)
            seen.add(case)
            case_rows = validate_curve(case, read(path))
            rows.extend(case_rows)
            sources.append(dict(case=case, origin=label, summary_sha256=sha(path), points=51))
            if label != 'cloud':
                for r in case_rows:
                    point = path.parent/f'vg_{r["index"]:03d}'
                    if read(point/'comparison.json') != r: raise ValueError('Point summary differs')
                    f = read(point/'fields.json')
                    if f['reconstruction_passed'] is not True or not bounded(f['density_reconstruction'],1e-12) or not bounded(f['source_reconstruction'],1e-7):
                        raise ValueError('Density/SRH reconstruction failed')
                    fields.extend(dict(case=case,index=r['index'],**v) for v in f['fields'])
                    fields.append(dict(case=case,index=r['index'],field='SRH',unit='weighted_relative_L1',max_abs=f['SRH_weighted_L1']))
    if seen != CASES or len(rows) != 816: raise ValueError('Full matrix coverage missing')
    with (cloud/'supervisor/physical_fields.csv').open() as f: fields.extend(csv.DictReader(f))
    keys = {(r['case'],r['index']) for r in rows}
    validate_fields(fields,keys)
    worst = max(rows,key=lambda r:max(abs(v) for v in r['Id_error_percent'].values()))
    summary = dict(numerical_passed=True,points=816,curves=16,field_rows=len(fields),
        verified_files=verified,core_files=len(expected_core),source_identity='raw SHA256 identical',
        max_abs_Id_error_percent=max(abs(v) for r in rows for v in r['Id_error_percent'].values()),
        worst_current=dict(case=worst['case'],vg=worst['vg'],errors=worst['Id_error_percent']),
        dual_max={k:max(r['dual'][k] for r in rows) for k in LIMITS},
        original_retry_exit='missing: preserved empty file, not repaired or inferred',
        retry_readonly_driver_replay_exit=0,replay=replay,full_regression_passed=False,
        native_fields='descriptive comparison; no new acceptance threshold',
        remaining=['Full CTest','Code review and commit','Mobility field definition coverage','Controlled performance work'])
    output.mkdir(parents=True,exist_ok=True)
    csvout(output/'comparison.csv',[dict(case=r['case'],index=r['index'],vg=r['vg'],cold_Id_error_percent=r['Id_error_percent']['cold'],native_Id_error_percent=r['Id_error_percent']['native'],**r['dual']) for r in rows])
    csvout(output/'physical_fields.csv',fields)
    for name,data in [('summary.json',summary),('sources.json',sources)]:
        (output/name).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    (output/'seal.json').write_text(json.dumps({p.name:sha(p) for p in output.iterdir() if p.name in ['comparison.csv','physical_fields.csv','summary.json','sources.json']},indent=2)+'\n')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('cloud','low','retry','core-manifest','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();print(json.dumps(join(a.cloud,a.low,a.retry,a.core_manifest,a.output),indent=2))
