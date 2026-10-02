"""First restoration stage: paired Masetti/no-BGN and Masetti/OldSlotboom DC.

Remote execution is explicit. Inputs, acceptance and failed outputs are immutable.
"""
import argparse
import math
import shutil
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import simplemos_masetti_curve_native_20260908 as prior

a, d, exporter = prior.a, prior.d, prior.exporter
REPO = prior.REPO
LOCAL = REPO / 'build-release/simplemos_bgn_restore_20260908'
OUT = REPO / 'reference_tcad/simplemos_sentaurus2022/bgn_restore_20260908'
REMOTE = '/tmp/vela_simplemos_bgn_restore_20260908'
INDICES = (0, 10, 40, 50)
MODELS = ('no_bgn', 'old_slotboom')


def deck(header, model):
    assert model in MODELS
    assert header.count('EffectiveIntrinsicDensity(NoBandGapNarrowing)') == 1
    assert 'Mobility(DopingDependence)' in header and prior.prior.MATH in header
    if model == 'old_slotboom':
        header = header.replace('EffectiveIntrinsicDensity(NoBandGapNarrowing)',
                                'EffectiveIntrinsicDensity(BandGapNarrowing(OldSlotboom))')
    result = header + 'Solve {\n'
    for index in INDICES:
        result += f''' Load(FilePrefix="vg_{index:03d}")
 NewCurrentPrefix="check_{index:03d}_"
 Coupled {{ Poisson Electron Hole }}
 Plot(FilePrefix="result_{index:03d}")
 Save(FilePrefix="result_{index:03d}")
'''
    return result + '}\n'


def prepare():
    assert not (OUT / 'native_freeze.json').exists()
    a.verify(prior.OUT / 'native_evidence.json')
    cases = a.read(prior.OUT / 'native_contract.json')['cases']
    files = [Path(__file__).resolve(), Path(prior.__file__), Path(exporter.__file__),
             prior.OUT / 'native_evidence.json']
    for case in cases:
        src = prior.LOCAL / 'native_raw/bundle' / case['case']
        header = (src / 'native_des.cmd').read_text().split('Solve {', 1)[0]
        for model in MODELS:
            dest = LOCAL / 'bundle' / model / case['case']
            dest.mkdir(parents=True, exist_ok=False)
            names = ['input_fps.tdr'] + [f'vg_{i:03d}{suffix}.sav' for i in INDICES
                                       for suffix in ('_des', '_circuit_des')]
            for name in names:
                shutil.copyfile(src / name, dest / name)
                files += [src / name, dest / name]
            (dest / 'native_des.cmd').write_text(deck(header, model), newline='\n')
            files += [src / 'native_des.cmd', dest / 'native_des.cmd']
    # Four independent simulator processes; each deck has four saved-bias targets.
    shell = '''#!/bin/bash
set -u
cd "$(dirname "$0")"
count=0
for path in bundle/*/*; do
 (
  cd "$path"
  /atctools/Synopsys/tcad/T-2022.03/bin/sdevice native_des.cmd > console.log 2>&1
  code=$?
  printf '%s\\n' "$code" > exit_code.txt
 ) &
 count=$((count+1))
 if [ "$((count % 4))" -eq 0 ]; then wait; fi
done
wait
tar -czf results.tgz bundle
printf 'complete\\n' > complete.txt
'''
    (LOCAL / 'run.sh').write_text(shell, newline='\n')
    files.append(LOCAL / 'run.sh')
    a.write(OUT / 'native_contract.json', dict(
        cases=[dict(case=c['case'], device=c['device'], vd=c['vd']) for c in cases],
        models=MODELS, indices=INDICES, vg=[prior.GRID[i] for i in INDICES],
        target_states=32, remote_root=REMOTE, math=prior.prior.MATH,
        single_axis='Only EffectiveIntrinsicDensity: explicit NoBandGapNarrowing versus BandGapNarrowing(OldSlotboom).',
        fixed='Masetti, SRH(DopingDependence), mesh, saved no-BGN input states, temperature, contacts and strengthened native Math.',
        gate=dict(native_exit_zero=True, runtime='T-2022.03-SP2', bias_error_V=1e-10,
                  kcl_over_Id=1e-8, bgn_max_abs_eV=1e-12),
        path='At each bias independently load qualified no-BGN native state and reclose. No model ramp or tolerance relaxation.',
        semantic_audit='Record actual runtime model, deltaEg(total impurity), effective ni, base ni inferred from exported Boltzmann identity; never fit material parameters.',
        acceptance_changed=False, defaults_changed=False))
    d.matrix.freeze(OUT / 'native_freeze.json', files + [OUT / 'native_contract.json'])
    with tarfile.open(LOCAL / 'input.tgz', 'w:gz') as tar:
        tar.add(LOCAL / 'bundle', arcname='bundle')
        tar.add(LOCAL / 'run.sh', arcname='run.sh')
    print('Prepared 8 decks / 32 native target states:', LOCAL / 'input.tgz', flush=True)


def unpack():
    a.verify(OUT / 'native_freeze.json')
    root = LOCAL / 'native_raw'
    assert not root.exists()
    with tarfile.open(LOCAL / 'results.tgz', 'r:gz') as tar:
        for item in tar.getmembers():
            if not (root/item.name).resolve().is_relative_to(root.resolve()) or item.issym() or item.islnk():
                raise ValueError('Unsafe archive member')
        tar.extractall(root, filter='data')
    for path in (LOCAL / 'bundle').rglob('*'):
        if path.is_file():
            assert a.sha(path) == a.sha(root/'bundle'/path.relative_to(LOCAL/'bundle'))
    points, jobs = [], []
    for model in MODELS:
        for c in a.read(OUT / 'native_contract.json')['cases']:
            raw = root / 'bundle' / model / c['case']
            code = int((raw / 'exit_code.txt').read_text())
            log = (raw / 'console.log').read_text(errors='replace')
            good = code == 0 and 'T-2022.03-SP2' in log and 'Good Bye' in log
            for index in INDICES:
                path = raw / f'check_{index:03d}_native_des.plt'
                row = dict(**c, model=model, index=index, vg=prior.GRID[index],
                           exit_code=code, native_qualified=False)
                if path.exists():
                    values = exporter.pltrows(path)
                    assert len(values) == 1, (path, len(values))
                    p = values[0]
                    assert abs(p['gate OuterVoltage']-row['vg']) <= 1e-10
                    assert abs(p['drain OuterVoltage']-c['vd']) <= 1e-10
                    currents = [p[t+' TotalCurrent'] for t in ('drain', 'source', 'gate', 'substrate')]
                    kcl = abs(math.fsum(currents)) / max(abs(currents[0]), 1e-300)
                    row.update(Id_A_per_um=currents[0], electron_Id_A_per_um=p['drain eCurrent'],
                               hole_Id_A_per_um=p['drain hCurrent'], kcl_over_Id=kcl,
                               native_qualified=good and kcl <= 1e-8)
                    tdr = raw / f'result_{index:03d}_des.tdr'
                    assert tdr.exists()
                    jobs.append(dict(case=c['case'], model=model, index=index, tdr=str(tdr),
                                     export=str(LOCAL/'native_exports'/model/c['case']/f'vg_{index:03d}')))
                points.append(row)
    keys = list(dict.fromkeys(k for r in points for k in r))
    a.write_csv(OUT/'native_points.csv', [{k:r.get(k, '') for k in keys} for r in points])
    a.write(OUT/'export_contract.json', dict(jobs=jobs))
    d.matrix.freeze(OUT/'native_evidence.json', [OUT/'native_freeze.json', OUT/'native_points.csv',
        OUT/'export_contract.json', LOCAL/'results.tgz'] + [p for p in root.rglob('*') if p.is_file()])
    print('Native DC qualified', sum(r['native_qualified'] for r in points), '/', len(points), flush=True)


def export():
    a.verify(OUT/'native_evidence.json')
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(exporter.export_one, a.read(OUT/'export_contract.json')['jobs']))
    d.matrix.freeze(OUT/'export_evidence.json', [OUT/'native_evidence.json'] +
                    [p for p in (LOCAL/'native_exports').rglob('*') if p.is_file()])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('prepare', 'unpack', 'export'))
    globals()[parser.parse_args().action]()
