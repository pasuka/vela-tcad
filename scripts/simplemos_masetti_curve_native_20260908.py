"""Immutable 51-point Masetti native curves; remote execution is explicit."""
import argparse
import math
import shutil
import tarfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import prepare_simplemos_reference_subset_20260906 as prior
import export_simplemos_fullfield_native_20260906 as exporter

a = prior.a
d = prior.d
REPO = Path(__file__).resolve().parents[1]
LOCAL = REPO / 'build-release/simplemos_masetti_curves_20260908'
OUT = REPO / 'reference_tcad/simplemos_sentaurus2022/masetti_curves_20260908'
REMOTE = '/tmp/vela_simplemos_masetti_curves_20260908'
GRID = [round(i * .02, 12) for i in range(51)]


def deck(header):
    assert 'Mobility(DopingDependence)' in header and prior.MATH in header
    assert 'NoBandGapNarrowing' in header
    text = header + '''Solve {
 Load(FilePrefix="result_016")
 Coupled { Poisson Electron Hole }
 Quasistationary(InitialStep=0.05 Increment=1.4 MinStep=1e-6 MaxStep=0.1 Goal { Name="gate" Voltage=0 }) {
  Coupled { Poisson Electron Hole }
 }
'''
    for i, vg in enumerate(GRID):
        if i:
            text += f''' NewCurrentPrefix="ramp_{i:03d}_"
 Quasistationary(InitialStep=0.5 Increment=1.4 MinStep=1e-6 MaxStep=1 Goal {{ Name="gate" Voltage={vg:.12g} }}) {{
  Coupled {{ Poisson Electron Hole }}
 }}
'''
        text += f''' NewCurrentPrefix="point_{i:03d}_"
 Coupled {{ Poisson Electron Hole }}
 Plot(FilePrefix="vg_{i:03d}")
 Save(FilePrefix="vg_{i:03d}")
'''
    return text + '}\n'


def prepare():
    assert not (OUT / 'native_freeze.json').exists()
    files = [Path(__file__).resolve(), Path(prior.__file__), Path(exporter.__file__)]
    cases = []
    for device in ('n19', 'n23'):
        for vd in (.05, 1.):
            tag = f"m65_{device}_vd_{vd:.6f}_endpoint".replace('.', 'p')
            src = prior.LOCAL / 'native_raw/bundle/masetti' / tag
            dest = LOCAL / 'bundle' / tag
            dest.mkdir(parents=True, exist_ok=False)
            header = (src / 'native_des.cmd').read_text().split('Solve {', 1)[0]
            (dest / 'native_des.cmd').write_text(deck(header), newline='\n')
            for name in ('input_fps.tdr', 'result_016_des.sav', 'result_016_circuit_des.sav'):
                shutil.copyfile(src / name, dest / name)
                files += [src / name, dest / name]
            files += [src / 'native_des.cmd', dest / 'native_des.cmd']
            cases.append(dict(case=tag, device=device, vd=vd, vg=GRID))
    shell = '''#!/bin/bash
set -u
cd "$(dirname "$0")"
for path in bundle/*; do
 (
  cd "$path"
  sdevice native_des.cmd > console.log 2>&1
  code=$?
  printf '%s\\n' "$code" > exit_code.txt
  printf '%s %s\\n' "$path" "$code"
 )
done
tar -czf results.tgz bundle
printf 'complete\\n' > complete.txt
'''
    (LOCAL / 'run.sh').write_text(shell, newline='\n')
    files.append(LOCAL / 'run.sh')
    a.write(OUT / 'native_contract.json', dict(
        cases=cases, remote_root=REMOTE, curves=4, target_states=204,
        scope='Current accepted 300 K Masetti/no-BGN/SRH subset, n19/n23 and two drain biases; no full-model promotion.',
        physics='Header, native mesh and saved Masetti seed copied from accepted reference subset; no fitted parameters.',
        math=prior.MATH, gate_step_V=.02,
        path='Reclose accepted native Vg=.8 seed, ramp down to zero, then ascend to one. Export every target separately.',
        gates=dict(native_exit_zero=True, runtime='T-2022.03-SP2', bias_error_V=1e-10, kcl_over_Id=1e-8),
        acceptance_changed=False, production_changed=False))
    files.append(OUT / 'native_contract.json')
    d.matrix.freeze(OUT / 'native_freeze.json', files)
    with tarfile.open(LOCAL / 'input.tgz', 'w:gz') as tar:
        tar.add(LOCAL / 'bundle', arcname='bundle')
        tar.add(LOCAL / 'run.sh', arcname='run.sh')
    print('Frozen 4 native curves / 204 points:', LOCAL / 'input.tgz', flush=True)


def unpack():
    a.verify(OUT / 'native_freeze.json')
    root = LOCAL / 'native_raw'
    assert not root.exists()
    with tarfile.open(LOCAL / 'results.tgz', 'r:gz') as tar:
        for item in tar.getmembers():
            if not (root / item.name).resolve().is_relative_to(root.resolve()) or item.issym() or item.islnk():
                raise ValueError('Unsafe archive member')
        tar.extractall(root, filter='data')
    for path in (LOCAL / 'bundle').rglob('*'):
        if path.is_file():
            assert a.sha(path) == a.sha(root / 'bundle' / path.relative_to(LOCAL / 'bundle'))
    points, jobs = [], []
    for case in a.read(OUT / 'native_contract.json')['cases']:
        raw = root / 'bundle' / case['case']
        code = int((raw / 'exit_code.txt').read_text())
        log = (raw / 'console.log').read_text(errors='replace')
        good = code == 0 and 'T-2022.03-SP2' in log and 'Good Bye' in log
        for i, vg in enumerate(GRID):
            path = raw / f'point_{i:03d}_native_des.plt'
            row = dict(case=case['case'], device=case['device'], vd=case['vd'], vg=vg, index=i,
                       exit_code=code, native_qualified=False)
            if path.exists():
                values = exporter.pltrows(path)
                assert len(values) == 1, (path, len(values))
                v = values[0]
                assert abs(v['gate OuterVoltage']-vg) <= 1e-10
                assert abs(v['drain OuterVoltage']-case['vd']) <= 1e-10
                currents = [v[c+' TotalCurrent'] for c in ('drain', 'source', 'gate', 'substrate')]
                kcl = abs(math.fsum(currents)) / max(abs(currents[0]), 1e-300)
                row.update(Id_A_per_um=currents[0], kcl_over_Id=kcl, native_qualified=good and kcl <= 1e-8)
                tdr = raw / f'vg_{i:03d}_des.tdr'
                assert tdr.exists()
                jobs.append(dict(case=case['case'], index=i, tdr=str(tdr), export=str(LOCAL/'native_exports'/case['case']/f'vg_{i:03d}')))
            points.append(row)
    keys = list(dict.fromkeys(k for r in points for k in r))
    a.write_csv(OUT/'native_points.csv', [{k:r.get(k, '') for k in keys} for r in points])
    a.write(OUT/'export_contract.json', dict(jobs=jobs))
    d.matrix.freeze(OUT/'native_evidence.json', [OUT/'native_freeze.json', OUT/'native_points.csv',
                     OUT/'export_contract.json', LOCAL/'results.tgz'] + [x for x in root.rglob('*') if x.is_file()])
    print('Native qualified', sum(x['native_qualified'] for x in points), '/', len(points), flush=True)


def export():
    a.verify(OUT/'native_evidence.json')
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(exporter.export_one, a.read(OUT/'export_contract.json')['jobs']))
    d.matrix.freeze(OUT/'export_evidence.json', [OUT/'native_evidence.json'] +
                    [x for x in (LOCAL/'native_exports').rglob('*') if x.is_file()])


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('prepare', 'unpack', 'export'))
    globals()[p.parse_args().action]()
