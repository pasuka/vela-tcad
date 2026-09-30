"""Extend the explicit PhuMob box candidate to eight low-Vg controls.

All prior failures remain immutable. Remote execution is an explicit external
step; this script prepares and validates its exact input/output bundle.
"""
import argparse
import copy
import math
import shutil
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import validate_simplemos_phumob_box_qualified_20260909 as previous
import analyze_simplemos_phumob_box_20260909 as comparison
import simplemos_phumob_supported_export_20260908 as native

a, d, q = previous.a, previous.d, previous.q
REPO = previous.v.REPO
LOCAL = REPO / 'build-release/phumob_lowvg_20260909'
OUT = REPO / 'reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909'
NLOCAL = LOCAL / 'native'
NOUT = OUT / 'supported_export'
REMOTE = '/tmp/vela_simplemos_phumob_lowvg_20260909'
INDICES = (0, 10)


def prepare():
    a.verify(previous.OUT / 'completion_evidence.json')
    prior = a.read(previous.OUT / 'contract.json')
    files = [Path(__file__).resolve(), previous.OUT / 'completion_evidence.json',
             q.run.RUNNER, REPO / 'build-release/libvela_core.a']
    jobs = []
    for old in a.read(q.OUT / 'validation_contract.json')['jobs']:
        if old['model'] != 'old_slotboom' or old['index'] not in INDICES:
            continue
        job = dict(old, model='phumob')
        cfg = a.read(Path(old['original_config']))
        cfg['solver']['mobility'].update(model='phumob', edge_averaging='element_box_phumob')
        dest = LOCAL / 'inputs' / f"{old['case']}_{old['index']}" / old['arm'] / 'config.json'
        a.write(dest, cfg)
        job['original_config'] = str(dest)
        jobs.append(job)
        files += [Path(old['original_config']), dest, Path(job['seed'])]
        files += [Path(cfg[k]) for k in ('mesh_file', 'materials_file', 'node_doping_file')]
    assert len(jobs) == 16
    contract = dict(jobs=jobs, gates=prior['gates'], low_indices=INDICES,
                    physics='Plain PhuMob, OldSlotboom, matched base ni, 300 K; legacy versus explicit element_box_phumob only.',
                    initialization='The same independent Vela and native-potential seed paths used by the prior BGN controls. Saved guesses need not be accepted final PhuMob states; final outputs must pass the unchanged gates.',
                    recovery='At most one same-bias saved-state reload per failed attempt; preserve both.',
                    fixed_jacobian_qualification=str(previous.OUT / 'fixed_analysis_evidence.json'),
                    final_jvp_gate=1e-4, source_relative_floor=1e-10,
                    no_default_or_acceptance_change=True, native_G_floor_algorithm_qualified=False)
    a.write(OUT / 'contract.json', contract)
    d.matrix.freeze(OUT / 'input_freeze.json', files + [OUT / 'contract.json'])

    # Reuse the supported runtime callback; physical fields come from TDR Plot.
    oldnative = a.read(native.OUT / 'native_contract.json')
    one = oldnative['jobs'][0]
    callback = native.LOCAL / 'bundle' / one['model'] / one['key'] / 'runtime.tcl'
    tcl = callback.read_text()
    assert 'ReadVector' not in tcl and 'Potential' not in tcl
    nfiles = [callback, Path(native.__file__), OUT / 'input_freeze.json']
    njobs = []
    for c in a.read(native.original.b.OUT / 'native_contract.json')['cases']:
        source = native.original.b.LOCAL / 'native_raw/bundle/old_slotboom' / c['case']
        header = (source / 'native_des.cmd').read_text().split('Plot {', 1)[0]
        assert 'Mobility(DopingDependence)' in header
        for model in ('masetti_control', 'phumob'):
            for index in INDICES:
                key = f"{c['case']}_vg_{index:03d}"
                dest = NLOCAL / 'bundle' / model / key
                dest.mkdir(parents=True, exist_ok=False)
                for name in ('input_fps.tdr', f'result_{index:03d}_des.sav', f'result_{index:03d}_circuit_des.sav'):
                    shutil.copyfile(source / name, dest / name)
                    nfiles += [source / name, dest / name]
                text = header if model == 'masetti_control' else header.replace('Mobility(DopingDependence)', 'Mobility(PhuMob)')
                text += '''Plot { eDensity hDensity Potential eQuasiFermi hQuasiFermi
 eMobility/Element hMobility/Element eCurrent/Vector/Element hCurrent/Vector/Element
 Doping DonorConcentration AcceptorConcentration SRHRecombination
 BandGapNarrowing EffectiveIntrinsicDensity }
CurrentPlot { Tcl(tcl="source runtime.tcl") }
''' + native.original.b.prior.prior.MATH + f'''
Solve {{ Load(FilePrefix="result_{index:03d}") Coupled {{ Poisson Electron Hole }}
 Plot(FilePrefix="final") Save(FilePrefix="final") }}
'''
                (dest / 'native_des.cmd').write_text(text, newline='\n')
                (dest / 'runtime.tcl').write_text(tcl, newline='\n')
                nfiles += [source / 'native_des.cmd', dest / 'native_des.cmd', dest / 'runtime.tcl']
                njobs.append(dict(**c, key=key, model=model, index=index, vg=index * .02))
    shutil.copyfile(native.LOCAL / 'run.sh', NLOCAL / 'run.sh')
    ncontract = dict(oldnative, jobs=njobs, remote=REMOTE,
                     scope='Eight low-Vg PhuMob states and eight same-bias Masetti controls; qualified supported export callback; no other physics change.')
    a.write(NOUT / 'native_contract.json', ncontract)
    d.matrix.freeze(NOUT / 'native_freeze.json', nfiles + [NOUT / 'native_contract.json', NLOCAL / 'run.sh'])
    with tarfile.open(NLOCAL / 'input.tgz', 'w:gz') as tar:
        tar.add(NLOCAL / 'bundle', arcname='bundle')
        tar.add(NLOCAL / 'run.sh', arcname='run.sh')
    print('Prepared 16 native jobs and 32 local first attempts, unchanged gates.', flush=True)


def unpack():
    native.original.LOCAL, native.original.OUT = NLOCAL, NOUT
    native.original.unpack()
    assert all(r['native_qualified'] == 'True' for r in a.rows(NOUT / 'native_points.csv'))


def export():
    native.LOCAL, native.OUT = NLOCAL, NOUT
    native.export()


def solve(stage):
    a.verify(OUT / 'input_freeze.json')
    q.run.LOCAL, q.run.OUT = LOCAL / stage, OUT / stage
    jobs = []
    for old in a.read(OUT / 'contract.json')['jobs']:
        cfg = a.read(Path(old['original_config']))
        if stage == 'legacy':
            cfg['solver']['mobility']['edge_averaging'] = 'legacy'
        dest = LOCAL / stage / 'inputs' / f"{old['case']}_{old['index']}" / old['arm'] / 'config.json'
        a.write(dest, cfg)
        jobs.append(dict(old, original_config=str(dest)))
    a.write(q.run.OUT / 'validation_contract.json', dict(jobs=jobs))
    d.matrix.freeze(q.run.OUT / 'validation_freeze.json', [OUT / 'input_freeze.json', q.run.RUNNER, Path(q.run.__file__)] + [Path(j['original_config']) for j in jobs])
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = [r for rr in pool.map(q.run.solve_job, jobs) for r in rr]
    q.run.v.csv_union(q.run.OUT / 'attempts.csv', rows)
    d.matrix.freeze(q.run.OUT / 'dc_evidence.json', [q.run.OUT / 'validation_freeze.json', q.run.OUT / 'attempts.csv'] + [p for p in (q.run.LOCAL / 'dc').rglob('*') if p.is_file()])
    print(stage, 'accepted', sum(r['qualified'] for r in rows), '/', len(rows), flush=True)


def compare():
    a.verify(NOUT / 'native_evidence.json')
    assert all(r['native_qualified'] == 'True' for r in a.rows(NOUT / 'native_points.csv'))
    comparison.OUT, comparison.LOCAL = OUT, LOCAL
    comparison.v.previous.p.OUT = OUT
    comparison.compare()


def post():
    comparison.OUT, comparison.LOCAL = OUT, LOCAL
    comparison.post_jvp()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('prepare', 'unpack', 'export', 'legacy', 'candidate', 'compare', 'post'))
    action = parser.parse_args().action
    solve(action) if action in ('legacy', 'candidate') else globals()[action]()
