"""Package existing qualified native HFS data for an isolated portable replay.

Read-only historical adapters are used here, never on the compute host. No old
runner, freeze gate, result, or scientific setting is overwritten.
"""
import argparse
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import tarfile

import simplemos_hfs_curves_20260914 as old

R = Path(__file__).resolve().parents[1]
NAMES = ('ElectrostaticPotential', 'eQuasiFermiPotential',
         'hQuasiFermiPotential', 'eDensity', 'hDensity', 'srhRecombination')


def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8'))


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def main(dest):
    assert not dest.exists(), dest
    dest.mkdir(parents=True)
    manifest = read(R/'build/local_main_merge_20260925/codespaces_manifest.json')
    manifest.update(read(R/'build/local_main_merge_20260925/cloud_delta_manifest.json'))
    sources = {p: h for p, h in manifest.items()
               if p.startswith(('src/', 'include/', 'cmake/')) or p == 'CMakeLists.txt'}
    io_files={'src/io/StateArchive.cpp','src/io/DDSolutionState.cpp',
              'src/tools/vela_example_runner.cpp','tests/regression/test_hdf5_state_worker.py',
              'scripts/state_archive.py','scripts/migrate_state_seed_to_hdf5.py',
              'tests/test_state_archive.cpp','tests/regression/test_state_archive.py'}
    io_changes={}
    for p, h in sources.items():
        actual=sha(R/p)
        if actual != h:
            assert p in io_files, ('source_changed_since_merge_validation', p)
            io_changes[p]=dict(before=h,after=actual)
            sources[p]=actual
    for p in io_files:
        sources[p]=sha(R/p)
        if sources[p]!=manifest[p]:io_changes[p]=dict(before=manifest[p],after=sources[p])
    write(dest/'source_hashes.json', sources)
    write(dest/'source_io_compatibility.json',dict(changes=io_changes,
          reason='Preserve all packed high/low DD coordinates through the merged HDF5 production restart interface. Physics and solver equations unchanged.'))
    provenance = {}
    frozen = read(old.O/'export_evidence.json')['input_hashes']
    frozen.update(read(old.O/'native_evidence.json')['input_hashes'])

    def record(path, verify_native=False):
        path = Path(path)
        rel = path.relative_to(R).as_posix()
        actual = sha(path)
        if verify_native:
            assert frozen[rel] == actual, rel
        provenance[rel] = actual

    high = read(old.q.O/'inputs.json')['cases']
    low = old.a.rows(old.ENORMAL/'continuation_attempts.csv')
    points = old.a.rows(old.O/'native_points.csv')
    assert len(points) == 204 and all(p['native_qualified'] == 'True' for p in points)
    record(old.O/'native_points.csv', True)
    shutil.copyfile(old.O/'native_points.csv', dest/'native_points.csv')
    cases = []
    for cc in read(old.O/'vela_contract.json')['cases']:
        root = dest/cc['case']
        root.mkdir()
        cfg = read(cc['template'])
        record(cc['template'])
        for key in ('mesh_file', 'node_doping_file', 'materials_file'):
            src = Path(cfg[key])
            record(src)
            target = root/(key + src.suffix)
            if key == 'materials_file':
                material = read(src)
                assert material['schema'] == 'vela.transportmodels.sentaurus2022.materials.v1'
                # Current main rejects unknown schema tags. Use its supported
                # legacy object loader with every material entry unchanged.
                schema = material.pop('schema')
                write(target, material)
                assert read(target)['materials'] == read(src)['materials']
                write(root/'material_format_migration.json', dict(
                    original_schema=schema, original_sha256=sha(src),
                    migrated_sha256=sha(target), changed_keys=['schema'],
                    all_material_entries_identical=True,
                    unit_scaling_unchanged=cfg['scaling']['mode']))
            else:
                shutil.copyfile(src, target)
            cfg[key] = target.relative_to(dest).as_posix()
        write(root/'template.json', cfg)
        geo, mask = old.q.V.m.previous.prior.support(cc)
        all_si = old.d.matrix.spatial.old.m78.supports(cc['device'], geo, .05)[0]['all_si']
        write(root/'geometry.json', dict(count=geo.count,
              free_si=[i for i, keep in enumerate(mask) if keep],
              all_si=[i for i, keep in enumerate(all_si) if keep],
              contacts=list(map(int, geo.contact_nodes)),
              barycentric_si=list(map(float, geo.volumes['barycentric_si'])),
              all_cell=list(map(float, geo.volumes['all_cell']))))
        assert sum(mask) == 907
        original_cfg = read(cc['template'])
        doping = {int(r['node_id']): float(r['donors_cm3']) + float(r['acceptors_cm3'])
                  for r in old.a.rows(Path(original_cfg['node_doping_file']))}
        ni = next(m['ni'] for m in read(original_cfg['materials_file'])['materials'] if m['name'] == 'Si') * 1e6
        for index in range(51):
            src = old.L/'native_exports'/cc['case']/f'vg_{index:03d}'
            # RAW_READER uses geometry mapping plus all-region native fields.
            for f in src.rglob('*'):
                if f.is_file():
                    record(f, True)
            psi, en, hp, spread = old.c.bgn.mapping.RAW_READER(src, geo)
            assert spread < 1e-12
            scalar = old.c.bgn.mapping.original.m73.scalar
            native = {name: scalar(src/'fields'/(name+'_region0.csv')) for name in NAMES}
            assert all(set(values) == set(map(int, __import__('numpy').flatnonzero(all_si))) for values in native.values())
            fn, fp = native['eQuasiFermiPotential'], native['hQuasiFermiPotential']
            for i in fn:
                effective = ni * math.exp(old.c.bgn.delta_eg(doping[i], 'old_slotboom')/(2*old.c.bgn.VT))
                en[i] = effective * math.exp((psi[i]-fn[i])/old.c.bgn.VT)
                hp[i] = effective * math.exp((fp[i]-psi[i])/old.c.bgn.VT)
            target = root/'native'/f'vg_{index:03d}'
            target.mkdir(parents=True)
            with (target/'state.csv').open('w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=['node_id','psi','phin','phip','electrons_m3','holes_m3'])
                writer.writeheader()
                writer.writerows(dict(node_id=i, psi=psi[i], phin=fn.get(i,0.), phip=fp.get(i,0.),
                                      electrons_m3=en[i], holes_m3=hp[i]) for i in range(geo.count))
            write(target/'fields.json', native)
        for index in (0,10,40,50):
            if index >= 40:
                src = Path(next(r['enormal_seed'] for r in high if r['case']==cc['case'] and r['index']==index))
            else:
                src = Path(next(r['dest'] for r in low if r['case']==cc['case'] and int(r['index'])==index and r['qualified']=='True'))/'state.csv'
            record(src)
            shutil.copyfile(src, root/f'enormal_{index:03d}.csv')
        cases.append({k:cc[k] for k in ('case','device','vd')})
        print('Packaged', cc['case'], flush=True)
    # Historical acceptance replay fixture: numerical gate equivalence, not a new result.
    historical = next(r for r in old.a.rows(old.O/'low_control_attempts.csv') if r['qualified']=='True')
    fixture = dest/'acceptance_fixture'
    fixture.mkdir()
    for name in ('config.status.json','all_row.status.json','acceptance_edges.status.json',
                 'functional.status.json','all_row.csv','acceptance_edges.csv','result.json'):
        src=Path(historical['dest'])/name
        record(src)
        shutil.copyfile(src, fixture/name)
    write(fixture/'case.json', dict(case=historical['case']))
    write(dest/'provenance.json', provenance)
    write(dest/'contract.json', dict(cases=cases, gates=old.v.GATES, VT=old.c.bgn.VT,
          git_head='06d0acb425b052f0df424341ab86a5994ed45bf5', source_includes_working_changes=True,
          controls=[0,10,40,50], curve_indices=list(range(51)), source_volume='all_cell',
          recovery='One same-bias reload; retain every attempt; only qualified states propagate.',
          scope='Merged source HFS baseline, four 0..1 V curves. Controls must pass first.',
          current_metric='100*(Id_Vela/Id_Sentaurus-1); descriptive 1% comparison separate from numerical gates.',
          source_unchanged_since_merge_validation=not io_changes,
          io_compatibility_patch='source_io_compatibility.json', acceptance_changed=False))
    write(dest/'manifest.json', {p.relative_to(dest).as_posix():sha(p)
          for p in sorted(dest.rglob('*')) if p.is_file()})
    with tarfile.open(dest.with_suffix('.tgz'), 'w:gz') as archive:
        archive.add(dest, arcname='inputs')
    print(dest.with_suffix('.tgz'), sha(dest.with_suffix('.tgz')), flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('destination', type=Path)
    main(parser.parse_args().destination.resolve())
