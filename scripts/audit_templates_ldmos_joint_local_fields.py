"""Score qualified D0 curves against the user-approved local-field contract.

Native TDR region ordering is verified against an existing importer export and
exact non-state geometry hashes. This does not perform Newton or alter states.
"""
import argparse
import csv
import hashlib
import json
import os
from electrothermal_state import read_bound_record
from pathlib import Path
import subprocess

import h5py
import numpy as np
from extract_templates_ldmos_d0_temperature import geometry_hash, values
from evidence_paths import candidate_path


def read(path):
    return read_bound_record(path)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def carrier_metric(reference, candidate, areas, selected, floor_fraction):
    reference, candidate, areas = map(np.asarray, (reference, candidate, areas))
    selected = np.asarray(selected, dtype=bool)
    if not np.isfinite(floor_fraction) or not 0 <= floor_fraction < 1:
        raise ValueError('Density floor fraction must be finite and in [0, 1)')
    if not (reference.shape == candidate.shape == areas.shape == selected.shape):
        raise ValueError('Local-field arrays must have identical shapes')
    if not all(np.all(np.isfinite(a)) for a in (reference, candidate, areas)):
        raise ValueError('Nonfinite local fields')
    if np.any(reference <= 0) or np.any(candidate <= 0) or np.any(areas < 0):
        raise ValueError('Positive carrier densities and nonnegative areas required')
    floor = reference.max()*floor_fraction
    active = selected & (reference > floor) & (areas > 0)
    if not np.any(active):
        raise ValueError('No resolved positive-volume nodes in local-field group')
    relative = (candidate[active]-reference[active])/reference[active]
    rms = np.sqrt(np.sum(areas[active]*relative**2)/np.sum(areas[active]))
    return dict(selected_nodes=int(selected.sum()), resolved_nodes=int(active.sum()),
                reference_floor_m3=float(floor), weighted_relative_rms=float(rms),
                maximum_absolute_relative_error=float(np.max(np.abs(relative))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--curves', type=Path, required=True)
    parser.add_argument('--prefix', default='p3_full_vg')
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--native', type=Path, default=Path('reference_staging/templates_ldmos_d0_electrothermal_20260912/native_fields_r1'))
    parser.add_argument('--native-raw', type=Path, default=Path('reference_staging/templates_ldmos_d0_electrothermal_20260912/native_full_r1/raw'))
    parser.add_argument('--contract', type=Path, default=Path('reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_joint_acceptance_20260914_v2.json'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--indices', type=int, nargs='+', default=list(range(31)))
    parser.add_argument('--candidate-path-map', nargs=2, metavar=('SOURCE_ROOT','COPIED_ROOT'))
    args = parser.parse_args()
    if any(i not in range(31) for i in args.indices) or len(set(args.indices)) != len(args.indices):
        raise ValueError('Indices must be distinct exact reference point indices')
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = read(args.native/'manifest.json')
    paths = {Path(p): h for p, h in manifest['sources_sha256'].items()}
    baseline = next(p for p in paths if p.suffix == '.tdr')
    map_file = next(p for p in paths if p.name == 'LatticeTemperature_region0.csv')
    assert sha(baseline) == paths[baseline] and sha(map_file) == paths[map_file]
    with map_file.open(encoding='utf-8') as f:
        mapped = list(csv.DictReader(f))
    node_ids = np.array([int(r['node_id']) for r in mapped])
    assert len(node_ids) == len(set(node_ids))
    with h5py.File(baseline, 'r') as f:
        assert geometry_hash(f) == manifest['geometry_sha256']
        assert np.array_equal(values(f, 'LatticeTemperature', 0), [float(r['component0']) for r in mapped])
    contract = read(args.contract)['local_fields']
    inputs = {str(p.resolve()): sha(p) for p in (args.probe, args.contract, args.native/'manifest.json', baseline, map_file, Path(__file__))}
    report = dict(status='running', scope=__doc__, points=[], sources_sha256=inputs,
                  candidate_path_map=args.candidate_path_map,
                  contract=contract, all_62_points=(sorted(args.indices) == list(range(31))))
    def save():
        (args.output/'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    env = dict(os.environ)
    if os.name == 'nt':
        env['PATH'] = 'D:/msys64/ucrt64/bin;'+env['PATH']
    for gate in (4, 8):
        ledger_path = args.curves/f'{args.prefix}{gate}'/'ledger.json'
        ledger = read(ledger_path)
        assert ledger['status'] == 'complete' and len(ledger['exact_points']) == 31
        inputs[str(ledger_path.resolve())] = sha(ledger_path)
        zero_input = candidate_path(ledger['runs'][0]['directory'], args.candidate_path_map)/'input.json'
        cfg = read(zero_input)
        inputs[str(zero_input.resolve())] = sha(zero_input)
        mesh_path = candidate_path(cfg['mesh_file'], args.candidate_path_map)
        mesh = read(mesh_path)
        inputs[str(mesh_path.resolve())] = sha(mesh_path)
        silicon, oxide = set(), set()
        for cell in mesh['triangles']:
            (silicon if cell['region_id'] == 0 else oxide).update(cell['node_ids'])
        assert set(node_ids) == silicon
        interface = np.array([i in silicon & oxide for i in node_ids])
        areas = np.array(cfg['silicon_area_m2'])[node_ids]
        # The manifest baseline may live outside the native full-run directory.
        candidates = list(args.native_raw.glob(f'field_vg{gate}_*_des.tdr'))
        assert len(candidates) == 31
        candidates.sort()
        for index in args.indices:
            field = next(x for x in manifest['fields'] if x['gate_V'] == gate and x['point_index'] == index)
            native_path = candidates[index]
            assert sha(native_path) == field['tdr_sha256']
            with h5py.File(native_path, 'r') as f:
                assert geometry_hash(f) == manifest['geometry_sha256']
                native = {k: values(f, name, 0)*factor for k, name, factor in (
                    ('electrons_m3','eDensity',1e6), ('holes_m3','hDensity',1e6),
                    ('conduction_band_eV','ConductionBandEnergy',1.), ('valence_band_eV','ValenceBandEnergy',1.))}
            point = ledger['exact_points'][index]
            assert abs(point['bias_V']-index*40/30) < 1e-9
            result_path = candidate_path(point['result'], args.candidate_path_map)
            result = read(result_path)
            x = result['referenced_state_interleaved']
            origin = result.get('potential_origin_V', 0.)
            states = [dict(id=int(i), potential_V=x[4*i]+origin,
                electron_qf_V=x[4*i+1], hole_qf_V=x[4*i+2],
                electron_qf_reference_V=result['electron_qf_reference_V'][i]+origin,
                hole_qf_reference_V=result['hole_qf_reference_V'][i]+origin,
                temperature_K=x[4*i+3], donors_m3=cfg['donors_m3'][i], acceptors_m3=cfg['acceptors_m3'][i]) for i in node_ids]
            directory = args.output/f'vg{gate}_{index:02d}'
            directory.mkdir()
            inp, out = directory/'input.json', directory/'output.json'
            inp.write_text(json.dumps(dict(states_SI=states, auger_with_generation=cfg['auger_with_generation'])), encoding='utf-8')
            with out.open('w', encoding='utf-8') as stream:
                subprocess.run([str(args.probe.resolve()), str(inp.resolve())], stdout=stream, env=env, check=True)
            local = read(out)['results']
            assert [r['id'] for r in local] == list(node_ids)
            row = dict(gate_V=gate, point_index=index, bias_V=point['bias_V'], carriers={}, bands={}, pass_gate=True)
            for key in contract['carriers']:
                candidate = np.array([r[key]['value'] for r in local])
                row['carriers'][key] = {}
                for group, selection in (('all_silicon', np.ones(len(node_ids),dtype=bool)), ('silicon_oxide_interface',interface)):
                    metric = carrier_metric(native[key], candidate, areas, selection, contract['reference_density_floor_fraction_of_carrier_peak'])
                    metric['pass'] = metric['weighted_relative_rms'] <= contract['volume_weighted_relative_rms_limit']
                    row['pass_gate'] &= metric['pass']
                    row['carriers'][key][group] = metric
            for key in contract['band_edges']:
                candidate = np.array([r[key]['value'] for r in local])
                delta = np.abs(candidate-native[key]); worst=int(np.argmax(delta))
                metric = dict(maximum_absolute_error_eV=float(delta[worst]), node_id=int(node_ids[worst]),
                              pass_gate=bool(delta[worst] <= contract['maximum_absolute_band_edge_error_eV']))
                row['pass_gate'] &= metric['pass_gate']; row['bands'][key] = metric
            report['points'].append(row)
            inputs[str(native_path.resolve())] = sha(native_path)
            inputs[str(result_path.resolve())] = sha(result_path)
            inputs[str(out.resolve())] = sha(out)
            save()
            print(json.dumps(dict(gate=gate,index=index,pass_gate=row['pass_gate'])),flush=True)
    report['status'] = 'pass' if all(p['pass_gate'] for p in report['points']) else 'fail'
    save()


if __name__ == '__main__':
    main()
