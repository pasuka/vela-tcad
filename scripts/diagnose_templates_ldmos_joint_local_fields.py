"""Fixed-state counterfactual diagnosis of failed local-field gates.

Uses the complete audit's verified native/node mapping. Replacing one state
component is diagnostic only: it neither changes the solution nor qualifies a
modified model. Counterfactual reductions are not additive causal shares.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess

import h5py
import numpy as np
from audit_templates_ldmos_joint_local_fields import read, sha, carrier_metric
from extract_templates_ldmos_d0_temperature import values


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audit', type=Path, required=True)
    p.add_argument('--curves', type=Path, required=True)
    p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    audit = read(a.audit/'summary.json')
    # Recheck native, mapping, state and output hashes; the scorer source can
    # legitimately have acquired new CLI options since this immutable run.
    for name, digest in audit['sources_sha256'].items():
        if Path(name).suffix != '.py' and sha(Path(name)) != digest:
            raise ValueError(f'Audit source changed: {name}')
    if sha(a.probe) != audit['sources_sha256'][str(a.probe.resolve())]:
        raise ValueError('Use the same local physics executable')
    a.output.mkdir(parents=True, exist_ok=False)
    contract = audit['contract']
    sources = {str(a.audit/'summary.json'): sha(a.audit/'summary.json'),
               str(Path(__file__).resolve()): sha(Path(__file__))}
    report = dict(scope=__doc__, points=[], sources_sha256=sources)
    env = dict(os.environ)
    if os.name == 'nt':
        env['PATH'] = 'D:/msys64/ucrt64/bin;'+env['PATH']
    for point in audit['points']:
        if point['pass_gate']:
            continue
        gate, index = point['gate_V'], point['point_index']
        directory = a.output/f'vg{gate}_{index:02d}'
        directory.mkdir()
        ledger = read(a.curves/f'p3_full_vg{gate}'/'ledger.json')
        cfg = read(Path(ledger['runs'][0]['directory'])/'input.json')
        mesh = read(Path(cfg['mesh_file']))
        silicon, oxide = set(), set()
        for cell in mesh['triangles']:
            (silicon if cell['region_id'] == 0 else oxide).update(cell['node_ids'])
        original = read(a.audit/f'vg{gate}_{index:02d}'/'input.json')
        states = original['states_SI']
        ids = np.array([r['id'] for r in states])
        native_paths = sorted(Path(name) for name in audit['sources_sha256']
                              if Path(name).name.startswith(f'field_vg{gate}_') and name.endswith('.tdr'))
        if len(native_paths) != 31:
            raise ValueError('Diagnosis requires a complete 62-point audit')
        with h5py.File(native_paths[index], 'r') as f:
            native = {key: values(f, field, 0) for key, field in (
                ('potential_V', 'ElectrostaticPotential'),
                ('electron_qf_V', 'eQuasiFermiPotential'),
                ('hole_qf_V', 'hQuasiFermiPotential'),
                ('temperature_K', 'LatticeTemperature'),
                ('electrons_m3', 'eDensity'), ('holes_m3', 'hDensity'))}
        areas = np.array(cfg['silicon_area_m2'])[ids]
        interface = np.array([n in silicon & oxide for n in ids])
        candidates = {'full': read(a.audit/f'vg{gate}_{index:02d}'/'output.json')['results']}
        groups = {'native_temperature': ('temperature_K',),
                  'native_potential': ('potential_V',),
                  'native_quasi_fermi': ('electron_qf_V', 'hole_qf_V'),
                  'all_native_state': ('potential_V', 'electron_qf_V', 'hole_qf_V', 'temperature_K')}
        for mode, components in groups.items():
            modified = [dict(row) for row in states]
            for j, row in enumerate(modified):
                for key in components:
                    row[key] = float(native[key][j])
                    if key in ('electron_qf_V', 'hole_qf_V'):
                        row[key.replace('_V', '_reference_V')] = 0.
            inp, out = directory/f'{mode}_input.json', directory/f'{mode}_output.json'
            inp.write_text(json.dumps(dict(states_SI=modified, auger_with_generation=False)), encoding='utf-8')
            with out.open('w', encoding='utf-8') as stream:
                subprocess.run([str(a.probe.resolve()), str(inp.resolve())], stdout=stream, env=env, check=True)
            candidates[mode] = read(out)['results']
            sources[str(out)] = sha(out)
        entry = dict(gate_V=gate, point_index=index, bias_V=point['bias_V'], carriers={})
        for carrier in contract['carriers']:
            ref = native[carrier]*1e6
            entry['carriers'][carrier] = {}
            rows = []
            for mode, results in candidates.items():
                candidate = np.array([r[carrier]['value'] for r in results])
                entry['carriers'][carrier][mode] = carrier_metric(
                    ref, candidate, areas, interface, contract['reference_density_floor_fraction_of_carrier_peak'])
            full = np.array([r[carrier]['value'] for r in candidates['full']])
            active = interface & (ref > ref.max()*contract['reference_density_floor_fraction_of_carrier_peak']) & (areas > 0)
            contributions = areas*((full-ref)/ref)**2
            total = contributions[active].sum()
            for j in np.where(active)[0]:
                node = mesh['nodes'][int(ids[j])]
                rows.append(dict(node_id=int(ids[j]), coordinates=json.dumps(node), reference_m3=ref[j],
                    relative_error=(full[j]-ref[j])/ref[j], weighted_squared_error_share=contributions[j]/total if total else 0.,
                    temperature_difference_K=states[j]['temperature_K']-native['temperature_K'][j],
                    potential_difference_V=states[j]['potential_V']-native['potential_V'][j],
                    electron_qf_difference_V=states[j]['electron_qf_V']+states[j]['electron_qf_reference_V']-native['electron_qf_V'][j],
                    hole_qf_difference_V=states[j]['hole_qf_V']+states[j]['hole_qf_reference_V']-native['hole_qf_V'][j]))
            rows.sort(key=lambda r: r['weighted_squared_error_share'], reverse=True)
            with (directory/f'{carrier}_interface.csv').open('w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader(); writer.writerows(rows)
        report['points'].append(entry)
        (a.output/'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print(json.dumps(dict(gate=gate, index=index, completed=True)), flush=True)


if __name__ == '__main__':
    main()
