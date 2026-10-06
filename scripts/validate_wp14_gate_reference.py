#!/usr/bin/env python3
"""Compare the five WP14 native MOS-capacitor controls with emitted Vela contacts.

Raw native files and converted fields stay in the supplied ignored directories.
The gate is fixed at 1 mV; the independently refined meshes share geometry,
doping and material parameters. Native Plot Time is affine within scan stage 2.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import statistics
import subprocess
from pathlib import Path

import check_sdevice_gate_contact_mapping as gate
from vela_sdevice.cmd_parser import parse_cmd_text
from vela_sdevice.deck_emit import EmitOptions, emit_stage_plan
from vela_sdevice.materials import active_models_from_ir, emit_materials
from vela_sdevice.mesh_import import import_grid
from vela_sdevice.plt.document import read_plt_document


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(archive: Path, output: Path, importer: Path, runner: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    inputs_path = archive/'source/wp14_gate_inputs.json'
    inputs = json.loads(inputs_path.read_text())
    if set(inputs['variants']) != {'barrier_zero', 'barrier_positive', 'barrier_negative',
                                  'workfunction', 'workfunction_zero'}:
        raise ValueError('the five native controls are required')
    qualification = import_grid(archive/'mesh/gate_msh.tdr', output/'native_mesh', importer=str(importer))
    fixture = gate.build_mesh()
    work = output/'vela_fixture'
    work.mkdir(exist_ok=True)
    (work/'mesh.json').write_text(json.dumps(fixture['mesh']))
    results = {}
    for name, variant in inputs['variants'].items():
        cmd = archive/'source'/f'{name}_sdevice.cmd'
        ir = parse_cmd_text(cmd.read_text(), path_name=cmd.name)
        materials = emit_materials(archive/'source/models.par', active_models=active_models_from_ir(ir),
            temperature_K=300., mesh_materials=('Si', 'SiO2'), allow_lossy=True,
            effective_models=ir['physics']['effective'])
        emitted = emit_stage_plan(ir, mesh_qualification=qualification, materials_report=materials.report,
            work_dir='gate_vela', cmd_basename=cmd.name, cmd_sha256=ir['source']['sha256'],
            options=EmitOptions(allow_lossy=True))
        contact = next(c for c in emitted.decks[-1]['contacts'] if c['name'] == 'gate')
        phi = contact['work_function_eV']
        _, current_rows = read_plt_document((archive/'plt'/f'{name}.plt').read_text())
        low, high = variant['scan_biases_V']
        if abs(current_rows[-1]['gate']['outer_voltage_V']-high) > 1e-9 or not any(
                abs(row['gate']['outer_voltage_V']-low) <= 1e-9 for row in current_rows):
            raise ValueError(f'{name}: native scan endpoints are missing')
        samples = {}
        fingerprints = []
        for tdr in sorted((archive/'tdr').glob(f'{name}_scan_*_des.tdr')):
            export = output/'native_fields'/tdr.stem
            export.mkdir(parents=True, exist_ok=True)
            with (export/'import.log').open('w') as log:
                subprocess.run([str(importer), '--tdr', str(tdr), '--export-dir', str(export)],
                               stdout=log, stderr=subprocess.STDOUT, check=True)
            nodes = {int(r['id']): float(r['y_um']) for r in csv.DictReader((export/'nodes.csv').open())}
            fields = json.loads((export/'field_manifest.json').read_text())['fields']
            potential, = [f for f in fields if f['name'] == 'ElectrostaticPotential' and f['region_name'] == 'R.Si']
            if potential['unit'] != 'V':
                raise ValueError('native potential is not in V')
            values = {int(r['node_id']): float(r['component0']) for r in csv.DictReader(
                (export/'fields'/potential['csv_file']).open())}
            surface = [v for i, v in values.items() if abs(nodes[i]-1.) < 1e-10]
            bulk = [v for i, v in values.items() if abs(nodes[i]) < 1e-10]
            index = int(tdr.stem.split('_')[-2])
            if not surface or not bulk or index not in range(41) or index in samples:
                raise ValueError(f'{name}: malformed native sample {index}')
            samples[index] = (low+(high-low)*index/40., statistics.mean(surface)-statistics.mean(bulk))
            fingerprints.append({'file': tdr.name, 'sha256': sha(tdr)})
        if set(samples) != set(range(41)):
            raise ValueError(f'{name}: missing native plot samples')
        native = gate.flatband_from_samples(dict(samples.values()))
        candidate = gate.check_barrier(runner, work, phi, fixture['probes'], 0., 1e17)
        error = candidate['measured_flatband_V']-native
        results[name] = {'status': 'pass' if abs(error) <= .001 else 'fail',
            'emitted_phi_ms_V': phi, 'native_flatband_V': native,
            'vela_flatband_V': candidate['measured_flatband_V'], 'difference_V': error,
            'source_sha256': sha(cmd), 'native_tdr_sha256': fingerprints, 'vela_extraction': candidate}
    return {'schema': 'vela.wp14.gate_native_comparison.v1', 'limit_V': .001,
        'status': 'pass' if all(r['status'] == 'pass' for r in results.values()) else 'fail',
        'inputs_sha256': sha(inputs_path), 'native_mesh': qualification,
        'native_mesh_sha256': sha(archive/'mesh/gate_msh.tdr'), 'controls': results}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--importer', type=Path, required=True)
    parser.add_argument('--runner', type=Path, required=True)
    args = parser.parse_args()
    result = validate(args.archive.absolute(), args.output.absolute(), args.importer.absolute(), args.runner.absolute())
    (args.output/'gate_native_comparison.json').write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], {name: r['difference_V'] for name, r in result['controls'].items()})
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
