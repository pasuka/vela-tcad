#!/usr/bin/env python3
"""Materialize WP14 BVmethods sources without changing the declared models.

The ablation template declares four SRH/Enormal branches; the other two decks
are copied byte-for-byte. This prepares local inputs only, never launches a VM.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
from check_sdevice_gate_contact_mapping import (
    build_mesh, intrinsic_work_function_eV, predicted_flatband_V,
    sentaurus_gate_reference_eV,
)

REPO = Path(__file__).resolve().parents[1]
SOURCE = REPO / 'reference_tcad/bvmethods_sentaurus2018/source'
VARIANTS = {
    'a_constant_no_enormal': (False, False),
    'b_doping_no_enormal': (True, False),
    'c_constant_enormal': (False, True),
    'd_doping_enormal': (True, True),
}


def numerical_precision_controls(text: str, name: str) -> str:
    """Apply the independently validated D precision repair; other decks stay exact."""
    if name == 'd_doping_enormal':
        return text.replace('Iterations=20', 'Digits=8\n  Iterations=60')
    raise ValueError('only the D precision repair has independently passed')


def materialize_bvmethods(output: Path, *, numerical_precision: bool = False) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    template = (SOURCE / 'full_physics_ablation_sdevice.cmd').read_text(encoding='utf-8')
    inputs = {}
    for name in ('models.par', 'full_physics_constant_srh.par',
                 'aba_poisson_sdevice.cmd', 'external_resistor_sdevice.cmd'):
        shutil.copyfile(SOURCE / name, output / name)
        inputs[name] = hashlib.sha256((SOURCE / name).read_bytes()).hexdigest()
    for name, (doping_srh, enormal) in VARIANTS.items():
        text = template.replace('SRH_PARAMETER_FILE',
                                'models.par' if doping_srh else 'full_physics_constant_srh.par')
        text = text.replace('SRH_MODEL', 'DopingDep' if doping_srh else '')
        text = text.replace('ENORMAL_MODEL', 'Enormal' if enormal else '')
        # Empty SRH parentheses are removed; model selection is otherwise exact.
        text = text.replace('SRH()', 'SRH')
        text = text.replace('full_physics_ablation_des', name + '_des')
        (output / (name + '_sdevice.cmd')).write_text(text, encoding='utf-8', newline='\n')
    inputs['full_physics_ablation_sdevice.cmd'] = hashlib.sha256(
        (SOURCE / 'full_physics_ablation_sdevice.cmd').read_bytes()).hexdigest()
    manifest = {'schema': 'vela.wp14.materialized_sources.v1',
                'template_sha256': inputs['full_physics_ablation_sdevice.cmd'],
                'source_inputs': inputs, 'variants': {k: {'srh_doping_dep': v[0], 'enormal': v[1]}
                                                      for k, v in VARIANTS.items()},
                'mesh_origin': 'bvmethods_sentaurus2018/v1/mesh/n1_msh.tdr',
                'mesh_alias': 'bvmethods_nmos_msh.tdr',
                'parameter_policy': 'preserve repository SRH overrides; archive native default dump separately'}
    if numerical_precision:
        for name in ('d_doping_enormal',):
            path = output / (name + '_sdevice.cmd')
            path.write_text(numerical_precision_controls(path.read_text(), name), encoding='utf-8', newline='\n')
        manifest['numerical_precision_controls'] = {
            name: hashlib.sha256((output / (name + '_sdevice.cmd')).read_bytes()).hexdigest()
            for name in ('d_doping_enormal',)}
    (output / 'wp14_inputs.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8', newline='\n')
    return manifest


def materialize_gate(output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(REPO / 'share/vela/materials/models.par', output / 'models.par')
    # Geometry and doping match the WP4a fixture. The native grid is refined
    # independently; all five variants share this one grid.
    sde = '''(sde:clear)
(sdegeo:create-rectangle (position 0 0 0) (position 0.05 1 0) "Silicon" "R.Si")
(sdegeo:create-rectangle (position 0 1 0) (position 0.05 1.01 0) "Oxide" "R.Ox")
(sdegeo:define-contact-set "substrate" 4 (color:rgb 1 0 0) "##")
(sdegeo:define-contact-set "gate" 4 (color:rgb 0 0 1) "##")
(sdegeo:define-2d-contact (find-edge-id (position 0.025 0 0)) "substrate")
(sdegeo:define-2d-contact (find-edge-id (position 0.025 1.01 0)) "gate")
(sdedr:define-constant-profile "ND" "PhosphorusActiveConcentration" 1e17)
(sdedr:define-constant-profile-region "ND.place" "ND" "R.Si")
(sdedr:define-refinement-size "Si.mesh" 0.025 0.02 0.025 0.005)
(sdedr:define-refinement-region "Si.place" "Si.mesh" "R.Si")
(sdedr:define-refinement-size "Ox.mesh" 0.025 0.005 0.025 0.005)
(sdedr:define-refinement-region "Ox.place" "Ox.mesh" "R.Ox")
(sde:build-mesh "snmesh" "" "gate")
'''
    (output / 'gate_sde.cmd').write_text(sde, encoding='utf-8', newline='\n')
    branches = {'barrier_zero': ('Barrier=0', 0.0),
                'barrier_positive': ('Barrier=0.55', 0.55),
                'barrier_negative': ('Barrier=-0.3', -0.3),
                'workfunction': (f'WorkFunction={intrinsic_work_function_eV():.17g}',
                                 intrinsic_work_function_eV() - sentaurus_gate_reference_eV()),
                'workfunction_zero': (f'WorkFunction={sentaurus_gate_reference_eV():.17g}', 0.0)}
    manifest = {'schema': 'vela.wp14.gate_reference_sources.v1',
                'geometry_um': {'width': 0.05, 'silicon_thickness': 1.0, 'oxide_thickness': 0.01},
                'donors_cm3': 1e17, 'temperature_K': 300.0,
                'reference_intrinsic_work_function_eV': intrinsic_work_function_eV(),
                'sentaurus_gate_reference_eV': sentaurus_gate_reference_eV(), 'variants': {}}
    for name, (attribute, barrier) in branches.items():
        expected = predicted_flatband_V(barrier, 1e17, 0.0)
        low, high = expected - 0.02, expected + 0.02
        text = f'''File {{
 Grid="gate_msh.tdr" Parameter="models.par"
 Plot="{name}_des.tdr" Current="{name}.plt" Output="{name}.log"
}}
Electrode {{
 {{ Name="substrate" Voltage=0 }}
 {{ Name="gate" Voltage=0 {attribute} }}
}}
Physics {{ Temperature=300 }}
Plot {{ Potential eDensity hDensity eQuasiFermi hQuasiFermi ElectricField/Vector Doping }}
Math {{ RelErrControl Digits=8 Iterations=100 NotDamped=100 Method=Blocked }}
Solve {{
 Coupled(Iterations=100) {{ Poisson }}
 Coupled {{ Poisson Electron Hole }}
 Quasistationary(InitialStep=0.01 MaxStep=0.1 MinStep=1e-8
   Goal {{ Name="gate" Voltage={low:.17g} }}) {{ Coupled {{ Poisson Electron Hole }} }}
 Quasistationary(InitialStep=0.025 MaxStep=0.025 MinStep=1e-8
   Goal {{ Name="gate" Voltage={high:.17g} }}) {{
   Coupled {{ Poisson Electron Hole }}
   Plot(FilePrefix="{name}_scan" Time=(Range=(0 1) Intervals=40) NoOverWrite)
 }}
}}
'''
        (output / (name + '_sdevice.cmd')).write_text(text, encoding='utf-8', newline='\n')
        manifest['variants'][name] = {'attribute': attribute, 'expected_flatband_V': expected,
                                      'scan_biases_V': [low, high], 'samples': 41}
    # Keep the Vela fixture as an inspectable geometric/doping reference.
    (output / 'vela_fixture.json').write_text(json.dumps(build_mesh(), indent=2)+'\n', encoding='utf-8', newline='\n')
    (output / 'wp14_gate_inputs.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8', newline='\n')
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bvmethods-output', type=Path)
    parser.add_argument('--gate-output', type=Path)
    parser.add_argument('--numerical-precision', action='store_true',
                        help='Apply validated D precision controls; keep external resistor exact and unqualified')
    args = parser.parse_args()
    if not (args.bvmethods_output or args.gate_output):
        parser.error('at least one output directory is required')
    if args.bvmethods_output: materialize_bvmethods(args.bvmethods_output, numerical_precision=args.numerical_precision)
    if args.gate_output: materialize_gate(args.gate_output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
