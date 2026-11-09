"""Freeze the qualified 2026-11-09 BVMethods text snapshot; never replace v1-v4."""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile

VERSION = 'u2022_12_20261109'
INTEGRATION_COMMIT = '77a889a95b87f6a17c7a5943c135a0b9d04cf6b0'
EVIDENCE_SHA256 = '6c6ae6e0ecb8c2467413c911e26a699e21fc390c2364c547627694d8ae1cdd3e'
NATIVE_ACCEPTED_SHA256 = '867fd250fb9522ee4a4c910b80c3cd08cdc5ffa37d7b93d46b85faeb2689339c'
NATIVE_SUMMARY_SHA256 = '5ecb51d94be5478acf496b58091ccbbf6c9adf8cf8b0bed2f22b11aa33a27c07'
CMD_SHA256 = '47ed4a2330f20bcb190535175aa5dc6fb180cf49315093659c79d8e450c2a518'
sha = lambda data: hashlib.sha256(data).hexdigest()

def encode(value):
    return (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n').encode('utf-8')

def csv_bytes(rows, columns):
    out = io.StringIO(newline='')
    writer = csv.DictWriter(out, fieldnames=columns, lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    return out.getvalue().encode('utf-8')

def freeze(evidence: Path, native_summary: Path, native_accepted: Path, output: Path):
    if output.exists():
        raise ValueError('version output already exists; archived snapshots are immutable')
    if output.name != VERSION:
        raise ValueError('this freezer writes only the named qualified snapshot')
    if sha(evidence.read_bytes()) != EVIDENCE_SHA256:
        raise ValueError('unqualified Vela evidence archive')
    if sha(native_summary.read_bytes()) != NATIVE_SUMMARY_SHA256:
        raise ValueError('unqualified native summary')
    native = json.loads(native_summary.read_text(encoding='utf-8'))
    native_bytes = native_accepted.read_bytes()
    if native['version'] != 'U-2022.12' or sha(native_bytes) != NATIVE_ACCEPTED_SHA256 or sha(native_bytes) != native['analysis_hashes']['accepted_iv.csv']:
        raise ValueError('native version or accepted-curve SHA mismatch')
    raw_native_rows = list(csv.DictReader(io.StringIO(native_bytes.decode('utf-8-sig'))))
    # The archived PLT includes three initialization records at 0 V. The
    # qualified scan ledger retains the last state per outer bias.
    native_rows = list({float(r['drain OuterVoltage']): r for r in raw_native_rows}.values())
    if len(native_rows) != native['accepted'] or not native['reached_100V']:
        raise ValueError('incomplete native bounded scan')
    files = {'.gitattributes': b'* -text\nsource/*.cmd -whitespace\nsource/*.par -whitespace\nsource/*.json whitespace=cr-at-eol\n'}
    candidates = {}
    with zipfile.ZipFile(evidence) as z:
        members = json.loads(z.read('terminal_SHA256.json'))
        if set(z.namelist()) != set(members) | {'terminal_SHA256.json'} or not all(sha(z.read(n)) == h for n, h in members.items()):
            raise ValueError('Vela evidence member SHA mismatch')
        identity = json.loads(z.read('final_identity.json'))
        if not identity['no_owned'] or not identity['source_unchanged']:
            raise ValueError('terminal identity qualification missing')
        prefix = 'cmd_original_r3/'
        for name in ('pp5_des.cmd', 'pp5_des.par', 'pp5_des.vela.json'):
            files['source/' + name] = z.read(prefix + name)
        if sha(files['source/pp5_des.cmd']) != CMD_SHA256:
            raise ValueError('original CMD changed')
        for name in ('materials.json', 'simulation_stage00.json', 'simulation_stage01.json'):
            data = z.read(prefix + 'pp5_des_vela/' + name)
            if re.search(rb'/home/|[A-Za-z]:\\', data):
                raise ValueError('nonportable configuration path')
            files['vela/' + name] = data
        external = {}
        for name, member in [('n1_msh.tdr', prefix + 'n1_msh.tdr'),
                             ('mesh.json', prefix + 'pp5_des_vela/mesh.json'),
                             ('doping.csv', prefix + 'pp5_des_vela/doping.csv')]:
            data = z.read(member)
            external[name] = {'sha256': sha(data), 'bytes': len(data), 'checked_in': False}
        files['source/inputs.json'] = encode({'schema': 'vela.bvmethods.external_inputs.v1', 'inputs': external})
        ordinary_status = json.loads(z.read('scan_status.json'))
        for case, base in [('ordinary_full_physical', 'ordinary_full_physical/'),
                           ('ordinary_full_unit', 'ordinary_full_unit/'), ('cmd_full_bv', prefix + 'pp5_des_vela/')]:
            is_cmd = case == 'cmd_full_bv'
            iv = list(csv.DictReader(io.StringIO(z.read(base + ('stage01.csv' if is_cmd else 'iv.csv')).decode('utf-8'))))
            iv = [row for row in iv if row['converged'] == '1']
            qualification = json.loads(z.read('cmd_r3_qualification.json')) if is_cmd else ordinary_status['groups'][case]
            if not qualification['qualified'] or qualification['timed_out']:
                raise ValueError('unqualified candidate case')
            curve = [{'outer_voltage_V': r['outer_voltage_V'], 'inner_voltage_V': r['inner_voltage_V'],
                      'drain_total_current_A_per_um': r['current_total_A_per_um'], 'newton_iterations': r['iterations']} for r in iv]
            curve_path = 'vela/curves/' + case + '.csv'
            files[curve_path] = csv_bytes(curve, list(curve[0]))
            reports = [json.loads(r['validation_diagnostics']) for r in iv]
            candidates[case] = {'role': 'Vela_candidate_not_native_reference', 'curve': curve_path,
                'entry': 'original CMD Poisson/DD + HDF stage chain' if is_cmd else 'formal cold-Poisson JSON entry',
                'accepted': len(iv), 'attempts': None, 'failed': None,
                'successful_newton_iterations': sum(int(r['iterations']) for r in iv),
                'failed_newton_iterations': None, 'seconds': qualification['seconds'],
                'solver_exit': 0, 'timed_out': False, 'stop_current_reached': True,
                'BV_inner_V': qualification['BV_inner_V'] if is_cmd else qualification['BV']['inner_V'],
                'F_pass': sum(d['coupled_device_residual_inf'] <= 1e-8 for d in reports),
                'E_pass': sum(d['density_update_valid'] and 0 <= d['density_update_rms'] < 1 for d in reports),
                'load_KCL_finite_all': all(d['actual_load_line_satisfied'] and d['terminal_kcl_satisfied'] and d['physical_residuals_finite'] for d in reports),
                'carrier_report_failed': qualification['carrier_report_failed'],
                'continuity_report_failed': qualification['continuity_report_failed']}
        runtime = {'core_sha256': identity['binary']['libvela_core.a'],
                   'runner_sha256': identity['binary']['vela_example_runner'],
                   'frontend_archive_sha256': identity['formal_frontend_SHA256'],
                   'MKL_runtime_sha256': identity['final_MKL_runtime_sha256'],
                   'backend': 'pardiso_mkl', 'threads_each': 1}
    ncurve = [{'outer_voltage_V': r['drain OuterVoltage'], 'inner_voltage_V': r['drain InnerVoltage'],
               'drain_total_current_A_per_um': r['drain TotalCurrent']} for r in native_rows]
    native_curve = 'reference_curves/native_u2022_12_adaptive_0_100.csv'
    files[native_curve] = csv_bytes(ncurve, list(ncurve[0]))
    native_public = {'role': 'Sentaurus_native_reference', 'version': 'U-2022.12', 'curve': native_curve,
        'scope_V': [0, 100], 'accepted': native['accepted'],
        'raw_accepted_csv_rows': len(raw_native_rows),
        'zero_initialization_records': sum(float(r['drain OuterVoltage']) == 0 for r in raw_native_rows),
        'selection': 'last record per outer voltage; matches the qualified native ledger',
        'attempts_including_zero_anchor': native['attempts_including_zero_anchor'], 'failed': native['failed'],
        'nonzero_goal_iterations': native['nonzero_goal_iterations'], 'failed_iterations': native['failed_iterations'],
        'zero_anchor_iterations': None, 'full_BV_qualified': False,
        'source_accepted_csv_sha256': sha(native_bytes), 'source_summary_sha256': sha(native_summary.read_bytes()),
        'regions': native['regions'], 'parameters': native['native_parameters']}
    summary = {'schema': 'vela.bvmethods.validation_summary.v1', 'integration_commit': INTEGRATION_COMMIT,
        'physical_constants_profile': 'sentaurus_u2022_12', 'native': native_public, 'candidates': candidates,
        'contract': {'nodes': 2600, 'triangles': 5020, 'temperature_K': 300, 'resistance_ohm_um': 1e7,
                     'BTBT': 'off', 'statistics': 'fermi_dirac_joyce_dixon',
                     'avalanche_drive': 'quasi_fermi_gradient_with_contact_electric_field_fallback',
                     'primary_coordinates': 'unit', 'protection_retained_fraction': 0.01,
                     'post_protection_derivatives': 'omit_mobility_and_avalanche_only',
                     'stopping': 'approximate E<1 OR mapped RHS infinity norm<1e-5',
                     'existing_BV_acceptance_anchor_V': 6.673833, 'BV_relative_tolerance': 0.005},
        'qualifications': {'ordinary_observer_full_closure': True, 'CMD_full_BV_chain': True,
                           'strict_native_stopping_equivalence': False, 'speed_noninferiority': False,
                           'all_report_only_conservation_checks_passed': False},
        'runtime_identity': runtime,
        'evidence': {'terminal_archive_sha256': EVIDENCE_SHA256, 'terminal_archive_members': len(members),
                     'report': '../../../docs/validation/bvmethods_production_closure_2026-11-09.md'},
        'limits': ['Native curve covers only 0-100 V; it is not a new native full-BV measurement.',
                   'Original CMD and JSON use different initialization/stage chains; compare their workload with that distinction.',
                   'Ordinary/CMD attempts and failures were not exported; null does not mean zero.',
                   'Native MaxStep is 100 V, Vela MaxStep 2500 V; endpoints and initialization differ.',
                   'Vela candidate curves never replace historical or native reference curves.']}
    files['comparison.json'] = encode(summary)
    files['README.md'] = b'''# Qualified U-2022.12 compatibility snapshot (2026-11-09)

This text-only snapshot accompanies integration commit `77a889a9` and the
[closure report](../../../docs/validation/bvmethods_production_closure_2026-11-09.md).
The parent directory's 2018 name is historical. Its old sources, curves,
materials and artifact manifests remain frozen. Vela candidates are not
native references; this snapshot is not a new remote-store registration.

## Contents and qualification

- `source/` retains the exact validated CMD/par bytes and explicit basic-state
  Plot sidecar. Unsupported current/source Plot exports are not qualified.
- `vela/materials.json` and `simulation_stage00/01.json` are the actual matching
  material and portable stage configurations, without numerical retuning.
- `vela/curves/` contains ordinary JSON full-BV candidates (62 points each)
  and the original-CMD candidate (106 points). Successful sweep Newton counts
  are 351, 351 and 732. Initialization/stage chains differ.
- `reference_curves/` contains the native U-2022.12 adaptive 0-100 V curve:
  45 accepted / 48 attempts / 3 failures, 308 nonzero-target Newton iterations.
  The raw PLT has three initialization records at zero; the neutral curve
  retains the last record per outer bias, matching the qualified ledger.
  Native zero-bias iterations are unknown in this nonzero-goal accounting.
- `comparison.json` records scope, units, report deficits, identities,
  unknown counters, qualification limits and evidence hashes.

Candidates reach the 1.3e-3 A/um stop and pass applicable F/E, load-line,
KCL/finite checks and the existing BV interpolation gate. Carrier/continuity
report failures remain 3/7, 3/8 and 70/13. Performance noninferiority and strict
native RHSMin equivalence are unqualified. The BV anchor is an existing
acceptance contract, not a BV measured from the bounded native curve.

## Reproduction inputs

Binary mesh and generated bulk mesh/doping remain external. Verify the mesh
against `source/inputs.json`, placing `n1_msh.tdr` beside `pp5_des.cmd`.
Use a matching `sentaurus_u2022_12` build with HDF5 import and PARDISO/MKL.
Set `PYTHONPATH` to the absolute repository `scripts/` directory (or use an
installed frontend package). Select `VELA_LINEAR_SOLVER=pardiso_mkl`, set
`VELA_MKL_RT_LIBRARY` to your local runtime and
`MKL_THREADING_LAYER=SEQUENTIAL`. From this snapshot directory:

```text
python -m vela_sdevice device source/pp5_des.cmd --allow-lossy --threads 1
```

The frontend generates mesh/doping/materials and stages in `source/pp5_des_vela`.
For direct use of the saved `vela/` configurations, first stage the verified
`mesh.json` and `doping.csv` in that folder and change into `vela/`. Run `simulation_stage00.json`, then
`simulation_stage01.json`, which loads `stage00_state.h5`. Backend selection
uses the documented runtime configuration. Shared-VM wall time is not an
exclusive speed benchmark. Output paths are relative to the run directory.
Use an external 12-hour supervisor as in the qualification; these decks do
not install a wall-clock supervisor by themselves.

No numerical run was repeated for this synchronization. Regenerate with
`scripts/freeze_bvmethods_u2022_validation.py` from the repository root,
explicit validated evidence/native inputs and a fresh version output folder.
'''
    manifest = {'schema': 'vela.bvmethods.text_snapshot.v1', 'version': VERSION,
                'kind': 'local_text_qualification_not_native_store_registration',
                'integration_commit': INTEGRATION_COMMIT,
                'files': {name: {'sha256': sha(data), 'bytes': len(data)} for name, data in sorted(files.items())}}
    output.mkdir(parents=True)
    for name, data in files.items():
        p = output / name; p.parent.mkdir(parents=True, exist_ok=True); p.write_bytes(data)
    (output / 'manifest.json').write_bytes(encode(manifest))
    return manifest

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-zip', type=Path, required=True)
    parser.add_argument('--native-summary', type=Path, required=True)
    parser.add_argument('--native-accepted', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = freeze(args.evidence_zip, args.native_summary, args.native_accepted, args.output)
    print(json.dumps({'version': manifest['version'], 'files': len(manifest['files'])}))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
