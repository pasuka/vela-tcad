"""Freeze the engineering profile and audit existing evidence, without solving."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def index_points(records, expected):
    result = {}
    for row in records:
        vg = float(row['vg'])
        vd = float(row['vd'])
        if not math.isfinite(vg) or not math.isfinite(vd):
            raise ValueError('Nonfinite bias')
        i = round(vg / .05)
        if abs(vg - i*.05) > 1e-10:
            raise ValueError('Off-grid bias')
        key = (row['device'], vd, i)
        if key in result:
            raise ValueError(f'Duplicate point {key}')
        result[key] = row
    if set(result) != set(expected):
        raise ValueError(f'Coverage mismatch: missing {set(expected)-set(result)}, extra {set(result)-set(expected)}')
    return result


def audit(root=ROOT):
    base = root/'build/outlier_analysis_20260928'
    inputs = base/'inputs'
    input_hashes = read(inputs/'hashes.json')
    for rel, expected in input_hashes.items():
        path = inputs/rel
        if not path.exists() and rel.endswith('_reference.csv'):
            path = base/'references'/Path(rel).name
        if not path.exists():
            path = root/'build/hfs_srh_20260926/original_inputs'/rel
        if sha(path) != expected:
            raise ValueError(f'Input hash mismatch: {rel}')
    contract = read(inputs/'contract.json')
    keys = {(c['device'], float(c['vd']), i) for c in contract['cases'] for i in range(51)}
    if len(keys) != 816:
        raise ValueError('Expected original 816-point matrix')
    current = index_points(rows(base/'all_points.csv'), keys)
    strict = index_points(rows(root/'build/outlier_m60_review_20260928/joined_816.csv'), keys)
    profile = None
    for device in sorted({k[0] for k in keys}):
        template = inputs/device/'template.json'
        if not template.exists():
            template = root/'build/hfs_srh_20260926/original_inputs'/device/'template.json'
        cfg = read(template)
        selected = {k:cfg[k] for k in ('solver','scaling','mesh_geometry')}
        if profile is not None and selected != profile:
            raise ValueError(f'Nonuniform engineering profile: {device}')
        profile = selected
    errors = []
    for key, row in current.items():
        if row['qualified'] != 'True':
            raise ValueError(f'Unqualified frozen point: {key}')
        v = float(row['current_A_per_um'])
        s = strict[key]
        if v != float(s['vela_Id_A_per_um']):
            raise ValueError(f'Frozen current identity mismatch: {key}')
        error = 100*(v/float(s['tight_default_Id_A_per_um'])-1)
        if not math.isfinite(error) or abs(error-float(s['tight_error_percent']))>1e-10:
            raise ValueError(f'Error reconstruction mismatch: {key}')
        errors.append(error)
    attempts, last_rows = [], []
    for case in contract['cases']:
        folder = base/'sweeps'/case['case']
        cfg = read(folder/'config.json')
        for name in profile:
            if cfg[name] != profile[name]:
                raise ValueError(f'Sweep/profile mismatch: {case["case"]}/{name}')
        trace = rows(folder/'iterations.csv')
        last = {}
        for r in trace:
            last[(r['run_id'],r['segment_id'],r['attempt_id'])] = r
        for r in rows(folder/'attempts.csv'):
            key = (r['run_id'],r['segment_id'],r['attempt_id'])
            if key not in last:
                raise ValueError(f'Missing trace: {case["case"]}/{key}')
            attempts.append(dict(case=case['case'],**r))
            last_rows.append(dict(case=case['case'],status=r['status'],reason=r['reason'],**last[key]))
    accepted = [r for r in attempts if r['status']=='accepted']
    for r in last_rows:
        if r['status']!='accepted':continue
        if int(r['carrier_row_violations'])!=0 or not float(r['carrier_row_max_ratio'])<=1e-6:
            raise ValueError(f'Accepted trace violates carrier gate: {r["case"]}')
        if r['reason']=='abstol' and not float(r['residual_norm'])<=profile['solver']['abstol']:
            raise ValueError('Absolute acceptance mismatch')
        if r['reason']=='reltol' and not float(r['relative_residual_norm'])<=profile['solver']['reltol']:
            raise ValueError('Relative acceptance mismatch')
    overlap = read(root/'build/simplemos_ep_20260928/probe_results/complete_cold_overlap.json')
    evidence = dict(original_points=816,qualified_points=len(current),
        m60_current_passed=sum(abs(e)<=2 for e in errors),max_abs_m60_error_percent=max(map(abs,errors)),
        cold_warm_control_comparisons=overlap,full_816_dual_initialization_complete=False,
        gate_sweep_attempts=len(attempts),accepted_gate_attempts=len(accepted),
        accepted_reasons=dict(Counter(r['reason'] for r in accepted)),
        failed_reasons=dict(Counter(r['reason'] for r in attempts if r['status']!='accepted')),
        observed_max_final_residual=max(float(r['final_residual_norm']) for r in accepted),
        observed_max_final_carrier_row_ratio=max(float(r['carrier_row_max_ratio']) for r in last_rows if r['status']=='accepted'),
        trace_scope='Gate sweeps only; equilibrium and drain-ramp traces are not present in this local collection',
        input_files_verified=len(input_hashes),reference_changed=False,acceptance_changed=False)
    policy = dict(schema='vela.simplemos.engineering-contract.v1',
        scope=dict(devices=sorted({k[0] for k in keys}),vd_V=[.05,1.],vg_V=[i*.05 for i in range(51)],temperature_K=300),
        engineering_reference=dict(name='M60',algorithm='default',Digits=8,ErrRef_e_cm3=100,ErrRef_h_cm3=100,
            RhsMin_default=1e-5,precision='binary64',current_error_percent=2),
        calibration_reference=dict(name='EP128',Digits=15,RhsMin=1e-15,precision='double-double',full_matrix=False),
        numerical_gates=dict(carrier_row_ratio=1e-6,kcl_over_Id=1e-8,port_relative=1e-8),
        dual_gates=dict(potential_max_V=1e-6,density_relative=1e-4,Id_relative=1e-6),
        native_fields=dict(mode='report_only',new_thresholds_frozen=False),
        unsupported_scope='No global defaults; no additional temperatures, quantum, avalanche, Fermi or arbitrary mesh qualification',
        frozen_input_hashes=input_hashes,
        profile=profile)
    return policy,evidence,attempts,last_rows


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    policy,evidence,attempts,last_rows=audit()
    args.out.mkdir(parents=True,exist_ok=True)
    for name,obj in [('engineering_contract.json',policy),('evidence_audit.json',evidence)]:
        (args.out/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for name,data in [('gate_attempts.csv',attempts),('gate_final_iterations.csv',last_rows)]:
        with (args.out/name).open('w',encoding='utf-8',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    print(json.dumps(evidence,indent=2))


if __name__=='__main__':
    main()
