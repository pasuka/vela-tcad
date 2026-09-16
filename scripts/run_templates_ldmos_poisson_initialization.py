"""Frozen two-point Poisson preparation controls; no solver/default/gate changes."""
import argparse
import copy
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

from analyze_templates_ldmos_predictor_study import state_delta
from run_templates_ldmos_density_projection import read, sha
from run_templates_ldmos_electrothermal_curve import state_gate


STATE_KEYS = ('state_interleaved', 'referenced_state_interleaved',
              'electron_qf_reference_V', 'hole_qf_reference_V', 'potential_origin_V')
VARIANTS = ('baseline', 'ptc_v1_s1', 'ptc_defect_ser_s1', 'ptc_defect_model_s1')


def restored_config(original, prepared):
    """Copy state only; original contact laws and all physical settings survive."""
    state_delta(original, prepared)  # Also checks finite values, sizes and origins.
    result = copy.deepcopy(original)
    for key in STATE_KEYS:
        result[key] = copy.deepcopy(prepared[key])
    result.update(initialization='provided_state', solve_mode='coupled')
    return result


def preparation_gate(original, prepared, audit):
    delta = state_delta(original, prepared)
    audit_delta = state_delta(prepared, audit)
    reasons = []
    if prepared['diagnostic_stop'] != 'diagnostic_scaled_residual':
        reasons.append('poisson_solve_not_closed')
    for name, data in (('prepared', prepared), ('restored_coupled', audit)):
        block = data['electrical_block_gates'][0]
        if not (block['satisfied'] and math.isfinite(block['weighted_l2'])
                and block['weighted_l2'] <= block['limit']):
            reasons.append(name+'_poisson_gate')
    if any(d > t for d, t in zip(delta[1:], (1e-14, 1e-14, 1e-12))):
        reasons.append('qf_or_temperature_changed')
    if audit['newton_updates'] != 0 or any(d > 1e-14 for d in audit_delta):
        reasons.append('audit_changed_prepared_state')
    return dict(pass_gate=not reasons, reasons=reasons,
                state_change=delta, audit_state_change=audit_delta,
                poisson_block=audit['electrical_block_gates'][0])


def variant_config(original, variant):
    cfg = copy.deepcopy(original)
    if variant not in VARIANTS:
        raise ValueError('Unknown initialization control')
    if variant != 'baseline':
        cfg.update(diagnostic_pseudo_transient=True, diagnostic_pseudo_time_scale=1.,
                   diagnostic_density_projection='v1')
    if variant.startswith('ptc_defect_'):
        cfg['diagnostic_pseudo_acceptance'] = 'defect_'+variant.split('_')[2]
    return cfg


def run_point(binary, directory, cfg):
    directory.mkdir()
    cfg = copy.deepcopy(cfg)
    cfg.update(performance_profiling=True, diagnostic_iteration_trace=True)
    inp, output = directory/'input.json', directory/'output.json'
    inp.write_text(json.dumps(cfg, allow_nan=False), encoding='utf-8')
    start = time.perf_counter()
    with (directory/'run.log').open('x') as log:
        proc = subprocess.Popen([str(binary), str(inp), str(output)], stdout=log,
            stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
        while True:
            try:
                code = proc.wait(timeout=30)
                break
            except subprocess.TimeoutExpired:
                print(directory.name+' running', flush=True)
    row = dict(directory=directory.name, exit_code=code, wall_seconds=time.perf_counter()-start,
               input_sha256=sha(inp))
    data = read(output) if output.exists() else None
    if data is not None:
        row.update(output_sha256=sha(output), newton_updates=data['newton_updates'],
            attempts=len(data['history']), stop=data['diagnostic_stop'],
            trials=sum(h['line_search_trials'] for h in data['history']),
            assemblies=data['performance']['assembly_calls'],
            factorizations=data['performance']['factorizations'])
    print(json.dumps(row), flush=True)
    return row, data


def total_cost(*rows):
    return {key: sum(r[key] for r in rows) for key in
            ('newton_updates', 'attempts', 'trials', 'assemblies', 'factorizations', 'wall_seconds')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--control', type=Path, required=True, help='Frozen previous acceptance matrix')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    # Reuse the exact frozen Release executable and DLLs, not an unverified rebuild.
    from analyze_templates_ldmos_pseudo_transient import checked_matrix
    previous, controls = checked_matrix(args.control)
    if {(r['variant'], r['vg']) for r, _ in controls} != {(v,g) for v in VARIANTS for g in (4,8)} or len(controls)!=8:
        raise ValueError('Expected four frozen variants at both gates')
    if sha(args.manifest) != previous['manifest_sha256']:
        raise ValueError('Frozen manifest differs from previous controls')
    cases = [c for c in read(args.manifest) if c['variant']=='linear_baseline'
             and (c['gate_V'], c['baseline_index']) in ((4,11),(8,9))]
    if len(cases)!=2:
        raise ValueError('Expected both frozen 5.333333 V cases')
    for case in cases:
        for row, _ in controls:
            if row['vg']==case['gate_V'] and (sha(Path(case['input']))!=row['source_input_sha256']
                    or sha(Path(case['expected']))!=row['source_output_sha256']):
                raise ValueError('Source state or reference changed')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    binary=out/'electrothermal_probe.exe'
    shutil.copy2(args.control/binary.name,binary)
    for dll in args.control.glob('*.dll'):
        shutil.copy2(dll,out/dll.name)
    if os.name=='nt':
        import ctypes
        ctypes.windll.kernel32.SetErrorMode(0x0001|0x0002|0x8000)
    report=dict(schema='vela.poisson_initialization.v1',status='running',
        control_directory=str(args.control.resolve()),control_summary_sha256=sha(args.control/'summary.json'),
        manifest_sha256=sha(args.manifest),binary_sha256=sha(binary),
        runtime_sha256={p.name:sha(p) for p in out.glob('*.dll')},preparations=[],runs=[],
        cost_policy='Full preparation and zero-update coupled audit charged to each variant; preparation executed once per gate.')
    def save():
        temp=out/'summary.tmp';temp.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
        os.replace(temp,out/'summary.json')
    save()
    for case in cases:
        name='vg%d_p%d'%(case['gate_V'],case['baseline_index'])
        original=read(Path(case['input']))
        cfg=copy.deepcopy(original);cfg.update(solve_mode='poisson',initialization='provided_state')
        prep, prepared=run_point(binary,out/('prepare_'+name),cfg)
        item=dict(vg=case['gate_V'],source_input_sha256=sha(Path(case['input'])),preparation=prep,ready=False)
        report['preparations'].append(item);save()
        if prep['exit_code']!=0 or prepared is None:
            continue
        restored=restored_config(original,prepared)
        cfg=copy.deepcopy(restored);cfg['diagnostic_newton_max_iterations']=0
        audit,audited=run_point(binary,out/('audit_'+name),cfg)
        item['audit']=audit
        if audit['exit_code']==0 and audited is not None:
            item['gate']=preparation_gate(original,prepared,audited)
            item['ready']=item['gate']['pass_gate']
        save()
        if not item['ready']:
            continue
        for variant in VARIANTS:
            row,data=run_point(binary,out/(variant+'_'+name),variant_config(restored,variant))
            row.update(variant=variant,vg=case['gate_V'],qualified=False,
                       charged_cost=total_cost(prep,audit,row) if data is not None else None)
            if data is not None:
                gate=state_gate(data,case['bias_V']);delta=state_delta(data,read(Path(case['expected'])))
                row.update(gate=gate,max_state_delta=delta,
                    qualified=row['exit_code']==0 and gate['pass_gate']
                    and all(d<=t for d,t in zip(delta,(1e-8,1e-8,1e-8,1e-7))))
            report['runs'].append(row);save()
    report['status']='completed_diagnostic_matrix';save()


if __name__=='__main__':
    main()
