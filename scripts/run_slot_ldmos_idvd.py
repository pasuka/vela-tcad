"""Bounded cold IdVd validation with original full-device gates.

Uses an explicit, separately reviewed Fermi material/model configuration.
No imported native state, Boltzmann poly reset, or acceptance relaxation.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
from run_slot_ldmos_idl import qualify, save, sha

ROOT = Path(__file__).resolve().parents[1]


def midpoint(parent, target):
    return parent + (target - parent) / 2


def adaptive_growth(iterations, increment=1.5, reference_iterations=25):
    """Steady ramp rule; reference iterations do not alter Newton's own limit."""
    n = int(iterations)
    if n < 0 or n != float(iterations):
        raise ValueError('Invalid accepted Newton iteration count')
    return 1 + max(0., (increment-1)*(1-(max(1,n)-1)/(.75*reference_iterations)))


def advance_adaptive_gate(attempt, observe, *, target=5., initial=.005,
                          maximum=.25, minimum=5e-8, max_attempts=100,
                          max_failures=8, predictor_maximum=None, start=0.):
    """Only qualified updates move the parent; abandon failed large targets."""
    if not all(math.isfinite(v) and v>0 for v in (target,initial,maximum,minimum)):
        raise ValueError('Invalid adaptive gate voltage limits')
    if predictor_maximum is not None and (not math.isfinite(predictor_maximum) or predictor_maximum<=0):
        raise ValueError('Invalid gate predictor capability')
    effective_maximum=min(maximum,predictor_maximum) if predictor_maximum is not None else maximum
    if minimum>effective_maximum:raise ValueError('Gate minimum exceeds effective maximum')
    if not math.isfinite(start) or not 0<=start<target:raise ValueError('Invalid gate start')
    current=start;step=initial;failures=0
    for index in range(max_attempts):
        requested_step=step
        step=min(step,effective_maximum)
        proposed=min(target,current+step)
        span=proposed-current
        if span<=0:raise RuntimeError('gate_step_not_representable')
        result=attempt(proposed)
        decision=dict(parent_V=current,target_V=proposed,span_V=span,
                      qualified=bool(result['qualified']),attempt_index=index,
                      requested_step_V=requested_step,effective_maximum_V=effective_maximum,
                      limit_reason='predictor_capability' if requested_step>effective_maximum and
                          predictor_maximum is not None and predictor_maximum<maximum else
                          ('policy_maximum' if requested_step>effective_maximum else None))
        if result['qualified']:
            n=int(result['iv']['newton_iterations'])
            growth=adaptive_growth(n)
            current=proposed;failures=0;step=span*growth
            decision.update(newton_iterations=n,growth=growth,
                            next_requested_step_V=step,next_step_V=min(effective_maximum,step))
        else:
            failures+=1;step=span/2
            decision.update(consecutive_failures=failures,next_step_V=step)
        observe(decision)
        if current==target:return current
        if failures>=max_failures:raise RuntimeError('gate_adaptive_failure_limit')
        if step<minimum:raise RuntimeError('gate_adaptive_minimum_step')
    raise RuntimeError('gate_adaptive_attempt_limit')


def advance_progressive(current, target, last_span, attempt, max_attempts=32):
    """Bounded native interval: grow accepted spans, halve failures, one warm retry.

    attempt(voltage, allow_secant) returns an evidence entry. Failed trials never
    update the parent. This policy is independent of the device qualification.
    """
    if not (all(math.isfinite(v) for v in (last_span,current,target)) and last_span > 0
            and target > current and type(max_attempts) is int and max_attempts>0):
        raise ValueError('Progressive drain interval requires positive history')
    step=min(1.,1.25*last_span,target-current)
    failures=0; used=0; warm_used=False
    while current < target:
        proposed=min(target,current+step)
        if proposed<=current:raise RuntimeError('drain_step_not_representable')
        if used>=max_attempts:raise RuntimeError('drain_interval_attempt_limit')
        result=attempt(proposed,True);used+=1
        # A timeout is already expensive; do not repeat it from a worse initial guess.
        if (not result['qualified'] and result.get('predictor') is not None
                and result.get('returncode')!=124
                and (result.get('iv') or {}).get('failure_reason')!='soft_time_limit'
                and not warm_used and used<max_attempts):
            warm_used=True;result=attempt(proposed,False);used+=1
        if result['qualified']:
            accepted_span=proposed-current;current=proposed;failures=0
            step=min(1.,1.25*accepted_span,target-current)
        else:
            failures+=1
            if failures>3:raise RuntimeError('drain_strict_failure_after_progressive_halvings')
            step=(proposed-current)/2
    return current


def make_gate_seed(config, parent, output, target):
    from slot_fermi_poly_predictor import predict
    return predict(config, parent, output, ROOT, target)


def make_drain_seed(config, older, parent, output, target):
    from slot_idvd_secant_predictor import predict
    return predict(config, older, parent, output, ROOT, target)


def physics_identity(config):
    """Exclude output plumbing and the gate preparation bias, not physics/gates."""
    c=copy.deepcopy(config)
    c.pop('sweep',None);c.pop('output_csv',None)
    c['solver'].pop('performance_profiling',None)
    c['solver'].pop('progress_directory',None)
    c['solver'].pop('soft_time_limit_seconds',None)
    c['solver']['carrier_row_convergence'].pop('diagnostic_csv',None)
    for contact in c['contacts']:
        if contact['name']=='gate':contact['bias']=5.
    return c


def load_drain_history(path, base, resources):
    """Import two saved qualified states; this is not an independent re-solve."""
    import slot_state_archive as states
    source=json.loads(path.read_text())
    if not source.get('completed'):raise ValueError('Import batch is not complete')
    records=source.get('records',[])[-2:]
    if len(records)!=2:raise ValueError('Two qualified drain states required')
    for r in records:
        state=Path(r['state']);d=state.parent
        cfg=json.loads((d/'config.json').read_text())
        evidence=json.loads((d/'exit.json').read_text())
        with (d/'iv.csv').open() as f:row=list(csv.DictReader(f))[-1]
        if (not evidence.get('qualified') or evidence.get('returncode')!=0
            or sha(d/'config.json')!=evidence.get('config_sha256')
            or sha(state)!=r['state_sha256'] or row!=r['iv']
            or row!=evidence.get('iv') or not qualify(row,r['voltage_V'])):
            raise ValueError('Imported evidence/qualification mismatch')
        if physics_identity(cfg)!=physics_identity(base):raise ValueError('Imported physics/gates differ')
        if cfg['sweep']['contact']!='drain':raise ValueError('Not a drain history state')
        _,meta=states.read(state,cfg,ROOT)
        if meta.get('state_role')!='accepted' or meta.get('qualified') is False:
            raise ValueError('Imported state not accepted')
        biases=meta['contact_biases_V']
        if (abs(biases['gate']-5)>1e-12 or abs(biases['drain']-r['voltage_V'])>1e-12
            or any(abs(v)>1e-12 for k,v in biases.items() if k not in ('gate','drain'))):
            raise ValueError('Imported physical contact biases differ')
        for key in ('mesh_file','materials_file','node_doping_file'):
            if meta['input_file_sha256'][key]!=sha(base[key]):raise ValueError('Imported physical input hash differs')
        for p in (state,d/'config.json',d/'exit.json',d/'iv.csv'):resources[str(p)]=sha(p)
    if not 0<=records[0]['voltage_V']<records[1]['voltage_V']:
        raise ValueError('Import drain history must increase')
    resources[str(path)]=sha(path)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--runner', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--budget', type=float, default=2400.)
    parser.add_argument('--child-cap', type=float, default=120.)
    parser.add_argument('--gate-initialization', choices=['none', 'fermi_poly_poisson'],
                        default='none', help='Explicit predictor/Poisson handoff, gate stage only')
    parser.add_argument('--drain-initialization', choices=['warm', 'secant'], default='warm')
    parser.add_argument('--gate-step-policy', choices=['fixed','adaptive'], default='fixed')
    parser.add_argument('--gate-only', action='store_true', help='Stop after full qualification at Vg=5 V, Vd=0')
    parser.add_argument('--drain-step-policy', choices=['native_bisection','progressive'], default='native_bisection')
    parser.add_argument('--newton-progress', action='store_true', help='Diagnostic snapshots and cooperative timeout; requires compatible runner')
    parser.add_argument('--drain-import-summary', type=Path,
                        help='Import last two qualified drain states from a completed evidence batch')
    parser.add_argument('--drain-stop', type=float, default=60.,
                        help='Last native drain sample to request (bounded integration runs)')
    args = parser.parse_args()
    if not all(math.isfinite(x) and x>0 for x in (args.budget,args.child_cap,args.drain_stop)):
        raise ValueError('Budgets must be positive')
    if args.drain_step_policy=='progressive' and args.drain_initialization!='secant':
        raise ValueError('Progressive drain stepping requires secant initialization')
    if args.gate_only and args.drain_import_summary:
        raise ValueError('Gate-only cannot import drain states')
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    base = json.loads(args.config.read_text())
    def absolute(v):
        if isinstance(v, dict): return {k:absolute(x) for k,x in v.items()}
        if isinstance(v, list): return [absolute(x) for x in v]
        if isinstance(v, str) and v.startswith(('build/', 'configs/')): return str(ROOT / v)
        return v
    base = absolute(base)
    resources = {str(p.resolve()):sha(p) for p in
        [args.config, args.native, Path(__file__), ROOT/'scripts/run_slot_ldmos_idl.py']}
    if args.gate_initialization != 'none':
        for name in ('slot_fermi_poly_predictor.py', 'slot_poly_equilibrium.py',
                     'slot_state_archive.py', 'state_archive.py'):
            p=ROOT/'scripts'/name
            resources[str(p)]=sha(p)
            shutil.copy2(p, out/name)
    if args.drain_initialization=='secant' or args.drain_import_summary:
        for name in ('slot_idvd_secant_predictor.py','slot_state_archive.py','state_archive.py'):
            p=ROOT/'scripts'/name;resources[str(p)]=sha(p);shutil.copy2(p,out/name)
    for key in ['mesh_file', 'node_doping_file', 'materials_file']:
        resources[base[key]] = sha(base[key])
    external = base.get('mesh_geometry', {}).get('external_averagebox_couples_file')
    if external: resources[external] = sha(external)
    shutil.copy2(args.runner, out/'runner.exe')
    resources[str(out/'runner.exe')] = sha(out/'runner.exe')
    shutil.copy2(__file__, out/'driver.py')
    with args.native.open() as f: native = list(csv.DictReader(f))
    if any(abs(float(r['gate_V'])-5)>1e-12 for r in native):
        raise ValueError('Native reference must have Vg=5 V')
    native_by_v = {float(r['drain_V']):r for r in native}
    if args.drain_stop not in native_by_v:raise ValueError('Drain stop must be an exact native sample')
    history=load_drain_history(args.drain_import_summary.resolve(),base,resources) if args.drain_import_summary else []
    if history and args.drain_stop<=history[-1]['voltage_V']:raise ValueError('Drain stop must exceed imported parent')
    summary = dict(completed=False, goal_reached=False, cold_start=not bool(history),
        imported_records=copy.deepcopy(history),requested_drain_stop_V=args.drain_stop,
        records=[], attempts=[], resources=resources, budget_s=args.budget,
        child_cap_s=args.child_cap, gates=dict(eps_row=.001, closure=.01),
        policy=dict(gate_step=.1, drain_max_step=1., max_halvings=3,
                    gate_initialization=args.gate_initialization,
                    drain_initialization=args.drain_initialization,secant_ratio_max=3.,
                    secant_protocol='fermi_split_qf_log_density_v1',
                    failed_secant_retry='one_warm_same_target_then_bisection'))
    t0 = time.monotonic()
    summary['policy'].update(drain_step_policy=args.drain_step_policy,
        gate_step_policy=args.gate_step_policy,gate_only=args.gate_only,
        progressive_growth=1.25,progressive_max_attempts_per_native_interval=32,
        progressive_max_halvings=3,progressive_warm_retries_per_interval=1,
        newton_progress=args.newton_progress,soft_timeout_fraction=.8)
    if args.gate_step_policy=='adaptive':
        predictor_warning=None
        if args.gate_initialization=='fermi_poly_poisson':
            from slot_fermi_poly_predictor import GATE_BIAS_WARNING_THRESHOLD_V
            predictor_warning=GATE_BIAS_WARNING_THRESHOLD_V
        summary['policy'].update(gate_step=None,gate_initial_step_V=.005,
            gate_max_step_V=.25,gate_min_step_V=5e-8,gate_increment=1.5,
            gate_decrement=2.,gate_iteration_reference=25,gate_max_attempts=100,
            gate_max_consecutive_failures=8,gate_iteration_metric='full_newton_only_minimum_one',
            gate_predictor_warning_threshold_V=predictor_warning,
            gate_predictor_warning_only=True,gate_effective_max_step_V=.25)
        summary['step_decisions']=[]
    if args.drain_step_policy=='progressive':
        summary['policy']['failed_secant_retry']='one_warm_per_native_interval_except_timeout_then_bounded_halving'
    env = os.environ.copy()
    env['PATH'] = 'D:/msys64/ucrt64/bin;D:/msys64/usr/bin;' + env['PATH']
    if env.get('VELA_LINEAR_SOLVER'): raise ValueError('Unexpected linear solver override')
    save(out/'base.json', base)
    parent = Path(history[-1]['state']) if history else None
    def persist():
        summary['elapsed_s'] = time.monotonic()-t0
        save(out/'summary.json', summary)
    def run_child(directory, config, cap):
        remaining=args.budget-(time.monotonic()-t0)
        if remaining < 2: raise TimeoutError('bounded_budget_exhausted')
        actual_cap=min(cap, remaining-1.)
        if args.newton_progress and config['solver']['method']=='newton':
            config['solver']['progress_directory']=str(directory/'progress')
            config['solver']['soft_time_limit_seconds']=.8*actual_cap
        save(directory/'config.json',config)
        tick=time.monotonic()
        result=dict(child_cap_s=actual_cap, config_sha256=sha(directory/'config.json'))
        summary['active']=dict(directory=str(directory), method=config['solver']['method'],
                               target_V=config['sweep']['bias_points'][0], child_cap_s=actual_cap)
        persist()
        with (directory/'stdout.log').open('wb') as stdout, (directory/'stderr.log').open('wb') as stderr:
            try:
                result['returncode']=subprocess.run([str(out/'runner.exe'),'--config',str(directory/'config.json')],
                    cwd=ROOT,env=env,stdout=stdout,stderr=stderr,timeout=actual_cap).returncode
            except subprocess.TimeoutExpired: result['returncode']=124
        result['elapsed_s']=time.monotonic()-tick
        csv_path=Path(config['output_csv'])
        if csv_path.exists():
            with csv_path.open() as f: rows=list(csv.DictReader(f))
            if rows: result['iv']=rows[-1]
        summary.pop('active',None)
        return result
    def execute(stage, target, allow_secant=True):
        nonlocal parent
        remaining=args.budget-(time.monotonic()-t0)
        if remaining < 2: raise TimeoutError('bounded_budget_exhausted')
        child_cap = min(args.child_cap, remaining-1.)
        d=out/f'{len(summary["attempts"]):03d}_{stage}'
        d.mkdir()
        c=copy.deepcopy(base)
        c['output_csv']=str(d/'iv.csv')
        c['solver']['carrier_row_convergence']['diagnostic_csv']=str(d/'rows.csv')
        c['solver']['performance_profiling']=dict(enabled=True,json_file=str(d/'performance.json'))
        for contact in c['contacts']:
            contact['bias']=5. if stage=='drain' and contact['name']=='gate' else 0.
        c['sweep']=dict(mode='iv',contact='gate' if stage=='gate' else 'drain',
            current_contact='drain',bias_points=[target],start=target,stop=target,step=.005,
            write_vtk=False,write_state_file=str(d/'state.h5'),
            diagnostics={'newton_history':dict(enabled=True,attempts_csv_file=str(d/'attempts.csv'),
                iterations_csv_file=str(d/'iterations.csv'),rejected_state_directory=str(d/'rejected_states'))})
        if parent: c['sweep']['initial_state_file']=str(parent)
        else: c['sweep']['initialization']=dict(mode='poisson_block',
            diagnostic_csv=str(d/'initialization.csv'),write_state_file=str(d/'poisson.h5'))
        entry=dict(stage=stage,target_V=target,directory=str(d),qualified=False,
            parent_sha256=sha(parent) if parent else None,child_cap_s=child_cap)
        tick=time.monotonic()
        initialization_ok=True
        if stage=='drain' and args.drain_initialization=='secant' and allow_secant and len(history)>=2:
            old,latest=history[-2:]
            if sha(old['state'])!=old['state_sha256'] or sha(latest['state'])!=latest['state_sha256']:
                raise ValueError('Drain history hash changed')
            span=latest['voltage_V']-old['voltage_V']
            ratio=(target-latest['voltage_V'])/span if span>0 else float('inf')
            if 0<ratio<=3.:
                prep=time.monotonic()
                entry['predictor']=make_drain_seed(c,old['state'],latest['state'],d/'predictor.h5',target)
                entry['prediction_s']=time.monotonic()-prep
                entry['prediction_history_sha256']=[old['state_sha256'],latest['state_sha256']]
                c['sweep']['initial_state_file']=str(d/'predictor.h5')
            else:entry['prediction_skip_reason']='ratio_outside_(0,3]'
        if stage=='gate' and args.gate_initialization=='fermi_poly_poisson':
            seed=d/'predictor.h5'
            entry['predictor']=make_gate_seed(c,parent,seed,target)
            pd=d/'poisson';pd.mkdir()
            pc=copy.deepcopy(c)
            pc['solver']['method']='poisson_only'
            pc['solver']['poisson_only_quasi_fermi_initialization']='preserve_initial'
            pc['solver']['carrier_row_convergence']['diagnostic_csv']=str(pd/'rows.csv')
            pc['solver']['performance_profiling']['json_file']=str(pd/'performance.json')
            pc['output_csv']=str(pd/'iv.csv')
            pc['sweep'].update(initial_state_file=str(seed),write_state_file=str(pd/'state.h5'),
                diagnostics={'newton_history':dict(enabled=True,attempts_csv_file=str(pd/'attempts.csv'),
                    iterations_csv_file=str(pd/'iterations.csv'),rejected_state_directory=str(pd/'rejected_states'))})
            prep=run_child(pd,pc,min(args.child_cap,240.))
            prep['qualified']=False  # A Poisson solve cannot qualify a device.
            save(pd/'exit.json',prep);entry['initialization']=prep
            initialization_ok=(prep['returncode']==0 and prep.get('iv',{}).get('converged')=='1'
                               and (pd/'state.h5').exists())
            c['sweep']['initial_state_file']=str(pd/'state.h5')
        if initialization_ok:
            if parent:entry['initial_state_sha256']=sha(c['sweep']['initial_state_file'])
            entry.update(run_child(d,c,args.child_cap))
            entry['qualified']=(entry['returncode']==0 and (d/'state.h5').exists()
                                and qualify(entry.get('iv',{}),target))
        else:
            entry['failure_stage']='gate_poisson_initialization'
            entry['returncode']=entry['initialization']['returncode']
        entry['elapsed_s']=time.monotonic()-tick
        summary['attempts'].append(entry)
        if entry['qualified']:
            parent=d/'state.h5'
            record=dict(stage=stage,voltage_V=target,state=str(parent),state_sha256=sha(parent),iv=entry['iv'])
            if stage=='drain' and target in native_by_v:
                ni=float(native_by_v[target]['Id_A_per_um'])
                vi=float(entry['iv']['current_total_A_per_um'])
                record['native_Id_A_per_um']=ni
                record['current_relative_difference']=(vi-ni)/ni if abs(ni)>1e-18 else None
            summary['records'].append(record)
            if stage=='drain':history.append(record)
        save(d/'exit.json',entry);persist()
        print(json.dumps({k:v for k,v in entry.items() if k!='iv'}),flush=True)
        return entry['qualified']
    def advance(stage, current, target, depth=0):
        if execute(stage,target): return
        if stage=='drain' and 'predictor' in summary['attempts'][-1]:
            if execute(stage,target,allow_secant=False):return
        if depth>=3 or target==current: raise RuntimeError(f'{stage}_strict_failure_after_bounded_retries')
        mid=midpoint(current,target)
        advance(stage,current,mid,depth+1)
        advance(stage,mid,target,depth+1)
    try:
        persist()
        if not history:
            if not execute('equilibrium',0.): raise RuntimeError('cold_equilibrium_not_qualified')
            current=0.
            if args.gate_step_policy=='adaptive':
                def gate_attempt(v):
                    execute('gate',v)
                    return summary['attempts'][-1]
                def observe(decision):
                    summary['step_decisions'].append(decision);persist()
                current=advance_adaptive_gate(gate_attempt,observe)
            else:
                for k in range(1,51):
                    target=k/10;advance('gate',current,target);current=target
        if args.gate_only:
            summary['goal_reached']=True
            summary['stop_reason']='gate_5V_zero_drain_qualified'
            return 0
        current=history[-1]['voltage_V'] if history else 0.
        for target in sorted(native_by_v):
            if target>args.drain_stop or (history and target<=current):continue
            if args.drain_step_policy=='progressive' and len(history)>=2:
                def attempt(v,allow):
                    execute('drain',v,allow_secant=allow)
                    return summary['attempts'][-1]
                current=advance_progressive(current,target,
                    history[-1]['voltage_V']-history[-2]['voltage_V'],attempt)
                continue
            while target-current>1.+1e-12:
                advance('drain',current,current+1.);current+=1.
            advance('drain',current,target);current=target
        summary['goal_reached']=True
        summary['stop_reason']='requested_drain_segment_qualified' if args.drain_import_summary or args.drain_stop<60 else 'full_native_matched_idvd_path_qualified'
    except Exception as error:
        summary['stop_reason']=type(error).__name__+': '+str(error)
    finally:
        summary['completed']=True
        summary['resources_unchanged']=all(sha(p)==h for p,h in resources.items())
        persist()
    return 0 if summary['goal_reached'] else 2


if __name__=='__main__': raise SystemExit(main())
