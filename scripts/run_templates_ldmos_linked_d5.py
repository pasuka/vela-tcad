#!/usr/bin/env python3
"""Run linked D5 or explicit experimental D4 sweeps from hashed external inputs.

No historical Python adapters are imported. Run --preflight to check inputs
without creating output. The template uses the repository unit_scaling schema.
"""
from __future__ import annotations
import argparse
import csv
import ctypes
from ctypes import wintypes
from datetime import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import shutil
import time
from copy import deepcopy

from run_templates_ldmos_stage4_d5 import assert_d5_contract, read_points
from analyze_templates_ldmos_stage4_d5 import curve_error, kcl_audit, LIMITS
from summarize_templates_ldmos_idvd_ablation import read_curve
from translate_dd_state import translate_state_csv
from dc_worker_client import DCWorker

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUNDLE = ROOT/'reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_units_inputs.json'
INITIAL, MINIMUM, MAXIMUM, GROWTH, BUDGET = .0025, .0025, .2, 1.35, 160
PHYSICS_PROFILE = 'D5'
FRAME_PIVOT, FRAME = 26.6666666666667, 28.
GATES = dict(mode='enforce', psi_residual_ceiling=5e-8,
             electron_residual_ceiling=1e-11, hole_residual_ceiling=3e-10)

def assert_linked_physics_contract(config):
    if PHYSICS_PROFILE == 'D5':
        return assert_d5_contract(config)
    if PHYSICS_PROFILE != 'D4':
        raise ValueError('Unknown linked physics profile')
    solver=config['solver'];mobility=solver['mobility']
    if mobility.get('model')!='ialmob' or not mobility.get('jacobian_field_derivatives',True):
        raise ValueError('D4 requires coupled IALMob derivatives')
    if solver.get('impact_ionization',{}).get('model')!='none':
        raise ValueError('D4 requires avalanche disabled')
    if any('quantum' in str(key).lower() for key in solver):
        raise ValueError('D4 classical profile excludes quantum potential')
    if not mobility['ialmob'].get('high_field',True) or mobility['ialmob'].get('reference_density_m3',1e18)!=1e18:
        raise ValueError('D4 requires the original HFS and RefDens settings')
    if config.get('discretization',{}).get('poisson_charge_volume_policy')!='material_local':
        raise ValueError('D4 requires qualified material-local charge volume')
    if config.get('mesh_geometry',{}).get('carrier_transport_couple_profile')!='templates_ldmos_external_averagebox':
        raise ValueError('D4 requires qualified external AverageBox transport')
    if 'predictor' in config['sweep'] or 'continuation' in config['sweep']:
        raise ValueError('D4 rejects unqualified internal predictor/continuation')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frame_offset_for(bundle,gate,override=None):
    value=float(override if override is not None else bundle.get('frame_offsets_V',{}).get(str(gate),28.))
    if not math.isfinite(value) or not 0.<value<=40.:
        raise ValueError('Invalid frame offset')
    return value

def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def stamp():
    return datetime.now().astimezone().isoformat()

def write(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    for retry in range(40):
        try:
            temp.replace(path)
            return
        except PermissionError:
            if retry == 39:
                raise
            time.sleep(.1)

def rows(path):
    with path.open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))

def csv_out(path, data):
    with path.open('x', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)

def close_bias(a,b):return abs(a-b)<=32*math.ulp(max(abs(a),abs(b),1.))

def state_difference(a,b):
    left,right=rows(a),rows(b)
    assert len(left)==len(right)==10241
    out={}
    for key in ('psi','phin','phip','electrons_m3','holes_m3'):
        pairs=[(float(x[key]),float(y[key])) for x,y in zip(left,right)]
        assert all(x['node_id']==y['node_id'] for x,y in zip(left,right))
        out[key]={'max_absolute':max(abs(x-y) for x,y in pairs),
            'max_relative_floor_1':max(abs(x-y)/max(abs(y),1.) for x,y in pairs)}
    return out

def cpu_times(proc):
    creation,exit_time,kernel,user=(wintypes.FILETIME() for _ in range(4))
    get_times=ctypes.WinDLL('kernel32',use_last_error=True).GetProcessTimes
    get_times.argtypes=[wintypes.HANDLE]+[ctypes.POINTER(wintypes.FILETIME)]*4
    get_times.restype=wintypes.BOOL
    if not get_times(wintypes.HANDLE(int(proc._handle)),ctypes.byref(creation),
                     ctypes.byref(exit_time),ctypes.byref(kernel),ctypes.byref(user)):
        raise ctypes.WinError(ctypes.get_last_error())
    seconds=lambda x:((x.dwHighDateTime<<32)|x.dwLowDateTime)/1e7
    return {'kernel':seconds(kernel),'user':seconds(user),'total':seconds(kernel)+seconds(user)}

def next_step(mode, proposed, actual, updates, density_used=False):
    if not math.isfinite(actual) or actual <= 0 or updates < 0:
        raise ValueError('Invalid accepted step or work count')
    if mode == 'fixed':
        return min(MAXIMUM, proposed * GROWTH), GROWTH
    if mode != 'newton_iterations':
        raise ValueError(mode)
    factor = 1. if density_used else 1. + (GROWTH - 1.) * max(
        0., 1. - max(updates - 1, 0) / (.75 * BUDGET))
    return min(MAXIMUM, max(MINIMUM, actual * factor)), factor

def ports(dest, target):
    terminals = rows(dest / 'terminal_balance.csv')
    assert len(terminals) == 4
    assert {r['contact'] for r in terminals} == {'source', 'drain', 'gate', 'substrate'}
    assert all(close_bias(float(r['bias_V']), target) for r in terminals)
    values = {r['contact']: float(r['current_total_A_per_um']) for r in terminals}
    assert all(math.isfinite(v) for v in values.values())
    return values, abs(math.fsum(values.values())) / max(max(map(abs, values.values())), 1e-30)

def blocks_pass(row):
    return all(math.isfinite(float(row[key])) and float(row[key]) <= limit for key, limit in (
        ('final_psi_residual_norm', 5e-8), ('final_electron_continuity_residual_norm', 1e-11),
        ('final_hole_continuity_residual_norm', 3e-10)))

def raw_good(result, dest, target):
    if result['returncode'] != 0 or len(result['curve']) != 1:
        return False
    row = result['curve'][0]
    return (row['converged'] == '1' and close_bias(float(row['bias_V']), target)
            and blocks_pass(row) and int(row['carrier_row_violations']) == 0
            and float(row['carrier_row_max_ratio']) <= 1e-8 and ports(dest, target)[1] <= 1e-8)

def density_eligible(result):
    return bool(result['curve'] and blocks_pass(result['curve'][-1])
                and int(result['curve'][-1]['carrier_row_violations']) > 0)

def frame_of(dest):
    return -next(c['bias'] for c in read(dest/'control.json')['contacts'] if c['name']=='source')

def good(result, dest, physical_target):
    return raw_good(result, dest, physical_target-frame_of(dest))

def prepare_config(base, dest, seed, target, cap, gate, frame=0.):
    """Relocate outputs explicitly; keep input physics and original gates intact."""
    def expand(value):
        if isinstance(value,dict): return {k:expand(v) for k,v in value.items()}
        if isinstance(value,list): return [expand(v) for v in value]
        if isinstance(value,str):
            return value.replace('@workspace',str(ROOT)).replace('@output',str(dest))
        return value
    cfg=expand(deepcopy(base))
    for contact in cfg['contacts']:
        contact['bias']=(gate if contact['name']=='gate' else target if contact['name']=='drain' else 0.)-frame
    raw=target-frame
    cfg['sweep'].update(start=raw, stop=raw, bias_points=[raw],
        step=cap, initial_step=cap, min_step=cap, max_step=cap, initial_state_file=str(seed))
    cfg['solver']['quasi_fermi_update_limit_V']=cap
    assert_linked_physics_contract(cfg)
    if cfg['solver']['block_absolute_convergence'] != GATES:
        raise ValueError('Changed original block gates')
    return cfg


def run_child(name, parent, parent_bias, target, cap, frame):
    global WORKER
    dest=HERE/name
    dest.mkdir(parents=True,exist_ok=False)
    cfg=prepare_config(BASE,dest,parent,target,cap,args.gate,frame)
    write(dest/'control.json',cfg)
    if digest(RUNNER)!=EXPECTED: raise ValueError('Runner changed')
    status=dict(case=name,seed=str(parent),seed_sha256=digest(parent),
        runner_sha256=EXPECTED,config_sha256=digest(dest/'control.json'),
        status='running',started_at=stamp())
    started=time.perf_counter()
    with (dest/'stdout.log').open('w',encoding='utf-8') as log:
        if args.worker:
            if WORKER is None: WORKER=DCWorker(RUNNER,HERE,ENV,cpu_times,reuse_linear_analysis=args.reuse_linear_analysis)
            status['pid']=WORKER.proc.pid;write(dest/'status.json',status)
            response,cpu=WORKER.run(dest/'control.json')
            log.write(json.dumps(response)+'\n')
            status['worker_request_id']=response['id']
            code=response['returncode']
            if 'error' in response: raise RuntimeError(response['error'])
        else:
            proc=subprocess.Popen([str(RUNNER),'--config',str(dest/'control.json')],
                cwd=dest,env=ENV,stdout=log,stderr=subprocess.STDOUT)
            status['pid']=proc.pid;write(dest/'status.json',status)
            try:
                code=proc.wait(timeout=600)
            except subprocess.TimeoutExpired:
                proc.kill();proc.wait()
                status.update(status='timeout',returncode=proc.returncode)
                write(dest/'status.json',status)
                raise
            cpu=cpu_times(proc)
        status.update(status='completed',returncode=code,wall_seconds=time.perf_counter()-started,
            cpu_seconds=cpu,finished_at=stamp(),
            attempts=rows(dest/'newton_attempts.csv'),curve=rows(dest/'curve.csv'))
    trace=rows(dest/'newton_iterations.csv')
    accepted=sum(r['event']=='accepted_iteration' for r in trace)
    updates=sum(int(r['newton_iterations']) for r in status['attempts'])
    if accepted!=updates: raise ValueError('Newton update accounting mismatch')
    status['trace_updates']=accepted
    write(dest/'status.json',status)
    write(dest/'lineage.json',dict(physical_parent_V=parent_bias,target_V=target,
        parent_state=str(parent),parent_sha256=digest(parent),native_sweep=False,budget=BUDGET,cap_V=cap))
    return status


def translate(source,destination,subtract):
    mesh=read(Path(BASE['mesh_file'].replace('@workspace',str(ROOT))))
    si={int(r['id']) for r in mesh['regions'] if r['material'].lower() in ('si','silicon')}
    nodes={int(n) for t in mesh['triangles'] if int(t['region_id']) in si for n in t['node_ids']}
    translate_state_csv(source,destination,subtract,nodes)

class CheckpointSweep:

    def __init__(self, mode):
        self.mode, self.out = mode, HERE / mode
        self.out.mkdir(exist_ok=False)
        self.ledger = dict(mode=mode, status='running', started_at=stamp(), accepted_bias_V=0.,
                           accepted_state=None, next_step_V=INITIAL, runs=[], transfers=[],
                           rollbacks=[], exact_points=[], active_child=None)
        self.save()

    def save(self):
        done = [r for r in self.ledger['runs'] if r.get('status') == 'completed']
        self.ledger['total_child_cpu_seconds'] = sum(r['cpu_seconds'] for r in done)
        self.ledger['total_child_wall_seconds'] = sum(r['wall_seconds'] for r in done)
        self.ledger['total_Newton_updates'] = sum(r['Newton_updates'] for r in done)
        self.ledger['updated_at'] = stamp()
        write(self.out / 'ledger.json', self.ledger)

    def point(self, dest, bias):
        state = rows(dest / 'state.csv')
        assert len(state) == len({int(r['node_id']) for r in state}) == 10241
        for row in state:
            assert all(math.isfinite(float(row[k])) for k in ('psi', 'phin', 'phip', 'electrons_m3', 'holes_m3'))
            assert min(float(row['electrons_m3']), float(row['holes_m3'])) >= 0.
        return dict(bias_V=bias, case=dest.relative_to(HERE).as_posix(), state=str(dest / 'state.csv'),
                    sha256=digest(dest / 'state.csv'))

class FrameSweep(CheckpointSweep):

    def __init__(self):
        super().__init__('fixed')
        self.ledger.update(frame_V=0., frame_events=[])
        self.save()

    def child(self, parent, parent_bias, target, cap, stage, density=False, frame=None):
        frame = self.ledger['frame_V'] if frame is None else frame
        name = f'fixed/child_{len(self.ledger["runs"]):05d}_{stage}'
        record = dict(case=name, stage=stage, status='running', parent_state=str(parent),
            parent_sha256=digest(parent), parent_V=parent_bias, target_V=target, cap_V=cap,
            frame_V=frame, raw_parent_V=parent_bias-frame, raw_target_V=target-frame,
            started_at=stamp())
        self.ledger['runs'].append(record)
        self.ledger['active_child'] = name
        self.save()
        result = run_child(name, parent, parent_bias, target, cap, frame)
        record.update(status='completed', returncode=result['returncode'],
            cpu_seconds=result['cpu_seconds']['total'], wall_seconds=result['wall_seconds'],
            Newton_updates=sum(int(r['newton_iterations']) for r in result['attempts']),finished_at=stamp())
        self.ledger['active_child'] = None
        self.save()
        return result, HERE/name

    def point(self, dest, bias):
        point = super().point(dest,bias)
        frame = frame_of(dest)
        point.update(frame_V=frame, raw_bias_V=bias-frame)
        return point

    def change_frame(self):
        assert self.ledger['frame_V']==0. and self.ledger['accepted_bias_V']==FRAME_PIVOT
        parent = Path(self.ledger['accepted_state'])
        suffix='' if not self.ledger['frame_events'] else f'_{len(self.ledger["frame_events"])}'
        inputs = self.out/f'frame_inputs{suffix}'
        inputs.mkdir(exist_ok=False)
        translated = inputs/f'frame{FRAME:g}_seed.csv'
        translate(parent,translated,FRAME)
        event = dict(status='running',physical_bias_V=FRAME_PIVOT,old_frame_V=0.,new_frame_V=FRAME,
            parent_state=str(parent),parent_sha256=digest(parent),translated_seed=str(translated),
            translated_seed_sha256=digest(translated),child_start=len(self.ledger['runs']))
        self.ledger['frame_events'].append(event)
        self.save()
        result, dest = self.child(translated,FRAME_PIVOT,FRAME_PIVOT,.1,'frame_baseline',frame=FRAME)
        for stage in ('frame_reclose','frame_density'):
            if good(result,dest,FRAME_PIVOT): break
            density = stage=='frame_density'
            if density and not density_eligible(result): break
            rejected = [r for r in result['attempts'] if r['status']=='rejected']
            if not rejected: break
            result,dest = self.child(Path(rejected[-1]['rejected_final_state_file']),FRAME_PIVOT,
                FRAME_PIVOT,.1,stage,density,frame=FRAME)
        event['child_end']=len(self.ledger['runs'])
        if not good(result,dest,FRAME_PIVOT):
            event['status']='failed_closure';self.save()
            raise RuntimeError('Planned frame change did not qualify under original gates')
        canonical = inputs/'canonical_for_comparison.csv'
        translate(dest/'state.csv',canonical,-FRAME)
        delta = state_difference(canonical,parent)
        new_ports,new_kcl = ports(dest,FRAME_PIVOT-FRAME)
        old_ports,old_kcl = ports(parent.parent,FRAME_PIVOT)
        current = max(abs(new_ports[k]-old_ports[k]) for k in new_ports)/max(max(map(abs,old_ports.values())),1e-30)
        checks = dict(current=current<=1e-8,KCL=max(new_kcl,old_kcl)<=1e-8,
            potentials=max(delta[k]['max_absolute'] for k in ('psi','phin','phip'))<=1e-7,
            densities=max(delta[k]['max_relative_floor_1'] for k in ('electrons_m3','holes_m3'))<=1e-6)
        event.update(status='qualified' if all(checks.values()) else 'failed_equivalence',
            state_difference=delta,current_difference_over_max_terminal=current,KCL=[old_kcl,new_kcl],checks=checks,
            canonical_for_comparison=str(canonical),canonical_sha256=digest(canonical),
            accepted=self.point(dest,FRAME_PIVOT))
        self.save()
        if not all(checks.values()): raise RuntimeError('Frame-equivalence guard failed')
        self.ledger.update(frame_V=FRAME,accepted_state=str(dest/'state.csv'))
        self.save()
        print('FRAME_QUALIFIED',FRAME_PIVOT,'frame',FRAME,flush=True)

class Sweep(FrameSweep):

    def prediction_guard(self, stage, parent_bias, target):
        transfers = self.ledger['transfers']
        retried = any(abs(r['parent_V']-parent_bias)<1e-12 for r in self.ledger['rollbacks'])
        guard_reason = 'non_direct' if stage != 'direct' else 'no_history' if not transfers else 'retry' if retried else None
        if guard_reason is None and transfers[-1].get('frame_V', 0.) != self.ledger['frame_V']:
            guard_reason = 'frame_changed'
        if guard_reason is None:
            previous_step = parent_bias - transfers[-1]['parent_V']
            ratio_raw = (target-parent_bias)/previous_step
            clipped = transfers[-1]['accepted_step_V'] < transfers[-1]['proposed_step_V'] - 1e-12
            if clipped:
                guard_reason = 'previous_step_clipped'
            elif ratio_raw > 2.:
                guard_reason = 'extrapolation_ratio_exceeded'
        return guard_reason

    def advance(self, reference):
        rollback_count = 0
        while not close_bias(self.ledger['accepted_bias_V'], reference):
            if (HERE / 'STOP').exists():
                raise RuntimeError('Requested administrative stop')
            start = self.ledger['accepted_bias_V']
            parent = Path(self.ledger['accepted_state'])
            requested = self.ledger['next_step_V']
            planned_target = min(start + requested, reference)
            guard_reason = self.prediction_guard('direct', start, planned_target)
            proposal = min(requested, .1) if guard_reason else requested
            target = min(start + proposal, reference)
            cap = target - start
            first_run = len(self.ledger['runs'])
            self.ledger.setdefault('step_predictions', []).append(dict(
                child_start=first_run, parent_V=start, reference_V=reference,
                requested_step_V=requested, effective_proposal_V=proposal,
                actual_step_V=cap, predictor_guard_reason=guard_reason,
                limited=proposal < requested))
            density_used = False
            result, dest = self.child(parent, start, target, cap, 'direct')
            for stage in ('reclose', 'density'):
                if good(result, dest, target):
                    break
                if stage == 'density' and not density_eligible(result):
                    break
                rejected = [r for r in result['attempts'] if r['status'] == 'rejected']
                if not rejected:
                    break
                recovery_seed = Path(rejected[-1]['rejected_final_state_file'])
                density_used = stage == 'density'
                result, dest = self.child(recovery_seed, target, target, cap, stage, density_used)
            if not good(result, dest, target):
                smaller = cap * .5
                self.ledger['rollbacks'].append(dict(parent_V=start, rejected_target_V=target,
                    next_step_V=smaller, child_start=first_run, child_end=len(self.ledger['runs'])))
                self.ledger['next_step_V'] = smaller
                self.save()
                rollback_count += 1
                if smaller < MINIMUM or rollback_count > 12:
                    raise RuntimeError(f'{self.mode}: cutback exhausted after {start} V; failed {target} V')
                continue
            work = sum(r['Newton_updates'] for r in self.ledger['runs'][first_run:])
            next_proposal, factor = next_step(self.mode, proposal, cap, work, density_used)
            point = self.point(dest, target)
            self.ledger['transfers'].append(dict(**point, parent_V=start, parent_state=str(parent),
                proposed_step_V=proposal, accepted_step_V=cap, Newton_updates=work,
                child_start=first_run, child_end=len(self.ledger['runs']), density_recovery=density_used,
                growth_factor=factor, next_step_V=next_proposal))
            self.ledger.update(accepted_bias_V=target, accepted_state=point['state'], next_step_V=next_proposal)
            self.save()
            rollback_count = 0
            print('ACCEPT', self.mode, repr(target), 'N', work, 'next', repr(next_proposal), flush=True)
        self.ledger['exact_points'].append(self.point(Path(self.ledger['accepted_state']).parent, reference))
        self.save()
        print('EXACT', self.mode, repr(reference), flush=True)

    def child(self, parent, parent_bias, target, cap, stage, density=False, frame=None):
        metadata = None
        transfers = self.ledger['transfers']
        guard_reason = self.prediction_guard(stage, parent_bias, target)
        if guard_reason is None:
            last = transfers[-1]
            assert abs(last['bias_V']-parent_bias)<1e-12
            previous = Path(last['parent_state'])
            ratio = min(2., (target-parent_bias)/(parent_bias-last['parent_V']))
            current_rows, previous_rows = rows(parent), rows(previous)
            assert len(current_rows)==len(previous_rows)==10241
            for current, old in zip(current_rows, previous_rows):
                assert current['node_id']==old['node_id']
                for field in ('psi','phin','phip'):
                    current[field] = format(float(current[field])+ratio*(float(current[field])-float(old[field])), '.17g')
                for field, reference, increment in [('phin','electron_qf_reference_V','electron_qf_increment_V'),('phip','hole_qf_reference_V','hole_qf_increment_V')]:
                    current[reference]=current[field]
                    current[increment]='0'
            prediction=HERE/'predictor_inputs'/f"seed_{len(self.ledger['runs']):05d}.csv"
            prediction.parent.mkdir(exist_ok=True)
            csv_out(prediction,current_rows)
            metadata=dict(mode='outer_secant_guarded_v1', ratio=ratio, previous_state=str(previous),
                previous_sha256=digest(previous), current_state=str(parent), current_sha256=digest(parent),
                predicted_state=str(prediction), predicted_sha256=digest(prediction))
            parent=prediction
        result, folder=super().child(parent,parent_bias,target,cap,stage,density,frame)
        self.ledger['runs'][-1]['outer_predictor']=metadata
        self.ledger['runs'][-1]['predictor_guard_reason']=guard_reason
        self.save()
        return result, folder

def initialize(sweep):
    result,dest=sweep.child(SEED,0.,0.,.1,'initial')
    for stage in ('initial_reclose','initial_density'):
        if good(result,dest,0.):break
        density=stage=='initial_density'
        if density and not density_eligible(result):break
        rejected=[r for r in result['attempts'] if r['status']=='rejected']
        if not rejected:break
        result,dest=sweep.child(Path(rejected[-1]['rejected_final_state_file']),0.,0.,.1,stage,density)
    if not good(result,dest,0.):raise RuntimeError('Initial candidate state did not qualify')
    point=sweep.point(dest,0.)
    sweep.ledger['exact_points'].append(point)
    sweep.ledger.update(accepted_bias_V=0.,accepted_state=point['state']);sweep.save()

def score(sweep):
    target_points=[(v,i) for v,i in read_curve(REFERENCE) if 0.-1e-12<=v<=stop+1e-12]
    assert [p['bias_V'] for p in sweep.ledger['exact_points']]==[v for v,i in target_points]
    out=HERE/'score';out.mkdir(exist_ok=False)
    curve=[];balance=[];provenance=[]
    for index,point in enumerate(sweep.ledger['exact_points']):
        dest=HERE/point['case'];status=read(dest/'status.json');bias=point['bias_V']
        assert good(status,dest,bias)
        assert digest(Path(point['state']))==point['sha256']
        curve.append(dict(bias_V=bias,current_total_A_per_um=float(status['curve'][0]['current_total_A_per_um'])))
        for row in rows(dest/'terminal_balance.csv'):
            row.update(point_index=str(index),bias_V=repr(bias));balance.append(row)
        provenance.append(point)
    csv_out(out/'curve.csv',curve);csv_out(out/'terminal_balance.csv',balance)
    report=dict(exact_points=len(curve),full_curve=stop==40.,provenance=provenance,KCL=kcl_audit(out/'terminal_balance.csv'))
    if stop==40.:
        metrics=curve_error(target_points,[(x['bias_V'],x['current_total_A_per_um']) for x in curve]);report['metrics']=metrics;report['verdicts']={}
        for level,limits in LIMITS.items():
            observed=dict(median=metrics['relative_error_percent']['median'],p95=metrics['relative_error_percent']['p95'],ron=metrics['low_vd_differential_resistance_error_percent'],endpoint=metrics['endpoint_current_error_percent'],kcl=report['KCL']['max_normalized_kcl_percent'])
            checks={k:observed[k]<=limits[k] for k in observed};report['verdicts'][level]=dict(checks=checks,pass_all=all(checks.values()))
    write(out/'summary.json',report)
    return report.get('verdicts',{'partial_curve':f'not a full {PHYSICS_PROFILE} score'})

def audit(sweep,plan,*,write_output=True):
    for path,expected in plan['frozen_files'].items():assert digest(Path(path))==expected,path
    accepted=sweep.ledger['exact_points']+sweep.ledger['transfers']+[e['accepted'] for e in sweep.ledger['frame_events'] if e['status']=='qualified']
    for run in sweep.ledger['runs']:
        dest=HERE/run['case'];status=read(dest/'status.json');cfg=read(dest/'control.json')
        assert status['config_sha256']==digest(dest/'control.json')
        assert status['runner_sha256']==EXPECTED
        assert status['seed_sha256']==run['parent_sha256']==digest(Path(run['parent_state']))
        expected=prepare_config(BASE,dest,Path(run['parent_state']),run['target_V'],
                                run['cap_V'],args.gate,run['frame_V'])
        assert cfg['solver']['mobility']==expected['solver']['mobility']
        assert cfg['solver']['block_absolute_convergence']==GATES
        assert cfg['solver']['carrier_row_convergence']['mode']=='enforce' and cfg['solver']['carrier_row_convergence']['eps_row']==1e-8
        assert cfg['solver']['quasi_fermi_update_limit_V']==run['cap_V']
        assert_linked_physics_contract(cfg)
        assert all(r['predictor_state_hash']=='' for r in status['attempts'])
    seen=set()
    for point in accepted:
        if point['state'] in seen:continue
        seen.add(point['state']);dest=HERE/point['case']
        assert good(read(dest/'status.json'),dest,point['bias_V'])
        assert sweep.point(dest,point['bias_V'])['sha256']==point['sha256']
    summary=dict(integrity_pass=True,accepted_states=len(seen),exact_points=len(sweep.ledger['exact_points']),transfers=len(sweep.ledger['transfers']),children=len(sweep.ledger['runs']),rollbacks=len(sweep.ledger['rollbacks']),Newton_updates=sweep.ledger['total_Newton_updates'],child_wall_seconds=sweep.ledger['total_child_wall_seconds'],child_cpu_seconds=sweep.ledger['total_child_cpu_seconds'])
    summary['solver_requests']=len(sweep.ledger['runs'])
    summary['process_count']=1 if args.worker else len(sweep.ledger['runs'])
    if write_output:write(HERE/'audit_summary.json',summary)
    return summary


def validate_resume_contract(plan,ledger,targets,runner_sha,gate,max_step,frame,worker,linear_solver='sparselu',reuse_linear_analysis=False):
    """Only resume an immutable, inactive prefix of the requested experiment."""
    if ledger.get('active_child') is not None or ledger.get('status')=='running':
        raise ValueError('Cannot resume an active solver run')
    if plan.get('reuse_linear_analysis',False)!=reuse_linear_analysis:
        raise ValueError('Cannot change cross-request linear analysis reuse on resume')
    if plan['runner_sha256']!=runner_sha or plan['gate_V']!=gate:
        raise ValueError('Resume runner/gate mismatch')
    if plan.get('linear_solver','sparselu')!=linear_solver:
        raise ValueError('Resume linear backend mismatch')
    if plan['max_step_V']!=max_step or plan['original_blocks']!=GATES:
        raise ValueError('Resume step or numerical gate mismatch')
    if plan.get('execution_mode','subprocess')!=('dc_worker' if worker else 'subprocess'):
        raise ValueError('Resume execution mode mismatch')
    previous=[p['bias_V'] for p in ledger['exact_points']]
    if not previous or previous!=targets[:len(previous)] or len(previous)>=len(targets):
        raise ValueError('Resume must be a proper exact-point prefix')
    accepted=ledger['accepted_bias_V']
    if not previous[-1]<=accepted<targets[len(previous)]:
        raise ValueError('Resume accepted bias is outside the next exact interval')
    if ledger['frame_V']!=0. and ledger['frame_V']!=frame:
        raise ValueError('Cannot change an already accepted reference frame')
    if any(r['status']!='completed' for r in ledger['runs']):
        raise ValueError('Resume contains an incomplete solver request')


def load_resume(source,targets,bundle,manifest):
    global HERE
    source=source.resolve()
    if HERE.is_relative_to(source):
        raise ValueError('Resume output cannot be inside its source')
    plan=read(source/'plan.json');ledger=read(source/'fixed/ledger.json')
    validate_resume_contract(plan,ledger,targets,manifest['runner_sha256'],args.gate,
                             MAXIMUM,FRAME,args.worker,args.linear_solver,args.reuse_linear_analysis)
    for relative,expected in bundle['files'].items():
        if plan['frozen_files'].get(str(ROOT/relative))!=expected:
            raise ValueError(f'Resume physical input differs: {relative}')
    previous_here=HERE
    try:
        HERE=source
        sweep=object.__new__(Sweep);sweep.mode='fixed';sweep.out=source/'fixed';sweep.ledger=ledger
        audit(sweep,plan,write_output=False)
    finally:HERE=previous_here
    return source,plan,ledger

def validate_backend_manifest(manifest, backend):
    # Older qualified manifests describe SparseLU and predate the explicit key.
    if manifest.get('linear_solver', 'sparselu') != backend:
        raise ValueError('Runner manifest linear backend mismatch')


def preflight(bundle_path, workspace, runner_manifest, gate):
    bundle=read(bundle_path)
    for relative,expected in bundle['files'].items():
        path=workspace/relative
        if not path.is_file() or digest(path)!=expected:
            raise ValueError(f'Missing or changed dependency: {path}')
    manifest=read(runner_manifest)
    runner=Path(manifest['runner'])
    if digest(runner)!=manifest['runner_sha256']: raise ValueError('Runner hash mismatch')
    for path,expected in manifest['frozen_sources'].items():
        if digest(Path(path))!=expected: raise ValueError(f'Frozen source changed: {path}')
    return bundle,manifest


def main():
    global args, ROOT, HERE, BASE, RUNNER, EXPECTED, ENV, SEED, REFERENCE, stop, MAXIMUM, WORKER, PHYSICS_PROFILE, FRAME
    WORKER=None
    if not __debug__: raise RuntimeError("Run without Python -O; acceptance assertions are required")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--physics-profile',choices=['D5','D4'],default='D5')
    parser.add_argument('--workspace',type=Path,default=ROOT)
    parser.add_argument('--bundle',type=Path,default=DEFAULT_BUNDLE,
                        help='Hashed input bundle; default uses qualified cm^6/s Auger coefficients')
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--gate',type=int,choices=[4,8],required=True)
    parser.add_argument('--points',type=int,choices=[2,8,31],default=8)
    parser.add_argument('--linear-solver',choices=['sparselu','umfpack','sparselu_metis','umfpack_metis','mumps','mumps_metis','superlu_mt','superlu_mt_metis','strumpack'],default='sparselu')
    parser.add_argument('--max-step',type=float,default=.2)
    parser.add_argument('--frame-offset',type=float,help='Potential gauge offset; defaults to bundle policy or 28 V; original equivalence gates remain enforced')
    parser.add_argument('--resume-from',type=Path,help='Read-only accepted prefix; copy its evidence into a new output and continue')
    parser.add_argument('--preflight',action='store_true')
    parser.add_argument('--worker',action='store_true',help='Reuse one sequential DC worker with prepared input cache')
    parser.add_argument('--reuse-linear-analysis',action='store_true',help='Default-off sequential Newton linear context; requires --worker')
    args=parser.parse_args()
    if args.reuse_linear_analysis and not args.worker:parser.error('--reuse-linear-analysis requires --worker')
    PHYSICS_PROFILE=args.physics_profile
    ROOT=args.workspace.resolve(); HERE=args.output.resolve()
    if HERE.exists(): raise ValueError('Refusing to overwrite an experiment')
    if not math.isfinite(args.max_step) or not .1<=args.max_step<=2.: raise ValueError('Invalid maximum step')
    MAXIMUM=args.max_step
    bundle,manifest=preflight(args.bundle,ROOT,args.manifest,args.gate)
    validate_backend_manifest(manifest,args.linear_solver)
    FRAME=frame_offset_for(bundle,args.gate,args.frame_offset)
    if PHYSICS_PROFILE=='D4' and bundle.get('schema')!='vela.templates_ldmos.linked_d4_inputs.v1':
        raise ValueError('D4 requires its own hashed reference/config bundle')
    BASE=read(ROOT/bundle['template'])
    SEED=ROOT/bundle['seeds'][str(args.gate)]
    REFERENCE=ROOT/bundle['references'][str(args.gate)]
    RUNNER=Path(manifest['runner']);EXPECTED=manifest['runner_sha256']
    ENV=dict(os.environ,VELA_LINEAR_SOLVER=args.linear_solver)
    ENV['PATH']=r'D:\msys64\ucrt64\bin'+os.pathsep+ENV.get('PATH','')
    ENV.pop('GMON_OUT_PREFIX',None)
    targets=read_points(REFERENCE)[:args.points];stop=targets[-1]
    for target in targets: prepare_config(BASE,HERE,SEED,target,.1,args.gate,FRAME if target>=FRAME_PIVOT else 0.)
    resume=load_resume(args.resume_from,targets,bundle,manifest) if args.resume_from else None
    print(f'PREFLIGHT_PASS gate={args.gate} points={len(targets)}',flush=True)
    if args.preflight:return
    HERE.mkdir(parents=True,exist_ok=False)
    frozen={str(ROOT/p):h for p,h in bundle['files'].items()}
    frozen.update(manifest['frozen_sources'])
    if resume:
        origin,origin_plan,origin_ledger=resume
        frozen.update(origin_plan['frozen_files'])
        frozen[str(origin/'plan.json')]=digest(origin/'plan.json')
        frozen[str(origin/'fixed/ledger.json')]=digest(origin/'fixed/ledger.json')
        shutil.copytree(origin/'fixed',HERE/'fixed')
        if (origin/'predictor_inputs').exists():
            shutil.copytree(origin/'predictor_inputs',HERE/'predictor_inputs')
    for path in (Path(__file__),args.bundle.resolve(),args.manifest.resolve(),RUNNER,
                 ROOT/'scripts/dc_worker_client.py',ROOT/'scripts/translate_dd_state.py',ROOT/'scripts/analyze_templates_ldmos_stage4_d5.py',
                 ROOT/'scripts/run_templates_ldmos_stage4_d5.py',ROOT/'scripts/summarize_templates_ldmos_idvd_ablation.py'):
        snapshot=HERE/'source'/path.name
        snapshot.parent.mkdir(exist_ok=True)
        shutil.copy2(path,snapshot)
        frozen[str(snapshot)]=digest(snapshot)
    plan=dict(gate_V=args.gate,targets_V=targets,runner=str(RUNNER),runner_sha256=EXPECTED,
        frozen_files=frozen,seed=str(SEED),max_step_V=MAXIMUM,initial_step_V=INITIAL,
        original_blocks=GATES,local_rows_eps=1e-8,KCL_ratio=1e-8,
        backend=manifest['backend'],reference=str(REFERENCE),frame_pivot_V=FRAME_PIVOT,
        frame_offset_V=FRAME,created_at=stamp(),scope=f'{PHYSICS_PROFILE} linked sweep; original gates')
    plan['execution_mode']='dc_worker' if args.worker else 'subprocess'
    plan['reuse_linear_analysis']=args.reuse_linear_analysis
    plan['linear_solver']=args.linear_solver
    if resume:
        plan['resume_origin']=dict(directory=str(origin),plan_sha256=digest(origin/'plan.json'),
            accepted_bias_V=origin_ledger['accepted_bias_V'],accepted_state=origin_ledger['accepted_state'],
            original_status=origin_ledger['status'],prior_frame_offset_V=origin_plan['frame_offset_V'])
    write(HERE/'plan.json',plan)
    if resume:
        sweep=object.__new__(Sweep);sweep.mode='fixed';sweep.out=HERE/'fixed';sweep.ledger=deepcopy(origin_ledger)
        sweep.ledger.update(status='running',active_child=None)
        sweep.ledger.pop('error',None);sweep.ledger.pop('finished_at',None)
        sweep.ledger.setdefault('resume_history',[]).append(plan['resume_origin']);sweep.save()
    else:sweep=Sweep()
    progress=dict(status='running',pid=os.getpid(),started_at=stamp())
    started=time.perf_counter();write(HERE/'progress.json',progress)
    try:
        if not resume:initialize(sweep)
        elif sweep.ledger['accepted_bias_V']==FRAME_PIVOT and sweep.ledger['frame_V']==0.:
            sweep.change_frame()
        for target in targets[len(sweep.ledger['exact_points']):]:
            progress['target_reference_V']=target;write(HERE/'progress.json',progress)
            sweep.advance(target)
            if target==FRAME_PIVOT:sweep.change_frame()
        verdicts=score(sweep);audit(sweep,plan)
        passed=all(v['pass_all'] for v in verdicts.values()) if args.points==31 else True
        progress.update(status='completed' if passed else 'failed_accuracy',verdicts=verdicts)
        sweep.ledger.update(status=progress['status'],finished_at=stamp())
    except Exception as error:
        progress.update(status='failed',error=str(error));sweep.ledger.update(status='stopped_after_failure',error=str(error))
        raise
    finally:
        if WORKER is not None:
            progress['worker_requests']=WORKER.count
            WORKER.close()
        sweep.save();progress.update(finished_at=stamp(),controller_wall_seconds=time.perf_counter()-started)
        write(HERE/'progress.json',progress)
    if progress['status']!='completed':raise RuntimeError('Original curve accuracy gate failed')

if __name__=='__main__':
    main()
