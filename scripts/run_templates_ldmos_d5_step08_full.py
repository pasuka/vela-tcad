"""Isolated paired full curves: protected 0.8 V continuation / fixed 0.2 V QF cap.

Frozen D5 executable and input package are read-only. No production defaults change.
"""
import argparse
from copy import deepcopy
import importlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace
from run_templates_ldmos_d5_targeted_controls import secant_rows, relocate, RUNNER_HASH

LOW_END=28/3


def qf_limit_for(target,cap,policy,predicted=True):
    if policy=='uniform_qf':return .2
    if policy not in ('high_only_qf','protected_qf'):raise ValueError('Unknown QF policy')
    if target<=LOW_END+1e-10:return cap
    if policy=='protected_qf' and not predicted:return min(cap,.2)
    return .2


def select_history(points,parent,target,frame,barrier):
    """Prefer a recent exact state, else closest span; never use future/recovered-frame history."""
    step=target-parent
    if step<=0:return None
    valid=[p for p in points if p['frame_V']==frame and barrier-1e-10<=p['bias_V']<parent-1e-10
           and step/2-1e-10<=parent-p['bias_V']<=4/3+1e-10]
    if not valid:return None
    return min(valid,key=lambda p:(not p.get('exact',False),abs(parent-p['bias_V']-step),-p['bias_V']))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--policy',choices=['uniform_qf','high_only_qf','protected_qf'],default='protected_qf')
    a=p.parse_args();a.package=a.package.resolve();a.output=a.output.resolve()
    a.output.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(a.package/'scripts'))
    import run_templates_ldmos_linked_d5 as linked
    from analyze_templates_ldmos_stage4_d5 import ratio_error,read_vela_curve,LIMITS
    batch=dict(status='running',pid=os.getpid(),started_at=linked.stamp(),cases=[],planned_curves=4,
        candidate=dict(low_max_step_V=.2,high_max_step_V=.8,transition_V=LOW_END,qf_limit_V=.2,
                       qf_policy=a.policy,max_history_span_V=4/3,max_extrapolation_ratio=2),runner_sha256=RUNNER_HASH)
    def save():linked.write(a.output/'summary.json',batch)
    save()

    def run_curve(mode,gate):
        b=importlib.reload(linked)
        b.ROOT=a.package;b.HERE=a.output/f'{mode}_vg{gate}';b.HERE.mkdir()
        origin=a.package/f'results/t470p/d5_p31_r0_vg{gate}_U1'
        original_plan=b.read(origin/'plan.json');original_ledger=b.read(origin/'fixed/ledger.json')
        if b.digest(Path(original_plan['runner']))!=RUNNER_HASH:raise ValueError('Executable changed')
        for path,expected in original_plan['frozen_files'].items():
            if b.digest(Path(path))!=expected:raise ValueError('Frozen input changed: '+path)
        initial=origin/original_ledger['runs'][0]['case']
        b.BASE=relocate(b.read(initial/'control.json'),initial,Path('@output'))
        b.SEED=Path(original_plan['seed']);b.REFERENCE=Path(original_plan['reference'])
        b.RUNNER=Path(original_plan['runner']);b.EXPECTED=RUNNER_HASH
        b.FRAME=original_plan['frame_offset_V'];b.PHYSICS_PROFILE='D5';b.WORKER=None
        b.args=SimpleNamespace(gate=gate,worker=True,reuse_linear_analysis=True)
        b.MAXIMUM=.2;b.stop=40.
        b.ENV=dict(os.environ,VELA_LINEAR_SOLVER='umfpack',VELA_LINEAR_THREADS='1',VELA_BLAS_THREADS='1',
            OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_MAX_ACTIVE_LEVELS='1',
            VELA_LINEAR_FACTOR_STATS='0')
        b.ENV['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+b.ENV.get('PATH','')
        original_prepare=b.prepare_config
        if mode=='candidate':
            def prepare(*args,**kwargs):
                cfg=original_prepare(*args,**kwargs)
                predicted=Path(args[2]).parent==b.HERE/'predictor_inputs'
                cfg['solver']['quasi_fermi_update_limit_V']=qf_limit_for(args[3],args[4],a.policy,predicted)
                return cfg
            b.prepare_config=prepare

        class Candidate(b.Sweep):
            recovery_barrier=-math.inf
            def history(self,parent,target):
                frame=self.ledger['frame_V']
                points={p['state']:dict(p,exact=False) for p in self.ledger['transfers']}
                for event in self.ledger['frame_events']:
                    if event['status']=='qualified':
                        q=event['accepted'];points[q['state']]=dict(q,exact=True)
                for q in self.ledger['exact_points']:points[q['state']]=dict(q,exact=True)
                return select_history(list(points.values()),parent,target,frame,self.recovery_barrier)
            def prediction_guard(self,stage,parent,target):
                if parent<LOW_END-1e-10:return super().prediction_guard(stage,parent,target)
                if stage!='direct':return 'non_direct'
                if any(abs(r['parent_V']-parent)<1e-12 for r in self.ledger['rollbacks']):return 'retry'
                return None if self.history(parent,target) else 'no_safe_same_frame_history'
            def advance(self,reference):
                b.MAXIMUM=.2 if reference<=LOW_END+1e-10 else .8
                super().advance(reference)
            def child(self,parent,parent_bias,target,cap,stage,density=False,frame=None):
                if parent_bias<LOW_END-1e-10:
                    result,folder=super().child(parent,parent_bias,target,cap,stage,density,frame)
                else:
                    guard=self.prediction_guard(stage,parent_bias,target);metadata=None
                    if guard is None:
                        previous=self.history(parent_bias,target)
                        path=Path(previous['state'])
                        if b.digest(path)!=previous['sha256']:raise ValueError('History changed')
                        ratio=(target-parent_bias)/(parent_bias-previous['bias_V'])
                        seed=b.HERE/'predictor_inputs'/f"seed_{len(self.ledger['runs']):05d}.csv"
                        seed.parent.mkdir(exist_ok=True)
                        b.csv_out(seed,secant_rows(b.rows(parent),b.rows(path),ratio))
                        metadata=dict(mode='protected_history_step08_v1',ratio=ratio,
                            previous_V=previous['bias_V'],previous_state=str(path),previous_sha256=b.digest(path),
                            current_state=str(parent),current_sha256=b.digest(parent),
                            predicted_state=str(seed),predicted_sha256=b.digest(seed),
                            frame_V=self.ledger['frame_V'],recovery_barrier_V=None if not math.isfinite(self.recovery_barrier) else self.recovery_barrier)
                        parent=seed
                    result,folder=b.FrameSweep.child(self,parent,parent_bias,target,cap,stage,density,frame)
                    self.ledger['runs'][-1].update(outer_predictor=metadata,predictor_guard_reason=guard)
                if b.good(result,folder,target) and any(int(r.get('carrier_row_recovery_attempted',0)) for r in result['curve']):
                    self.recovery_barrier=target
                    self.ledger.setdefault('history_resets',[]).append(dict(target_V=target,case=str(folder),reason='internal_density_recovery'))
                self.save()
                return result,folder

        plan=dict(original_plan,mode=mode,candidate=batch['candidate'] if mode=='candidate' else None,
                  source_curve=str(origin),original_plan_sha256=b.digest(origin/'plan.json'),
                  driver_sha256=b.digest(Path(__file__)),max_step_V=.8 if mode=='candidate' else .2)
        b.write(b.HERE/'plan.json',plan)
        shutil.copy2(Path(__file__),b.HERE/'driver_source.py')
        case=dict(mode=mode,gate=gate,directory=str(b.HERE),status='running',started_at=b.stamp())
        batch['cases'].append(case);save()
        sweep=Candidate() if mode=='candidate' else b.Sweep()
        targets=b.read_points(b.REFERENCE)
        if len(targets)!=31 or targets[-1]!=40:raise ValueError('Expected full 31-point curve')
        comparisons=[];started=time.perf_counter()
        def compare_latest():
            index=len(sweep.ledger['exact_points'])-1
            current=sweep.ledger['exact_points'][index];ref=original_ledger['exact_points'][index]
            if current['bias_V']!=ref['bias_V']:raise ValueError('Reference voltage differs')
            delta=b.state_difference(Path(current['state']),Path(ref['state']))
            potential=max(delta[k]['max_absolute'] for k in ('psi','phin','phip'))
            row=dict(bias_V=current['bias_V'],max_potential_difference_V=potential,state_difference=delta,
                     qualified=potential<=1e-8)
            comparisons.append(row);b.write(b.HERE/'state_comparison.json',comparisons)
            if not row['qualified']:raise RuntimeError('Original 1e-8 V state equivalence failed')
        try:
            b.initialize(sweep);compare_latest()
            for target in targets[1:]:
                if (a.output/'STOP').exists():raise RuntimeError('Requested administrative stop')
                sweep.advance(target);compare_latest()
                case.update(exact_points=len(sweep.ledger['exact_points']),bias_V=target);save()
                if len(comparisons)==8:print('FIRST8_PASS',mode,gate,flush=True)
                if target==b.FRAME_PIVOT:sweep.change_frame()
            verdicts=b.score(sweep)
            if not all(x['pass_all'] for x in verdicts.values()):raise RuntimeError('Original curve gates failed')
            # Full configuration equality, lineage and original gates replace the old cap==step audit.
            recovery_points=0
            for r in sweep.ledger['runs']:
                dest=b.HERE/r['case'];status=b.read(dest/'status.json');cfg=b.read(dest/'control.json')
                expected=b.prepare_config(b.BASE,dest,Path(r['parent_state']),r['target_V'],r['cap_V'],gate,r['frame_V'])
                if cfg!=expected or status['config_sha256']!=b.digest(dest/'control.json'):raise ValueError('Configuration changed')
                if status['seed_sha256']!=r['parent_sha256'] or b.digest(Path(r['parent_state']))!=r['parent_sha256']:raise ValueError('Seed changed')
                b.assert_linked_physics_contract(cfg)
                recovery_points+=sum(int(x.get('carrier_row_recovery_attempted',0)) for x in status['curve'])
            for path,expected in original_plan['frozen_files'].items():
                if b.digest(Path(path))!=expected:raise ValueError('Frozen evidence changed')
            case.update(status='completed',verdicts=verdicts,advances=len(sweep.ledger['transfers']),
                updates=sweep.ledger['total_Newton_updates'],services=len(sweep.ledger['runs']),
                rollbacks=len(sweep.ledger['rollbacks']),internal_recovery_points=recovery_points,
                max_potential_difference_V=max(r['max_potential_difference_V'] for r in comparisons))
            sweep.ledger['status']='completed'
        except Exception as error:
            case.update(status='failed',error=str(error));sweep.ledger.update(status='stopped_after_failure',error=str(error))
            raise
        finally:
            if b.WORKER is not None:b.WORKER.close()
            case.update(wall_seconds=time.perf_counter()-started,finished_at=b.stamp())
            sweep.save();save()
        print('CURVE_DONE',mode,gate,case['updates'],case['wall_seconds'],flush=True)

    try:
        # Candidate gate4 first; independent from-zero controls in an alternating serial order.
        for mode,gate in (('candidate',4),('baseline',4),('baseline',8),('candidate',8)):
            run_curve(mode,gate)
        batch['joint']={}
        for mode in ('candidate','baseline'):
            references={};curves={}
            for gate in (4,8):
                plan=linked.read(a.output/f'{mode}_vg{gate}'/'plan.json')
                references[f'Vg{gate}']=linked.read_curve(Path(plan['reference']))
                curves[f'Vg{gate}']=read_vela_curve(a.output/f'{mode}_vg{gate}'/'score/curve.csv')
            ratio=ratio_error(references,curves)
            checks={level:ratio['endpoint']<=limits['gate_ratio'] for level,limits in LIMITS.items()}
            batch['joint'][mode]=dict(ratio_error_percent=ratio,checks=checks,pass_all=all(checks.values()))
            if not all(checks.values()):raise RuntimeError('Original dual-gate ratio gate failed')
        batch['status']='completed'
    except Exception as error:
        batch.update(status='failed',error=str(error));raise
    finally:
        batch['finished_at']=linked.stamp();save()


if __name__=='__main__':main()
