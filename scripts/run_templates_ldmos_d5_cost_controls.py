"""Isolated D5 cost-feedback continuation and nonoverlapping parent timers."""
import argparse
import importlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

from run_templates_ldmos_d5_gprof import relocate, trajectory_signature
from analyze_templates_ldmos_d5_gprof import aggregate
from ldmos_cost_continuation import Timings, work_units, next_proposal
from run_templates_ldmos_d5_kernel_controls import configure_candidate
from ldmos_file_efficiency import CompletedFileCache, EfficiencyPolicy
from ldmos_binary_adapter import install_binary
from ldmos_memory_seeds import MemorySeeds, replay_prediction
# Import this request-capable client before adding the frozen package to sys.path.
from dc_worker_client import DCWorker


def install_candidate(b,service_budget=600,efficiency=False):
    class Candidate(b.Sweep):
        def __init__(self):
            super().__init__();self.pending=[];self.barrier=-math.inf
            self.efficiency=EfficiencyPolicy(maximum=b.MAXIMUM) if efficiency else None
        def prediction_guard(self,stage,parent_bias,target):
            reason=super().prediction_guard(stage,parent_bias,target)
            if reason is None and parent_bias<=self.barrier+1e-12:
                return 'recovery_history_reset'
            return reason
        def child(self,parent,parent_bias,target,cap,stage,density=False,frame=None):
            if len(self.ledger['runs'])>=service_budget:raise RuntimeError('Experimental service budget exhausted')
            result,dest=super().child(parent,parent_bias,target,cap,stage,density,frame)
            if stage in ('direct','reclose','density'):
                data=b.read(dest/'performance_profile.json')
                trace=b.rows(dest/'newton_iterations.csv')
                alphas=[float(r['damping']) for r in trace if r['event']=='accepted_iteration']
                recovered=(stage!='direct' or any(int(r.get('carrier_row_recovery_attempted',0)) for r in result['curve']))
                sample=dict(work=work_units(data['counters']),min_alpha=min(alphas,default=1.),
                            recovery=recovered,case=str(dest))
                if self.efficiency is not None:
                    sample['cost_seconds']=result['cpu_seconds']['total']+.38
                self.pending.append(sample);self.ledger['runs'][-1]['cost_sample']=sample
                if recovered:self.barrier=target
                self.save()
            return result,dest
        def feedback(self,mode,proposed,actual,updates,density_used=False):
            work=sum(x['work'] for x in self.pending)
            alpha=min((x['min_alpha'] for x in self.pending),default=1.)
            recovered=density_used or any(x['recovery'] for x in self.pending)
            detail={}
            if self.efficiency is not None:
                cost=sum(x['cost_seconds'] for x in self.pending)
                value,detail=self.efficiency.next(proposed,actual,cost,alpha,recovered,self.ledger['frame_V'])
                reason=detail['reason'];detail['estimated_total_cost_seconds']=cost
            else:value,reason=next_proposal(proposed,actual,work,alpha,recovered,maximum=b.MAXIMUM)
            self.ledger.setdefault('cost_feedback',[]).append(dict(proposed_V=proposed,actual_V=actual,
                work_units=work,min_alpha=alpha,recovery=recovered,next_step_V=value,reason=reason,
                attempts=list(self.pending),efficiency=detail))
            self.pending=[]
            return value,value/proposed
    return Candidate


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('package','runner','baseline','output'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--points',type=int,choices=(8,31),default=8)
    p.add_argument('--rounds',type=int,choices=(1,2,3),default=1)
    p.add_argument('--candidate',choices=('cost','files','efficiency','binary','memory_seeds','hdf5','flatbuffers','npy'),default='cost')
    p.add_argument('--control',choices=('baseline','files','binary'),default='baseline')
    p.add_argument('--candidate-only',action='store_true',help='Full-curve qualification, not paired timing')
    a=p.parse_args()
    public_extensions={'hdf5':'.h5','flatbuffers':'.vfb','npy':'.npy'}
    if a.candidate in public_extensions:
        import dd_state_public as public_codec  # Dependency startup outside paired clocks.
    if not a.candidate_only and a.control==a.candidate:
        p.error('Control and candidate must differ')
    if not __debug__:raise RuntimeError('Assertions must stay enabled')
    a.output=a.output.resolve();a.package=a.package.resolve()
    a.output.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(a.package/'scripts'))
    import run_templates_ldmos_linked_d5 as linked
    from analyze_templates_ldmos_stage4_d5 import ratio_error,read_vela_curve,LIMITS
    runtime=a.output/'runtime';runtime.mkdir()
    runner=runtime/'vela_example_runner.exe';shutil.copy2(a.runner,runner)
    manifest=linked.read(a.package/'package_manifest.json')
    for dll in (a.package/'build-release').glob('*.dll'):
        if linked.digest(dll)!=manifest['files']['build-release/'+dll.name]:raise ValueError('Frozen DLL changed')
        shutil.copy2(dll,runtime/dll.name)
    hashes={str(f):linked.digest(f) for f in runtime.iterdir()}
    env=dict(os.environ,VELA_LINEAR_SOLVER='umfpack',VELA_LINEAR_THREADS='1',VELA_BLAS_THREADS='1',
             OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_MAX_ACTIVE_LEVELS='1',VELA_LINEAR_FACTOR_STATS='0')
    env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env.get('PATH','')
    batch=dict(status='running',pid=os.getpid(),started_at=linked.stamp(),points=a.points,
               planned_curves=(2 if a.candidate_only else 4)*a.rounds,cases=[],joint=[],runtime_sha256=hashes,
               candidate=a.candidate,control=a.control,paired_timing=not a.candidate_only,
               driver_sha256=linked.digest(Path(__file__)))
    def save():linked.write(a.output/'summary.json',batch)
    save()
    try:
        for repeat in range(a.rounds):
            order=[(a.control,4),(a.candidate,4),(a.candidate,8),(a.control,8)]
            if a.candidate_only:order=[(a.candidate,4),(a.candidate,8)]
            if repeat%2:order.reverse()
            for mode,gate in order:
                if (a.output/'STOP').exists():raise RuntimeError('Requested administrative stop')
                b=importlib.reload(linked)
                b.ROOT=a.package;b.HERE=a.output/f'r{repeat}_{mode}_vg{gate}';b.HERE.mkdir()
                origin=a.baseline/f'baseline_vg{gate}'
                oldplan=b.read(origin/'plan.json');oldledger=b.read(origin/'fixed/ledger.json')
                if oldledger['status']!='completed':raise ValueError('Unqualified baseline')
                for path,h in oldplan['frozen_files'].items():
                    if b.digest(Path(path))!=h:raise ValueError('Frozen input changed')
                initial=origin/oldledger['runs'][0]['case']
                b.BASE=configure_candidate(relocate(b.read(initial/'control.json'),initial,Path('@output')),'structure')
                b.BASE['solver']['diagnostic_fermi_node_cache']=False
                b.SEED=Path(oldplan['seed']);b.REFERENCE=Path(oldplan['reference'])
                b.RUNNER=runner;b.EXPECTED=b.digest(runner);b.FRAME=oldplan['frame_offset_V']
                targets=b.read_points(b.REFERENCE)[:a.points]
                adaptive=mode in ('cost','efficiency')
                b.PHYSICS_PROFILE='D5';b.WORKER=None;b.MAXIMUM=.8 if adaptive else .2;b.stop=targets[-1];b.ENV=env
                b.args=SimpleNamespace(gate=gate,worker=True,reuse_linear_analysis=True)
                store=MemorySeeds() if mode=='memory_seeds' else None
                state_actual=install_binary(b,store) if mode in ('binary','memory_seeds') else lambda p:Path(p)
                if mode in public_extensions:
                    state_actual=install_binary(b,extension=public_extensions[mode],codec=public_codec)
                original_prepare=b.prepare_config
                if adaptive:
                    def prepare(*args,**kwargs):
                        cfg=original_prepare(*args,**kwargs)
                        cfg['solver']['quasi_fermi_update_limit_V']=min(.2,args[4])
                        return cfg
                    b.prepare_config=prepare
                raw_rows,raw_digest=b.rows,b.digest
                cache=CompletedFileCache(raw_rows,raw_digest) if mode in ('files','efficiency','binary','memory_seeds',*public_extensions) else None
                if cache is not None:
                    b.rows=lambda p:raw_rows(state_actual(p)) if store is not None and store.contains(state_actual(p)) else cache.rows(state_actual(p))
                    b.digest=lambda p:raw_digest(state_actual(p)) if store is not None and store.contains(state_actual(p)) else cache.digest(state_actual(p))
                timers=Timings()
                for name in ('rows','read','write','digest','prepare_config','good','score'):
                    setattr(b,name,timers.wrap(name,getattr(b,name)))
                original_worker_run=b.DCWorker.run
                def worker_run(worker,config,*args,**kwargs):
                    if store is not None:
                        seed=Path(b.read(Path(config))['sweep']['initial_state_file'])
                        if store.contains(seed):
                            kwargs['initial_state_vds']=store.get(seed)
                            store.stats['worker_memory_requests']+=1
                    result=original_worker_run(worker,config,*args,**kwargs)
                    if cache is not None:cache.seal(Path(config).parent)
                    return result
                b.DCWorker.run=timers.wrap('worker_roundtrip',worker_run)
                plan=dict(oldplan,runner=str(runner),runner_sha256=b.EXPECTED,mode=mode,points=a.points,
                          baseline_directory=str(origin),runtime_sha256=hashes,driver_sha256=batch['driver_sha256'])
                plan.update(max_step_V=b.MAXIMUM,control_policy='efficiency_v2' if mode=='efficiency' else 'cost_feedback_v1' if mode=='cost' else 'frozen_v1',
                            qf_limit_policy='min(step,0.2V)' if adaptive else 'step',
                            cost_proxy='(factorizations+Jacobians+0.25*residuals)/2.25',nominal_work_budget=6.,service_budget=600)
                plan.update(file_cache=cache is not None,cold_final_audit=True,
                            restart_format=mode if mode in public_extensions else 'VDS1' if mode in ('binary','memory_seeds') else 'CSV',
                            memory_predictor_seeds=store is not None,solver_states_persisted=True,
                            efficiency_cost='child CPU seconds + 0.38 seconds per service',
                            efficiency_regression_ratio=1.2,efficiency_step_ratio=1.1)
                plan['frozen_files']=dict(oldplan['frozen_files'],**hashes)
                for name in ('ldmos_cost_continuation.py','ldmos_file_efficiency.py','run_templates_ldmos_d5_kernel_controls.py',
                             'run_templates_ldmos_d5_gprof.py','analyze_templates_ldmos_d5_gprof.py',
                             'dd_state_binary.py','ldmos_binary_adapter.py','ldmos_memory_seeds.py','dc_worker_client.py'):
                    helper=Path(__file__).parent/name;plan['frozen_files'][str(helper)]=b.digest(helper)
                if a.candidate in public_extensions:
                    for name in ('dd_state_public.py','generated/vela_state/DDState.py'):
                        helper=Path(__file__).parent/name;plan['frozen_files'][str(helper)]=b.digest(helper)
                b.write(b.HERE/'plan.json',plan)
                shutil.copy2(Path(__file__),b.HERE/'driver_source.py')
                case=dict(mode=mode,gate=gate,repeat=repeat,status='running',directory=str(b.HERE))
                batch['cases'].append(case);save()
                sweep=install_candidate(b,efficiency=mode=='efficiency')() if adaptive else b.Sweep();comparisons=[]
                if adaptive:b.next_step=sweep.feedback
                sweep.point=timers.wrap('point_validation',sweep.point)
                sweep.save=timers.wrap('ledger_save',sweep.save)
                def compare():
                    index=len(sweep.ledger['exact_points'])-1
                    current,old=sweep.ledger['exact_points'][index],oldledger['exact_points'][index]
                    if current['bias_V']!=old['bias_V']:raise ValueError('Changed exact target')
                    with timers.stage('state_comparison'):
                        delta=b.state_difference(Path(current['state']),Path(old['state']))
                    error=max(delta[k]['max_absolute'] for k in ('psi','phin','phip'))
                    comparisons.append(dict(bias_V=current['bias_V'],max_potential_V=error,state_difference=delta))
                    b.write(b.HERE/'state_comparison.json',comparisons)
                    if not math.isfinite(error) or error>1e-8:raise ValueError('State differs from frozen reference')
                timers.data.clear()
                started=time.perf_counter()
                timing_scope=timers.stage('controller');timing_scope.__enter__()
                try:
                    b.initialize(sweep);compare()
                    for target in targets[1:]:
                        if (a.output/'STOP').exists():raise RuntimeError('Requested administrative stop')
                        sweep.advance(target);compare()
                        case.update(exact_points=len(comparisons),bias_V=target);save()
                        if target==b.FRAME_PIVOT:sweep.change_frame()
                    if cache is not None:cache.disable()
                    # Recheck every exact state from disk, independent of cached parses.
                    with timers.stage('cold_state_audit'):
                        for current,old in zip(sweep.ledger['exact_points'],oldledger['exact_points']):
                            cold=b.state_difference(Path(current['state']),Path(old['state']))
                            if max(cold[k]['max_absolute'] for k in ('psi','phin','phip'))>1e-8:
                                raise ValueError('Cold state comparison failed')
                    verdicts=b.score(sweep)
                    if a.points==31 and not all(v['pass_all'] for v in verdicts.values()):raise ValueError('Original curve gates failed')
                    oldruns=oldledger['runs'][:len(sweep.ledger['runs'])]
                    same_trajectory=trajectory_signature(sweep.ledger)==trajectory_signature(dict(runs=oldruns))
                    if not adaptive and not same_trajectory:
                        raise ValueError('Baseline trajectory changed')
                    for run in sweep.ledger['runs']:
                        dest=b.HERE/run['case']
                        expected=b.prepare_config(b.BASE,dest,Path(run['parent_state']),run['target_V'],run['cap_V'],gate,run['frame_V'])
                        if b.read(dest/'control.json')!=expected:raise ValueError('Configuration changed')
                        status=b.read(dest/'status.json')
                        if status['config_sha256']!=b.digest(dest/'control.json') or status['runner_sha256']!=b.EXPECTED:
                            raise ValueError('Provenance changed')
                        if status['seed_sha256']!=run['parent_sha256'] or b.digest(Path(run['parent_state']))!=run['parent_sha256']:
                            raise ValueError('Parent changed')
                    for path,h in plan['frozen_files'].items():
                        if b.digest(Path(path))!=h:raise ValueError('Input/runtime changed')
                    points=sweep.ledger['exact_points']+sweep.ledger['transfers']+[e['accepted'] for e in sweep.ledger['frame_events'] if e['status']=='qualified']
                    for point in {x['state']:x for x in points}.values():
                        dest=b.HERE/point['case']
                        if not b.good(b.read(dest/'status.json'),dest,point['bias_V']):raise ValueError('Transferred state failed gates')
                        if sweep.point(dest,point['bias_V'])['sha256']!=point['sha256']:raise ValueError('State changed')
                    if store is not None:
                        # All predictor parents are persisted solver states. Rebuild
                        # each seed after disabling file caches; no disk-only claim
                        # is made for the intentionally memory-resident seeds.
                        for run in sweep.ledger['runs']:
                            metadata=run.get('outer_predictor')
                            if metadata:
                                replay_prediction(metadata,b.rows,b.digest)
                                store.stats['cold_parent_replays']+=1
                        evidence=store.audit()
                        b.write(b.HERE/'memory_seed_audit.json',evidence)
                        case['memory_seed_stats']={k:v for k,v in evidence.items() if k!='seeds'}
                    audit=dict(integrity_pass=True,exact_points=len(comparisons),Newton_updates=sweep.ledger['total_Newton_updates'],
                               services=len(sweep.ledger['runs']),advances=len(sweep.ledger['transfers']),rollbacks=len(sweep.ledger['rollbacks']))
                    b.write(b.HERE/'audit_summary.json',audit)
                    case.update(verdicts=verdicts,audit=audit,trajectory_exact=same_trajectory,
                                max_potential_difference_V=max(c['max_potential_V'] for c in comparisons),
                                file_cache_stats=dict(cache.stats) if cache is not None else {},cold_final_audit=True)
                    sweep.ledger['status']='completed'
                except BaseException as e:
                    if store is not None:store.flush()
                    case.update(status='failed',error=repr(e));sweep.ledger['status']='stopped_after_failure';raise
                finally:
                    if b.WORKER is not None:b.WORKER.close()
                    case.update(wall_seconds=time.perf_counter()-started)
                    sweep.save();save()
                    timing_scope.__exit__(None,None,None)
                    b.DCWorker.run=original_worker_run
                    case['parent_timing']=json.loads(json.dumps(timers.data))
                    b.write(b.HERE/'parent_timing.json',case['parent_timing'])
                case['performance']=aggregate(b.HERE)
                c=case['performance']['counters']
                if c.get('jacobian.structure_cache_hits',0)==0:raise ValueError('Structure reuse did not execute')
                case.update(status='completed',finished_at=b.stamp());save()
                print('CURVE_DONE',repeat,mode,gate,case['wall_seconds'],flush=True)
            if a.points==31:
                for mode in ([a.candidate] if a.candidate_only else [a.control,a.candidate]):
                    refs={};curves={}
                    for gate in (4,8):
                        directory=a.output/f'r{repeat}_{mode}_vg{gate}'
                        plan=linked.read(directory/'plan.json')
                        refs[f'Vg{gate}']=linked.read_curve(Path(plan['reference']))
                        curves[f'Vg{gate}']=read_vela_curve(directory/'score/curve.csv')
                    ratio=ratio_error(refs,curves)
                    checks={k:ratio['endpoint']<=v['gate_ratio'] for k,v in LIMITS.items()}
                    if not all(checks.values()):raise ValueError('Dual-gate ratio failed')
                    batch['joint'].append(dict(repeat=repeat,mode=mode,ratio=ratio,checks=checks))
            save()
        batch['status']='completed'
    except BaseException as e:
        batch.update(status='failed',error=repr(e));raise
    finally:
        batch['finished_at']=linked.stamp();save()


if __name__=='__main__':main()
