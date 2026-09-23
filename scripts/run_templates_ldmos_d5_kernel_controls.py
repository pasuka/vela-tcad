"""Serial, same-executable D5 kernel controls with frozen physics and trajectories."""
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


def configure_candidate(base, mode):
    base = json.loads(json.dumps(base))
    if mode not in ('baseline', 'structure', 'fermi', 'combined'):
        raise ValueError('Unknown candidate')
    base['solver']['diagnostic_reuse_jacobian_structure'] = mode in ('structure','combined')
    if mode in ('fermi','combined'):
        base['solver']['diagnostic_fermi_node_cache'] = True
    return base


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('package','runner','baseline','output'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--points',type=int,choices=(8,31),default=8)
    p.add_argument('--rounds',type=int,choices=(1,2,3),default=1)
    p.add_argument('--candidate',choices=('structure','fermi','combined'),default='structure')
    p.add_argument('--control',choices=('baseline','structure'),default='baseline')
    p.add_argument('--candidate-only',action='store_true',help='Full-curve qualification, not paired timing')
    a=p.parse_args()
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
                b.BASE=configure_candidate(relocate(b.read(initial/'control.json'),initial,Path('@output')),mode)
                b.SEED=Path(oldplan['seed']);b.REFERENCE=Path(oldplan['reference'])
                b.RUNNER=runner;b.EXPECTED=b.digest(runner);b.FRAME=oldplan['frame_offset_V']
                targets=b.read_points(b.REFERENCE)[:a.points]
                b.PHYSICS_PROFILE='D5';b.WORKER=None;b.MAXIMUM=.2;b.stop=targets[-1];b.ENV=env
                b.args=SimpleNamespace(gate=gate,worker=True,reuse_linear_analysis=True)
                plan=dict(oldplan,runner=str(runner),runner_sha256=b.EXPECTED,mode=mode,points=a.points,
                          baseline_directory=str(origin),runtime_sha256=hashes,driver_sha256=batch['driver_sha256'])
                plan['frozen_files']=dict(oldplan['frozen_files'],**hashes)
                b.write(b.HERE/'plan.json',plan)
                shutil.copy2(Path(__file__),b.HERE/'driver_source.py')
                case=dict(mode=mode,gate=gate,repeat=repeat,status='running',directory=str(b.HERE))
                batch['cases'].append(case);save()
                sweep=b.Sweep();comparisons=[]
                def compare():
                    index=len(sweep.ledger['exact_points'])-1
                    current,old=sweep.ledger['exact_points'][index],oldledger['exact_points'][index]
                    if current['bias_V']!=old['bias_V']:raise ValueError('Changed exact target')
                    delta=b.state_difference(Path(current['state']),Path(old['state']))
                    error=max(delta[k]['max_absolute'] for k in ('psi','phin','phip'))
                    comparisons.append(dict(bias_V=current['bias_V'],max_potential_V=error,state_difference=delta))
                    b.write(b.HERE/'state_comparison.json',comparisons)
                    if not math.isfinite(error) or error>1e-8:raise ValueError('State differs from frozen reference')
                started=time.perf_counter()
                try:
                    b.initialize(sweep);compare()
                    for target in targets[1:]:
                        if (a.output/'STOP').exists():raise RuntimeError('Requested administrative stop')
                        sweep.advance(target);compare()
                        case.update(exact_points=len(comparisons),bias_V=target);save()
                        if target==b.FRAME_PIVOT:sweep.change_frame()
                    verdicts=b.score(sweep);audit=b.audit(sweep,plan)
                    if a.points==31 and not all(v['pass_all'] for v in verdicts.values()):raise ValueError('Original curve gates failed')
                    oldruns=oldledger['runs'][:len(sweep.ledger['runs'])]
                    if trajectory_signature(sweep.ledger)!=trajectory_signature(dict(runs=oldruns)):
                        raise ValueError('Kernel candidate changed continuation/Newton trajectory')
                    for run in sweep.ledger['runs']:
                        dest=b.HERE/run['case']
                        expected=b.prepare_config(b.BASE,dest,Path(run['parent_state']),run['target_V'],run['cap_V'],gate,run['frame_V'])
                        if b.read(dest/'control.json')!=expected:raise ValueError('Configuration changed')
                    case.update(verdicts=verdicts,audit=audit,trajectory_exact=True,
                                max_potential_difference_V=max(c['max_potential_V'] for c in comparisons))
                    sweep.ledger['status']='completed'
                except BaseException as e:
                    case.update(status='failed',error=repr(e));sweep.ledger['status']='stopped_after_failure';raise
                finally:
                    if b.WORKER is not None:b.WORKER.close()
                    case.update(wall_seconds=time.perf_counter()-started)
                    sweep.save();save()
                case['performance']=aggregate(b.HERE)
                c=case['performance']['counters']
                if mode in ('structure','combined') and c.get('jacobian.structure_cache_hits',0)==0:raise ValueError('Structure reuse did not execute')
                if mode in ('fermi','combined') and c.get('jacobian.fermi_node_cache_assemblies',0)==0:raise ValueError('Fermi node reuse did not execute')
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
