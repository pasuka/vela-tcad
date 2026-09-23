"""Profile frozen D5 baseline full curves; reject changed trajectories or states.

Uses a separately built Release -pg runner. No candidate continuation policy is
installed. DLLs, input files, seeds and numerical gates stay frozen.
"""
import argparse
import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace


def relocate(value,source,dest):
    if isinstance(value,dict):return {k:relocate(v,source,dest) for k,v in value.items()}
    if isinstance(value,list):return [relocate(v,source,dest) for v in value]
    if isinstance(value,str):return value.replace(str(source),str(dest)).replace(source.as_posix(),dest.as_posix())
    return value


def trajectory_signature(ledger):
    return [(r['parent_V'],r['target_V'],r['cap_V'],r['stage'],r['frame_V'],r['Newton_updates'])
            for r in ledger['runs']]


def validate_completed_curve(case,ledger,baseline):
    if case['status']!='completed' or ledger['status']!='completed':
        raise ValueError('Cannot reuse incomplete curve')
    if not case.get('trajectory_exact') or not case['audit']['integrity_pass']:
        raise ValueError('Unqualified completed curve')
    if not all(v['pass_all'] for v in case['verdicts'].values()):
        raise ValueError('Original curve gates failed')
    if len(ledger['exact_points'])!=31 or ledger['exact_points'][-1]['bias_V']!=40:
        raise ValueError('Not a full curve')
    if trajectory_signature(ledger)!=trajectory_signature(baseline):
        raise ValueError('Changed completed trajectory')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('package','runner','baseline','output','gprof'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--reuse-completed-vg4',type=Path,
                   help='Read-only completed Vg4 directory from an earlier interrupted batch')
    a=p.parse_args()
    if not __debug__:raise RuntimeError('Python assertions must remain enabled')
    a.package=a.package.resolve();a.output=a.output.resolve();a.runner=a.runner.resolve()
    a.output.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(a.package/'scripts'))
    import run_templates_ldmos_linked_d5 as linked
    from analyze_templates_ldmos_stage4_d5 import ratio_error,read_vela_curve,LIMITS
    from run_bv_performance_benchmarks import generate_gprof
    runtime=a.output/'runtime';runtime.mkdir()
    runner=runtime/'vela_example_runner.exe';shutil.copy2(a.runner,runner)
    manifest=linked.read(a.package/'package_manifest.json')
    for dll in (a.package/'build-release').glob('*.dll'):
        if linked.digest(dll)!=manifest['files']['build-release/'+dll.name]:raise ValueError('Frozen DLL changed')
        shutil.copy2(dll,runtime/dll.name)
    runtime_hashes={str(f):linked.digest(f) for f in runtime.iterdir()}
    env=dict(os.environ,VELA_LINEAR_SOLVER='umfpack',VELA_LINEAR_THREADS='1',VELA_BLAS_THREADS='1',
             OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_MAX_ACTIVE_LEVELS='1',
             VELA_LINEAR_FACTOR_STATS='0')
    env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env.get('PATH','')
    env.pop('GMON_OUT_PREFIX',None)
    smoke=a.output/'smoke';smoke.mkdir()
    with (smoke/'help.log').open('x') as f:
        subprocess.run([str(runner),'--help'],cwd=smoke,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
    batch=dict(status='running',pid=os.getpid(),started_at=linked.stamp(),cases=[],
        runtime_sha256=runtime_hashes,driver_sha256=linked.digest(Path(__file__)),
        scope='D5 baseline full dual-gate gprof; no v3 policy; frozen UMFPACK DLLs, 1/1 threads')
    def save():linked.write(a.output/'summary.json',batch)
    save()
    try:
        for gate in (4,8):
            if gate==4 and a.reuse_completed_vg4:
                previous=a.reuse_completed_vg4.resolve()
                previous_summary=linked.read(previous.parent/'summary.json')
                case=dict(next(c for c in previous_summary['cases'] if c['gate']==4))
                if Path(case['directory']).resolve()!=previous:raise ValueError('Wrong completed directory')
                ledger=linked.read(previous/'fixed/ledger.json')
                origin=a.baseline/'baseline_vg4'
                baseline=linked.read(origin/'fixed/ledger.json')
                validate_completed_curve(case,ledger,baseline)
                plan=linked.read(previous/'plan.json')
                if plan['baseline_ledger_sha256']!=linked.digest(origin/'fixed/ledger.json'):
                    raise ValueError('Baseline changed')
                for path,expected in plan['frozen_files'].items():
                    if linked.digest(Path(path))!=expected:raise ValueError('Reused input/runtime changed')
                old_runtime={Path(k).name:v for k,v in plan['runtime_sha256'].items()}
                if old_runtime!={Path(k).name:v for k,v in runtime_hashes.items()}:
                    raise ValueError('Reused runtime differs')
                if linked.digest(previous/'gmon.out')!=case['gmon_sha256']:raise ValueError('Reused gmon changed')
                for current,old in zip(ledger['exact_points'],baseline['exact_points'],strict=True):
                    if current['bias_V']!=old['bias_V']:raise ValueError('Changed reused exact target')
                    delta=linked.state_difference(Path(current['state']),Path(old['state']))
                    if max(delta[k]['max_absolute'] for k in ('psi','phin','phip'))>1e-8:
                        raise ValueError('Reused state differs from baseline')
                case.update(reused_from=str(previous),rechecked_at=linked.stamp(),
                            source_summary_sha256=linked.digest(previous.parent/'summary.json'))
                batch['cases'].append(case);save()
                print('PROFILE_REUSED',gate,flush=True)
                continue
            b=importlib.reload(linked)
            b.ROOT=a.package;b.HERE=a.output/f'vg{gate}';b.HERE.mkdir()
            origin=a.baseline/f'baseline_vg{gate}'
            oldplan=b.read(origin/'plan.json');oldledger=b.read(origin/'fixed/ledger.json')
            if oldledger['status']!='completed':raise ValueError('Unqualified baseline')
            for path,expected in oldplan['frozen_files'].items():
                if b.digest(Path(path))!=expected:raise ValueError('Frozen baseline input changed')
            initial=origin/oldledger['runs'][0]['case']
            b.BASE=relocate(b.read(initial/'control.json'),initial,Path('@output'))
            b.SEED=Path(oldplan['seed']);b.REFERENCE=Path(oldplan['reference'])
            b.RUNNER=runner;b.EXPECTED=b.digest(runner);b.FRAME=oldplan['frame_offset_V']
            b.PHYSICS_PROFILE='D5';b.WORKER=None;b.MAXIMUM=.2;b.stop=40.;b.ENV=env
            b.args=SimpleNamespace(gate=gate,worker=True,reuse_linear_analysis=True)
            plan=dict(oldplan,runner=str(runner),runner_sha256=b.EXPECTED,baseline_directory=str(origin),
                      baseline_ledger_sha256=b.digest(origin/'fixed/ledger.json'),runtime_sha256=runtime_hashes,
                      driver_sha256=b.digest(Path(__file__)))
            plan['frozen_files']=dict(oldplan['frozen_files'],**runtime_hashes)
            b.write(b.HERE/'plan.json',plan);shutil.copy2(Path(__file__),b.HERE/'driver_source.py')
            case=dict(gate=gate,status='running',directory=str(b.HERE),started_at=b.stamp())
            batch['cases'].append(case);save()
            sweep=b.Sweep();targets=b.read_points(b.REFERENCE);comparisons=[]
            if len(targets)!=31 or targets[-1]!=40:raise ValueError('Expected full reference sequence')
            def compare():
                index=len(sweep.ledger['exact_points'])-1
                current,old=sweep.ledger['exact_points'][index],oldledger['exact_points'][index]
                if current['bias_V']!=old['bias_V']:raise ValueError('Changed exact target')
                delta=b.state_difference(Path(current['state']),Path(old['state']))
                error=max(delta[k]['max_absolute'] for k in ('psi','phin','phip'))
                comparisons.append(dict(bias_V=current['bias_V'],max_potential_V=error,state_difference=delta))
                b.write(b.HERE/'state_comparison.json',comparisons)
                if error>1e-8:raise ValueError('Instrumented state differs from baseline')
            started=time.perf_counter()
            try:
                b.initialize(sweep);compare()
                for target in targets[1:]:
                    if (a.output/'STOP').exists():raise RuntimeError('Requested administrative stop')
                    sweep.advance(target);compare()
                    case.update(exact_points=len(comparisons),bias_V=target);save()
                    if target==b.FRAME_PIVOT:sweep.change_frame()
                verdicts=b.score(sweep);audit=b.audit(sweep,plan)
                if not all(v['pass_all'] for v in verdicts.values()):raise ValueError('Original curve gates failed')
                if trajectory_signature(sweep.ledger)!=trajectory_signature(oldledger):
                    raise ValueError('Instrumented continuation/Newton trajectory differs from baseline')
                case.update(verdicts=verdicts,audit=audit,trajectory_exact=True,
                            max_potential_difference_V=max(c['max_potential_V'] for c in comparisons))
                sweep.ledger['status']='completed'
            except BaseException as error:
                case.update(status='failed',error=repr(error));sweep.ledger['status']='stopped_after_failure';raise
            finally:
                if b.WORKER is not None:b.WORKER.close()
                case.update(instrumented_wall_seconds=time.perf_counter()-started)
                sweep.save();save()
            case['status']='generating_profile';save()
            generate_gprof(b.HERE,runner,a.gprof.resolve())
            case.update(status='completed',gmon_sha256=b.digest(b.HERE/'gmon.out'),finished_at=b.stamp());save()
            print('PROFILE_DONE',gate,case['audit']['Newton_updates'],flush=True)
        refs={};curves={}
        for gate in (4,8):
            directory=Path(next(c['directory'] for c in batch['cases'] if c['gate']==gate))
            plan=linked.read(directory/'plan.json')
            refs[f'Vg{gate}']=linked.read_curve(Path(plan['reference']))
            curves[f'Vg{gate}']=read_vela_curve(directory/'score/curve.csv')
        ratio=ratio_error(refs,curves)
        checks={k:ratio['endpoint']<=v['gate_ratio'] for k,v in LIMITS.items()}
        if not all(checks.values()):raise ValueError('Dual-gate ratio failed')
        batch.update(status='completed',joint=dict(ratio=ratio,checks=checks))
        for path,expected in runtime_hashes.items():
            if linked.digest(Path(path))!=expected:raise ValueError('Runtime changed')
    except BaseException as error:
        batch.update(status='failed',error=repr(error));raise
    finally:
        batch['finished_at']=linked.stamp();save()


if __name__=='__main__':main()
