"""T470p 1/2/4-thread screening and bounded D5 curves; original numerical gates."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import time

from run_templates_ldmos_linked_d5 import read, write, digest, preflight, cpu_times
from run_templates_ldmos_host_pilot import verify_package, host_info
from run_templates_ldmos_analysis_reuse import aggregate, compare
from audit_templates_ldmos_linked_d5 import audit
from run_templates_ldmos_full_curve_timing import stop_tree
from windows_system_load import cpu_snapshot, cpu_interval

ROOT=Path(__file__).resolve().parents[1]
BACKENDS=('sparselu','umfpack','mumps','superlu_mt','strumpack')


def configurations():
    rows=[]
    for level in (1,2,4):
        order=BACKENDS if level!=2 else BACKENDS[::-1]
        for backend in order:
            if backend=='sparselu' and level!=1:continue
            rows.append(dict(key=f'{backend}_t{level}',backend=backend,level=level,
                solver_threads=level if backend in ('mumps','superlu_mt','strumpack') else 1,
                blas_threads=level if backend=='umfpack' else 1,
                omp_threads=level if backend!='sparselu' else 1))
    return rows


def environment(config):
    env=dict(os.environ,VELA_LINEAR_SOLVER=config['backend'],
             VELA_LINEAR_THREADS=str(config['solver_threads']),
             VELA_BLAS_THREADS=str(config['blas_threads']),OPENBLAS_NUM_THREADS=str(config['blas_threads']),
             OMP_NUM_THREADS=str(config['omp_threads']),OMP_DYNAMIC='FALSE',OMP_MAX_ACTIVE_LEVELS='1',
             VELA_LINEAR_FACTOR_STATISTICS='0')
    env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env.get('PATH','')
    for name in ('GMON_OUT_PREFIX','VELA_LINEAR_CAPTURE_DIR'):env.pop(name,None)
    return env


def audit_threads(profiles,config):
    active=[p for p in profiles if p.get('counters',{}).get('linear.solve_calls',0)>0]
    if not active:raise ValueError('No active solve profiles')
    expected={'linear.openblas_threads':config['blas_threads']}
    if config['backend'] in ('mumps','superlu_mt','strumpack'):
        prefix='linear.'+config['backend']+'_omp_'
        expected.update({prefix+'threads':config['solver_threads'],prefix+'max_active_levels':1})
    for profile in active:
        for name,value in expected.items():
            obs=profile.get('observations',{}).get(name,{})
            if obs.get('min')!=value or obs.get('max')!=value or obs.get('count',0)<1:
                raise ValueError('Missing or wrong thread observation: '+name)
    return dict(verified=True,profiles=len(active),expected=expected)


def execute(command,log,env,out,report,save,timeout):
    if (out/'PAUSE').exists():raise KeyboardInterrupt('PAUSE requested')
    samples=[];previous=cpu_snapshot();start=time.perf_counter();proc=None
    try:
        with log.open('x') as handle:
            proc=subprocess.Popen(list(map(str,command)),cwd=ROOT,env=env,stdout=handle,stderr=subprocess.STDOUT)
            report['active_pid']=proc.pid;save()
            while proc.poll() is None:
                try:proc.wait(timeout=10)
                except subprocess.TimeoutExpired:pass
                now=cpu_snapshot();samples.append(cpu_interval(previous,now));previous=now
                report['active_elapsed_s']=time.perf_counter()-start;save()
                if (out/'PAUSE').exists():raise KeyboardInterrupt('PAUSE requested')
                if time.perf_counter()-start>timeout and proc.poll() is None:
                    stop_tree(proc)
                    return dict(status='rejected',reason='timeout',wall_seconds=time.perf_counter()-start,
                                returncode=proc.returncode,cpu_seconds=cpu_times(proc),system_cpu_samples=samples)
        return dict(status='pass' if proc.returncode==0 else 'rejected',returncode=proc.returncode,
                    wall_seconds=time.perf_counter()-start,cpu_seconds=cpu_times(proc),system_cpu_samples=samples)
    except BaseException:
        if proc:stop_tree(proc)
        raise
    finally:
        report.pop('active_pid',None);report.pop('active_elapsed_s',None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if os.name!='nt' or not __debug__:raise RuntimeError('Windows with assertions required')
    package=verify_package(ROOT)
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'manifests').mkdir();(out/'screen').mkdir()
    runner=ROOT/'build-release/vela_example_runner.exe';replay=ROOT/'build-release/linear_solver_replay.exe'
    captures=sorted((ROOT/'captures').rglob('*.bin'))
    if len(captures)!=36:raise ValueError('Expected the frozen 36-system R11 set')
    deps={str(ROOT/p):sha for p,sha in package['files'].items()}
    bundle=ROOT/'reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_no_generation_inputs.json'
    configs=configurations()
    for cfg in configs:
        manifest=out/'manifests'/f"{cfg['key']}.json"
        write(manifest,dict(runner=str(runner),runner_sha256=digest(runner),backend=cfg['backend'],
                           linear_solver=cfg['backend'],frozen_sources=deps,build_type='UCRT64 Release -O3 -DNDEBUG'))
        preflight(bundle,ROOT,manifest,8)
    report=dict(status='running',configurations=configs,host=host_info(),screen=[],cases=[],comparisons={},
                package_sha256=digest(ROOT/'package_manifest.json'),runner_sha256=digest(runner),
                scope='Vg8 D5 first8, reuse, original gates; one curve per qualified configuration',
                screen_rounds=3,screen_systems=36,gate_V=8,points=8)
    lock=out/'batch.lock'
    with lock.open('x') as f:f.write(str(os.getpid()))
    def save():write(out/'summary.json',report)
    rejected=set();paths={}
    try:
        save()
        # Independent processes and serial rounds expose nondeterministic failures.
        for repeat in range(3):
            for cfg in (configs if repeat%2==0 else configs[::-1]):
                if cfg['key'] in rejected:continue
                name=f"r{repeat}_{cfg['key']}";dest=out/'screen'/f'{name}.json'
                report.update(stage='matrix_screen',active=name);save()
                result=execute([replay,cfg['backend'],dest,*captures],out/'screen'/f'{name}.log',
                               environment(cfg),out,report,save,180)
                row=dict(key=cfg['key'],round=repeat,**result)
                if row['status']=='pass':
                    try:
                        data=read(dest)
                        if data.get('status')!='pass' or len(data['systems'])!=36:raise ValueError('Incomplete replay')
                        row['thread_audit']=audit_threads([s['profiling'] for s in data['systems']],cfg)
                        row['worst_scaled_backward_error']=max(s['quality']['scaled']['normwise_backward_error'] for s in data['systems'])
                    except Exception as error:row.update(status='rejected',reason=repr(error))
                if row['status']!='pass':rejected.add(cfg['key'])
                report['screen'].append(row);save()
                print('SCREEN',name,row['status'],round(row['wall_seconds'],3),flush=True)
        for cfg in configs:
            if cfg['key'] in rejected:
                report['cases'].append(dict(**cfg,status='not_run',reason='fixed_matrix_screen_rejected'));save();continue
            if cfg['level']!=1 and cfg['backend']+'_t1' not in paths:
                report['cases'].append(dict(**cfg,status='not_run',reason='single_thread_control_not_qualified'));save();continue
            report.update(stage='curves',active=cfg['key']);save()
            dest=out/cfg['key']
            command=[sys.executable,ROOT/'scripts/run_templates_ldmos_linked_d5.py','--workspace',ROOT,
                     '--bundle',bundle,'--manifest',out/'manifests'/f"{cfg['key']}.json",'--output',dest,
                     '--gate','8','--points','8','--linear-solver',cfg['backend'],'--worker','--reuse-linear-analysis']
            row=dict(**cfg,**execute(command,out/f"{cfg['key']}.log",environment(cfg),out,report,save,900))
            if row['status']=='pass':
                try:
                    row['audit']=audit(dest);ledger=read(dest/'fixed/ledger.json')
                    if len(ledger['exact_points'])!=8:raise ValueError('Incomplete exact points')
                    row['thread_audit']=audit_threads([read(dest/r['case']/'performance_profile.json') for r in ledger['runs']],cfg)
                    row['totals']=aggregate(dest)
                    row['process_tree_cpu_seconds']=row['cpu_seconds']['total']+row['audit']['child_cpu_seconds']
                    if cfg['level']!=1:
                        comparison=compare(paths[cfg['backend']+'_t1'],dest)
                        report['comparisons'][cfg['key']]=comparison
                        if not comparison['equivalent']:raise ValueError('State differs from single-thread control')
                    paths[cfg['key']]=dest
                except Exception as error:row.update(status='rejected',reason=repr(error))
            report['cases'].append(row);save()
            print('CURVE',cfg['key'],row['status'],round(row['wall_seconds'],3),flush=True)
        verify_package(ROOT)
        report['status']='pass' if all(c['status']=='pass' for c in report['cases']) else 'completed_with_rejections'
    except KeyboardInterrupt as error:report.update(status='paused',error=repr(error))
    except BaseException as error:report.update(status='failed',error=repr(error));raise
    finally:
        report.pop('active',None);save();lock.unlink()
    print('THREAD_PILOT_END',report['status'],flush=True)


if __name__=='__main__':main()
