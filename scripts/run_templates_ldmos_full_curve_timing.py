"""Serial frozen D5/D4 full-curve timing plan, with explicit stage boundaries."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from run_templates_ldmos_linked_d5 import read, write, digest, preflight, cpu_times
from run_templates_ldmos_analysis_reuse import aggregate, compare
from run_templates_ldmos_isothermal_backends import audit_runtime_threads
from audit_templates_ldmos_linked_d5 import audit
from analyze_templates_ldmos_stage4_d5 import analyze
from windows_system_load import cpu_snapshot, cpu_interval

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = {
    'U0': dict(backend='umfpack', threads=1, reuse=False),
    'U1': dict(backend='umfpack', threads=1, reuse=True),
    'T1': dict(backend='strumpack', threads=1, reuse=True),
    'L1': dict(backend='sparselu', threads=1, reuse=True),
    'M1': dict(backend='mumps', threads=1, reuse=True),
    'S1': dict(backend='superlu_mt', threads=1, reuse=True),
}
STAGES = ['P1', 'P2', 'P3short', 'P3full']


def schedule():
    result = []
    def group(stage, profile, points, repeat, order):
        for gate in (4, 8):
            for config in (order if gate == 4 else order[::-1]):
                result.append(dict(stage=stage, profile=profile, points=points,
                    round=repeat, gate=gate, config=config,
                    key=f'{profile.lower()}_p{points}_r{repeat}_vg{gate}_{config}'))
    group('P1', 'D5', 31, 0, ['U0', 'U1', 'T1'])
    group('P1', 'D5', 31, 0, ['L1', 'M1', 'S1'])
    for repeat in (1, 2):
        order = ['U0', 'U1', 'T1']
        group('P2', 'D5', 31, repeat, order[repeat:] + order[:repeat])
    group('P3short', 'D4', 8, 0, ['U0', 'U1', 'T1'])
    for repeat in range(3):
        order = ['U0', 'U1', 'T1']
        group('P3full', 'D4', 31, repeat, order[repeat:] + order[:repeat])
    return result


def require_joint_pass(result):
    if any(result.get(level, {}).get('status') != 'pass' for level in ('engineering', 'final')):
        raise ValueError('Original dual-gate qualification failed')


def validate_resume(report):
    if report['status'] != 'paused' or report.get('active'):
        raise ValueError('Resume requires an explicitly paused, inactive batch')
    if report['configurations'] != CONFIGS or report['schedule'] != schedule():
        raise ValueError('Frozen configuration or schedule changed')


def bundle_path(profile):
    return ROOT/f'reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_{profile.lower()}_auger_no_generation_inputs.json'


def freeze(out):
    flags = (ROOT/'build-release/build.ninja').read_text()
    for token in ('-O3 -DNDEBUG', 'VELA_HAS_OPENBLAS_THREAD_CONTROL=1',
                  'VELA_HAS_UMFPACK=1', 'VELA_HAS_STRUMPACK=1',
                  'VELA_HAS_MUMPS=1', 'VELA_HAS_SUPERLU_MT=1'):
        if token not in flags: raise ValueError('Missing build feature: '+token)
    if '-pg ' in flags: raise ValueError('Profiling build is not a timing control')
    binary = out/'binary'; binary.mkdir()
    for source in [ROOT/'build-release/vela_example_runner.exe', ROOT/'build-release/build.ninja',
                   ROOT/'build-release/CMakeCache.txt', *(ROOT/'build-release').glob('*.dll')]:
        shutil.copy2(source, binary/source.name)
    frozen, live_scripts = {}, {}
    for base in ('src', 'include', 'scripts'):
        for source in (ROOT/base).rglob('*'):
            if source.is_file() and source.suffix in ('.cpp', '.h', '.py'):
                target = out/'source'/source.relative_to(ROOT)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target); frozen[str(target)] = digest(target)
                if base == 'scripts': live_scripts[str(source)] = digest(source)
    shutil.copy2(ROOT/'CMakeLists.txt', out/'source/CMakeLists.txt')
    frozen[str(out/'source/CMakeLists.txt')] = digest(out/'source/CMakeLists.txt')
    runtime = {str(p): digest(p) for p in binary.iterdir()}
    runner = binary/'vela_example_runner.exe'
    inputs = {}
    for profile in ('D5', 'D4'):
        bundle = bundle_path(profile); inputs[str(bundle)] = digest(bundle)
        for relative, sha in read(bundle)['files'].items(): inputs[str(ROOT/relative)] = sha
    for name, config in CONFIGS.items():
        write(binary/f'{name}.json', dict(runner=str(runner), runner_sha256=digest(runner),
            backend=config['backend'], linear_solver=config['backend'], frozen_sources=frozen,
            runtime_sha256=runtime, build_type='UCRT64 Release -O3 -DNDEBUG'))
    manifests = {str(binary/f'{name}.json'): digest(binary/f'{name}.json') for name in CONFIGS}
    return dict(status='prepared', configurations=CONFIGS, schedule=schedule(), cases=[],
        comparisons={}, repeat_comparisons={}, joint_qualifications={},
        hashes={**frozen, **runtime, **inputs, **manifests}, live_scripts=live_scripts,
        runner_sha256=digest(runner), logical_processors=os.cpu_count(),
        power_scheme=subprocess.check_output(['powercfg', '/getactivescheme']).decode(errors='replace'))


def verify_frozen(report):
    for path, expected in {**report['hashes'], **report['live_scripts']}.items():
        if digest(Path(path)) != expected: raise ValueError('Frozen dependency changed: '+path)


def completed(report):
    return {c['key']: c for c in report['cases'] if c['status'] == 'pass'}


def check_ready(out, report):
    done = completed(report)
    control = report.get('control_config', 'U0')
    for key, case in done.items():
        base = f"{case['profile'].lower()}_p{case['points']}_r{case['round']}_vg{case['gate']}_{control}"
        if case['config'] != control and base in done and key not in report['comparisons']:
            comparison = compare(out/done[base]['name'], out/case['name'])
            report['comparisons'][key] = dict(control=base, **comparison)
            if not comparison['equivalent']: raise ValueError('State comparison failed: '+key)
        first = f"{case['profile'].lower()}_p{case['points']}_r0_vg{case['gate']}_{case['config']}"
        if case['round'] and first in done and key not in report['repeat_comparisons']:
            comparison = compare(out/done[first]['name'], out/case['name'])
            report['repeat_comparisons'][key] = dict(control=first, **comparison)
            if not comparison['equivalent']: raise ValueError('Repeated state mismatch: '+key)
        if case['points'] != 31: continue
        group = f"{case['profile'].lower()}_r{case['round']}_{case['config']}"
        keys = {g: f"{case['profile'].lower()}_p31_r{case['round']}_vg{g}_{case['config']}" for g in (4,8)}
        if group in report['joint_qualifications'] or not all(k in done for k in keys.values()): continue
        bundle = read(bundle_path(case['profile']))
        dirs = {g: out/done[k]['name'] for g,k in keys.items()}
        result = analyze({f'Vg{g}': ROOT/bundle['references'][str(g)] for g in dirs},
            {f'Vg{g}': d/'score/curve.csv' for g,d in dirs.items()},
            {f'Vg{g}': d/'score/terminal_balance.csv' for g,d in dirs.items()},
            out/(group+'_joint'), physics_profile=case['profile'])
        report['joint_qualifications'][group] = result
        require_joint_pass(result)


def stop_tree(proc):
    if proc.poll() is None:
        subprocess.run(['taskkill', '/PID', str(proc.pid), '/T', '/F'], capture_output=True)
        proc.wait(timeout=30)


def run_case(out, report, item, save):
    config = CONFIGS[item['config']]
    attempt = 1 + sum(c['key'] == item['key'] for c in report['cases'])
    name = item['key'] + (f'_attempt{attempt}' if attempt > 1 else '')
    dest = out/name
    command = [sys.executable, out/'source/scripts/run_templates_ldmos_linked_d5.py',
        '--workspace', ROOT, '--physics-profile', item['profile'], '--bundle', bundle_path(item['profile']),
        '--manifest', out/'binary'/f"{item['config']}.json", '--output', dest,
        '--gate', item['gate'], '--points', item['points'], '--linear-solver', config['backend'], '--worker']
    if config['reuse']: command.append('--reuse-linear-analysis')
    env = dict(os.environ, VELA_LINEAR_SOLVER=config['backend'], VELA_LINEAR_THREADS=str(config['threads']),
        OMP_NUM_THREADS=str(config['threads']), VELA_BLAS_THREADS='1', OPENBLAS_NUM_THREADS='1',
        OMP_DYNAMIC='FALSE', OMP_MAX_ACTIVE_LEVELS='1', VELA_LINEAR_FACTOR_STATISTICS='0')
    env['PATH'] = 'D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env.get('PATH', '')
    for key in ('GMON_OUT_PREFIX', 'VELA_LINEAR_CAPTURE_DIR'): env.pop(key, None)
    case = dict(item, name=name, status='running', argv=list(map(str,command)), system_cpu_samples=[])
    report['cases'].append(case); report['active']=name;report['active_elapsed_s']=0.;save()
    proc = None
    start=time.perf_counter()
    try:
        previous=cpu_snapshot(); start=time.perf_counter()
        with (out/(name+'.log')).open('x') as log:
            proc=subprocess.Popen(list(map(str,command)),cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
            case['pid']=proc.pid;save()
            while proc.poll() is None:
                try: proc.wait(timeout=20)
                except subprocess.TimeoutExpired: pass
                current=cpu_snapshot();case['system_cpu_samples'].append(cpu_interval(previous,current));previous=current
                report['active_elapsed_s']=time.perf_counter()-start;save()
                print(name,round(report['active_elapsed_s'],1),flush=True)
                if (out/'PAUSE').exists() and proc.poll() is None: raise KeyboardInterrupt('PAUSE requested')
            case.update(wall_seconds=time.perf_counter()-start,cpu_seconds=cpu_times(proc),returncode=proc.returncode)
        if proc.returncode: raise RuntimeError('Curve process failed: '+name)
        case['audit']=audit(dest)
        if case['audit']['exact_points'] != item['points']: raise ValueError('Incomplete exact points')
        if item['points'] == 31:
            score=read(dest/'score/summary.json')
            if not score['verdicts']['final']['pass_all']: raise ValueError('Original curve score failed')
        case['thread_audit']=audit_runtime_threads(dest,config['backend'],config['threads'])
        case['totals']=aggregate(dest);case['status']='pass';save()
    except BaseException as error:
        if proc: stop_tree(proc)
        case.setdefault('wall_seconds',time.perf_counter()-start)
        if proc:
            case.setdefault('cpu_seconds',cpu_times(proc))
            case.setdefault('returncode',proc.returncode)
        case.update(status='interrupted' if isinstance(error,KeyboardInterrupt) else 'failed',error=repr(error))
        raise
    finally:
        report.pop('active',None);report.pop('active_elapsed_s',None);save()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--through',choices=STAGES,default='P1')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if not __debug__ or os.name != 'nt': raise RuntimeError('Requires Windows and assertions')
    out=args.output.resolve()
    if args.resume:
        report=read(out/'summary.json');validate_resume(report);verify_frozen(report)
        if (out/'PAUSE').exists(): raise ValueError('Remove the PAUSE request before explicit resume')
    else:
        out.mkdir(parents=True,exist_ok=False);report=freeze(out)
    # Exclusive batch ownership; a stale lock is evidence requiring inspection.
    lock=out/'batch.lock'
    with lock.open('x') as handle: handle.write(str(os.getpid()))
    def save(): write(out/'summary.json',report)
    try:
        verify_frozen(report)
        for profile in ('D5','D4'):
            for gate in (4,8): preflight(bundle_path(profile),ROOT,out/'binary/U0.json',gate)
        report.update(status='running',through=args.through);save()
        check_ready(out,report);save()
        for item in report['schedule']:
            if STAGES.index(item['stage']) > STAGES.index(args.through): break
            if item['key'] in completed(report): continue
            if (out/'PAUSE').exists(): raise KeyboardInterrupt('PAUSE requested')
            run_case(out,report,item,save)
            check_ready(out,report);save()
        verify_frozen(report)
        report['status']='pass' if len(completed(report))==len(report['schedule']) else 'paused'
        report['stop_reason']='complete' if report['status']=='pass' else 'stage_boundary'
    except KeyboardInterrupt as error:
        report.update(status='paused',stop_reason='user_interrupt',error=repr(error))
    except BaseException as error:
        report.update(status='failed',error=repr(error));raise
    finally:
        save();lock.unlink()
    print('FULL_CURVE_STAGE_COMPLETE',report['status'],len(completed(report)),flush=True)


if __name__ == '__main__': main()
