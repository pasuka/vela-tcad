"""Default-off D5 cross-request analysis reuse study; serial, frozen, no overwrite."""
import argparse
from collections import Counter
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from run_templates_ldmos_linked_d5 import read, write, digest, preflight, cpu_times, state_difference
from audit_templates_ldmos_linked_d5 import audit
from run_templates_ldmos_isothermal_backends import audit_runtime_threads
from windows_system_load import cpu_snapshot, cpu_interval

ROOT = Path(__file__).resolve().parents[1]
MODES = {'subprocess': [], 'worker': ['--worker'],
         'reuse': ['--worker', '--reuse-linear-analysis']}
BACKENDS = ('sparselu', 'umfpack', 'mumps', 'superlu_mt', 'strumpack')


def validate_options(backend, threads, modes):
    if len(modes)<2 or len(set(modes))!=len(modes) or 'reuse' not in modes:
        raise ValueError('Require distinct control and reuse modes')
    if backend in ('sparselu','umfpack') and threads!=1:
        raise ValueError('SparseLU/UMFPACK controls require one solver thread')
    return next(m for m in modes if m!='reuse')


def aggregate(directory):
    ledger = read(directory/'fixed/ledger.json')
    counters, stages = Counter(), Counter()
    peaks = []
    for run in ledger['runs']:
        profile = read(directory/run['case']/'performance_profile.json')
        counters.update(profile['counters'])
        for stage in profile['stages']:
            stages[stage['name']] += stage['total_ns']/1e9
        peak = profile.get('resources', {}).get('process_peak_working_set_bytes')
        if peak is not None:
            peaks.append(peak)
    return dict(counters=dict(counters), seconds=dict(stages),
                peak_solver_process_bytes=max(peaks) if peaks else None)


def compare(left, right):
    a, b = (read(p/'fixed/ledger.json') for p in (left, right))
    if len(a['exact_points']) != len(b['exact_points']):
        raise ValueError('Exact point count mismatch')
    points = []
    for x, y in zip(a['exact_points'], b['exact_points']):
        if x['bias_V'] != y['bias_V']:
            raise ValueError('Exact voltage mismatch')
        delta = state_difference(Path(x['state']), Path(y['state']))
        points.append(dict(bias_V=x['bias_V'], delta=delta))
    maximum = max(p['delta'][k]['max_absolute'] for p in points for k in ('psi', 'phin', 'phip'))
    return dict(points=points, max_state_difference_V=maximum, equivalent=maximum<=1e-8,
                same_target_sequence=[r['target_V'] for r in a['runs']]==[r['target_V'] for r in b['runs']],
                same_update_sequence=[r['Newton_updates'] for r in a['runs']]==[r['Newton_updates'] for r in b['runs']])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--points', type=int, choices=[2, 8], default=8)
    parser.add_argument('--gates', type=int, nargs='+', choices=[4, 8], default=[4, 8])
    parser.add_argument('--rounds', type=int, default=1)
    parser.add_argument('--backend', choices=BACKENDS, default='strumpack')
    parser.add_argument('--threads', type=int, choices=[1,2,4], default=None)
    parser.add_argument('--modes', choices=list(MODES), nargs='+', default=list(MODES))
    args = parser.parse_args()
    args.threads=args.threads if args.threads is not None else (2 if args.backend=='strumpack' else 1)
    try: control=validate_options(args.backend,args.threads,args.modes)
    except ValueError as error: parser.error(str(error))
    if not __debug__ or os.name != 'nt':
        raise RuntimeError('Requires assertions and Windows CPU accounting')
    if args.rounds < 1 or len(set(args.gates)) != len(args.gates):
        parser.error('Require positive rounds and unique gates')
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binary = out/'binary'; binary.mkdir()
    flags = (ROOT/'build-release/build.ninja').read_text()
    assert '-O3 -DNDEBUG' in flags and '-pg ' not in flags
    assert 'VELA_HAS_OPENBLAS_THREAD_CONTROL=1' in flags
    if args.backend!='sparselu':assert 'VELA_HAS_'+args.backend.upper()+'=1' in flags
    runner = binary/'vela_example_runner.exe'
    for source in [ROOT/'build-release/vela_example_runner.exe',
                   ROOT/'build-release/build.ninja', ROOT/'build-release/CMakeCache.txt',
                   *(ROOT/'build-release').glob('*.dll')]:
        shutil.copy2(source, binary/source.name)
    frozen = {}
    for base in ['src', 'include', 'scripts']:
        for source in (ROOT/base).rglob('*'):
            if source.is_file() and source.suffix in ('.cpp', '.h', '.py'):
                target = out/'source'/source.relative_to(ROOT)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target); frozen[str(target)] = digest(target)
    shutil.copy2(ROOT/'CMakeLists.txt', out/'source/CMakeLists.txt')
    frozen[str(out/'source/CMakeLists.txt')] = digest(out/'source/CMakeLists.txt')
    runtime = {str(p):digest(p) for p in binary.iterdir() if p.is_file()}
    manifest = binary/'manifest.json'
    write(manifest, dict(runner=str(runner), runner_sha256=digest(runner), backend=args.backend,
                         linear_solver=args.backend, frozen_sources=frozen,
                         runtime_sha256=runtime, build_type='UCRT64 Release -O3 -DNDEBUG'))
    bundle = ROOT/'reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_no_generation_inputs.json'
    for gate in args.gates:
        preflight(bundle, ROOT, manifest, gate)
    env = dict(os.environ, VELA_LINEAR_SOLVER=args.backend, VELA_LINEAR_THREADS=str(args.threads),
               OMP_NUM_THREADS=str(args.threads), VELA_BLAS_THREADS='1', OPENBLAS_NUM_THREADS='1',
               OMP_DYNAMIC='FALSE', OMP_MAX_ACTIVE_LEVELS='1', VELA_LINEAR_FACTOR_STATISTICS='0')
    env['PATH'] = 'D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env.get('PATH', '')
    for name in ('GMON_OUT_PREFIX', 'VELA_LINEAR_CAPTURE_DIR'):
        env.pop(name, None)
    report = dict(status='running', points=args.points, rounds=args.rounds,
                  backend=args.backend, solver_threads=args.threads, blas_threads=1,
                  modes=args.modes, control=control, cases=[], comparisons=[],
                  runner_sha256=digest(runner), logical_processors=os.cpu_count())
    def save(): write(out/'summary.json', report)
    save()
    try:
        for repeat in range(args.rounds):
            for gate in args.gates:
                paths = {}
                order = args.modes
                shift = (repeat + args.gates.index(gate)) % len(order)
                for mode in order[shift:]+order[:shift]:
                    name = f'd5_vg{gate}_r{repeat}_{mode}'; dest = out/name
                    command = [sys.executable, out/'source/scripts/run_templates_ldmos_linked_d5.py',
                               '--workspace', ROOT, '--bundle', bundle, '--manifest', manifest,
                               '--output', dest, '--gate', gate, '--points', args.points,
                               '--linear-solver', args.backend, *MODES[mode]]
                    record = dict(name=name, mode=mode, gate=gate, round=repeat, status='running',
                                  argv=list(map(str, command)), system_cpu_samples=[])
                    report['cases'].append(record); report['active']=name
                    report['active_elapsed_s']=0.0; save()
                    start = time.perf_counter(); previous = cpu_snapshot()
                    with (out/(name+'.log')).open('x') as log:
                        proc = subprocess.Popen(list(map(str, command)), cwd=ROOT, env=env,
                                                stdout=log, stderr=subprocess.STDOUT)
                        record['pid']=proc.pid; save()
                        while proc.poll() is None:
                            try: proc.wait(timeout=20)
                            except subprocess.TimeoutExpired: pass
                            current=cpu_snapshot()
                            record['system_cpu_samples'].append(cpu_interval(previous,current));previous=current
                            report['active_elapsed_s']=time.perf_counter()-start;save()
                            print(name, round(report['active_elapsed_s'],1), flush=True)
                        record.update(wall_seconds=time.perf_counter()-start, cpu_seconds=cpu_times(proc),
                                      returncode=proc.returncode)
                    if proc.returncode:
                        record['status']='failed';save()
                        raise RuntimeError(f'Curve failed: {name}')
                    record['audit']=audit(dest)
                    record['thread_audit']=audit_runtime_threads(dest,args.backend,args.threads)
                    record['totals']=aggregate(dest)
                    record['status']='pass';paths[mode]=dest;save()
                for mode in (m for m in args.modes if m!=control):
                    result=compare(paths[control],paths[mode])
                    report['comparisons'].append(dict(gate=gate,round=repeat,control=control,candidate=mode,**result));save()
                    if not result['equivalent']:raise RuntimeError('State equivalence failed: '+mode)
        for path, expected in {**frozen, **runtime}.items():
            if digest(Path(path)) != expected: raise RuntimeError('Frozen file changed: '+path)
        report['status']='pass'
    except BaseException as error:
        report.update(status='failed',error=repr(error));raise
    finally:
        report.pop('active',None);report.pop('active_elapsed_s',None);save()
    print('CROSS_POINT_ANALYSIS_COMPLETE',flush=True)


if __name__ == '__main__': main()
