"""Isolated serial matrix experiment; no production defaults or physics change."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from run_templates_ldmos_linked_d5 import digest, write, cpu_times
from windows_system_load import cpu_snapshot, cpu_interval

ROOT = Path(__file__).resolve().parents[1]


def configurations():
    rows = []
    for mode in ('sparselu', 'umfpack'):
        for blas in (1, 2, 4):
            rows.append(dict(name=f'{mode}_blas{blas}', mode=mode,
                             eigen_blas=mode == 'sparselu', solver_threads=1, blas_threads=blas))
    rows.append(dict(name='sparselu_control', mode='sparselu', eigen_blas=False,
                     solver_threads=1, blas_threads=1))
    for threads in (1, 2, 4):
        for mode in ('strumpack', 'blr6', 'blr8', 'hss6', 'hss8', 'amalg',
                     'strumpack_amd', 'strumpack_mmd', 'strumpack_and'):
            rows.append(dict(name=f'{mode}_t{threads}', mode=mode, eigen_blas=False,
                             solver_threads=threads, blas_threads=1))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--captures', type=Path, required=True)
    p.add_argument('--rounds', type=int, default=3)
    p.add_argument('--timeout', type=float, default=300)
    p.add_argument('--configs', nargs='+')
    p.add_argument('--cold-analysis', action='store_true', help='Rebuild analysis for every matrix, including identical patterns')
    p.add_argument('--reset-per-directory', action='store_true', help='Reset analysis between captured point directories')
    args = p.parse_args()
    if args.cold_analysis and args.reset_per_directory:p.error('Conflicting analysis policies')
    configs = configurations()
    if args.configs:
        if set(args.configs) - {x['name'] for x in configs}: p.error('Unknown configuration')
        if len(set(args.configs)) != len(args.configs): p.error('Duplicate configuration')
        by_name = {x['name']: x for x in configs}
        configs = [by_name[name] for name in args.configs]
    if args.rounds < 1 or args.timeout <= 0: p.error('Invalid repeat/timeout')
    files = sorted(args.captures.resolve().rglob('*.bin'))
    if not files: p.error('No captures')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    binary = out / 'binary'; binary.mkdir()
    for name in ('linear_solver_acceleration_study.exe', 'linear_solver_acceleration_blas_study.exe'):
        shutil.copy2(ROOT / 'build-release' / name, binary / name)
    for src in (ROOT / 'build-release').glob('*.dll'): shutil.copy2(src, binary / src.name)
    sources = [Path(__file__), ROOT/'src/tools/linear_solver_acceleration_study.cpp',
               ROOT/'src/tools/LinearReplaySystem.h', ROOT/'CMakeLists.txt',
               ROOT/'build-release/build.ninja', ROOT/'build-release/CMakeCache.txt']
    for src in sources:
        dst = out/'source'/src.relative_to(ROOT); dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    frozen = {str(x): digest(x) for x in [*files, *binary.iterdir(),
              *[s for s in (out/'source').rglob('*') if s.is_file()]]}
    report = dict(status='running', configurations=configs, cases=[], sha256=frozen,
                  cold_analysis=args.cold_analysis,
                  reset_per_directory=args.reset_per_directory,
                  scope='solver-input matrices; original error gates; background load recorded')
    rejected = set()
    try:
        for r in range(args.rounds):
            order = configs[r % len(configs):] + configs[:r % len(configs)]
            if r % 2: order = list(reversed(order))
            for c in order:
                if c['name'] in rejected: continue
                name = f"r{r}_{c['name']}"
                env = dict(os.environ, OMP_NUM_THREADS=str(max(c['solver_threads'],c['blas_threads'])),
                           VELA_LINEAR_THREADS=str(c['solver_threads']), VELA_BLAS_THREADS=str(c['blas_threads']),
                           OMP_DYNAMIC='FALSE', OMP_MAX_ACTIVE_LEVELS='1')
                env['VELA_STUDY_COLD_ANALYSIS']='1' if args.cold_analysis else '0'
                env['VELA_STUDY_RESET_PER_DIRECTORY']='1' if args.reset_per_directory else '0'
                exe = binary / ('linear_solver_acceleration_blas_study.exe' if c['eigen_blas'] else 'linear_solver_acceleration_study.exe')
                start = time.perf_counter(); snap = cpu_snapshot()
                with (out/(name+'.log')).open('x') as log:
                    proc = subprocess.Popen([str(exe), c['mode'], str(out/(name+'.json')), *map(str, files)],
                                            cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
                    try: rc = proc.wait(timeout=args.timeout)
                    except subprocess.TimeoutExpired: proc.kill(); proc.wait(); rc = 124
                row = dict(name=name, configuration=c, round=r, returncode=rc,
                           wall_seconds=time.perf_counter()-start, cpu_seconds=cpu_times(proc),
                           system_cpu=cpu_interval(snap, cpu_snapshot()))
                if rc: rejected.add(c['name'])
                report['cases'].append(row); write(out/'summary.json', report)
                print(json.dumps(row), flush=True)
        for path, sha in frozen.items():
            if digest(Path(path)) != sha: raise RuntimeError('Frozen evidence changed: '+path)
        report.update(status='completed_with_rejections' if rejected else 'pass', rejected=sorted(rejected))
    except BaseException as e:
        report.update(status='failed', error=repr(e)); raise
    finally: write(out/'summary.json', report)


if __name__ == '__main__': main()
