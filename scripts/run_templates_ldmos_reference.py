"""Run the checked-in D5 LDMOS fixture with production analysis reuse defaults."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import time

from run_templates_ldmos_linked_d5 import read, write, digest, preflight
from audit_templates_ldmos_linked_d5 import audit
from analyze_templates_ldmos_stage4_d5 import analyze
from run_templates_ldmos_full_curve_timing import stop_tree

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'reference_tcad/templates_ldmos_sentaurus2022/validation_d5/inputs.json'
BACKENDS=('umfpack','sparselu','strumpack','mumps','superlu_mt')


def verify_fixture(workspace=ROOT, bundle_path=FIXTURE):
    bundle=read(bundle_path)
    for name,expected in bundle['files'].items():
        path=(workspace/name).resolve()
        if not path.is_relative_to(bundle_path.parent.resolve()) or digest(path)!=expected:
            raise ValueError('Missing, external or changed reference input: '+name)
    return bundle


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runner',type=Path,default=ROOT/'build-release/vela_example_runner.exe')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--points',type=int,choices=(8,31),default=8)
    parser.add_argument('--linear-solver',choices=BACKENDS,default='umfpack')
    parser.add_argument('--preflight',action='store_true')
    args=parser.parse_args()
    if not __debug__:raise RuntimeError('Run without Python -O')
    bundle=verify_fixture()
    runner=args.runner.resolve()
    if not runner.is_file():raise ValueError('Build the Release runner first')
    if args.preflight:
        print('REFERENCE_INPUTS_PASS',len(bundle['files']),'files');return
    if os.name!='nt':raise RuntimeError('This linked driver uses Windows process CPU accounting')
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    # Freeze the source actually used, without requiring historical build folders.
    sources={str(p):digest(p) for p in (ROOT/'scripts').rglob('*.py')}
    manifest=out/'runner_manifest.json'
    write(manifest,dict(runner=str(runner),runner_sha256=digest(runner),
        backend=args.linear_solver,linear_solver=args.linear_solver,frozen_sources=sources))
    env=dict(os.environ,VELA_LINEAR_SOLVER=args.linear_solver,VELA_LINEAR_THREADS='1',
             OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_MAX_ACTIVE_LEVELS='1')
    for key in ('GMON_OUT_PREFIX','VELA_LINEAR_CAPTURE_DIR'):env.pop(key,None)
    # Only numerical qualification is claimed here. An opt-in OpenBLAS API
    # setting remains the caller's choice and must not inherit a stale value.
    env.pop('VELA_BLAS_THREADS',None)
    report=dict(status='running',profile='D5',points_per_gate=args.points,
        backend=args.linear_solver,reuse_linear_analysis=True,runner_sha256=digest(runner),cases=[])
    def save():write(out/'summary.json',report)
    save()
    try:
        for gate in (4,8):
            preflight(FIXTURE,ROOT,manifest,gate)
            dest=out/f'vg{gate}'
            command=[sys.executable,ROOT/'scripts/run_templates_ldmos_linked_d5.py',
                '--workspace',ROOT,'--bundle',FIXTURE,'--manifest',manifest,
                '--output',dest,'--gate',str(gate),'--points',str(args.points),
                '--linear-solver',args.linear_solver,'--worker','--reuse-linear-analysis']
            started=time.perf_counter()
            print('REFERENCE_RUNNING',gate,flush=True)
            with (out/f'vg{gate}.log').open('x') as log:
                proc=subprocess.Popen(list(map(str,command)),cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
                try:
                    if proc.wait():raise RuntimeError(f'Vg{gate} curve failed; see log')
                except BaseException:
                    stop_tree(proc);raise
            result=audit(dest)
            report['cases'].append(dict(gate_V=gate,wall_seconds=time.perf_counter()-started,audit=result))
            save();print('REFERENCE_PASS',gate,args.points,flush=True)
        if args.points==31:
            joint=analyze({f'Vg{g}':ROOT/bundle['references'][str(g)] for g in (4,8)},
                {f'Vg{g}':out/f'vg{g}/score/curve.csv' for g in (4,8)},
                {f'Vg{g}':out/f'vg{g}/score/terminal_balance.csv' for g in (4,8)},out/'joint',physics_profile='D5')
            report['joint']=joint
            if any(joint[level]['status']!='pass' for level in ('engineering','final')):
                raise RuntimeError('Original dual-gate curve gates failed')
        report.update(status='pass',qualification='full_dual_gate' if args.points==31 else 'prefix_numerical_only')
    except BaseException as error:
        report.update(status='failed',error=repr(error));raise
    finally:save()


if __name__=='__main__':main()
