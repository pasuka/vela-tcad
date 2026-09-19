"""Serial fixed-system screening through the production LinearSolver adapter."""
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

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--backends',nargs='+',required=True)
    p.add_argument('--rounds',type=int,default=3)
    p.add_argument('--threads',type=int,choices=[1,2,4],default=1)
    p.add_argument('--timeout',type=float,default=300)
    p.add_argument('--statistics',choices=['on','off'],default='off')
    p.add_argument('--captures',type=Path,help='Explicit input capture directory; default is historical 14-system set')
    a=p.parse_args()
    if a.rounds<1 or a.timeout<=0 or len(set(a.backends))!=len(a.backends):p.error('Invalid rounds, timeout or duplicate backend')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    tool=out/'linear_solver_replay.exe';shutil.copy2(ROOT/'build-release'/tool.name,tool)
    for dll in (ROOT/'build-release').glob('*.dll'):shutil.copy2(dll,out/dll.name)
    files=sorted((a.captures or ROOT/'reference_staging/templates_ldmos_hotspot_execution_20260910/stage4/captures').rglob('*.bin'))
    if not files or (a.captures is None and len(files)!=14):raise ValueError('Missing frozen systems')
    env=dict(os.environ,OMP_NUM_THREADS=str(a.threads),OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE',
             VELA_LINEAR_THREADS=str(a.threads),VELA_LINEAR_FACTOR_STATISTICS='1' if a.statistics=='on' else '0')
    env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env.get('PATH','')
    env.pop('VELA_LINEAR_CAPTURE_DIR',None)
    record=dict(status='running',threads=a.threads,blas_threads=1,factor_statistics=a.statistics,
                load_scope='recorded_background_load_no_idle_claim',cases=[],
                input_sha256={str(x):digest(x) for x in files},
                binary_sha256={str(x):digest(x) for x in out.iterdir() if x.suffix in ('.dll','.exe')})
    sources=[ROOT/'CMakeLists.txt',Path(__file__),ROOT/'include/vela/solver/LinearSolver.h',
             *sorted((ROOT/'src/solver').glob('*Backend*')),*sorted((ROOT/'src/solver').glob('LinearCapture.*')),ROOT/'src/solver/LinearSolver.cpp',
             ROOT/'src/tools/linear_solver_replay.cpp',ROOT/'src/tools/LinearReplaySystem.h']
    for src in sources:
        dest=out/'source'/src.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
    record['source_sha256']={str(x):digest(x) for x in (out/'source').rglob('*') if x.is_file()}
    rejected=set()
    try:
        for r in range(a.rounds):
            offset=(r//2)%len(a.backends)
            order=a.backends[offset:]+a.backends[:offset]
            if r%2:order.reverse()
            for backend in order:
                if backend in rejected:continue
                name=f'r{r}_{backend}';start=time.perf_counter();snap=cpu_snapshot()
                with (out/(name+'.log')).open('x') as log:
                    proc=subprocess.Popen([str(tool),backend,str(out/(name+'.json')),*map(str,files)],env=env,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                    timed_out=False
                    try:rc=proc.wait(timeout=a.timeout)
                    except subprocess.TimeoutExpired:
                        proc.kill();proc.wait();rc=124;timed_out=True
                row=dict(round=r,backend=backend,returncode=rc,timeout=timed_out,wall_seconds=time.perf_counter()-start,cpu_seconds=cpu_times(proc),
                         system_cpu=cpu_interval(snap,cpu_snapshot()),status='pass' if rc==0 else 'rejected')
                record['cases'].append(row)
                if rc:rejected.add(backend)
                write(out/'summary.json',record)
                print(json.dumps(row),flush=True)
        for key in ['input_sha256','binary_sha256','source_sha256']:
            for path,sha in record[key].items():
                if digest(Path(path))!=sha:raise RuntimeError('Evidence changed: '+path)
        record.update(status='completed_with_rejections' if rejected else 'pass',rejected_backends=sorted(rejected))
    except BaseException as e:
        record.update(status='failed',error=repr(e));raise
    finally:write(out/'summary.json',record)

if __name__=='__main__':main()
