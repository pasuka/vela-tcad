"""Serial independent matrix strategy controls; run separately from curve timing."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from run_templates_ldmos_linked_d5 import digest,write
from windows_system_load import wait_for_idle

ROOT=Path(__file__).resolve().parents[1]
MODES=['eigen_colamd','eigen_amd','default','unsym','noscale','pivot1',
       'given_colamd','given_amd','given_colamd_noscale','given_colamd_unsym','given_colamd_pivot1']

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--rounds',type=int,default=3)
    p.add_argument('--modes',nargs='+',choices=MODES,default=MODES)
    p.add_argument('--idle-cpu-percent',type=float,default=10.);a=p.parse_args()
    if a.rounds<2:p.error('At least two rounds required')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    tool=out/'linear_solver_strategy_study.exe';shutil.copy2(ROOT/'build-release'/tool.name,tool)
    files=sorted((ROOT/'reference_staging/templates_ldmos_hotspot_execution_20260910/stage4/captures').rglob('*.bin'))
    if len(files)!=14:raise ValueError('Expected 14 frozen systems')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env.get('PATH','')
    record=dict(status='running',tool_sha256=digest(tool),input_sha256={str(x):digest(x) for x in files},cases=[])
    sources=[ROOT/'src/tools/linear_solver_strategy_study.cpp',ROOT/'src/tools/LinearReplaySystem.h',Path(__file__)]
    for src in sources:shutil.copy2(src,out/src.name)
    record['source_sha256']={str(out/x.name):digest(out/x.name) for x in sources}
    record['runtime_sha256']={str(x):digest(x) for x in (Path('D:/msys64/ucrt64/bin')/n for n in ['libumfpack.dll','libopenblas.dll','libwinpthread-1.dll','libgcc_s_seh-1.dll','libstdc++-6.dll'])}
    write(out/'summary.json',record)
    try:
        record['idle_cpu_percent']=a.idle_cpu_percent
        record['idle_samples']=wait_for_idle(a.idle_cpu_percent)
        rejected=set()
        for r in range(a.rounds):
            # Reverse then rotate to avoid consistently favoring a late/early mode.
            selected=a.modes
            modes=selected if r==0 else list(reversed(selected)) if r==1 else selected[len(selected)//2:]+selected[:len(selected)//2]
            for mode in modes:
                if mode in rejected:continue
                name=f'r{r}_{mode}';started=time.perf_counter()
                try:run=subprocess.run([str(tool),mode,str(out/(name+'.json')),*map(str,files)],env=env,cwd=ROOT,capture_output=True,text=True,timeout=300)
                except subprocess.TimeoutExpired as e:
                    record['cases'].append(dict(round=r,mode=mode,wall_seconds=time.perf_counter()-started,returncode=124,status='rejected_timeout'))
                    rejected.add(mode);record['rejected_modes']=sorted(rejected)
                    (out/(name+'.log')).write_text('Timed out at 300 seconds; not a qualified timing sample.\n'+str(e.stdout or '')+str(e.stderr or ''))
                    write(out/'summary.json',record);print(name+' rejected: timeout',flush=True);continue
                (out/(name+'.log')).write_text(run.stdout+run.stderr)
                record['cases'].append(dict(round=r,mode=mode,wall_seconds=time.perf_counter()-started,returncode=run.returncode))
                if run.returncode:
                    rejected.add(mode);record['rejected_modes']=sorted(rejected)
                    record['cases'][-1]['status']='rejected_quality_or_solver_failure'
                write(out/'summary.json',record)
                print(name+(' passed' if run.returncode==0 else ' rejected: solver/quality failure'),flush=True)
        for group in ['input_sha256','source_sha256','runtime_sha256']:
            for path,sha in record[group].items():
                if digest(Path(path))!=sha:raise ValueError('Changed evidence: '+path)
        if digest(tool)!=record['tool_sha256']:raise ValueError('Changed tool')
        record['status']='completed_with_rejected_variants' if rejected else 'pass'
    except BaseException as e:
        record.update(status='failed',error=repr(e))
        if hasattr(e,'samples'):record['idle_failure_samples']=e.samples
        raise
    finally:write(out/'summary.json',record)

if __name__=='__main__':main()
