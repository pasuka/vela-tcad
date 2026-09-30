"""Run the affected built Catch2 executables and preserve every exit status."""
import concurrent.futures
import json
import os
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'build-release/split_prod_0912/regressions'
TARGETS=('test_split_dd_runtime','test_split_dd_state','test_production_numerics',
         'test_element_box_transport','test_sg_flux','test_newton_solver','test_dc_sweep')

def one(target):
    env=os.environ.copy()
    env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env['PATH']
    for key in list(env):
        if key.startswith(('VELA_DIAGNOSTIC_','VELA_VALIDATE_','VELA_MINORITY_',
                           'VELA_SIMPLEMOS_','VELA_TEST_','VELA_NATIVE_GEOMETRY_','VELA_CANDIDATE_')):
            env.pop(key)
    start=time.monotonic()
    with (OUT/(target+'.log')).open('w',encoding='utf-8') as log:
        result=subprocess.run([str(ROOT/'build-release'/(target+'.exe'))],cwd=ROOT,
                              env=env,stdout=log,stderr=subprocess.STDOUT)
    row=dict(target=target,exit_code=result.returncode,seconds=time.monotonic()-start)
    print(row,flush=True)
    return row

if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(one,TARGETS))
    (OUT/'summary.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    raise SystemExit(any(r['exit_code'] for r in results))
