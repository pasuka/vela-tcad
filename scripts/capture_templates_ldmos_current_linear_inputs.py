"""Replay qualified R11 target attempts and capture current solver-input matrices.

These are diagnostic reruns, excluded from timing comparisons. VELALU02 contains
the actual scaled input to LinearSolver, not an unscaled physical Jacobian.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time
from run_templates_ldmos_linked_d5 import read,write,digest,state_difference

ROOT=Path(__file__).resolve().parents[1]

def relocate(value, old, new):
    if isinstance(value,dict):return {k:relocate(v,old,new) for k,v in value.items()}
    if isinstance(value,list):return [relocate(v,old,new) for v in value]
    if isinstance(value,str) and value.replace('\\','/').startswith(old.as_posix()+'/'):
        return str(new/Path(value.replace('\\','/')).relative_to(old))
    return value

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    runner=ROOT/'build-release/vela_example_runner.exe'
    bundle_path=ROOT/'reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_no_generation_inputs.json'
    bundle=read(bundle_path)
    for path,sha in bundle['files'].items():
        if digest(ROOT/path)!=sha:raise RuntimeError('Physical input changed: '+path)
    report=dict(status='running',coordinate_scope='solver_input_only',runner_sha256=digest(runner),cases=[],
                bundle_sha256=digest(bundle_path),physical_files=bundle['files'])
    try:
        for gate in [4,8]:
            source=ROOT/f'reference_staging/templates_ldmos_r11_isothermal_20260917/d5_vg{gate}'
            ledger=read(source/'fixed/ledger.json')
            # Use the actual R11 trajectory; 0.375 V belongs to older sequences.
            for bias in [.34621190996608003,4.,16./3.,20.,40.]:
                run=min(ledger['runs'],key=lambda r:abs(r['target_V']-bias))
                if abs(run['target_V']-bias)>1e-7:raise RuntimeError('Missing frozen target')
                old=source/run['case'];dest=out/f'vg{gate}_vd{bias:.6f}';dest.mkdir()
                config=relocate(read(old/'control.json'),old,dest)
                if config['sweep']['initial_state_file']==config['sweep']['write_state_file']:raise RuntimeError('Seed would be overwritten')
                write(dest/'control.json',config);captures=dest/'captures';captures.mkdir()
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE',
                         VELA_LINEAR_SOLVER='sparselu',VELA_LINEAR_FACTOR_STATISTICS='0',VELA_LINEAR_CAPTURE_DIR=str(captures))
                env['PATH']='D:/msys64/ucrt64/bin;'+env.get('PATH','')
                start=time.perf_counter()
                with (dest/'stdout.log').open('x') as log:
                    proc=subprocess.run([str(runner),'--config',str(dest/'control.json')],cwd=dest,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600)
                row=dict(gate=gate,bias_V=bias,source_control=str(old/'control.json'),source_sha256=digest(old/'control.json'),
                         seed_sha256=digest(Path(config['sweep']['initial_state_file'])),returncode=proc.returncode,wall_seconds=time.perf_counter()-start,
                         captures={str(x):digest(x) for x in captures.glob('*.bin')})
                report['cases'].append(row);write(out/'summary.json',report)
                if proc.returncode or not row['captures']:raise RuntimeError('Capture run failed: '+str(dest))
                row['state_difference']=state_difference(Path(config['sweep']['write_state_file']),old/'state.csv')
                if any(row['state_difference'][k]['max_absolute']>1e-8 for k in ['psi','phin','phip']):raise RuntimeError('Frozen state gate failed')
                write(out/'summary.json',report);print(f'CAPTURE_PASS Vg={gate} Vd={bias} matrices={len(row["captures"])}',flush=True)
        if digest(runner)!=report['runner_sha256']:raise RuntimeError('Runner changed during capture')
        report['status']='pass'
    except BaseException as e:report.update(status='failed',error=repr(e));raise
    finally:write(out/'summary.json',report)

if __name__=='__main__':main()
