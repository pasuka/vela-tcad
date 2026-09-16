"""Replay original R7 attempt inputs with opt-in trace and exact trajectory checks."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from analyze_templates_ldmos_newton_phases import analyze


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence-root',type=Path,required=True)
    p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--limit',type=int,default=0)
    args=p.parse_args();root=args.evidence_root.resolve();out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    report=dict(schema='vela.newton_trace_replay.v1',status='running',probe_sha256=sha(args.probe),runs=[])
    def save():
        temp=out/'summary.tmp';temp.write_text(json.dumps(report,indent=2),encoding='utf-8');os.replace(temp,out/'summary.json')
    save()
    for vg in (4,8):
        ledger_path=root/('results/r7_full_vg%d/ledger.json'%vg)
        ledger=read(ledger_path)
        for index,row in enumerate(ledger['runs']):
            if args.limit and index>=args.limit:break
            source=Path(row['directory'])
            source.relative_to(root)  # reject any source outside frozen VM evidence
            cfg=read(source/'input.json');original=read(source/'output.json')
            cfg['diagnostic_iteration_trace']=True
            case=out/('vg%d_%03d'%(vg,index));case.mkdir()
            inp=case/'input.json';inp.write_text(json.dumps(cfg),encoding='utf-8')
            start=time.perf_counter()
            if sha(args.probe)!=report['probe_sha256']:raise ValueError('Probe changed')
            with (case/'run.log').open('x') as log:
                proc=subprocess.Popen([str(args.probe),str(inp),str(case/'output.json')],stdout=log,stderr=subprocess.STDOUT,
                                      env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
                while True:
                    try:code=proc.wait(timeout=30);break
                    except subprocess.TimeoutExpired:print('trace vg%d index%d running'%(vg,index),flush=True)
            result=read(case/'output.json');clean=[{k:v for k,v in h.items() if k!='iteration_trace'} for h in result['history']]
            exact=clean==original['history'] and all(result[k]==original[k] for k in (
                'referenced_state_interleaved','electron_qf_reference_V','hole_qf_reference_V',
                'residual','diagnostic_stop','newton_updates'))
            record=dict(vg=vg,index=index,bias_V=row['bias_V'],exit_code=code,exact_trajectory=exact,
                source_input_sha256=sha(source/'input.json'),source_output_sha256=sha(source/'output.json'),
                ledger_sha256=sha(ledger_path),wall_seconds=time.perf_counter()-start,
                trace_seconds=result.get('iteration_trace_seconds'),phases=analyze(result)['phases'])
            report['runs'].append(record);save();print(json.dumps(record),flush=True)
            if code or not exact:
                report['status']='failed_replay_check';save();raise RuntimeError('Trace changed trajectory or replay failed')
    report['status']='complete';save()


if __name__=='__main__':main()
