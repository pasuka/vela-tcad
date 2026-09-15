"""Serial frozen single-point NGMRES recovery controls; never a full-curve qualification."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import shutil
import time

from analyze_templates_ldmos_predictor_study import state_delta
from run_templates_ldmos_electrothermal_curve import state_gate


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-root',type=Path,required=True)
    parser.add_argument('--profile',type=Path,required=True)
    parser.add_argument('--probe',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=args.evidence_root.resolve();out=args.output.resolve();probe=args.probe.resolve()
    profile=read(args.profile)
    for name,digest in profile['files_sha256'].items():
        if sha(root/name)!=digest:raise ValueError('Frozen input changed: '+name)
    recorded=PurePosixPath(profile['dependencies'][0]['recorded_path']).parents[2]
    def mapped(name):
        relative=PurePosixPath(name).relative_to(recorded)
        if '..' in relative.parts:raise ValueError('Invalid recorded path')
        return root.joinpath(*relative.parts)
    cases=[('r7_vg4_4V',root/'results/r7_full_vg4/ledger.json',4.),
           ('r8_vg8_2p472V',root/'r8_density_low_20260915/vela_vg8/ledger.json',2.4723958333333336)]
    out.mkdir(parents=True,exist_ok=False)
    report=dict(schema='vela.ldmos.stagnation_study.v1',status='running',probe_sha256=sha(probe),profile_sha256=sha(args.profile),
        scope='Same frozen stalled input and original budget; Newton and NGMRES updates counted separately; no full-curve or native timing qualification',runs=[])
    def save():
        temp=out/'summary.tmp';temp.write_text(json.dumps(report,indent=2),encoding='utf-8');os.replace(temp,out/'summary.json')
    save();env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    for name,ledger_path,bias in cases:
        ledger=read(ledger_path)
        matches=[r for r in ledger['runs'] if not r['gate']['pass_gate'] and abs(r['bias_V']-bias)<1e-12]
        if len(matches)!=1:raise ValueError('Expected one frozen failure: '+name)
        source=mapped(matches[0]['directory']);original=read(source/'output.json')
        if state_gate(original,bias)['pass_gate']:raise ValueError('Sample is not a failed point')
        if original['diagnostic_stop']!='diagnostic_stagnation_reject':raise ValueError('Expected a stagnation-watch sample')
        snapshot=out/('source_'+name);snapshot.mkdir()
        for filename in ('input.json','output.json'):shutil.copy2(source/filename,snapshot/filename)
        shutil.copy2(ledger_path,snapshot/'ledger.json')
        for variant in ('baseline','ngmres'):
            case=out/(name+'_'+variant);case.mkdir()
            cfg=read(source/'input.json');cfg['mesh_file']=str(root/'cases/data/mesh.json')
            if 'ialmob' in cfg['mobility_SI']:cfg['mobility_SI']['ialmob']['geometry_file']=str(root/'cases/data/ialmob_geometry.json')
            cfg['performance_profiling']=True
            if variant=='ngmres':cfg['diagnostic_ngmres_recovery']=True
            inp=case/'input.json';inp.write_text(json.dumps(cfg),encoding='utf-8')
            if sha(probe)!=report['probe_sha256']:raise ValueError('Probe changed during study')
            before=None
            if os.name!='nt':
                import resource
                before=resource.getrusage(resource.RUSAGE_CHILDREN)
            with (case/'run.log').open('x') as log:
                start=time.perf_counter();process=subprocess.Popen([str(probe),str(inp),str(case/'output.json')],stdout=log,stderr=subprocess.STDOUT,env=env)
                while True:
                    try:code=process.wait(timeout=30);break
                    except subprocess.TimeoutExpired:print(json.dumps(dict(sample=name,variant=variant,elapsed_s=time.perf_counter()-start)),flush=True)
                wall=time.perf_counter()-start
            row=dict(sample=name,bias_V=bias,variant=variant,source_ledger_sha256=sha(ledger_path),source_input_sha256=sha(source/'input.json'),
                source_output_sha256=sha(source/'output.json'),input_sha256=sha(inp),exit_code=code,wall_seconds=wall)
            if before is not None:
                after=resource.getrusage(resource.RUSAGE_CHILDREN);row['child_cpu_seconds']=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
            if (case/'output.json').exists():
                result=read(case/'output.json');row.update(gate=state_gate(result,bias),newton_updates=result['newton_updates'],
                    ngmres_updates=result.get('ngmres_recovery',{}).get('updates',0),stop=result['diagnostic_stop'],
                    recovery=result.get('ngmres_recovery'),performance=result['performance'],
                    line_search_trials=sum(h.get('line_search_trials',0) for h in result['history']))
                row['total_state_updates']=row['newton_updates']+row['ngmres_updates']
                if variant=='baseline':
                    delta=state_delta(result,original);row['source_state_delta_V_V_V_K']=delta
                    if any(delta) or result['newton_updates']!=original['newton_updates']:raise ValueError('Frozen failure did not reproduce exactly')
            report['runs'].append(row);save();print(json.dumps(row),flush=True)
    report['status']='completed_diagnostic_matrix';save()


if __name__=='__main__':
    main()
