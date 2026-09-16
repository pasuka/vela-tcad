"""F1 G1a frozen representative controls; all failures and trial costs retained."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from analyze_templates_ldmos_newton_phases import analyze
from analyze_templates_ldmos_predictor_study import state_delta
from run_templates_ldmos_electrothermal_curve import state_gate


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--variants',nargs='+',choices=['baseline','r8_60','v1','v2','nleq_r7','nleq_v1','nleq_v2','adaptive_r7','ptc_v1_s01','ptc_v1_s1','ptc_v1_s10','ptc_qf_s1','ptc_defect_ser_s1','ptc_defect_model_s1'],default=['baseline','r8_60','v1','v2'])
    p.add_argument('--cases',nargs='+',help='Frozen case names, e.g. vg4_p11; empty means all eight')
    p.add_argument('--no-trace',action='store_true')
    p.add_argument('--direction-audit',action='store_true',help='One iteration, read-only coefficient-direction comparisons; no point qualification')
    args=p.parse_args()
    if args.direction_audit and any(not v.startswith('ptc_') for v in args.variants):p.error('Direction audit requires only ptc variants')
    cases=[r for r in read(args.manifest) if r['variant']=='linear_baseline']
    if len(cases)!=8:raise ValueError('Expected eight frozen representative inputs')
    if args.cases:
        names={'vg%d_p%d'%(c['gate_V'],c['baseline_index']) for c in cases}
        if not set(args.cases)<=names:raise ValueError('Unknown frozen case; available: '+', '.join(sorted(names)))
        cases=[c for c in cases if 'vg%d_p%d'%(c['gate_V'],c['baseline_index']) in args.cases]
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    binary=out/args.probe.name;shutil.copy2(args.probe,binary)
    if os.name=='nt':
        import ctypes
        ctypes.windll.kernel32.SetErrorMode(0x0001|0x0002|0x8000)
        for dll in args.probe.parent.glob('*.dll'):shutil.copy2(dll,out/dll.name)
    report=dict(schema='vela.density_projection_g1a.v1',status='running',binary_sha256=sha(binary),
                manifest_sha256=sha(args.manifest),runtime_sha256={p.name:sha(p) for p in out.glob('*.dll')},direction_audit_only=args.direction_audit,runs=[])
    def save():
        temp=out/'summary.tmp';temp.write_text(json.dumps(report,indent=2),encoding='utf-8');os.replace(temp,out/'summary.json')
    save()
    for variant in args.variants:
        for case in cases:
            directory=out/('%s_vg%d_p%d'%(variant,case['gate_V'],case['baseline_index']));directory.mkdir()
            cfg=read(Path(case['input']));cfg['performance_profiling']=True
            cfg['diagnostic_iteration_trace']=not args.no_trace
            if variant=='r8_60':cfg['diagnostic_density_update_iterations']=60
            if variant in ('v1','v2'):cfg['diagnostic_density_projection']=variant
            if variant.startswith('nleq_'):
                cfg['diagnostic_natural_damping']=True
                if variant!='nleq_r7':cfg['diagnostic_density_projection']=variant[-2:]
            if variant=='adaptive_r7':cfg['diagnostic_adaptive_jacobian']=True
            if variant.startswith('ptc_'):
                cfg['diagnostic_pseudo_transient']=True
                cfg['diagnostic_pseudo_time_scale']={'s01':.1,'s1':1.,'s10':10.}[variant.rsplit('_',1)[1]]
                if '_v1_' in variant:cfg['diagnostic_density_projection']='v1'
                if variant.startswith('ptc_defect_'):
                    cfg['diagnostic_density_projection']='v1'
                    cfg['diagnostic_pseudo_acceptance']='defect_'+variant.split('_')[2]
            if args.direction_audit:cfg.update(diagnostic_pseudo_direction_audit=True,diagnostic_newton_max_iterations=1)
            inp=directory/'input.json';inp.write_text(json.dumps(cfg),encoding='utf-8')
            start=time.perf_counter()
            with (directory/'run.log').open('x') as log:
                proc=subprocess.Popen([str(binary),str(inp),str(directory/'output.json')],stdout=log,stderr=subprocess.STDOUT,
                    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
                while True:
                    try:code=proc.wait(timeout=30);break
                    except subprocess.TimeoutExpired:print(directory.name+' running',flush=True)
            row=dict(variant=variant,vg=case['gate_V'],bias_V=case['bias_V'],exit_code=code,
                wall_seconds=time.perf_counter()-start,input_sha256=sha(inp),
                source_input_sha256=sha(Path(case['input'])),source_output_sha256=sha(Path(case['expected'])),qualified=False)
            if (directory/'output.json').exists():
                result=read(directory/'output.json');gate=state_gate(result,case['bias_V'])
                delta=state_delta(result,read(Path(case['expected'])))
                row.update(output_sha256=sha(directory/'output.json'),natural_damping=result.get('natural_damping'),adaptive_jacobian=result.get('adaptive_jacobian'),pseudo_transient=result.get('pseudo_transient'),pseudo_direction_audit=result.get('pseudo_direction_audit'),newton_updates=result['newton_updates'],stop=result['diagnostic_stop'],gate=gate,
                    max_state_delta=delta,qualified=not args.direction_audit and code==0 and gate['pass_gate'] and all(v<=t for v,t in zip(delta,(1e-8,1e-8,1e-8,1e-7))),
                    phases=analyze(result)['phases'],performance=result['performance'],
                    trace_seconds=result.get('iteration_trace_seconds'),
                    trials=sum(h['line_search_trials'] for h in result['history']),
                    fallbacks=sum(h.get('density_update_fallback',False) for h in result['history']))
            report['runs'].append(row);save();print(json.dumps(row),flush=True)
            if variant=='baseline' and not row['qualified']:raise RuntimeError('Frozen baseline did not qualify')
    report['status']='completed_diagnostic_matrix';save()


if __name__=='__main__':main()
