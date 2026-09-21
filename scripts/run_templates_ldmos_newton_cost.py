"""Isolated fixed-target R9 controls; retain failed attempts and total costs."""
import argparse, copy, json, os, shutil, subprocess, time
from pathlib import Path
from run_templates_ldmos_contact_sweep import read, sha
from analyze_templates_ldmos_contact_sweep import numerical
from analyze_templates_ldmos_predictor_study import state_delta
from run_templates_ldmos_electrothermal_curve import state_gate

REMOTE='/home/tcad/sentaurus_runs/vela_oracle_2022/vela_build_20260914'
FULL='extrapolation_20260915/contact_sweep_20260916/full'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('evidence','probe','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--variants',nargs='+',choices=('baseline','roots','local_qf','neutral_newton'),default=['baseline','roots','local_qf'])
    p.add_argument('--repeats',type=int,default=1)
    a=p.parse_args();root=a.evidence.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    if a.repeats<1:raise ValueError('Positive repeats required')
    binary=out/a.probe.name;shutil.copy2(str(a.probe),str(binary))
    if os.name=='nt':
        import ctypes
        ctypes.windll.kernel32.SetErrorMode(0x0001|0x0002|0x8000)
        for dll in a.probe.parent.glob('*.dll'):shutil.copy2(str(dll),str(out/dll.name))
    resolve=lambda s:root/s[len(REMOTE)+1:] if s.startswith(REMOTE+'/') else Path(s)
    def remap(value):
        if isinstance(value,dict):return {k:remap(v) for k,v in value.items()}
        if isinstance(value,list):return [remap(v) for v in value]
        return str(resolve(value)) if isinstance(value,str) and value.startswith(REMOTE+'/') else value
    cases=[]
    for gate in (4,8):
        ledger=read(root/FULL/('r0_vg%d_contact'%gate)/'results/ledger.json')
        assert ledger['status']=='complete'
        for bias in (.375,40/30*4,16.,40/30*23):
            rows=[r for r in ledger['runs'] if abs(r['bias_V']-bias)<1e-10 and r['gate']['pass_gate']]
            assert len(rows)==1
            row=rows[0];directory=resolve(row['directory']);cfg=remap(read(directory/'input.json'))
            assert cfg['diagnostic_near_steady_contact_consistency']
            cases.append((gate,row['bias_V'],cfg,read(directory/'output.json'),str(directory)))
    report=dict(status='running',scope='Eight actual accepted R9 input states; fixed targets, original gates and budgets. No curve or complete-sweep timing qualification.',
        runtime_sha256={p.name:sha(p) for p in [binary,*out.glob('*.dll')]},runs=[])
    def save():
        tmp=out/'summary.tmp';tmp.write_text(json.dumps(report,indent=2));os.replace(str(tmp),str(out/'summary.json'))
    save();env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    for repeat in range(a.repeats):
        for index,(gate,bias,original,reference,source) in enumerate(cases):
            results={}
            for variant in (a.variants if repeat%2==0 else list(reversed(a.variants))):
                cfg=copy.deepcopy(original)
                if variant=='roots':cfg['reuse_neutral_contact_roots']=True
                if variant=='local_qf':cfg['diagnostic_local_qf_limiter']=True
                if variant=='neutral_newton':cfg['diagnostic_neutral_root_newton']=True
                directory=out/('r%d_case%d_%s'%(repeat,index,variant));directory.mkdir()
                inp=directory/'input.json';inp.write_text(json.dumps(cfg))
                start=time.perf_counter()
                with (directory/'run.log').open('w') as log:
                    proc=subprocess.Popen([str(binary),str(inp),str(directory/'output.json')],stdout=log,stderr=subprocess.STDOUT,env=env)
                    while True:
                        try:code=proc.wait(timeout=30);break
                        except subprocess.TimeoutExpired:print(directory.name+' running',flush=True)
                row=dict(repeat=repeat,case=index,variant=variant,gate_V=gate,bias_V=bias,source=source,
                    input_sha256=sha(inp),exit_code=code,wall_seconds=time.perf_counter()-start)
                if not (directory/'output.json').exists():raise RuntimeError('Missing output: '+str(directory))
                result=read(directory/'output.json');results[variant]=result
                delta=state_delta(result,reference);gate_check=state_gate(result,bias)
                row.update(updates=result['newton_updates'],trials=sum(h['line_search_trials'] for h in result['history']),
                    floor_updates=sum(h['scaled_l2_before']<1e-9 for h in result['history']),
                    max_state_delta=delta,gate=gate_check,performance=result['performance'],
                    qualified=code==0 and gate_check['pass_gate'] and all(d<=t for d,t in zip(delta,(1e-8,1e-8,1e-8,1e-7))),
                    output_sha256=sha(directory/'output.json'))
                report['runs'].append(row);save();print(json.dumps(row),flush=True)
                if variant=='baseline' and not row['qualified']:raise RuntimeError('Baseline failed qualification')
            if 'baseline' in results:
                for row in report['runs']:
                    if row['repeat']==repeat and row['case']==index:
                        row['numerical_exact_to_baseline']=numerical(results[row['variant']])==numerical(results['baseline'])
                save()
            if 'baseline' in results and 'roots' in results:
                assert numerical(results['baseline'])==numerical(results['roots']),'Root cache changed numerical trajectory'
    report['status']='complete';save()

if __name__=='__main__':main()
