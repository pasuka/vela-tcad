"""Frozen R10 seeds, explicit HFS and generated IALMob controls; unchanged gates."""
import argparse,copy,json,os,shutil,subprocess,time
from pathlib import Path
from run_templates_ldmos_contact_sweep import read,sha
from run_templates_ldmos_electrothermal_curve import state_gate
from analyze_templates_ldmos_predictor_study import state_delta
from analyze_templates_ldmos_contact_sweep import numerical

REMOTE='/home/tcad/sentaurus_runs/vela_oracle_2022/vela_build_20260914'
FULL='extrapolation_20260915/preparation_20260917/full'
FLAGS={'baseline':(False,False),'high':(True,False),'low':(False,True),'combined':(True,True)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('evidence','probe','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--variants',nargs='+',choices=tuple(FLAGS),default=['baseline','high','combined'])
    p.add_argument('--repeats',type=int,default=2)
    a=p.parse_args();root=a.evidence.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    assert a.repeats>0
    runtime=out/'runtime';runtime.mkdir();binary=runtime/a.probe.name;shutil.copy2(a.probe,binary)
    if os.name=='nt':
        import ctypes
        ctypes.windll.kernel32.SetErrorMode(0x0001|0x0002|0x8000)
        for dll in a.probe.parent.glob('*.dll'):shutil.copy2(dll,runtime/dll.name)
    frozen={p.name:sha(p) for p in runtime.iterdir()}
    def remap(x):
        if isinstance(x,dict):return {k:remap(v) for k,v in x.items()}
        if isinstance(x,list):return [remap(v) for v in x]
        return str(root/x[len(REMOTE)+1:]) if isinstance(x,str) and x.startswith(REMOTE+'/') else x
    cases=[]
    for g in (4,8):
        ledger=read(root/FULL/('r0_vg%d_combined/results/ledger.json'%g));assert ledger['status']=='complete'
        for target in (.375,40/30*4,16.,40/30*23):
            row=next(r for r in ledger['runs'] if abs(r['bias_V']-target)<1e-10 and r['gate']['pass_gate'])
            directory=Path(remap(row['directory']));cfg=remap(read(directory/'input.json'))
            assert cfg['reuse_ialmob_thermal_high_field'] and cfg['diagnostic_near_steady_contact_consistency']
            cases.append((g,row['bias_V'],directory,cfg,read(directory/'output.json')))
    report=dict(status='running',runtime_sha256=frozen,runs=[],scope=__doc__)
    def save():
        tmp=out/'summary.tmp';tmp.write_text(json.dumps(report,indent=2));os.replace(tmp,out/'summary.json')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1');save()
    try:
        for repeat in range(a.repeats):
            for i,(g,bias,source,original,reference) in enumerate(cases):
                results={}
                for variant in (a.variants if repeat%2==0 else list(reversed(a.variants))):
                    cfg=copy.deepcopy(original)
                    cfg['diagnostic_ialmob_explicit_high_field'],cfg['diagnostic_ialmob_generated_low_field']=FLAGS[variant]
                    directory=out/('r%d_case%d_%s'%(repeat,i,variant));directory.mkdir()
                    inp=directory/'input.json';inp.write_text(json.dumps(cfg));started=time.perf_counter()
                    with (directory/'run.log').open('x') as log:
                        code=subprocess.call([str(binary),str(inp),str(directory/'output.json')],env=env,stdout=log,stderr=subprocess.STDOUT)
                    wall=time.perf_counter()-started
                    if not (directory/'output.json').exists():
                        report['runs'].append(dict(repeat=repeat,case=i,gate_V=g,bias_V=bias,variant=variant,
                            exit_code=code,wall_seconds=wall,qualified=False,error='Missing point output',
                            input_sha256=sha(inp),source_input_sha256=sha(source/'input.json')))
                        save()
                        if variant=='baseline':raise RuntimeError('Baseline produced no output')
                        continue
                    result=read(directory/'output.json');results[variant]=result
                    delta=state_delta(result,reference);gate=state_gate(result,bias)
                    row=dict(repeat=repeat,case=i,gate_V=g,bias_V=bias,variant=variant,exit_code=code,wall_seconds=wall,
                        qualified=code==0 and gate['pass_gate'] and all(d<=t for d,t in zip(delta,(1e-8,1e-8,1e-8,1e-7))),
                        gate=gate,max_state_delta=delta,updates=result['newton_updates'],
                        trials=sum(h['line_search_trials'] for h in result['history']),performance=result['performance'],
                        source_input_sha256=sha(source/'input.json'),input_sha256=sha(inp),output_sha256=sha(directory/'output.json'))
                    report['runs'].append(row);save()
                    print(json.dumps({k:v for k,v in row.items() if k not in ('performance','gate')}),flush=True)
                    if variant=='baseline' and not row['qualified']:raise RuntimeError('Baseline failed')
                if 'baseline' in results:
                    for row in report['runs']:
                        if row['repeat']==repeat and row['case']==i and row['variant'] in results:
                            row['exact_to_baseline']=numerical(results[row['variant']])==numerical(results['baseline'])
                    save()
        assert all(sha(runtime/n)==h for n,h in frozen.items())
        report['status']='complete';report['all_qualified']=all(r['qualified'] for r in report['runs']);save()
    except BaseException as e:report.update(status='failed',error=repr(e));save();raise

if __name__=='__main__':main()
