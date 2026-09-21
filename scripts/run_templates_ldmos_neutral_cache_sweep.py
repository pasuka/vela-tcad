"""Serial full R9/root-cache/native controls from independent neutral states."""
import argparse,json,os,resource,shutil,subprocess,time
from pathlib import Path
from run_templates_ldmos_contact_sweep import read,sha,score

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('environment','runner','profile','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();envroot=a.environment.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    profile=read(a.profile)
    for name,digest in profile['files_sha256'].items():assert sha(envroot/name)==digest,name
    binary=out/'vela_example_runner';shutil.copy2(str(a.runner),str(binary))
    report=dict(status='running',runner_sha256=sha(binary),profile_sha256=sha(a.profile),runs=[],
        scope='One serial independent neutral 300K full pair per gate, R9 versus exact-key neutral-root cache. Same native output scope and unchanged gates. All attempts, initialization and output included in process timing.')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1');native=read(envroot/'benchmark_r7.json')
    def save():
        tmp=out/'summary.tmp';tmp.write_text(json.dumps(report,indent=2));os.replace(str(tmp),str(out/'summary.json'))
    save()
    try:
        for gate in (4,8):
            case=next(c for c in profile['cases'] if c['gate_V']==gate)
            for variant in (('native','baseline','roots') if gate==4 else ('roots','baseline','native')):
                directory=out/('r0_vg%d_%s'%(gate,variant));directory.mkdir()
                if variant=='native':
                    for name in ('IdVd.cmd','n1_fps.tdr','sdevice.par','Siliconc100.par'):
                        source=envroot/'cases_r7'/('native_vg%d'%gate)/name
                        assert sha(source)==native['cases_sha256']['native_vg%d/%s'%(gate,name)]
                        shutil.copy2(str(source),str(directory/name))
                    argv=['/atctools/Synopsys/tcad/T-2022.03/bin/sdevice','IdVd.cmd'];cwd=directory
                else:
                    cfg=read(envroot/case['input']);deck=read(envroot/case['deck'])
                    assert cfg['diagnostic_near_steady_contact_consistency']
                    assert 'reuse_neutral_contact_roots' not in cfg and 'diagnostic_local_qf_limiter' not in cfg
                    if variant=='roots':cfg['reuse_neutral_contact_roots']=True
                    deck['input_file']=str(directory/'input.json');deck['output_directory']=str(directory/'results')
                    (directory/'input.json').write_text(json.dumps(cfg));(directory/'deck.json').write_text(json.dumps(deck,indent=2))
                    argv=[str(binary),'--config',str(directory/'deck.json')];cwd=out
                row=dict(gate_V=gate,variant=variant,directory=str(directory),start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
                before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.perf_counter()
                with (directory/'run.log').open('x') as log:
                    proc=subprocess.Popen(argv,cwd=str(cwd),env=env,stdout=log,stderr=subprocess.STDOUT)
                    while True:
                        try:code=proc.wait(timeout=30);break
                        except subprocess.TimeoutExpired:
                            progress=dict(run=directory.name,elapsed_seconds=time.perf_counter()-start)
                            ledger=directory/'results/ledger.json'
                            if ledger.exists():
                                data=read(ledger);progress.update(bias_V=data['accepted_bias_V'],exact_points=len(data['exact_points']))
                            print(json.dumps(progress),flush=True)
                wall=time.perf_counter()-start;after=resource.getrusage(resource.RUSAGE_CHILDREN)
                row.update(exit_code=code,external_wall_seconds=wall,child_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
                assert code==0,row
                if variant!='native':
                    reference=envroot/Path(case['deck']).parent/'results'
                    row['qualification']=score(directory/'results',reference,31)
                    row.update(input_sha256=sha(directory/'input.json'),ledger_sha256=sha(directory/'results/ledger.json'))
                report['runs'].append(row);save();print(json.dumps(row),flush=True)
        report['status']='complete';save()
    except BaseException as error:
        report['status']='failed';report['error']=repr(error);save();raise

if __name__=='__main__':main()
