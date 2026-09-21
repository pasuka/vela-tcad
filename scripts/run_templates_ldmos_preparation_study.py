"""Serial R9 first8 controls for immutable sweep preparation, unchanged gates."""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from analyze_templates_ldmos_contact_sweep import numerical
from run_templates_ldmos_contact_sweep import read, sha
from run_templates_ldmos_gprof_first8 import inspect


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('bundle','runner','reference','output'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--repeats',type=int,default=2)
    p.add_argument('--thermal-high-field',action='store_true',
                   help='Verify the combined candidate; enables within-element reuse in both controls')
    a=p.parse_args()
    if a.repeats<1:p.error('Positive repeat count required')
    bundle,out=a.bundle.resolve(),a.output.resolve()
    manifest=read(bundle/'manifest.json')
    for name,digest in manifest['files_sha256'].items():assert sha(bundle/name)==digest,name
    cfg,deck=read(bundle/'input_vg8.json'),read(bundle/'vg8.json')
    if a.thermal_high_field:cfg['reuse_ialmob_thermal_high_field']=True
    assert cfg['diagnostic_near_steady_contact_consistency'] and cfg['performance_profiling']
    for key in ('reuse_neutral_contact_roots','diagnostic_neutral_root_newton','diagnostic_local_qf_limiter'):
        assert not cfg.get(key,False)
    deck['sweep']['bias_points_V']=deck['sweep']['bias_points_V'][:8]
    reference,reference_drain,reference_init,_=inspect(a.reference.resolve(),deck['sweep']['bias_points_V'])
    expected=numerical(reference_drain),numerical(reference_init)
    out.mkdir(parents=True,exist_ok=False)
    runtime=out/'runtime';runtime.mkdir()
    runner=a.runner.resolve();exe=runtime/runner.name;shutil.copy2(runner,exe)
    for path in runner.parent.glob('*.dll'):shutil.copy2(path,runtime/path.name)
    frozen={path.name:sha(path) for path in runtime.iterdir()}
    report=dict(status='running',runtime_sha256=frozen,thermal_high_field=a.thermal_high_field,bundle_manifest_sha256=sha(bundle/'manifest.json'),
                reference_ledger_sha256=sha(a.reference/'results/ledger.json'),runs=[])
    def save():
        (out/'summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    save();env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    try:
        for repeat in range(a.repeats):
            for enabled in ((False,True) if repeat%2==0 else (True,False)):
                name='r%d_%s'%(repeat,'static' if enabled else 'baseline')
                directory=out/name;directory.mkdir()
                run_deck=copy.deepcopy(deck)
                run_deck.update(reuse_static_preparation=enabled,input_file=str(directory/'input.json'),
                                output_directory=str(directory/'results'))
                (directory/'input.json').write_text(json.dumps(cfg),encoding='utf-8')
                (directory/'deck.json').write_text(json.dumps(run_deck,indent=2),encoding='utf-8')
                start=time.perf_counter()
                with (directory/'run.log').open('x',encoding='utf-8') as log:
                    proc=subprocess.Popen([str(exe),'--config',str(directory/'deck.json')],cwd=directory,
                                          env=env,stdout=log,stderr=subprocess.STDOUT)
                    while True:
                        try:code=proc.wait(timeout=30);break
                        except subprocess.TimeoutExpired:
                            progress=dict(run=name,elapsed_seconds=time.perf_counter()-start)
                            path=directory/'results/ledger.json'
                            if path.exists():
                                try:
                                    ledger=read(path);progress.update(bias_V=ledger['accepted_bias_V'],
                                        exact_points=len(ledger['exact_points']))
                                except (ValueError,OSError):pass
                            print(json.dumps(progress),flush=True)
                wall=time.perf_counter()-start;assert code==0,(name,code)
                ledger,drain,init,row=inspect(directory,deck['sweep']['bias_points_V'])
                assert (numerical(drain),numerical(init))==expected,'Numerical trajectory changed'
                for current,old in zip(ledger['runs'],reference['runs']):
                    assert current['bias_V']==old['bias_V'] and current['parent_bias_V']==old['parent_bias_V']
                    assert numerical(current['prediction'])==numerical(old['prediction'])
                results=init+drain
                hits=sum(bool(x['performance']['static_preparation_reused']) for x in results)
                assert hits==(len(results)-1 if enabled else 0),(name,hits)
                row.update(name=name,repeat=repeat,reuse_static_preparation=enabled,external_wall_seconds=wall,
                           numerical_exact_to_frozen=True,static_preparation_hits=hits,
                           input_sha256=sha(directory/'input.json'),deck_sha256=sha(directory/'deck.json'))
                report['runs'].append(row);save();print(json.dumps(row),flush=True)
        for name,digest in frozen.items():assert sha(runtime/name)==digest
        report['status']='complete';save()
    except BaseException as error:
        report.update(status='failed',error=repr(error));save();raise


if __name__=='__main__':main()
