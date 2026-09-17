"""Same-VM full R9/immutable-preparation + thermal-HFS/native repeat controls.

No physical defaults, gates, continuation or output scope changes. Linux only.
"""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from analyze_templates_ldmos_neutral_cache import inspect
from run_templates_ldmos_contact_sweep import read, sha, score


def prepare(cfg, deck, directory, enabled):
    cfg, deck = copy.deepcopy(cfg), copy.deepcopy(deck)
    assert cfg['diagnostic_near_steady_contact_consistency']
    for key in ('reuse_neutral_contact_roots', 'diagnostic_neutral_root_newton',
                'diagnostic_local_qf_limiter', 'diagnostic_ialmob_kernel_timing',
                'diagnostic_iteration_trace'):
        assert not cfg.get(key, False), key
    assert deck['initialization']['mode']=='neutral_300K'
    points=deck['sweep']['bias_points_V']
    assert len(points)==31 and points[0]==0 and points[-1]==40
    cfg['reuse_ialmob_thermal_high_field']=enabled
    deck.update(reuse_static_preparation=enabled, input_file=str(directory/'input.json'),
                output_directory=str(directory/'results'))
    return cfg,deck


def main():
    import resource
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('environment','runner','profile','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--repeats',type=int,default=2)
    args=parser.parse_args()
    if args.repeats<2:parser.error('At least two rounds required')
    envroot=args.environment.resolve();out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    profile=read(args.profile)
    for name,digest in profile['files_sha256'].items():
        assert sha(envroot/name)==digest,name
    binary=out/'vela_example_runner';shutil.copy2(str(args.runner),str(binary))
    native=read(envroot/'benchmark_r7.json')
    report=dict(status='running',repeats=args.repeats,runner_sha256=sha(binary),
        profile_sha256=sha(args.profile),runs=[],
        scope='Serial fresh neutral 300K full 31-point curves per gate; both preparation controls together; identical native output scope; process wall/CPU include initialization, all attempts and output. Audit time excluded.')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    def save():
        tmp=out/'summary.tmp';tmp.write_text(json.dumps(report,indent=2))
        os.replace(str(tmp),str(out/'summary.json'))
    save()
    references={}
    for case in profile['cases']:
        references[case['gate_V']]=inspect(envroot/Path(case['deck']).parent)[0]
    try:
        for repeat in range(args.repeats):
            for gate in (4,8):
                case=next(c for c in profile['cases'] if c['gate_V']==gate)
                order=('native','baseline','combined') if (repeat+(gate==8))%2==0 else ('combined','baseline','native')
                for variant in order:
                    directory=out/('r%d_vg%d_%s'%(repeat,gate,variant));directory.mkdir()
                    row=dict(repeat=repeat,gate_V=gate,variant=variant,directory=str(directory),
                             start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
                    if variant=='native':
                        for name in ('IdVd.cmd','n1_fps.tdr','sdevice.par','Siliconc100.par'):
                            source=envroot/'cases_r7'/('native_vg%d'%gate)/name
                            assert sha(source)==native['cases_sha256']['native_vg%d/%s'%(gate,name)]
                            shutil.copy2(str(source),str(directory/name))
                        argv=['/atctools/Synopsys/tcad/T-2022.03/bin/sdevice','IdVd.cmd'];cwd=directory
                    else:
                        cfg,deck=prepare(read(envroot/case['input']),read(envroot/case['deck']),directory,variant=='combined')
                        (directory/'input.json').write_text(json.dumps(cfg))
                        (directory/'deck.json').write_text(json.dumps(deck,indent=2))
                        argv=[str(binary),'--config',str(directory/'deck.json')];cwd=out
                    before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.perf_counter()
                    with (directory/'run.log').open('x') as log:
                        proc=subprocess.Popen(argv,cwd=str(cwd),env=env,stdout=log,stderr=subprocess.STDOUT)
                        while True:
                            try:code=proc.wait(timeout=30);break
                            except subprocess.TimeoutExpired:
                                progress=dict(run=directory.name,elapsed_seconds=time.perf_counter()-start)
                                path=directory/'results/ledger.json'
                                if path.exists():
                                    data=read(path);progress.update(bias_V=data['accepted_bias_V'],exact_points=len(data['exact_points']))
                                print(json.dumps(progress),flush=True)
                    wall=time.perf_counter()-start;after=resource.getrusage(resource.RUSAGE_CHILDREN)
                    row.update(exit_code=code,external_wall_seconds=wall,
                               child_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
                    assert code==0,row
                    if variant!='native':
                        reference=envroot/Path(case['deck']).parent/'results'
                        row['qualification']=score(directory/'results',reference,31)
                        signatures,costs=inspect(directory)
                        assert signatures==references[gate],('Changed complete numerical trajectory',directory)
                        ledger=read(directory/'results/ledger.json')
                        paths=[Path(r['result']) for r in ledger['initialization_runs']]+[Path(r['directory'])/'output.json' for r in ledger['runs']]
                        perf={}
                        for path in paths:
                            for key,value in read(path)['performance'].items():
                                if isinstance(value,(int,float)):perf[key]=perf.get(key,0)+value
                        assert perf['static_preparation_reused']==(len(paths)-1 if variant=='combined' else 0)
                        row.update(exact_to_frozen_R9=True,costs=costs,performance=perf,
                            input_sha256=sha(directory/'input.json'),ledger_sha256=sha(directory/'results/ledger.json'))
                    report['runs'].append(row);save();print(json.dumps(row),flush=True)
        assert sha(binary)==report['runner_sha256']
        report['status']='complete';save()
    except BaseException as error:
        report.update(status='failed',error=repr(error));save();raise


if __name__=='__main__':main()
