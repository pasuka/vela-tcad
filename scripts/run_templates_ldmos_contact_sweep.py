"""Serial original-seed contact-repair sweeps and same-VM timing controls.

Linux execution, frozen R7 profile, explicit output-only native decks. First8
includes an actual candidate pause/resume; full rounds always start afresh.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from analyze_templates_ldmos_predictor_study import state_delta
from run_templates_ldmos_electrothermal_curve import state_gate

FLAG = 'diagnostic_near_steady_contact_consistency'
STATE_KEYS = ('state_interleaved', 'referenced_state_interleaved',
              'electron_qf_reference_V', 'hole_qf_reference_V', 'residual', 'history')


def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8'))


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def prepare(original_input, original_deck, directory, variant, stage):
    if variant not in ('baseline', 'contact') or stage not in ('first8', 'full'):
        raise ValueError('Unknown variant or stage')
    cfg, deck = copy.deepcopy(original_input), copy.deepcopy(original_deck)
    if FLAG in cfg:
        raise ValueError('Expected unmodified frozen R7 input')
    if variant == 'contact':
        cfg[FLAG] = True
    if stage == 'first8':
        deck['sweep']['bias_points_V'] = deck['sweep']['bias_points_V'][:8]
    deck['input_file'] = str(directory/'input.json')
    deck['output_directory'] = str(directory/'results')
    return cfg, deck


def score(directory, reference, count):
    ledger = read(directory/'ledger.json')
    old = read(reference/'ledger.json')
    if ledger['status'] != 'complete' or len(ledger['exact_points']) != count:
        raise ValueError('Incomplete sweep')
    delta = [0.]*4
    for i, point in enumerate(ledger['exact_points']):
        target = old['exact_points'][i]
        if point['bias_V'] != target['bias_V']:
            raise ValueError('Exact voltage changed')
        result = read(point['result']); prior = read(target['result'])
        change = state_delta(result, prior)
        delta = [max(a,b) for a,b in zip(delta,change)]
        if not state_gate(result,point['bias_V'])['pass_gate']:
            raise ValueError('Original terminal gates failed')
    if any(v>t for v,t in zip(delta,(1e-8,1e-8,1e-8,1e-7))):
        raise ValueError('State differs from qualified reference: '+str(delta))
    costs = dict(updates=0,attempts=len(ledger['runs']),trials=0,failed_attempts=0,
                 contact_attempts=0,contact_accepted=0,contact_reassemblies=0,contact_seconds=0.)
    for row in ledger['runs']:
        result = read(Path(row['directory'])/'output.json')
        costs['updates'] += result['newton_updates']
        costs['trials'] += sum(h['line_search_trials'] for h in result['history'])
        costs['failed_attempts'] += not row['gate']['pass_gate']
        if row['gate']['pass_gate'] and not state_gate(result,row['bias_V'])['pass_gate']:
            raise ValueError('Continuation gate changed')
        for event in result.get('near_steady_contact_consistency',{}).get('events',[]):
            costs['contact_attempts'] += 1; costs['contact_accepted'] += event['accepted']
            costs['contact_reassemblies'] += event['reassembled']; costs['contact_seconds'] += event['seconds']
    if len(ledger['initialization_runs']) != len(old['initialization_runs']):
        raise ValueError('Initialization stage count changed')
    for a,b in zip(ledger['initialization_runs'],old['initialization_runs']):
        x,y=read(a['result']),read(b['result'])
        if any(x[k]!=y[k] for k in STATE_KEYS):
            raise ValueError('Initialization changed')
        if read(Path(a['result']).parent/'input.json').get(FLAG,False):
            raise ValueError('Contact candidate leaked into initialization')
    costs['initialization_updates']=sum(r['newton_updates'] for r in ledger['initialization_runs'])
    return dict(pass_gate=True,max_state_delta=delta,initialization_exact=True,costs=costs)


def main():
    import resource
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('environment','runner','profile','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--stage',choices=('first8','full'),required=True)
    p.add_argument('--repeats',type=int,default=2)
    args=p.parse_args()
    if args.repeats<1: p.error('Positive repeats required')
    envroot=args.environment.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    profile=read(args.profile)
    for name,digest in profile['files_sha256'].items():
        if sha(envroot/name)!=digest:raise ValueError('Frozen dependency changed: '+name)
    binary=out/'vela_example_runner';shutil.copy2(str(args.runner),str(binary))
    original=read(envroot/'benchmark_r7.json')
    report=dict(schema='vela.contact_sweep.v1',status='running',stage=args.stage,
        runner_sha256=sha(binary),profile_sha256=sha(args.profile),runs=[],
        scope='Serial independent neutral 300K initialization; original output scope; wall and child CPU include all runner work. First8 candidate includes pause/resume and is not a performance benchmark.')
    def save():
        tmp=out/'summary.tmp';tmp.write_text(json.dumps(report,indent=2),encoding='utf-8');os.replace(str(tmp),str(out/'summary.json'))
    save()
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    def execute(argv,cwd,log,progress=None,expected=0):
        before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.perf_counter()
        with log.open('x') as stream:
            process=subprocess.Popen(argv,cwd=str(cwd),env=env,stdout=stream,stderr=subprocess.STDOUT)
            while True:
                try:code=process.wait(timeout=30);break
                except subprocess.TimeoutExpired:
                    status=dict(run=log.parent.name,elapsed_seconds=time.perf_counter()-start)
                    if progress and progress.exists():
                        try:
                            ledger=read(progress);status.update(bias_V=ledger['accepted_bias_V'],exact_points=len(ledger['exact_points']),initialization_stages=len(ledger.get('initialization_runs',[])))
                        except (ValueError,OSError):pass
                    print(json.dumps(status),flush=True)
        wall=time.perf_counter()-start;after=resource.getrusage(resource.RUSAGE_CHILDREN)
        result=dict(exit_code=code,external_wall_seconds=wall,child_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
        if code!=expected:raise RuntimeError('Unexpected process exit: '+str(log)+' '+str(code))
        return result
    count=8 if args.stage=='first8' else 31
    repeats=1 if args.stage=='first8' else args.repeats
    for repeat in range(repeats):
        for gate in (4,8):
            order=('baseline','contact') if args.stage=='first8' else (('native','baseline','contact') if repeat%2==0 else ('contact','baseline','native'))
            for variant in order:
                directory=out/('r%d_vg%d_%s'%(repeat,gate,variant));directory.mkdir()
                row=dict(repeat=repeat,gate_V=gate,variant=variant,directory=str(directory),start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
                if variant=='native':
                    for name in ('IdVd.cmd','n1_fps.tdr','sdevice.par','Siliconc100.par'):
                        src=envroot/'cases_r7'/('native_vg%d'%gate)/name
                        if sha(src)!=original['cases_sha256']['native_vg%d/%s'%(gate,name)]:raise ValueError('Native input changed')
                        shutil.copy2(str(src),str(directory/name))
                    row.update(execute(['/atctools/Synopsys/tcad/T-2022.03/bin/sdevice','IdVd.cmd'],directory,directory/'run.log'))
                else:
                    cfg,deck=prepare(read(envroot/'cases_r7/data'/('input_vg%d.json'%gate)),read(envroot/'cases_r7'/('vela_vg%d.json'%gate)),directory,variant,args.stage)
                    (directory/'input.json').write_text(json.dumps(cfg),encoding='utf-8')
                    if args.stage=='first8' and variant=='contact':deck['pause_after_attempts']=3
                    path=directory/'deck.json';path.write_text(json.dumps(deck,indent=2))
                    argv=[str(binary),'--config',str(path)];ledger=directory/'results/ledger.json'
                    row.update(execute(argv,out,directory/'run.log',ledger,1 if 'pause_after_attempts' in deck else 0))
                    if 'pause_after_attempts' in deck:
                        checkpoint=read(ledger)
                        if checkpoint['status']!='stopped_at_checkpoint':raise ValueError('Not a controlled pause')
                        shutil.copy2(str(ledger),str(directory/'paused_ledger.json'))
                        deck.pop('pause_after_attempts');deck['resume']=True
                        path=directory/'resume.json';path.write_text(json.dumps(deck,indent=2))
                        resumed=execute([str(binary),'--config',str(path)],out,directory/'resume.log',ledger)
                        row['pause_resume']=True
                        for key in ('external_wall_seconds','child_cpu_seconds'):row[key]+=resumed[key]
                        row['exit_code']=resumed['exit_code']
                    row['qualification']=score(directory/'results',envroot/'results'/('r7_full_vg%d'%gate),count)
                    row['input_sha256']=sha(directory/'input.json');row['ledger_sha256']=sha(ledger)
                report['runs'].append(row);save();print(json.dumps(row),flush=True)
    report['status']='complete';save()


if __name__=='__main__':main()
