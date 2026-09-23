"""Isolated state/2 vs frozen VDS1 qualification; never a production format switch.

The plan names immutable binaries, input decks and historical reference curves.
Run each curve in a fresh Python process to isolate the old driver's module state.
Only the historical control loads its VDS1 adapter. Final acceptance is cold-read.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time
import traceback
from types import SimpleNamespace


def read(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p, value):
    p=Path(p); temporary=p.with_suffix(p.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temporary.replace(p)


def one_curve(plan, mode, gate, directory, points):
    # Same dependency initialization, outside the measured controller boundary.
    import h5py, numpy
    sys.path.insert(0,str(Path(plan['sources'][mode])/'scripts'))
    import run_templates_ldmos_linked_d5 as b
    from ldmos_file_efficiency import CompletedFileCache
    from analyze_templates_ldmos_d5_gprof import aggregate
    from run_templates_ldmos_d5_gprof import trajectory_signature
    b.ROOT=Path(plan['workspace']);b.HERE=directory;directory.mkdir(exist_ok=False)
    case=plan['gates'][str(gate)]
    b.BASE=read(case['templates'][mode]);b.SEED=Path(case['seeds'][mode])
    b.RUNNER=Path(plan['runners'][mode]);b.EXPECTED=plan['frozen'][str(b.RUNNER)]
    b.REFERENCE=Path(case['reference']);b.FRAME=case['frame_offset_V']
    b.MAXIMUM=.2;b.PHYSICS_PROFILE='D5';b.WORKER=None
    b.args=SimpleNamespace(gate=gate,worker=True,reuse_linear_analysis=True)
    b.ENV=dict(os.environ,VELA_LINEAR_SOLVER='umfpack',VELA_LINEAR_THREADS='1',
        VELA_BLAS_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',
        OMP_DYNAMIC='FALSE',OMP_MAX_ACTIVE_LEVELS='1',VELA_LINEAR_FACTOR_STATS='0')
    targets=b.read_points(b.REFERENCE)[:points];b.stop=targets[-1]
    actual=lambda p:Path(p)
    if mode=='vds':
        from ldmos_binary_adapter import install_binary
        actual=install_binary(b)
    raw_rows,raw_digest=b.rows,b.digest
    cache=CompletedFileCache(raw_rows,raw_digest)
    b.rows=lambda p:cache.rows(actual(p));b.digest=lambda p:cache.digest(actual(p))
    original_run=b.DCWorker.run
    def worker_run(worker, config, *args, **kwargs):
        result=original_run(worker,config,*args,**kwargs)
        cache.seal(Path(config).parent);return result
    b.DCWorker.run=worker_run
    frozen=dict(plan['frozen'])
    curve_plan=dict(frozen_files=frozen,runner=str(b.RUNNER),runner_sha256=b.EXPECTED,
        gate_V=gate,reference=str(b.REFERENCE),seed=str(b.SEED),mode=mode,
        targets_V=targets,max_step_V=.2,frame_offset_V=b.FRAME,initial_step_V=b.INITIAL,
        backend='umfpack',original_blocks=b.GATES,local_rows_eps=1e-8)
    write(directory/'plan.json',curve_plan)
    progress=dict(status='running',pid=os.getpid(),mode=mode,gate=gate,exact_points=0)
    write(directory/'progress.json',progress)
    sweep=b.Sweep();start=time.perf_counter()
    try:
        b.initialize(sweep)
        for target in targets[1:]:
            if (directory.parent/'STOP').exists(): raise RuntimeError('Requested boundary stop')
            sweep.advance(target)
            if target==b.FRAME_PIVOT:sweep.change_frame()
            progress.update(exact_points=len(sweep.ledger['exact_points']),bias_V=target)
            write(directory/'progress.json',progress)
        cache.disable()
        verdicts=b.score(sweep);audit=b.audit(sweep,curve_plan)
        if points==31 and not all(v['pass_all'] for v in verdicts.values()):
            raise ValueError('Original curve accuracy gates failed')
        reference=read(Path(case['historical_curve'])/'fixed/ledger.json')
        expected=trajectory_signature(reference)[:len(sweep.ledger['runs'])]
        if trajectory_signature(sweep.ledger)!=expected:
            raise ValueError('Frozen nonlinear trajectory changed')
        sweep.ledger['status']='completed'
        progress.update(status='completed',verdicts=verdicts,audit=audit,trajectory_exact=True,
            exact_points=len(sweep.ledger['exact_points']),bias_V=targets[-1])
    except BaseException as error:
        sweep.ledger.update(status='failed',error=repr(error));progress.update(status='failed',error=repr(error));raise
    finally:
        if b.WORKER is not None:b.WORKER.close()
        sweep.save();progress['wall_seconds']=time.perf_counter()-start
        write(directory/'progress.json',progress)
    progress['performance']=aggregate(directory)
    progress['file_cache_statistics']=dict(cache.stats)
    write(directory/'progress.json',progress)


def cold_fields(path, mesh, mode):
    import numpy as np
    import state_archive
    if mode=='hdf5':
        fields,_=state_archive.read(path,len(mesh['nodes']),state_archive.mesh_identity(mesh,1e-6))
    else:
        import dd_state_binary
        fields=state_archive.rows_to_fields(dd_state_binary.read(path))
    # No serializer normalization: signed zero and all split low bits participate.
    return {name:hashlib.sha256(np.asarray(v,dtype='<f8').tobytes()).hexdigest() for name,v in fields.items()}


def audit_batch(plan, output, cases, rounds, points):
    sys.path.insert(0,str(Path(plan['sources']['hdf5'])/'scripts'))
    from analyze_templates_ldmos_stage4_d5 import ratio_error,read_vela_curve,LIMITS
    from summarize_templates_ldmos_idvd_ablation import read_curve
    from run_templates_ldmos_d5_gprof import trajectory_signature
    signatures={};joint=[];sizes=[]
    for case in cases:
        r,m,g=case['repeat'],case['mode'],case['gate'];d=Path(case['directory'])
        ledger=read(d/'fixed/ledger.json');report=read(d/'progress.json')
        if report['status']!='completed' or ledger['status']!='completed':raise ValueError('Incomplete curve')
        if len(ledger['exact_points'])!=points:raise ValueError('Incomplete exact points')
        mesh=read(plan['gates'][str(g)]['mesh'])
        values=[(p['bias_V'],cold_fields(Path(p['state']),mesh,m)) for p in ledger['exact_points']]
        signatures[r,m,g]=(values,trajectory_signature(ledger),report['performance']['counters'])
        files=[p for p in d.rglob('*') if p.is_file() and p.suffix in ('.vds','.h5')]
        sizes.append(dict(repeat=r,mode=m,gate=g,files=len(files),bytes=sum(p.stat().st_size for p in files)))
    for r in range(rounds):
        for g in (4,8):
            if signatures[r,'vds',g]!=signatures[r,'hdf5',g]:raise ValueError('Cross-format fields, trajectory or counters differ')
        for m in ('vds','hdf5'):
            if points==31:
                refs={f'Vg{g}':read_curve(Path(plan['gates'][str(g)]['reference'])) for g in (4,8)}
                curves={f'Vg{g}':read_vela_curve(output/f'r{r}_{m}_vg{g}/score/curve.csv') for g in (4,8)}
                ratio=ratio_error(refs,curves)
                checks={k:ratio['endpoint']<=v['gate_ratio'] for k,v in LIMITS.items()}
                if not all(checks.values()):raise ValueError('Joint gate failed')
                joint.append(dict(repeat=r,mode=m,ratio=ratio,checks=checks))
            for g in (4,8):
                if signatures[r,m,g]!=signatures[0,m,g]:raise ValueError('Repeated fields, trajectory or counters differ')
    timing=[]
    for g in (4,8):
        v={m:[c['wall_seconds'] for c in cases if c['mode']==m and c['gate']==g] for m in ('vds','hdf5')}
        timing.append(dict(gate=g,seconds=v,overhead_percent=100*(sum(v['hdf5'])/sum(v['vds'])-1),
            paired_overheads_percent=[100*(h/b-1) for b,h in zip(v['vds'],v['hdf5'])]))
    total={m:sum(c['wall_seconds'] for c in cases if c['mode']==m) for m in ('vds','hdf5')}
    overhead=100*(total['hdf5']/total['vds']-1)
    return dict(numerical_integrity_pass=True,points=points,rounds=rounds,joint=joint,timing=timing,
        state_sizes=sizes,total_wall_seconds=total,overhead_percent=overhead,
        performance_qualified=rounds==3 and points==31 and overhead<=3.,limit_percent=3.)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--rounds',type=int,choices=(1,3),default=3)
    p.add_argument('--points',type=int,choices=(8,31),default=31)
    p.add_argument('--curve',choices=('vds','hdf5'));p.add_argument('--gate',type=int,choices=(4,8))
    a=p.parse_args();plan=read(a.plan)
    for path,digest in plan['frozen'].items():
        if sha(path)!=digest:raise ValueError('Frozen dependency changed: '+path)
    if a.curve:return one_curve(plan,a.curve,a.gate,a.output,a.points)
    a.output.mkdir(exist_ok=False)
    state=dict(status='running',pid=os.getpid(),planned_curves=4*a.rounds,points=a.points,cases=[],started=time.time())
    write(a.output/'batch.lock',dict(pid=os.getpid(),started=state['started']))
    try:
        for r in range(a.rounds):
            order=[('vds',4),('hdf5',4),('hdf5',8),('vds',8)]
            if r%2:order.reverse()
            for mode,gate in order:
                if (a.output/'STOP').exists():raise RuntimeError('Requested boundary stop')
                directory=a.output/f'r{r}_{mode}_vg{gate}'
                row=dict(repeat=r,mode=mode,gate=gate,directory=str(directory),status='running')
                state['cases'].append(row);write(a.output/'summary.json',state)
                cmd=[sys.executable,str(Path(__file__).resolve()),'--plan',str(a.plan.resolve()),'--output',str(directory.resolve()),
                     '--curve',mode,'--gate',str(gate),'--points',str(a.points)]
                with (a.output/f'r{r}_{mode}_vg{gate}.log').open('x') as log:
                    subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
                row.update(read(directory/'progress.json'));write(a.output/'summary.json',state)
                print('CURVE_DONE',r,mode,gate,row['wall_seconds'],flush=True)
        result=audit_batch(plan,a.output,state['cases'],a.rounds,a.points)
        write(a.output/'analysis.json',result)
        state.update(status='completed',numerical_integrity_pass=True,performance_qualified=result['performance_qualified'])
    except BaseException as error:
        state.update(status='failed',error=repr(error),traceback=traceback.format_exc());raise
    finally:
        state['finished']=time.time();write(a.output/'summary.json',state)


if __name__=='__main__':main()
