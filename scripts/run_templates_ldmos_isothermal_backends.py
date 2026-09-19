"""Serial isothermal backend controls; frozen inputs, original gates, no overwrite."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from run_templates_ldmos_linked_d5 import read,write,digest,preflight,cpu_times,rows,state_difference
from analyze_templates_ldmos_stage4_d5 import analyze
from windows_system_load import cpu_snapshot,cpu_interval,wait_for_idle
ROOT=Path(__file__).resolve().parents[1]
PROFILES=ROOT/'reference_tcad/templates_ldmos_sentaurus2022/profiles'

def backend_threads(backend,requested):
 return requested if backend in ['mumps','mumps_metis','superlu_mt','superlu_mt_metis','strumpack'] else 1

def validate_resume_settings(report,args):
 for key in ['backends','threads','points','rounds','profiles','factor_statistics','idle_cpu_percent']:
  if report.get(key)!=getattr(args,key):raise ValueError('Resume settings changed: '+key)

def archive_interrupted_artifacts(out,name,stamp):
 """Preserve original bytes and record relocation without rewriting old ledgers."""
 moves=[]
 for item in [out/name,out/(name+'.log'),out/(name+'.json')]:
  if not item.exists():continue
  if not item.resolve().is_relative_to(out.resolve()):raise ValueError('Artifact outside output root')
  backup=item.with_name(item.name+'.interrupted_'+stamp)
  if backup.exists():raise ValueError('Interrupted backup exists')
  moves.append((item,backup))
 for item,backup in moves:item.rename(backup)
 return [dict(original=str(item),archived=str(backup)) for item,backup in moves]

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--output',type=Path,required=True)
 p.add_argument('--resume',action='store_true',help='Reuse verified completed cases and restart interrupted curves with frozen binaries')
 p.add_argument('--runner',type=Path,default=ROOT/'build-release/vela_example_runner.exe')
 p.add_argument('--replay',type=Path,default=ROOT/'build-release/linear_solver_replay.exe')
 p.add_argument('--backends',nargs='+',default=['sparselu','umfpack'])
 p.add_argument('--threads',type=int,choices=[1,2,4],default=1)
 p.add_argument('--points',type=int,choices=[8,31],default=8)
 p.add_argument('--rounds',type=int,default=1)
 p.add_argument('--profiles',nargs='+',choices=['D5','D4'],default=['D5','D4'])
 p.add_argument('--factor-statistics',choices=['on','off'],default='on')
 p.add_argument('--idle-cpu-percent',type=float,default=None)
 args=p.parse_args()
 if len(set(args.backends))!=len(args.backends) or len(args.backends)<2:p.error('Need at least two distinct backends')
 if args.rounds<1:p.error("--rounds must be positive")
 if args.idle_cpu_percent is not None and not 0<args.idle_cpu_percent<=100:p.error('Invalid CPU gate')
 if len(set(args.profiles))!=len(args.profiles):p.error('Duplicate profiles')
 if not __debug__ or os.name!='nt':raise RuntimeError('Requires assertions and Windows CPU accounting')
 out=args.output.resolve()
 if not args.resume:out.mkdir(parents=True,exist_ok=False)
 env=dict(os.environ,OMP_NUM_THREADS=str(args.threads),OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE',VELA_LINEAR_THREADS=str(args.threads))
 env['VELA_LINEAR_FACTOR_STATISTICS']='1' if args.factor_statistics=='on' else '0'
 env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env.get('PATH','')
 env.pop('GMON_OUT_PREFIX',None)
 env.pop('VELA_LINEAR_CAPTURE_DIR',None)
 flags=((out/'binary' if args.resume else ROOT/'build-release')/'build.ninja').read_text()
 assert '-O3 -DNDEBUG' in flags and '-pg ' not in flags and 'VELA_HAS_UMFPACK=1' in flags
 frozen={}
 binary=out/'binary'
 runner=binary/'vela_example_runner.exe';replay=binary/'linear_solver_replay.exe'
 bundles={x:PROFILES/f'linked_{x.lower()}_auger_no_generation_inputs.json' for x in args.profiles}
 if args.resume:
  report=read(out/'summary.json')
  validate_resume_settings(report,args)
  if digest(runner)!=report['runner_sha256']:raise ValueError('Frozen runner changed')
  frozen=read(binary/f'{args.backends[0]}.json')['frozen_sources']
  for path,sha in frozen.items():
   if digest(Path(path))!=sha:raise ValueError('Frozen source changed: '+path)
  for backend in args.backends:
   for path,sha in read(binary/f'{backend}.json')['runtime_sha256'].items():
    if digest(Path(path))!=sha:raise ValueError('Frozen runtime changed: '+path)
  for path,sha in read(out/'matrix_inputs.json').items():
   if digest(Path(path))!=sha:raise ValueError('Frozen matrix input changed: '+path)
  # Keep the pre-resume record and every interrupted artifact; never append
  # resumed elapsed time to an earlier partial performance measurement.
  stamp=time.strftime('%Y%m%d_%H%M%S')
  shutil.copy2(out/'summary.json',out/f'summary_before_resume_{stamp}.json')
  report.setdefault('resumptions',[]).append(dict(at=stamp,driver_sha256=digest(Path(__file__)),prior_status=report['status']))
  shutil.copy2(__file__,out/f'resume_driver_{stamp}.py')
  report.update(status='running');report.pop('error',None)
 else:
  names=subprocess.check_output(['git','-c','core.fsmonitor=false','ls-files','src','include','CMakeLists.txt'],cwd=ROOT,text=True).splitlines()
  names += [str(x.relative_to(ROOT)) for base in ['src','include'] for x in (ROOT/base).rglob('*') if x.is_file() and x.suffix in ('.h','.cpp')]
  names+=['src/tools/linear_solver_replay.cpp']+[str(x.relative_to(ROOT)) for x in (ROOT/'scripts').glob('*.py')]
  names += [str(x.relative_to(ROOT)) for x in (ROOT/'src/tools').iterdir() if x.suffix in ('.h','.cpp')]
  for name in sorted(set(names)):
   src=ROOT/name;dst=out/'source'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst);frozen[str(dst)]=digest(dst)
  binary.mkdir()
  shutil.copy2(args.runner,runner);shutil.copy2(args.replay,replay)
  for dll in args.runner.parent.glob('*.dll'):shutil.copy2(dll,binary/dll.name)
  for name in ['CMakeCache.txt','build.ninja']:shutil.copy2(ROOT/'build-release'/name,binary/name)
  report=dict(status='running',backends=args.backends,threads=args.threads,points=args.points,rounds=args.rounds,profiles=args.profiles,
   factor_statistics=args.factor_statistics,idle_cpu_percent=args.idle_cpu_percent,
   logical_processors=os.cpu_count(),power_scheme=subprocess.check_output(['powercfg','/getactivescheme']).decode(errors='replace'),
   runner_sha256=digest(runner),cases=[],matrix_cases=[],comparisons=[])
 for backend in args.backends:
  if not args.resume:
   manifest=dict(runner=str(runner),runner_sha256=digest(runner),backend=backend,linear_solver=backend,frozen_sources=frozen,build_type='UCRT64 Release -O3 -DNDEBUG',threads=backend_threads(backend,args.threads),factor_statistics=args.factor_statistics,
    runtime_sha256={str(x):digest(x) for x in binary.glob('*.dll')})
   write(binary/f'{backend}.json',manifest)
  for bundle in bundles.values():
   for gate in [4,8]:preflight(bundle,ROOT,binary/f'{backend}.json',gate)
 def save():write(out/'summary.json',report)
 def execute(name,argv,backend):
  if args.resume:
   archived=archive_interrupted_artifacts(out,name,stamp)
   report['resumptions'][-1].setdefault('archived_artifacts',[]).extend(archived)
  effective_threads=backend_threads(backend,args.threads)
  row=dict(name=name,argv=list(map(str,argv)),backend=backend,status='running',solver_threads=effective_threads,blas_threads=1);report['active']=name;report['active_elapsed_s']=0.;save()
  if args.idle_cpu_percent is not None:
   print('Waiting for CPU idle gate: '+name,flush=True)
   try:row['pre_run_load']=wait_for_idle(args.idle_cpu_percent)
   except TimeoutError as e:
    write(out/(name+'_load_gate_failure.json'),dict(samples=getattr(e,'samples',[]),error=str(e)))
    raise
  previous=cpu_snapshot();row['system_cpu_samples']=[];report['active_run']=row
  started=time.perf_counter()
  with (out/(name+'.log')).open('x') as log:
   proc=subprocess.Popen(list(map(str,argv)),cwd=ROOT,env=dict(env,VELA_LINEAR_SOLVER=backend,OMP_NUM_THREADS=str(effective_threads),VELA_LINEAR_THREADS=str(effective_threads)),stdout=log,stderr=subprocess.STDOUT)
   row['pid']=proc.pid
   report['active_pid']=proc.pid;save()
   while proc.poll() is None:
    try:proc.wait(timeout=20)
    except subprocess.TimeoutExpired:
     current=cpu_snapshot();row['system_cpu_samples'].append(cpu_interval(previous,current));previous=current
     report['active_elapsed_s']=time.perf_counter()-started;save()
     print(json.dumps(dict(active=name,elapsed_s=round(report['active_elapsed_s'],1))),flush=True)
   current=cpu_snapshot();row['system_cpu_samples'].append(cpu_interval(previous,current))
   row.update(returncode=proc.returncode,wall_seconds=time.perf_counter()-started,cpu_seconds=cpu_times(proc),status='pass' if proc.returncode==0 else 'failed')
  report.pop('active_run',None)
  return row
 try:
  captures=sorted((ROOT/'reference_staging/templates_ldmos_hotspot_execution_20260910/stage4/captures').rglob('*.bin'))
  assert len(captures)==14
  if not args.resume:write(out/'matrix_inputs.json',{str(x):digest(x) for x in captures})
  for r in range(2):
   for backend in (args.backends if r==0 else list(reversed(args.backends))):
    name=f'matrix_r{r}_{backend}'
    if any(x['name']==name and x['status']=='pass' for x in report['matrix_cases']):continue
    row=execute(name,[replay,backend,out/(name+'.json'),*captures],backend)
    report['matrix_cases'].append(row);save();assert row['returncode']==0,name
  for r in range(args.rounds):
   for profile in args.profiles:
    for gate in [4,8]:
     pair={}
     for backend in (args.backends[r%len(args.backends):]+args.backends[:r%len(args.backends)]):
      name=f'{profile.lower()}_vg{gate}_r{r}_{backend}';dest=out/name
      row=next((x for x in report['cases'] if x['name']==name and x['status']=='pass' and 'audit' in x),None)
      if row is None:
       old=[x for x in report['cases'] if x['name']==name]
       report.setdefault('interrupted_cases',[]).extend(old)
       report['cases']=[x for x in report['cases'] if x['name']!=name]
       row=execute(name,[sys.executable,out/'source/scripts/run_templates_ldmos_linked_d5.py','--workspace',ROOT,'--physics-profile',profile,'--bundle',bundles[profile],'--manifest',binary/f'{backend}.json','--output',dest,'--gate',gate,'--points',args.points,'--linear-solver',backend],backend)
       report['cases'].append(row);save();assert row['returncode']==0,name
      row['audit']=read(dest/'audit_summary.json');assert row['audit']['integrity_pass']
      assert row['audit']['exact_points']==args.points and read(dest/'progress.json')['status']=='completed'
      pair[backend]=dest;save()
     for candidate in args.backends[1:]:
      left,right=(read(pair[b]/'fixed/ledger.json') for b in [args.backends[0],candidate])
      assert len(left['exact_points'])==len(right['exact_points'])==args.points
      comps=[]
      for a,b in zip(left['exact_points'],right['exact_points']):
       assert a['bias_V']==b['bias_V'];delta=state_difference(Path(a['state']),Path(b['state']))
       assert all(delta[k]['max_absolute']<=1e-8 for k in ['psi','phin','phip']), (name,a['bias_V'],delta)
       comps.append(dict(bias_V=a['bias_V'],state_difference=delta))
      result=dict(profile=profile,gate=gate,round=r,baseline=args.backends[0],candidate=candidate,points=comps)
      report['comparisons']=[x for x in report['comparisons'] if not all(x.get(k)==result[k] for k in ['profile','gate','round','baseline','candidate'])]
      report['comparisons'].append(result);save()
    if args.points==31:
     for backend in args.backends:
      bundle=read(bundles[profile])
      group=f'{profile.lower()}_r{r}_{backend}'
      result=analyze(
       {'Vg'+str(g):ROOT/bundle['references'][str(g)] for g in [4,8]},
       {'Vg'+str(g):out/f'{profile.lower()}_vg{g}_r{r}_{backend}'/'score/curve.csv' for g in [4,8]},
       {'Vg'+str(g):out/f'{profile.lower()}_vg{g}_r{r}_{backend}'/'score/terminal_balance.csv' for g in [4,8]},
       out/(group+'_joint'),physics_profile=profile)
      report.setdefault('joint_qualifications',{})[group]=result;save()
      assert result['engineering']['status']==result['final']['status']=='pass',group
  for path,sha in frozen.items():assert digest(Path(path))==sha,path
  assert digest(runner)==report['runner_sha256']
  for backend in args.backends:
   for path,sha in read(binary/f'{backend}.json')['runtime_sha256'].items():assert digest(Path(path))==sha,path
  report['status']='pass'
 except BaseException as e:report.update(status='failed',error=repr(e));raise
 finally:
  report.pop('active',None);report.pop('active_pid',None);report.pop('active_elapsed_s',None);report.pop('active_run',None);save()
 print('BACKEND_CONTROLS_COMPLETE',flush=True)
if __name__=='__main__':main()
