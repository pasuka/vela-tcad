"""Summarize actual shared-backend replays and completed linked controls."""
import argparse
import json
import statistics
from pathlib import Path
from run_templates_ldmos_linked_d5 import read,write,rows
from summarize_templates_ldmos_idvd_ablation import read_curve
from analyze_templates_ldmos_stage4_d5 import curve_error

def factor_statistics(profiles):
 """Weight observations by fresh factor count; never average per-process means."""
 names=set().union(*(d.get('observations',{}).keys() for d in profiles)) if profiles else set()
 result={}
 for name in sorted(names):
  if not name.startswith(('linear.numeric_factor_', 'linear.umfpack_', 'linear.sparselu_structural_', 'linear.mumps_', 'linear.strumpack_')):continue
  records=[d['observations'][name] for d in profiles if name in d.get('observations',{})]
  count=sum(x['count'] for x in records);total=sum(x['count']*x['average'] for x in records)
  result[name]=dict(count=count,total=total,average=total/count,min=min(x['min'] for x in records),max=max(x['max'] for x in records))
 peaks=[d.get('resources',{}).get('process_peak_working_set_bytes') for d in profiles]
 peaks=[x for x in peaks if x is not None]
 return dict(observations=result,max_solver_process_peak_working_set_bytes=max(peaks) if peaks else None,
  scope='OS process-lifetime high-water mark; includes inputs and assembly; maximum across child processes, not simultaneous process-tree memory',
  flop_semantics='SparseLU: scalar LU estimate from stored factor structure including supernode padding; UMFPACK: library-reported numerical operation count; neither is a hardware instruction counter')

def repeat_statistics(root,cases):
 groups={}
 for c in cases:
  parts=c['name'].split('_');key='_'.join(parts[:2]+parts[3:])
  groups.setdefault(key,[]).append(c)
 result=[]
 for key,items in groups.items():
  times=[x['wall_seconds'] for x in items]
  ledgers=[read(root/x['name']/'fixed/ledger.json') for x in items]
  baseline=ledgers[0]['exact_points']
  same=all([(p['bias_V'],p['sha256']) for p in d['exact_points']]==[(p['bias_V'],p['sha256']) for p in baseline] for d in ledgers)
  result.append(dict(case=key,samples=len(items),wall_seconds=times,median_wall_seconds=statistics.median(times),
   range_over_mean_percent=100*(max(times)-min(times))/statistics.mean(times),
   exact_state_hashes_identical=same,Newton_updates=[x['audit']['Newton_updates'] for x in items],
   same_target_and_update_history=all([(r['target_V'],r['Newton_updates']) for r in d['runs']]==[(r['target_V'],r['Newton_updates']) for r in ledgers[0]['runs']] for d in ledgers)))
 return result

def summarize(root):
 summary=read(root/'summary.json')
 matrix=[]
 for r in range(2):
  group={}
  for backend in ['sparselu','umfpack']:
   systems=read(root/f'matrix_r{r}_{backend}.json')['systems']
   total={stage:sum(next(t['total_ns'] for t in s['profiling']['stages'] if t['name']==stage)/1e9 for s in systems) for stage in ['linear.analyze','linear.factorize','linear.solve','linear.total']}
   group[backend]=dict(seconds=total,analyses=systems[-1]['cumulative_analyses'],
    factor_statistics=factor_statistics([s['profiling'] for s in systems]),
    worst_raw_normwise_backward_error=max(s['quality']['raw']['normwise_backward_error'] for s in systems),
    worst_raw_componentwise_backward_error=max(s['quality']['raw']['componentwise_backward_error'] for s in systems))
  matrix.append(dict(round=r,backends=group,total_linear_saving=1-group['umfpack']['seconds']['linear.total']/group['sparselu']['seconds']['linear.total'],numeric_saving=1-group['umfpack']['seconds']['linear.factorize']/group['sparselu']['seconds']['linear.factorize']))
 cases=[]
 for c in summary['cases']:
  if c['status']!='pass' or 'audit' not in c:continue
  dest=root/c['name'];plan=read(dest/'plan.json');ledger=read(dest/'fixed/ledger.json')
  profiles=[read(dest/r['case']/'performance_profile.json') for r in ledger['runs']]
  stages={k:sum(next((t['total_ns'] for t in d['stages'] if t['name']==k),0)/1e9 for d in profiles) for k in ['linear.total','linear.factorize','linear.solve','linear.analyze','linear.factor_statistics','dd.jacobian','dd.residual']}
  counters={k:sum(d['counters'].get(k,0) for d in profiles) for k in ['linear.solve_calls','linear.factorize_calls','linear.analyze_calls','linear.factorize_cache_hits','linear.analyze_cache_hits','linear.numeric_ordering_colamd','linear.numeric_ordering_amd','dd.jacobian_calls','dd.residual_calls']}
  curve=read_curve(dest/'score/curve.csv');ref=read_curve(Path(plan['reference']))
  errors=curve_error(ref[:len(curve)],curve)
  cases.append(dict(name=c['name'],wall_seconds=c['wall_seconds'],parent_cpu_seconds=c['cpu_seconds']['total'],audit=c['audit'],stages_seconds=stages,counters=counters,reference_error=errors,factor_statistics=factor_statistics(profiles)))
  samples=c.get('system_cpu_samples',[])
  if samples:
   duration=sum(x['interval_seconds'] for x in samples)
   total_cpu=c['cpu_seconds']['total']+c['audit']['child_cpu_seconds']
   cases[-1]['load']=dict(samples=len(samples),
    mean_busy_percent=sum(x['busy_percent']*x['interval_seconds'] for x in samples)/duration,
    max_busy_percent=max(x['busy_percent'] for x in samples),
    system_busy_cpu_seconds=sum(x['busy_cpu_seconds'] for x in samples),
    background_cpu_seconds_estimate=max(0.,sum(x['busy_cpu_seconds'] for x in samples)-total_cpu),
    process_tree_cpu_to_wall=total_cpu/c['wall_seconds'])
 pairs=[]
 for c in summary['comparisons']:
  stem=f"{c['profile'].lower()}_vg{c['gate']}_r{c['round']}"
  names=[stem+'_'+b for b in ['sparselu','umfpack']]
  left,right=[next(x for x in cases if x['name']==name) for name in names]
  ledgers=[read(root/name/'fixed/ledger.json') for name in names]
  target_sequences=[[r['target_V'] for r in d['runs']] for d in ledgers]
  curves=[read_curve(root/name/'score/curve.csv') for name in names]
  assert [v for v,i in curves[0]]==[v for v,i in curves[1]]
  old=root.parent/'templates_ldmos_r11_isothermal_20260917'/stem.split('_r')[0]/'fixed/ledger.json'
  old_points=read(old)['exact_points'] if old.exists() else []
  pairs.append(dict(profile=c['profile'],gate=c['gate'],round=c['round'],
   wall_saving_percent=100*(1-right['wall_seconds']/left['wall_seconds']),
   wall_saving_excluding_statistics_estimate_percent=100*(1-(right['wall_seconds']-right['stages_seconds']['linear.factor_statistics'])/(left['wall_seconds']-left['stages_seconds']['linear.factor_statistics'])),
   statistics_adjustment_scope='Arithmetic subtraction of measured diagnostic time, not a separate uninstrumented run',
   total_cpu_seconds=[x['parent_cpu_seconds']+x['audit']['child_cpu_seconds'] for x in [left,right]],
   max_potential_difference_V=max(p['state_difference'][k]['max_absolute'] for p in c['points'] for k in ['psi','phin','phip']),
   max_current_relative_difference=max(abs(b/a-1.) for (v,a),(w,b) in zip(*curves) if v>0 and a!=0),
   same_target_sequence=target_sequences[0]==target_sequences[1],
   sparselu_exact_hashes_match_prior_R11=(all(a['sha256']==b['sha256'] for a,b in zip(ledgers[0]['exact_points'],old_points)) if len(old_points)>=len(ledgers[0]['exact_points']) else None),
   changed_update_counts=[dict(target_V=a['target_V'],sparselu=a['Newton_updates'],umfpack=b['Newton_updates']) for a,b in zip(*[d['runs'] for d in ledgers]) if a['Newton_updates']!=b['Newton_updates']] if target_sequences[0]==target_sequences[1] else None))
 return dict(status=summary['status'],factor_statistics_mode=summary.get('factor_statistics','on'),
  matrix=matrix,cases=cases,pairs=pairs,state_comparisons=summary['comparisons'],repeats=repeat_statistics(root,cases),joint_qualifications=summary.get('joint_qualifications',{}))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);args=p.parse_args()
 result=summarize(args.root);write(args.root/'analysis.json',result)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
