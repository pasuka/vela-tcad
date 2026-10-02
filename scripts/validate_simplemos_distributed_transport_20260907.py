"""Freeze and calibrate a native-defined electron transport difference source."""
import argparse
import copy
import json
import math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess
import numpy as np
import audit_simplemos_masetti_box_mobility_20260907 as mobility
import build_simplemos_distributed_flux_20260907 as b
p=b.p;a=p.a;d=p.d;LOCAL=b.LOCAL;OUT=p.OUT/'distributed_response'

def prepare():
 a.verify(p.prior.OUT/'validation_evidence.json');a.verify(p.OUT/'box_mobility_scope.json')
 assert b.RUNNER.exists()
 cases=[];edge_rows=[];files=[Path(__file__).resolve(),Path(b.__file__).resolve(),b.RUNNER,p.OUT/'box_mobility_scope.json']
 for c in a.read(p.prior.OUT/'vela_contract.json')['cases']:
  if c['vg']!=1.:continue
  geo,weights,metrics=mobility.geometry(c['device']);root=LOCAL/c['key'];base=Path(c['base']);cfg=a.read(base/'config.json')
  edges=a.rows(p.prior.LOCAL/'vela'/c['key']/'strict/edges.csv');mesh=a.read(Path(cfg['mesh_file']))
  drain=set(next(x['node_ids'] for x in mesh['contacts'] if x['name']=='drain'))
  adj=d.ordered(p.prior.LOCAL/'vela'/c['key']/'adjoint.csv',geo.count);lam=np.array([float(x['lambda_electron']) for x in adj])
  rows={'combined':np.zeros(geo.count)};direct={'combined':[]};lines=[]
  for edge in edges:
   i,j=int(edge['node0']),int(edge['node1']);w=weights[(i,j)]['e']*1e-4
   g0=float(edge['couple_m'])/float(edge['length_m']);mu0=float(edge['electron_mobility_m2_V_s'])
   # Combined native per-element mobility/coefficient field, without inventing
   # a separately exported native edge mobility.
   if mu0<=0 or g0<=0:
    assert w==0;ratio=1.
   else:ratio=w/(g0*mu0)
   flux=float(edge['electron_flux'])*(ratio-1)
   physical=float(edge['electron_particle_line_flux_per_m_s'])*(ratio-1)
   native_units=physical/.01
   lines.append(f"{edge['edge_id']} {native_units:.17g}\n")
   rows['combined'][i]+=flux;rows['combined'][j]-=flux
   direct['combined'].append(-d.fixed.Q*1e-6*(int(i in drain)-int(j in drain))*physical)
   edge_rows.append(dict(key=c['key'],edge=int(edge['edge_id']),node0=i,node1=j,native_over_vela_mu_geometry=ratio,
     source_normalized=flux,source_particle_per_m_s=physical,source_injection_native_units=native_units))
  vector=rows['combined'].copy();vector[geo.contact_nodes]=0.
  feedback=-math.fsum(float(x*y) for x,y in zip(lam,vector));directI=math.fsum(direct['combined']);prediction=feedback+directI
  root.mkdir(parents=True,exist_ok=False);source=root/'source.txt';source.write_text(''.join(lines),newline='\n');files.append(source)
  jobs=[];current=a.read(base/'result.json')['Id_A_per_um']
  for label,mult in [('zero',0.),('plus_full',.001),('minus_full',-.001),('plus_half',.0005),('minus_half',-.0005)]:
   dest=root/label;deck=copy.deepcopy(cfg);deck['state_file']=str(base/'state.csv');deck['output_state_file']=str(dest/'state.csv')
   deck['solver']['local_update_diagnostics']['csv_file']=str(dest/'updates.csv');a.write(dest/'config.json',deck)
   q=copy.deepcopy(deck);q.pop('output_state_file');q['solver'].pop('local_update_diagnostics');q.update(simulation_type='newton_carrier_term_probe',state_file=str(dest/'state.csv'),output_csv=str(dest/'all_row.csv'),carrier_term_probe={'solved_equation_terms':True})
   q['solver']['carrier_row_convergence']['mode']='report';q['solver']['global_continuity_closure']=dict(mode='enforce',tolerance=1e-6,source_floor=1e-10);a.write(dest/'all_row.json',q)
   jobs.append(dict(label=label,scale=mult,config=str(dest/'config.json')));files += [dest/'config.json',dest/'all_row.json']
  for label in ('zero','unit'):
   q=copy.deepcopy(cfg);q.pop('output_state_file');q['solver'].pop('local_update_diagnostics');q['state_file']=str(base/'state.csv')
   q.update(simulation_type='terminal_current_functional_probe',contact='drain',residual_output_csv=str(root/'preflight'/label/'residual.csv'))
   a.write(root/'preflight'/label/'config.json',q);files.append(root/'preflight'/label/'config.json')
  cases.append(dict(key=c['key'],device=c['device'],base=str(base),vd=c['vd'],vg=1.,source=str(source),base_Id_A_per_um=current,
   unit_feedback_A_per_um=feedback,unit_direct_A_per_um=directI,unit_prediction_A_per_um=prediction,jobs=jobs))
  files += [base/'config.json',base/'state.csv',base/'result.json',p.prior.LOCAL/'vela'/c['key']/'strict/edges.csv',p.prior.LOCAL/'vela'/c['key']/'adjoint.csv']
 assert len(cases)==4
 a.write_csv(OUT/'source_edges.csv',edge_rows)
 a.write(OUT/'contract.json',dict(status='frozen_before_execution',cases=cases,new_DC=20,
  definition='For every edge at each qualified Vela state, delta electron SG flux = existing flux * (sum(native Si element coefficient * native box-weighted element mobility)/(Vela edge geometry * Vela mobility)-1). The difference is held state-independent during each solve.',
  direct='Same signed edge source in residual, term diagnostics and contact current. Explicit direct drain contribution is included; replaced contact rows are excluded from feedback.',
  amplitudes='plus/minus .001 and .0005 of the defined coefficient-difference source, with zero control.',
  frozen_gates=a.read(p.OUT/'native_contract.json')['eventual_response_gates'],
  exclusions='This is an electron distributed source calibration. It is not a self-consistent mobility-model replacement, a hole source calibration, or a Poisson correction.',
  production_changes=False,acceptance_changes=False))
 files += [OUT/'contract.json',OUT/'source_edges.csv']+[Path(x) for x in a.read(LOCAL/'build_command.json') if Path(x).is_file() and Path(x).is_relative_to(p.REPO)]
 d.matrix.freeze(OUT/'freeze.json',files);print('Frozen four distributed profiles and 20 DC solves',flush=True)

def execute(path,c,scale):
 target=path.with_suffix('.status.json')
 if target.exists():return a.read(target)
 r=subprocess.run([str(b.RUNNER),'--config',str(path),'--log','off'],env=b.env(path.parent,c['source'],scale),capture_output=True,text=True)
 path.with_suffix('.stdout.txt').write_text(r.stdout);path.with_suffix('.stderr.txt').write_text(r.stderr)
 s=json.loads(r.stdout.strip().splitlines()[-1]);s['exit_code']=r.returncode;a.write(target,s);return s

def run():
 a.verify(OUT/'freeze.json');source_rows=a.rows(OUT/'source_edges.csv')
 def one(c):
  root=LOCAL/c['key'];geo=d.matrix.spatial.m73.Geometry(c['device']);mask=d.matrix.spatial.old.m78.supports(c['device'],geo,.05)[0]['all_si'].copy();mask[geo.contact_nodes]=False
  zz=execute(root/'preflight/zero/config.json',c,0.);ff=execute(root/'preflight/unit/config.json',c,1.);assert zz['exit_code']==ff['exit_code']==0
  z=d.ordered(root/'preflight/zero/residual.csv',geo.count);f=d.ordered(root/'preflight/unit/residual.csv',geo.count);old=d.ordered(p.prior.LOCAL/'vela'/c['key']/'strict/residual.csv',geo.count)
  expect=np.zeros(geo.count);absrow=np.zeros(geo.count)
  for e in source_rows:
   if e['key']!=c['key']:continue
   i,j=int(e['node0']),int(e['node1']);flux=float(e['source_normalized']);expect[i]+=flux;expect[j]-=flux;absrow[i]+=abs(flux);absrow[j]+=abs(flux)
  expect[geo.contact_nodes]=0.
  for i,(x,y,o) in enumerate(zip(z,f,old)):
   for k in ('psi_residual','phin_residual','phip_residual'):assert x[k]==o[k]
   for k in ('psi_residual','phip_residual'):assert x[k]==y[k]
   actual=float(y['phin_residual'])-float(x['phin_residual'])
   assert abs(actual-expect[i])<=max(1e-8*absrow[i],1e-25)
  assert abs((ff['current_A_per_um']-zz['current_A_per_um'])-c['unit_direct_A_per_um'])<=max(1e-8*abs(c['unit_direct_A_per_um']),1e-14*abs(c['base_Id_A_per_um']))
  results=[]
  for job in c['jobs']:
   path=Path(job['config']);s=execute(path,c,job['scale']);q=execute(path.parent/'all_row.json',c,job['scale'])
   cc=s['contact_currents_A_per_um'];kcl=abs(math.fsum(cc.values()))/abs(cc['drain']);g=q['carrier_row_convergence']
   terms=d.ordered(path.parent/'all_row.csv',geo.count);ratios=[];zero=0
   for row,keep in zip(terms,mask):
    if not keep:continue
    for car in ('electron','hole'):
     scale=max(float(row[car+'_flux_abs_sum']),abs(float(row[car+'_recombination'])),abs(float(row[car+'_impact'])));zero+=scale==0
     ratios.append(abs(float(row[car+'_residual']))/scale if scale else math.inf)
   assert len(ratios)==g['qualified_row_count']==1814
   bad=sum(r>1e-6 for r in ratios);assert bad==g['violation_count']
   qualified=s['exit_code']==q['exit_code']==0 and s['converged'] and not bad and not zero and q['global_continuity_closure']['satisfied'] and kcl<=1e-8
   results.append(dict(key=c['key'],label=job['label'],current_A_per_um=cc['drain'],qualified=qualified,iterations=s['iterations'],row_violations=bad,zero_scale_rows=zero,max_row_ratio=max(ratios),kcl_over_Id=kcl,failure=s['failure_reason']))
  print(c['key'],'qualified',sum(r['qualified'] for r in results),'/5',flush=True);return results
 with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
 a.write_csv(OUT/'dc.csv',rows)

def analyze():
 a.verify(OUT/'freeze.json');rows=a.rows(OUT/'dc.csv');result=[]
 for c in a.read(OUT/'contract.json')['cases']:
  r={x['label']:x for x in rows if x['key']==c['key']};curr={n:float(x['current_A_per_um']) for n,x in r.items()};zero=curr['zero']
  odd={n:(curr['plus_'+n]-curr['minus_'+n])/2 for n in ('full','half')};linearity=abs(odd['full']/(2*odd['half'])-1)
  for name,amp in [('full',.001),('half',.0005)]:
   pred=amp*c['unit_prediction_A_per_um'];err=abs(odd[name]/pred-1);even=abs(((curr['plus_'+name]+curr['minus_'+name])/2-zero)/odd[name]);drift=abs(zero-c['base_Id_A_per_um']);snr=abs(odd[name])/max(drift,1e-300);dd=abs(math.log10(zero/c['base_Id_A_per_um']))
   signs=(curr['plus_'+name]-zero)*pred>0 and (curr['minus_'+name]-zero)*pred<0
   qualified=all(x['qualified']=='True' for x in r.values()) and err<=.001 and linearity<=.001 and even<=.01 and snr>=100 and dd<=1e-5 and signs
   result.append(dict(key=c['key'],amplitude=name,prediction_A_per_um=pred,odd_A_per_um=odd[name],prediction_relative_error=err,two_amplitude_relative=linearity,even_over_odd=even,signal_over_zero_drift=snr,zero_drift_dex=dd,signs_correct=signs,qualified=qualified))
 a.write_csv(OUT/'calibration.csv',result);print('Distributed response amplitudes qualified',sum(r['qualified'] for r in result),'/8',flush=True)

if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('action',choices=('prepare','run','analyze'));globals()[q.parse_args().action]()
