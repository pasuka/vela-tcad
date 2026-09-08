"""Separate native semiconductor geometry and subsequent mobility differences."""
import argparse
import copy
import math
from pathlib import Path
import numpy as np
import validate_simplemos_distributed_transport_20260907 as v
p=v.p;a=p.a;d=p.d
PARENT_OUT=v.OUT;PARENT_LOCAL=v.LOCAL
LOCAL=p.LOCAL/'factor_response';OUT=p.OUT/'factor_response'

def prepare():
 a.verify(PARENT_OUT/'freeze.json')
 original=a.read(PARENT_OUT/'contract.json')
 geometry_cache={dev:v.mobility.geom.geometry(dev) for dev in ('n19','n23')}
 for kind in ('geometry','mobility'):
  cases=[];edge_rows=[];files=[Path(__file__).resolve(),Path(v.__file__).resolve(),PARENT_OUT/'freeze.json',v.b.RUNNER]
  for c0 in original['cases']:
   c=copy.deepcopy(c0);base=Path(c['base']);root=LOCAL/kind/c['key'];root.mkdir(parents=True,exist_ok=True);assert not any(root.iterdir())
   geo,xy,cells,vol,K,parts,info,mod=geometry_cache[c['device']]
   mesh=a.read(Path(a.read(base/'config.json')['mesh_file']));drain=set(next(x['node_ids'] for x in mesh['contacts'] if x['name']=='drain'))
   adj=d.ordered(p.prior.LOCAL/'vela'/c['key']/'adjoint.csv',geo.count);lam=np.array([float(x['lambda_electron']) for x in adj]);lam[geo.contact_nodes]=0.
   full={int(x['edge']):x for x in a.rows(PARENT_OUT/'source_edges.csv') if x['key']==c['key']};rhs=np.zeros(geo.count);direct=[];lines=[]
   for edge in a.rows(p.prior.LOCAL/'vela'/c['key']/'strict/edges.csv'):
    i,j=int(edge['node0']),int(edge['node1']);g0=float(edge['couple_m'])/float(edge['length_m'])
    gn=math.fsum(x['coefficient'] for x in parts[(i,j)] if x['material']=='Si')
    if float(edge['electron_mobility_m2_V_s'])<=0 or g0<=0:
     assert gn==0;gr=1.
    else:assert g0>0;gr=gn/g0
    fr=float(full[int(edge['edge_id'])]['native_over_vela_mu_geometry'])
    delta=(gr-1) if kind=='geometry' else (fr-gr)
    flux=float(edge['electron_flux'])*delta;physical=float(edge['electron_particle_line_flux_per_m_s'])*delta;nu=physical/.01
    rhs[i]+=flux;rhs[j]-=flux;direct.append(-d.fixed.Q*1e-6*(int(i in drain)-int(j in drain))*physical)
    lines.append(f"{edge['edge_id']} {nu:.17g}\n")
    edge_rows.append(dict(key=c['key'],edge=int(edge['edge_id']),node0=i,node1=j,native_over_vela_mu_geometry=1+delta,
      source_normalized=flux,source_particle_per_m_s=physical,source_injection_native_units=nu))
   rhs[geo.contact_nodes]=0.;feedback=-math.fsum(float(x*y) for x,y in zip(lam,rhs));di=math.fsum(direct)
   c.update(source=str(root/'source.txt'),unit_feedback_A_per_um=feedback,unit_direct_A_per_um=di,unit_prediction_A_per_um=feedback+di)
   (root/'source.txt').write_text(''.join(lines),newline='\n');files.append(root/'source.txt')
   for job in c['jobs']:
    old=Path(job['config']);dest=root/job['label'];cfg=a.read(old);cfg['output_state_file']=str(dest/'state.csv');cfg['solver']['local_update_diagnostics']['csv_file']=str(dest/'updates.csv');a.write(dest/'config.json',cfg)
    q=a.read(old.parent/'all_row.json');q['state_file']=str(dest/'state.csv');q['output_csv']=str(dest/'all_row.csv');a.write(dest/'all_row.json',q)
    job['config']=str(dest/'config.json');files += [dest/'config.json',dest/'all_row.json']
   for label in ('zero','unit'):
    old=PARENT_LOCAL/c['key']/'preflight'/label/'config.json';cfg=a.read(old);cfg['residual_output_csv']=str(root/'preflight'/label/'residual.csv');a.write(root/'preflight'/label/'config.json',cfg);files.append(root/'preflight'/label/'config.json')
   cases.append(c)
  a.write_csv(OUT/kind/'source_edges.csv',edge_rows)
  contract=copy.deepcopy(original);contract.update(cases=cases,factor=kind,
   definition=('Change g_all to native summed Si g, while holding Vela edge mobility.' if kind=='geometry' else 'After the native Si geometry change, change Vela edge mobility to the native coefficient-weighted element mobility. This source plus geometry source exactly equals the previously calibrated combined source.'))
  a.write(OUT/kind/'contract.json',contract);files += [OUT/kind/'contract.json',OUT/kind/'source_edges.csv'];d.matrix.freeze(OUT/kind/'freeze.json',files)
  print('Frozen',kind,'20 states',flush=True)

def run():
 for kind in ('geometry','mobility'):
  v.OUT=OUT/kind;v.LOCAL=LOCAL/kind;v.run();v.analyze()

if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('action',choices=('prepare','run'));globals()[q.parse_args().action]()
