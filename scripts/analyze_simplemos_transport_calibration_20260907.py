"""Native continuity replay and attribution within the calibrated source direction."""
from pathlib import Path
import math
import numpy as np
import audit_simplemos_masetti_box_mobility_20260907 as m
import validate_simplemos_distributed_transport_20260907 as v
p=m.p;a=p.a;d=p.d;OUT=p.OUT

def main():
 a.verify(OUT/'box_mobility_scope.json');a.verify(v.OUT/'freeze.json')
 edges=a.rows(OUT/'native_sg_edges.csv');source=a.rows(v.OUT/'source_edges.csv');balance=[];contributions=[];summary=[]
 for c in a.read(v.OUT/'contract.json')['cases']:
  geo,xy,cells,vol,K,parts,info,mod=m.geom.geometry(c['device'])
  prior=next(x for x in a.read(p.prior.OUT/'vela_contract.json')['cases'] if x['key']==c['key'])
  fields=p.prior.prev.LOCAL/'native_exports/masetti'/prior['case']/'vg_020/fields'
  rec=m.scalar(fields/'SRHRecombination_region0.csv');mask=d.matrix.spatial.old.m78.supports(c['device'],geo,.05)[0]['all_si'].copy();mask[geo.contact_nodes]=False
  for mode in ('density','qf'):
   for carrier,field in [('electron','electron_particle_flux_per_cm_s'),('hole','hole_particle_flux_per_cm_s')]:
    terms=[[] for _ in range(geo.count)]
    for e in edges:
     if e['key']!=c['key'] or e['mode']!=mode:continue
     i,j=int(e['node0']),int(e['node1']);f=float(e[field]);terms[i].append(f);terms[j].append(-f)
    div=np.array([math.fsum(t) for t in terms]);R=np.array([rec.get(i,0)*vol['Si'][i]*1e4 for i in range(geo.count)])
    denom=prior['native_Id_A_per_um']/(m.Q*1e-4)
    balance.append(dict(key=c['key'],mode=mode,carrier=carrier,Si_free_nodes=int(sum(mask)),
      continuity_L1_over_terminal_particle_flux=float(sum(abs((div+R)[mask]))/denom),
      wrong_SRH_sign_L1_over_terminal_particle_flux=float(sum(abs((div-R)[mask]))/denom),
      interpretation='Independent reconstructed native continuity expression on 907 free Si nodes; not a replacement for the original nonlinear per-row acceptance.'))
  cfg=a.read(Path(c['base'])/'config.json');mesh=a.read(Path(cfg['mesh_file']));drain=set(next(x['node_ids'] for x in mesh['contacts'] if x['name']=='drain'))
  adj=d.ordered(p.prior.LOCAL/'vela'/c['key']/'adjoint.csv',geo.count);lam=np.array([float(r['lambda_electron']) for r in adj]);lam[geo.contact_nodes]=0.
  local=[]
  for e in source:
   if e['key']!=c['key']:continue
   i,j=int(e['node0']),int(e['node1']);feedback=-(lam[i]-lam[j])*float(e['source_normalized'])
   direct=-d.fixed.Q*1e-6*(int(i in drain)-int(j in drain))*float(e['source_particle_per_m_s'])
   local.append(dict(key=c['key'],edge=int(e['edge']),node0=i,node1=j,x_mid_um=(xy[i][0]+xy[j][0])/2,y_mid_um=(xy[i][1]+xy[j][1])/2,
     coefficient_ratio=float(e['native_over_vela_mu_geometry']),unit_feedback_A_per_um=feedback,unit_direct_A_per_um=direct,
     unit_total_A_per_um=feedback+direct,individual_edge_FD_qualified=False))
  total=math.fsum(r['unit_total_A_per_um'] for r in local)
  assert math.isclose(total,c['unit_prediction_A_per_um'],rel_tol=1e-12,abs_tol=1e-23)
  local.sort(key=lambda r:abs(r['unit_total_A_per_um']),reverse=True)
  for rank,r in enumerate(local,1):r['rank_absolute_contribution']=rank
  contributions.extend(local)
  summary.append(dict(key=c['key'],unit_dId_dalpha_A_per_um=total,unit_relative_Id_derivative=total/c['base_Id_A_per_um'],
    direct_fraction=c['unit_direct_A_per_um']/total,largest_edge=local[0]['edge'],largest_edge_fraction=local[0]['unit_total_A_per_um']/total,
    unit_alpha_recalculation_performed=False))
 a.write_csv(OUT/'native_continuity_replay.csv',balance);a.write_csv(OUT/'distributed_edge_contributions.csv',contributions);a.write_csv(OUT/'response_summary.csv',summary)
 a.write(OUT/'response_analysis_scope.json',dict(status='completed_continuity_and_calibrated_direction_analysis',
  accepted_amplitudes=[.001,.0005],unit_response='Derivative with respect to alpha at the frozen state; alpha=1 has not been simulated.',
  local_attribution='Edge contributions are an additive adjoint decomposition of a calibrated combined direction; individual edge experiments have not all been performed.',
  input_hashes={a.rel(f):a.sha(f) for f in [Path(__file__).resolve(),OUT/'box_mobility_scope.json',v.OUT/'freeze.json',v.OUT/'calibration.csv']}))
 for c in summary:print(c,flush=True)

if __name__=='__main__':main()
