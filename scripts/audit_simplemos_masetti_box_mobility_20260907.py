"""Calibrate box-share mobility and native SG flux; no fitted parameters."""
from pathlib import Path
import math
import numpy as np
from decimal import Decimal,localcontext
import prepare_simplemos_masetti_runtime_20260907 as p
import audit_simplemos_masetti_native_geometry_20260907 as geom
a=p.a;d=p.d;prior=p.prior;OUT=p.OUT;LOCAL=p.LOCAL
Q=1.602192e-19;VT=1.380662e-23*300/Q

def scalar(path):return geom.sc(path)

def geometry(device):
 geo,xy,cells,vol,K,edgeparts,info,mod=geom.geometry(device)
 debug=(prior.LOCAL/'native_raw/bundle'/device/'MeasureCoefficients.debug').read_text()
 measure=geom.box.parse_debug_block(debug,'Measure');vp=info['measure_permutation']
 src=prior.LOCAL/'native_exports'/device;nd=scalar(src/'fields/DonorConcentration_region0.csv');na=scalar(src/'fields/AcceptorConcentration_region0.csv')
 mus={};metrics=[]
 for carrier,pars in geom.mobility.PARAMETERS.items():
  native={int(x['cell_id']):float(x['component0']) for x in a.rows(src/f'fields/{carrier}Mobility_region0_cells.csv')};mm={};errors=[]
  for i,cell in cells.items():
   if cell['material']!='Si':continue
   values=[geom.mobility.formula(nd[n]+na[n],pars) for n in cell['nodes']]
   shares=[measure[i]['values'][vp[j]] for j in range(3)]
   value=math.fsum(w*mu for w,mu in zip(shares,values))/math.fsum(shares)
   mm[i]=value;errors.append(abs(value/native[i]-1))
  mus[carrier]=mm;metrics.append(dict(device=device,carrier=carrier,cells=len(mm),maximum_relative_error=max(errors),qualified=max(errors)<=1e-8))
 weights={edge:{car:math.fsum(x['coefficient']*mus[car][x['cell']] for x in parts if x['material']=='Si') for car in mus} for edge,parts in edgeparts.items()}
 return geo,weights,metrics

def B(x):
 if abs(x)<Decimal('1e-20'):return 1-x/2+x*x/12
 return x/(x.exp()-1)

def sg_flux(psi0,psi1,n0,n1,p0,p1,fn0,fn1,fp0,fp1,mode):
 with localcontext() as ctx:
  ctx.prec=60;D=lambda x:Decimal(str(x));vt=D(VT);z=(D(psi1)-D(psi0))/vt
  if mode=='density':return vt*(D(n0)*B(-z)-D(n1)*B(z)),vt*(D(p0)*B(z)-D(p1)*B(-z))
  # Exact constant-ni electrochemical SG identity, evaluated without subtracting
  # nearly equal drift/diffusion terms. Native qf sign is checked by port sums.
  return vt*D(n0)*B(-z)*(1-((D(fn0)-D(fn1))/vt).exp()),vt*D(p0)*B(z)*(1-((D(fp1)-D(fp0))/vt).exp())

def main():
 a.verify(prior.OUT/'validation_evidence.json');allmu=[];checks=[];rows=[]
 for device in ('n19','n23'):
  geo,weights,metrics=geometry(device);allmu.extend(metrics)
  assert all(x['qualified'] for x in metrics)
  for c in a.read(prior.OUT/'vela_contract.json')['cases']:
   if c['device']!=device or c['vg']!=1.:continue
   fields=prior.prev.LOCAL/'native_exports/masetti'/c['case']/'vg_020/fields'
   raw={n:scalar(fields/(field+'_region0.csv')) for n,field in [('psi','ElectrostaticPotential'),('n','eDensity'),('p','hDensity'),('fn','eQuasiFermiPotential'),('fp','hQuasiFermiPotential')]}
   cfg=a.read(Path(c['config']));mesh=a.read(Path(cfg['mesh_file']));edges=a.rows(prior.LOCAL/'vela'/c['key']/'strict/edges.csv')
   contacts={x['name']:set(map(int,x['node_ids'])) for x in mesh['contacts']}
   for mode in ('density','qf'):
    fluxes={};current={contact:[] for contact in contacts}
    for edge in edges:
     i,j=int(edge['node0']),int(edge['node1']);w=weights[(i,j)]
     if i not in raw['n'] or j not in raw['n']:en=hp=0.
     else:
      en,hp=sg_flux(raw['psi'][i],raw['psi'][j],raw['n'][i],raw['n'][j],raw['p'][i],raw['p'][j],raw['fn'][i],raw['fn'][j],raw['fp'][i],raw['fp'][j],mode)
      en=float(en)*w['e'];hp=float(hp)*w['h']
     fluxes[int(edge['edge_id'])]=(en,hp)
     for name,ids in contacts.items():current[name].append(-Q*1e-4*(int(i in ids)-int(j in ids))*(en-hp))
     rows.append(dict(key=c['key'],mode=mode,edge=int(edge['edge_id']),node0=i,node1=j,
      native_mu_geometry_e_cm2_Vs=w['e'],native_mu_geometry_h_cm2_Vs=w['h'],electron_particle_flux_per_cm_s=en,hole_particle_flux_per_cm_s=hp))
    cc={name:math.fsum(v) for name,v in current.items()};error=cc['drain']/c['native_Id_A_per_um']-1
    checks.append(dict(key=c['key'],mode=mode,Id_A_per_um=cc['drain'],target_A_per_um=c['native_Id_A_per_um'],signed_Id_relative_error=error,
      kcl_over_Id=abs(math.fsum(cc.values()))/abs(cc['drain']),port_qualified=abs(error)<=1e-6))
 a.write_csv(OUT/'box_mobility_identity.csv',allmu);a.write_csv(OUT/'native_sg_edges.csv',rows);a.write_csv(OUT/'native_sg_port_checks.csv',checks)
 a.write(OUT/'box_mobility_scope.json',dict(status='completed_box_weighted_mobility_and_SG_check',
  mobility='sum(native element-vertex box measures * Masetti(node total impurity)) / sum(native element-vertex box measures). Never replace denominator with coordinate triangle area.',
  physical_parameters_fitted=False,q_C=Q,thermal_voltage_V=VT,
  flux='Sum per-element native coefficient times element mobility; SG in cm units, conventional terminal current with A/cm to A/um factor 1e-4.',
  input_hashes={a.rel(x):a.sha(x) for x in [Path(__file__).resolve(),prior.OUT/'validation_evidence.json',p.OUT/'native_contract.json']}))
 for r in checks:print(r['key'],r['mode'],'Id error',r['signed_Id_relative_error'],'KCL',r['kcl_over_Id'],flush=True)

if __name__=='__main__':main()
