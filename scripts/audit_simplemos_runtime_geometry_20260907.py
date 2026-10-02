"""Map explicit native runtime topology and mobility to independently exported TDR."""
from pathlib import Path
import math
import tarfile
import numpy as np
from scipy.spatial import cKDTree
import prepare_simplemos_runtime_geometry_only_20260907 as r
import audit_simplemos_masetti_box_mobility_20260907 as m
p=r.p;a=p.a;LOCAL=r.LOCAL;OUT=r.OUT

def scalar(path):return {int(x['index']):float(x['value']) for x in a.rows(path)}

def main():
 a.verify(OUT/'freeze.json');dest=LOCAL/'native_raw'
 if not dest.exists():
  with tarfile.open(LOCAL/'results.tgz') as t:
   for item in t.getmembers():assert (dest/item.name).resolve().is_relative_to(dest.resolve()) and not item.issym() and not item.islnk()
   t.extractall(dest,filter='data')
 for f in (LOCAL/'bundle').rglob('*'):
  if f.is_file():assert a.sha(f)==a.sha(dest/'bundle'/f.relative_to(LOCAL/'bundle'))
 checks=[];mappings=[]
 for c in a.read(OUT/'contract.json')['jobs']:
  root=dest/'bundle'/c['key'];log=(root/'console.log').read_text(errors='replace')
  assert (root/'exit_code.txt').read_text().strip()=='0' and 'Good Bye' in log and 'T-2022.03-SP2' in log
  points=p.prior.prev.p.exporter.pltrows(root/'native_des.plt');assert points
  endpoint=points[-1];cc={n:endpoint[n+' TotalCurrent'] for n in ('drain','source','gate','substrate')}
  assert abs(endpoint['drain OuterVoltage']-c['vd'])<=1e-10 and abs(endpoint['gate OuterVoltage']-1.)<=1e-10
  drift=abs(cc['drain']/c['native_Id_A_per_um']-1);kcl=abs(math.fsum(cc.values()))/abs(cc['drain']);assert drift<=1e-8 and kcl<=1e-8
  snapshots=sorted(root.glob('runtime_*_meta.txt'));assert snapshots
  prefix=snapshots[-1].name.removesuffix('_meta.txt')
  geo,xy,cells,vol,K,parts,info,modified=m.geom.geometry(c['device']);coords=np.array([xy[i] for i in range(geo.count)])
  vertices=a.rows(root/(prefix+'_vertices.csv'));dist,idx=cKDTree(coords).query(np.array([[float(x['x_um']),float(x['y_um'])] for x in vertices]));assert max(dist)<=1e-12
  node_map={int(x['index']):int(j) for x,j in zip(vertices,idx)}
  edge_map={int(x['index']):tuple(sorted((node_map[int(x['start'])],node_map[int(x['end'])]))) for x in a.rows(root/(prefix+'_edges.csv'))}
  cell_lookup={tuple(sorted(x['nodes'])):i for i,x in cells.items() if x['material']=='Si'}
  ev=a.rows(root/(prefix+'_element_vertices.csv'));si={}
  for x in ev:
   if x['region']=='Silicon_1':si.setdefault(int(x['element']),[]).append(x)
  assert len(si)==1750
  map_cells={i:cell_lookup[tuple(sorted(node_map[int(x['vertex'])] for x in vs))] for i,vs in si.items()}
  raw=p.prior.LOCAL/'native_exports'/c['device'];nd=m.scalar(raw/'fields/DonorConcentration_region0.csv');na=m.scalar(raw/'fields/AcceptorConcentration_region0.csv')
  mu_error={};plot_error={}
  for carrier,pars in m.geom.mobility.PARAMETERS.items():
   mu=scalar(root/(prefix+'_element_'+carrier+'Mobility.csv'))
   plotted={int(x['cell_id']):float(x['component0']) for x in a.rows(raw/f'fields/{carrier}Mobility_region0_cells.csv')}
   ee=[];pp=[]
   for ci,vs in si.items():
    weights=[float(x['measure_um2']) for x in vs]
    values=[m.geom.mobility.formula(nd[node_map[int(x['vertex'])]]+na[node_map[int(x['vertex'])]],pars) for x in vs]
    calc=math.fsum(w*u for w,u in zip(weights,values))/math.fsum(weights)
    ee.append(abs(calc/mu[ci]-1));pp.append(abs(mu[ci]/plotted[map_cells[ci]]-1))
   mu_error[carrier]=max(ee);plot_error[carrier]=max(pp)
   assert max(ee)<=1e-8 and max(pp)<=1e-8
  expected={(x['cell'],edge):x['coefficient'] for edge,ps in parts.items() for x in ps if x['material']=='Si'}
  ce=[];native_shares=np.zeros(geo.count)
  for e in a.rows(root/(prefix+'_element_edges.csv')):
   ci=int(e['element'])
   if ci not in map_cells:continue
   pair=edge_map[int(e['edge'])];ce.append(abs(float(e['coefficient'])-expected[(map_cells[ci],pair)]))
  for vs in si.values():
   for row in vs:native_shares[node_map[int(row['vertex'])]]+=float(row['measure_um2'])*1e-12
  volume_error=max(abs(native_shares-vol['Si']))*1e12
  assert max(ce)<=1e-10 and volume_error<=1e-12
  checks.append(dict(key=c['key'],native_return_code=0,Id_relative_drift=drift,kcl_over_Id=kcl,snapshot=prefix,Si_cells=len(si),runtime_vertices=len(vertices),
   coordinate_error_um=float(max(dist)),runtime_to_debug_coefficient_max_absolute=max(ce),runtime_to_debug_node_volume_max_um2=float(volume_error),
   electron_box_mobility_max_relative=mu_error['e'],hole_box_mobility_max_relative=mu_error['h'],electron_runtime_to_plot_max_relative=plot_error['e'],hole_runtime_to_plot_max_relative=plot_error['h'],qualified=True))
  mappings.extend(dict(key=c['key'],runtime_cell=i,tdr_cell=j) for i,j in map_cells.items())
  print(c['key'],'runtime topology and mobility qualified',flush=True)
 a.write_csv(OUT/'checks.csv',checks);a.write_csv(OUT/'cell_mapping.csv',mappings)
 a.write(OUT/'evidence.json',dict(status='completed_native_runtime_geometry_mobility_identity',input_hashes={a.rel(f):a.sha(f) for f in [Path(__file__).resolve(),OUT/'freeze.json',LOCAL/'results.tgz']},
  edge_definition='No unsupported edge-mobility dataset used. Actual element-edge coefficients and element mobility are available; their sum defines the calibrated conservative SG coefficient.',
  vector_plots='No claim that a plotted current vector is an individual conservative edge flux.'))

if __name__=='__main__':main()
