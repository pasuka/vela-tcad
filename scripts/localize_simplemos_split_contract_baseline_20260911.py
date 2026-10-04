"""Attribute the failed baseline reconstruction at hole node 1087 without relaxing its gate."""
from decimal import Decimal as Q,localcontext
from pathlib import Path
import validate_simplemos_split_contract_v2_20260911 as v
F=lambda x:Q.from_float(float(x))
d=v.LOCAL/'states/n23/1.0/0';a=v.read(d/'operator.json');s=v.read(d/'results.csv.base.json');node=1087;N=len(a['nodes']);rows=v.rows(d/'results.csv');nodes={int(r['index']):r for r in rows if r['case']=='base' and r['kind']=='node'};edges={int(r['index']):r for r in rows if r['case']=='base' and r['kind']=='edge'}
with localcontext() as ctx:
 ctx.prec=100;scale=F(s['potential_scale_V']);x=[F(h)+F(l) for h,l in s['coordinates']];psi=[x[i]*scale for i in range(N)];qf=[x[2*N+i]*scale+F(s['hole_reference_V'][i]) for i in range(N)]
 def B(x):return Q(1) if x==0 else x/(x.exp()-1)
 def flux(e,mode):
  i,j=e['i'],e['j'];vt=F(a['nodes'][i][10]);n0=F(a['nodes'][i][9]);n1=F(a['nodes'][j][9]);p0,p1=psi[i],psi[j];f0,f1=qf[i],qf[j]
  if mode=='rounded_inputs':
   anchor=s['hole_reference_V'][i]
   p0=F(float(float(x[i])*float(scale))-anchor);p1=F(float(float(x[j])*float(scale))-anchor)
   f0=F(float(x[2*N+i])*float(scale));f1=F(float(F(s['hole_reference_V'][j])-F(anchor)+F(float(x[2*N+j])*float(scale))))
  eta=(p1-p0)/vt+(n0/n1).ln();mu=sum(F(w)*Q(nodes[k]['d']) for k,w in e['weights']);geo=F(e['coefficient'])
  if mode=='baseline_mobility':mu=F(e['baseline_mu'][1])
  return mu*geo*B(-eta)*n1*((f1-p1)/vt).exp()*(((f0-f1)/vt).exp()-1)
 sums={m:Q(0) for m in ('exact','rounded_inputs','baseline_mobility')};terms=[]
 for ei,e in enumerate(a['transport_edges']):
  if not e['weights'] or node not in (e['i'],e['j']):continue
  sign=1 if e['i']==node else -1;z={m:sign*flux(e,m) for m in sums}
  for m in sums:sums[m]+=z[m]
  terms.append(dict(edge=e['id'],node0=e['i'],node1=e['j'],**{m:str(z[m]) for m in z},cpp=str(sign*Q(edges[ei]['b']))))
 source=Q(nodes[node]['e'])*F(a['node_physics'][node][4]);new=Q(nodes[node]['h']);old=F(v.rows(d/'baseline.csv')[node]['phip_residual'])
 result=dict(node=node,packed_row=2*N+node,old=str(old),unified=str(new),source=str(source),difference=str(new-old),variants={m:dict(residual=str(z+source),difference_to_old=str(z+source-old),fraction_of_old_new_gap_remaining=str((z+source-old)/(new-old))) for m,z in sums.items()})
v.write(v.OUT/'hotspot/summary.json',result);v.csvout(v.OUT/'hotspot/edges.csv',terms);v.freeze(v.OUT/'hotspot/evidence.json',[Path(__file__).resolve(),d/'operator.json',d/'results.csv.base.json',d/'results.csv',d/'baseline.csv',v.OUT/'hotspot/summary.json',v.OUT/'hotspot/edges.csv']);print(result)
