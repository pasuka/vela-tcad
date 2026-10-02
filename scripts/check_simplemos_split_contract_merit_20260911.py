"""Whole-DD squared-merit difference using frozen weight / scale**2."""
from decimal import Decimal as Q,localcontext
from pathlib import Path
import validate_simplemos_split_contract_v2_20260911 as v
from analyze_simplemos_split_contract_v2_20260911 import F
result=[];ports=[]
with localcontext() as ctx:
 ctx.prec=160
 for c in v.rows(v.OUT/'runs.csv'):
  d=Path(c['dest']);N=len(v.read(d/'operator.json')['nodes']);g={}
  for r in v.rows(d/'results.csv'):g.setdefault(r['case'],{})[r['kind'],r['index']]=r
  old=v.rows(Path(c['old_dest'])/'trial_0.csv')
  for label in ('plus_full','plus_small'):
   blocks=[]
   for b in range(3):
    total=Q(0)
    for i in range(N):
     r=Q(g['base']['node',str(i)]['fgh'[b]]);t=Q(g[label]['node',str(i)]['fgh'[b]]);row=old[b*N+i]
     total+=(t-r)*(t+r)*F(row['weight'])/F(row['scale'])**2
    blocks.append(total)
   result.append(dict(device=c['device'],vd=c['vd'],index=c['index'],label=label,poisson=str(blocks[0]),electron=str(blocks[1]),hole=str(blocks[2]),total=str(sum(blocks)),descent=sum(blocks)<0))
  previous=sum(F(r['total_current_A_per_um']) for r in v.rows(d/'contact.csv'));current=Q(g['base']['port','drain']['a'])
  ports.append(dict(device=c['device'],vd=c['vd'],index=c['index'],old_drain_A_per_um=str(previous),new_drain_A_per_um=str(current),relative_difference=str((current-previous)/previous)))
for name,rs in (('merit.csv',result),('ports.csv',ports)):v.csvout(v.OUT/'merit_correct'/name,rs)
v.freeze(v.OUT/'merit_correct/evidence.json',[Path(__file__).resolve(),v.OUT/'run_evidence.json',v.OUT/'merit_correct/merit.csv',v.OUT/'merit_correct/ports.csv'])
print('descent',sum(r['descent'] for r in result),len(result));print(ports)
