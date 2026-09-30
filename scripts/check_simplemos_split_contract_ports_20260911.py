"""Audit sub-ULP edge/port responses and whole-DD merit on frozen directions."""
from decimal import Decimal as Q,localcontext
from pathlib import Path
import validate_simplemos_split_contract_v2_20260911 as v
from analyze_simplemos_split_contract_v2_20260911 import reference,F
checks=[];merits=[]
with localcontext() as ctx:
 ctx.prec=160
 for c in v.rows(v.OUT/'runs.csv'):
  d=Path(c['dest']);data=v.read(d/'operator.json');N=len(data['nodes']);group={}
  for r in v.rows(d/'results.csv'):group.setdefault(r['case'],{})[r['kind'],r['index']]=r
  refs={}
  for label in ('base','micro_0','micro_1','micro_2'):
   h=reference(data,v.read(d/('results.csv.'+label+'.json')),160)
   refs[label]={('edge',str(i),k):z for i,e in enumerate(h['edges']) for k,z in enumerate(e)}
   refs[label].update({('port',name,0):z for name,z in h['ports'].items()})
   if label=='base':continue
   for kind in ('edge','port'):
    active=[];changed=0
    for (kd,index,k),z in refs[label].items():
     if kd!=kind:continue
     delta=z-refs['base'][kd,index,k];actual=Q(group[label][kd,index]['ab'[k]])-Q(group['base'][kd,index]['ab'[k]])
     changed+=actual!=0
     if abs(delta)>Q('1e-85')*max(abs(actual),Q(1)):active.append(abs(actual-delta)/abs(delta))
    error=max(active,default=Q(0));assert error<Q('1e-10')
    checks.append(dict(device=c['device'],vd=c['vd'],index=c['index'],label=label,kind=kind,changed=changed,qualified_nonzero=len(active),max_relative=str(error)))
  old=v.rows(Path(c['old_dest'])/'trial_0.csv')
  for label in ('plus_full','plus_small'):
   block=[]
   for b in range(3):
    total=Q(0)
    for i in range(N):
     r=Q(group['base']['node',str(i)]['fgh'[b]]);t=Q(group[label]['node',str(i)]['fgh'[b]]);w=F(old[b*N+i]['weight']);total+=(t-r)*(t+r)*w*w
    block.append(total)
   merits.append(dict(device=c['device'],vd=c['vd'],index=c['index'],label=label,poisson=str(block[0]),electron=str(block[1]),hole=str(block[2]),total=str(sum(block)),descent=sum(block)<0))
  print(c['device'],c['vd'],c['index'],'completed',flush=True)
for name,rs in (('responses.csv',checks),('merit.csv',merits)):v.csvout(v.OUT/'response_extension'/name,rs)
v.freeze(v.OUT/'response_extension/evidence.json',[Path(__file__).resolve(),v.ROOT/'scripts/analyze_simplemos_split_contract_v2_20260911.py',v.OUT/'run_evidence.json',v.OUT/'response_extension/responses.csv',v.OUT/'response_extension/merit.csv'])
