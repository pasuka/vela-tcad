"""Nonzero sub-ULP terminal responses at actual drain-adjacent free nodes."""
from decimal import Decimal as Q,localcontext
from pathlib import Path
import math,subprocess
import validate_simplemos_split_contract_v2_20260911 as v
from analyze_simplemos_split_contract_v2_20260911 import reference
L=v.LOCAL/'port_micro';O=v.OUT/'port_micro';records=[];files=[]
v.verify(v.OUT/'final_kernel/build.json')
for device,vd,index in [('n19','0.05','40'),('n23','1.0','0')]:
 d=v.LOCAL/'states'/device/vd/index;model=v.read(d/'operator.json');state=v.read(d/'results.csv.base.json');N=len(model['nodes']);bc={i for i,z in model['boundary'][1]};node=None
 for edge in model['transport_edges']:
  if 'drain' not in edge['ports']:continue
  for i in (edge['i'],edge['j']):
   if i not in bc and all(state['coordinates'][b*N+i][0]!=0 for b in range(3)):node=i;break
  if node is not None:break
 assert node is not None
 dest=L/device/vd/index;dest.mkdir(parents=True,exist_ok=False);jobs=[dict(label='base',state=state)]
 for b in range(3):
  step=[0.]*(3*N);step[b*N+node]=math.ulp(state['coordinates'][b*N+node][0])/64;assert step[b*N+node]!=0
  jobs.append(dict(label=f'port_{b}',state=state,step=step,alpha=1.))
 v.write(dest/'jobs.json',jobs)
 p=subprocess.run([str(v.LOCAL/'final_kernel/driver.exe'),str(d/'operator.json'),str(dest/'jobs.json'),str(dest/'results.csv')],env=v.env(),capture_output=True,text=True);(dest/'run.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
 currents={r['case']:Q(r['a']) for r in v.rows(dest/'results.csv') if r['kind']=='port' and r['index']=='drain'}
 with localcontext() as ctx:
  ctx.prec=160;base=-reference(model,state,160)['ports']['drain']
  for b in range(3):
   label=f'port_{b}';actual=currents[label]-currents['base'];st=v.read(dest/('results.csv.'+label+'.json'));expected=-reference(model,st,160)['ports']['drain']-base
   assert actual!=0 and expected!=0
   error=abs(actual/expected-1);assert error<Q('1e-10')
   records.append(dict(device=device,vd=vd,index=index,node=node,block=b,step_ULP='1/64',delta_Id_A_per_um=str(actual),reference_delta=str(expected),relative=str(error)))
 files += [d/'operator.json']+list(dest.iterdir());print(device,vd,index,node,'passed',flush=True)
v.csvout(O/'responses.csv',records);v.freeze(O/'evidence.json',[Path(__file__).resolve(),v.ROOT/'scripts/analyze_simplemos_split_contract_v2_20260911.py',v.OUT/'final_kernel/build.json',O/'responses.csv']+files)
