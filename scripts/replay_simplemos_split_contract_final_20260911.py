"""Replay the signed-current repair; only port signs may change."""
from pathlib import Path
import subprocess
from decimal import Decimal as Q,getcontext
getcontext().prec=160
import validate_simplemos_split_contract_v2_20260911 as v
L=v.LOCAL/'final_kernel';O=v.OUT/'final_kernel'
L.mkdir(exist_ok=False)
src=L/'driver.cpp';src.write_text(v.DRIVER)
cmd=['D:/msys64/ucrt64/bin/c++.exe','-std=c++20','-O2','-I'+str(v.ROOT/'include'),str(src),'-o',str(L/'driver.exe')]
v.write(L/'compile.json',cmd)
p=subprocess.run(cmd,env=v.env(),capture_output=True,text=True);(L/'compile.log').write_text(p.stdout+p.stderr);assert p.returncode==0
v.freeze(O/'build.json',[Path(__file__).resolve(),src,L/'driver.exe',L/'compile.json',L/'compile.log',v.ROOT/'include/vela/numerics/SplitDDState.h',v.ROOT/'include/vela/equation/SplitDDOperator.h',v.ROOT/'include/vela/numerics/SplitCoordinate.h'])
records=[];files=[O/'build.json',v.OUT/'run_evidence.json']
for c in v.rows(v.OUT/'runs.csv'):
 d=Path(c['dest']);f=L/c['device']/c['vd']/c['index'];f.mkdir(parents=True)
 p=subprocess.run([str(L/'driver.exe'),str(d/'operator.json'),str(d/'jobs.json'),str(f/'results.csv')],env=v.env(),capture_output=True,text=True);(f/'run.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
 new=v.rows(f/'results.csv');old=v.rows(d/'results.csv');assert len(new)==len(old)
 for x,y in zip(new,old):
  if x['kind']=='port':
   assert (x['case'],x['kind'],x['index'])==(y['case'],y['kind'],y['index'])
   assert Q(x['a'])==-Q(y['a'])
  else:assert x==y
 current=next(Q(r['a']) for r in new if r['case']=='base' and r['kind']=='port' and r['index']=='drain')
 production=sum(Q.from_float(float(r['total_current_A_per_um'])) for r in v.rows(d/'contact.csv'))
 port_error=abs(current/production-1);assert port_error<Q('2e-12')
 for q in f.glob('results.csv.*.json'):assert v.prior.sha(q)==v.prior.sha(d/q.name)
 records.append(dict(device=c['device'],vd=c['vd'],index=c['index'],nonport_outputs_identical=True,ports_exactly_sign_reversed=True,jobs=9,drain_relative_to_production=str(port_error)));files += list(f.iterdir());print(records[-1],flush=True)
v.csvout(O/'replay.csv',records);v.freeze(O/'replay.json',files+[O/'replay.csv'])
