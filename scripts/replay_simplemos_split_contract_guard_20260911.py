"""Rebuild the final frame-guarded kernel and require exact replay of all 54 jobs."""
from pathlib import Path
import subprocess
import validate_simplemos_split_contract_v2_20260911 as v
L=v.LOCAL/'frame_guard';O=v.OUT/'frame_guard'
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
 assert v.prior.sha(f/'results.csv')==v.prior.sha(d/'results.csv')
 for q in f.glob('results.csv.*.json'):assert v.prior.sha(q)==v.prior.sha(d/q.name)
 records.append(dict(device=c['device'],vd=c['vd'],index=c['index'],all_outputs_identical=True,jobs=9));files += list(f.iterdir());print(records[-1],flush=True)
v.csvout(O/'replay.csv',records);v.freeze(O/'replay.json',files+[O/'replay.csv'])
