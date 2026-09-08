"""Retain failed unsupported runtime reads; use only vertex/element fields."""
from pathlib import Path
import tarfile
import shutil
import prepare_simplemos_masetti_runtime_20260907 as p
a=p.a;LOCAL=p.LOCAL/'supported';OUT=p.OUT/'supported_runtime'
REMOTE=p.REMOTE+'/supported'

def main():
 a.verify(p.OUT/'native_freeze.json');files=[Path(__file__).resolve(),p.OUT/'native_freeze.json'];jobs=a.read(p.OUT/'native_contract.json')['jobs']
 script=p.TCL.replace(' edge $::des_data_edge [$mesh size_edge]','').replace(' element_vertex $::des_data_element_vertex [$mesh size_element_vertex]','')
 script=script.replace('[list edge $::des_data_edge [$mesh size_edge] element $::des_data_element [$mesh size_element]]','[list element $::des_data_element [$mesh size_element]]')
 assert '$data ReadScalar $::des_data_edge' not in script
 for job in jobs:
  src=p.LOCAL/'bundle'/job['key'];dest=LOCAL/'bundle'/job['key'];dest.mkdir(parents=True,exist_ok=False)
  for name in ('input_fps.tdr','result_020_des.sav','result_020_circuit_des.sav','native_des.cmd'):
   shutil.copyfile(src/name,dest/name);files += [src/name,dest/name]
  (dest/'runtime.tcl').write_text(script,newline='\n');files.append(dest/'runtime.tcl')
 shell='''#!/bin/bash
set -u
cd "$(dirname "$0")"
count=0
for path in bundle/*; do
 (cd "$path"; sdevice native_des.cmd > console.log 2>&1; printf '%s\\n' "$?" > exit_code.txt) &
 count=$((count+1))
 if test "$count" -eq 2; then wait; count=0; fi
done
wait
tar -czf results.tgz bundle
printf 'complete\\n' > complete.txt
'''
 (LOCAL/'run.sh').write_text(shell,newline='\n');files.append(LOCAL/'run.sh')
 a.write(OUT/'contract.json',dict(status='frozen_before_execution',jobs=jobs,remote=REMOTE,
  prior_failure='Default Element mobility exposes no Edge-RegionWise eMobility runtime dataset. First diagnostic attempts are retained and do not qualify native DC exports.',
  changes='Restrict scalar mobility to vertex/element and vector current to element. Explicit topology still includes element-vertex indices, measures and element-edge coefficients.',
  physics_solver_and_gates_unchanged=True))
 files.append(OUT/'contract.json');p.d.matrix.freeze(OUT/'freeze.json',files)
 with tarfile.open(LOCAL/'input.tgz','w:gz') as t:t.add(LOCAL/'bundle',arcname='bundle');t.add(LOCAL/'run.sh',arcname='run.sh')
 print('Frozen supported-field exports',flush=True)

if __name__=='__main__':main()
