"""Use only already successful runtime mobility reads plus documented geometry."""
import tarfile
import shutil
from pathlib import Path
import prepare_simplemos_masetti_runtime_20260907 as p
import prepare_simplemos_runtime_supported_20260907 as s
a=p.a;LOCAL=p.LOCAL/'geometry_runtime';OUT=p.OUT/'geometry_runtime';REMOTE=p.REMOTE+'/geometry_runtime'

def main():
 a.verify(s.OUT/'freeze.json');files=[Path(__file__).resolve(),s.OUT/'freeze.json'];jobs=a.read(s.OUT/'contract.json')['jobs']
 example=(s.LOCAL/'bundle'/jobs[0]['key']/'runtime.tcl').read_text()
 start=example.index(' foreach name {eDensity')
 end=example.index(' set f [open "${prefix}_vertices.csv"')
 script=example[:start]+example[end:]
 assert 'Potential' not in script and 'ReadVector' not in script and 'ReadFlux' not in script
 for job in jobs:
  src=s.LOCAL/'bundle'/job['key'];dest=LOCAL/'bundle'/job['key'];dest.mkdir(parents=True,exist_ok=False)
  for name in ('input_fps.tdr','result_020_des.sav','result_020_circuit_des.sav','native_des.cmd'):
   shutil.copyfile(src/name,dest/name);files += [src/name,dest/name]
  (dest/'runtime.tcl').write_text(script,newline='\n');files.append(dest/'runtime.tcl')
 shutil.copyfile(s.LOCAL/'run.sh',LOCAL/'run.sh');files.append(LOCAL/'run.sh')
 a.write(OUT/'contract.json',dict(status='frozen_before_execution',jobs=jobs,remote=REMOTE,
  prior_read_failures=['Edge-RegionWise eMobility is undefined in this Element averaging configuration.','Plot alias Potential is not accepted by runtime ReadScalar.'],
  retained_scope='Only the already successful vertex/element e/h mobility reads, explicit vertices/edges/element connectivity, native coefficients and measures. Native field values are available in the standard TDR plots.',
  no_data_write_api=True,physics_solver_and_gates_unchanged=True))
 files.append(OUT/'contract.json');p.d.matrix.freeze(OUT/'freeze.json',files)
 with tarfile.open(LOCAL/'input.tgz','w:gz') as t:t.add(LOCAL/'bundle',arcname='bundle');t.add(LOCAL/'run.sh',arcname='run.sh')
 print('Frozen geometry and known mobility runtime reads',flush=True)

if __name__=='__main__':main()
