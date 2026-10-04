"""Corrected native inverse-scattering direction using B,C,delta,eta."""
import shutil,tarfile
from pathlib import Path
import simplemos_enormal_response_20260912 as old
a,d,p=old.a,old.d,old.p
L=old.L/'material_scale_v2';O=old.O/'material_scale_v2'
REMOTE='/tmp/vela_simplemos_enormal_material_scale_20260912'

def prepare():
    a.verify(old.O/'failure_evidence.json');a.verify(old.O/'native_freeze.json');assert not (O/'native_freeze.json').exists()
    jobs=a.read(old.O/'native_contract.json')['jobs'];files=[Path(__file__).resolve(),old.O/'failure_evidence.json']
    for job in jobs:
        source=old.L/'bundle'/job['name'];dest=L/'bundle'/job['name'];dest.mkdir(parents=True,exist_ok=False)
        for name in ('input_fps.tdr','anchor_des.sav','anchor_circuit_des.sav','native_des.cmd'):
            shutil.copyfile(source/name,dest/name);files.extend([source/name,dest/name])
        factor=1+job['delta'];en,hp=p.PARAMS
        text='Material="Silicon" {\n EnormalDependence {\n'
        for name,index in [('B',0),('C',1),('delta',3),('eta',4)]:text+=f'  {name} = {en[index]/factor:.17g}, {hp[index]/factor:.17g}\n'
        text+=' }\n}\n';(dest/'scale.par').write_text(text,newline='\n');files.append(dest/'scale.par')
    pilot_case=jobs[0]['case']
    for stage in ('pilot','rest'):
        selected=[j['name'] for j in jobs if (j['case']==pilot_case)==(stage=='pilot')]
        shell='#!/bin/bash\nset -u\ncd "$(dirname "$0")"\nprintf "%s\\n" "$$" > STAGE_launcher.pid\nrunning=0\nfor path in '+' '.join('bundle/'+x for x in selected)+'''; do
 (cd "$path" || exit 91
  test ! -f exit_code.txt || exit 92
  /atctools/Synopsys/tcad/T-2022.03/bin/sdevice native_des.cmd > console.log 2>&1 &
  pid=$!; printf '%s\\n' "$pid" > solver.pid
  wait "$pid"; code=$?; printf '%s\\n' "$code" > exit_code.txt
 ) &
 running=$((running+1))
 if [ "$running" -eq 2 ]; then wait; running=0; fi
done
wait
printf 'finished; inspect per-job exits\\n' > STAGE_complete.txt
'''
        (L/f'run_{stage}.sh').write_text(shell.replace('STAGE',stage),newline='\n');files.append(L/f'run_{stage}.sh')
    contract=a.read(old.O/'native_contract.json');contract.update(remote_root=REMOTE,pilot_case=pilot_case,stages='Five points for the first case must show a nonzero qualified two-amplitude response before launching the other fifteen.',parameter_direction='For f=1+delta divide B,C,delta,eta by f, both carriers; inverse acoustic and roughness scattering are then multiplied by f. Native a_ac/a_sr stay at defaults. No stress model is enabled.',corrects='Prior a_ac/a_sr experiment was a zero-signal failure, preserved separately.')
    a.write(O/'native_contract.json',contract);d.matrix.freeze(O/'native_freeze.json',files+[O/'native_contract.json'])
    with tarfile.open(L/'input.tgz','w:gz') as t:
        for name in ('bundle','run_pilot.sh','run_rest.sh'):t.add(L/name,arcname=name)
    print(L/'input.tgz',a.sha(L/'input.tgz'),flush=True)
if __name__=='__main__':prepare()
