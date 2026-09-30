"""Same-perturbation native controls for the Enormal inverse-scattering scale."""
import tarfile,shutil
from pathlib import Path
import calibrate_simplemos_enormal_20260912 as p
a,d=p.a,p.d;R=p.R
L=R/'build-release/enormal_response_20260912';O=R/'reference_tcad/simplemos_sentaurus2022/enormal_response_20260912'
REMOTE='/tmp/vela_simplemos_enormal_response_20260912'

def prepare():
    a.verify(p.O/'floor_position_evidence.json');s=a.read(p.O/'summary.json')
    assert s['identity_passed']==8 and s['observer_final_state_exact']
    assert a.read(p.O/'floor_position_summary.json')['left_max_relative']<=1e-7
    assert not (O/'native_freeze.json').exists()
    jobs=[];files=[Path(__file__).resolve(),p.O/'floor_position_evidence.json',p.L/'vendor_reference/enormal_native_models.par']
    for c in a.read(p.O/'native_contract.json')['jobs']:
        if c['arm']!='baseline' or c['index']!=40:continue
        source=p.L/'native_raw/bundle'/c['name']
        for label,delta in [('zero',0.),('minus_large',-.001),('minus_small',-.0005),('plus_small',.0005),('plus_large',.001)]:
            key=c['name'].removesuffix('_baseline')+'_'+label;dest=L/'bundle'/key;dest.mkdir(parents=True,exist_ok=False)
            header=(source/'native_des.cmd').read_text().split('Solve {')[0].replace('File {','File { Parameter="scale.par"',1)
            (dest/'native_des.cmd').write_text(header+'Solve {\n Load(FilePrefix="anchor")\n Coupled { Poisson Electron Hole }\n Plot(FilePrefix="final")\n Save(FilePrefix="final")\n}\n',newline='\n')
            factor=1+delta
            (dest/'scale.par').write_text('Material="Silicon" {\n EnormalDependence {\n'+f'  a_ac = {factor:.17g}, {factor:.17g}\n  a_sr = {factor:.17g}, {factor:.17g}\n'+' }\n}\n',newline='\n')
            for old,new in [('input_fps.tdr','input_fps.tdr'),('final_des.sav','anchor_des.sav'),('final_circuit_des.sav','anchor_circuit_des.sav')]:shutil.copyfile(source/old,dest/new);files.append(source/old)
            files.extend(dest.iterdir());jobs.append(dict(case=c['case'],device=c['device'],vd=c['vd'],vg=c['vg'],index=40,name=key,label=label,delta=delta,baseline=c['name']))
    for stage in ('zero','signed'):
        names=[j['name'] for j in jobs if (j['label']=='zero')==(stage=='zero')]
        shell='#!/bin/bash\nset -u\ncd "$(dirname "$0")"\nprintf "%s\\n" "$$" > STAGE_launcher.pid\nrunning=0\nfor path in '+ ' '.join('bundle/'+x for x in names)+'''; do
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
    a.write(O/'native_contract.json',dict(jobs=jobs,remote_root=REMOTE,parameter_direction='Same multiplier 1+delta on both a_ac and a_sr, both carriers. Keep bulk PhuMob, all other physical and numerical parameters fixed.',deltas=[-.001,-.0005,0,.0005,.001],stages='Four explicit zero controls first, then sixteen signed controls after zero identity.',gates=dict(native_exit_zero=True,kcl_over_Id=1e-8,bias_V=1e-10,zero_Id_relative=1e-12,zero_potential_V=1e-12,zero_density_relative=1e-10,two_amplitude_relative=1e-3,even_over_odd=.01,signal_over_zero_drift=100),original_acceptance_changed=False))
    d.matrix.freeze(O/'native_freeze.json',files+[O/'native_contract.json'])
    with tarfile.open(L/'input.tgz','w:gz') as t:
        for name in ('bundle','run_zero.sh','run_signed.sh'):t.add(L/name,arcname=name)
    print(L/'input.tgz',a.sha(L/'input.tgz'),flush=True)
if __name__=='__main__':prepare()
