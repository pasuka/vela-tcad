"""Freeze four native 0..1 V field sweeps; no implicit remote execution."""
from pathlib import Path
import shutil
import tarfile
import decompose_simplemos_calibrated_transport as d

a=d.a
LOCAL=d.REPO/'build-release/simplemos_fullfield_20260906'
OUT=d.ROOT/'fullfield_validation_20260906'
REMOTE='/tmp/vela_simplemos_fullfield_20260906'


def main():
    a.verify(d.FREEZE);a.verify(d.OUT/'evidence.json')
    cases=[];inputs=[Path(__file__).resolve(),d.FREEZE,d.OUT/'evidence.json']
    for c in a.read(d.CONTRACT)['cases']:
        tag=c['case'];device=c['device'];vd=float(c['vd'])
        old=d.REPO/'build-release/m69_stage_v2/sentaurus_bundle'/device/('vd_'+format(vd,'.6f').replace('.','p'))
        source=next(old.glob('*_des.cmd'))
        dest=LOCAL/'bundle'/tag;dest.mkdir(parents=True,exist_ok=False)
        shutil.copyfile(old/'input_fps.tdr',dest/'input_fps.tdr')
        original=source.read_text()
        physics=original[original.index('Physics {'):original.index('Solve {')]
        header='''File { Grid="input_fps.tdr" Plot="native_des.tdr" Current="native" Output="native.log" }
Electrode {
 { Name="source" Voltage=0 }
 { Name="drain" Voltage=0 }
 { Name="gate" Voltage=0 }
 { Name="substrate" Voltage=0 }
}
'''
        solve=f'''Solve {{
 Coupled(Iterations=100) {{ Poisson }}
 Coupled {{ Poisson Electron Hole }}
 Quasistationary(InitialStep=0.1 Increment=1.5 MinStep=1e-5 MaxStep=1 Goal {{ Name="drain" Voltage={vd} }}) {{ Coupled {{ Poisson Electron Hole }} }}
 NewCurrentPrefix="IdVg_"
 Coupled {{ Poisson Electron Hole }}
 Plot(FilePrefix="vg_000")
 Save(FilePrefix="vg_000")
'''
        for index in range(1,21):
            vg=index*.05
            solve+=f''' Quasistationary(InitialStep=0.5 Increment=1.5 MinStep=0.0001 MaxStep=1 Goal {{ Name="gate" Voltage={vg:.12g} }}) {{
  Coupled {{ Poisson Electron Hole }}
  CurrentPlot(Time=(1))
 }}
 Plot(FilePrefix="vg_{index:03d}")
 Save(FilePrefix="vg_{index:03d}")
'''
        solve+='}\n'
        (dest/'native_des.cmd').write_text(header+physics+solve,encoding='utf-8',newline='\n')
        inputs += [source,old/'input_fps.tdr',dest/'input_fps.tdr',dest/'native_des.cmd']
        cases.append(dict(case=tag,device=device,vd=vd,vg=[round(.05*i,12) for i in range(21)]))
    shell='''#!/bin/bash
set -u
cd "$(dirname "$0")"
for path in bundle/*; do
 (
  cd "$path"
  sdevice native_des.cmd > console.log 2>&1
  code=$?
  printf '%s\\n' "$code" > exit_code.txt
  printf '%s %s\\n' "$path" "$code"
 )
done
tar -czf results.tgz bundle
printf 'complete\\n' > complete.txt
'''
    (LOCAL/'run.sh').write_text(shell,encoding='utf-8',newline='\n');inputs.append(LOCAL/'run.sh')
    a.write(OUT/'native_contract.json',dict(status='frozen_before_execution',cases=cases,remote_root=REMOTE,
        purpose='Matched ni/no BGN 84 native target fields with ascending gate path; preserve source mesh and physics.',
        math='Inherited qualified native Digits=8, ErrRef electron/hole=100, Iterations=20, ExitOnFailure; no tolerance relaxation.',
        acceptance={'exit_zero':True,'runtime':'T-2022.03-SP2','bias_error_V':1e-10,'target_states':84,'preserve_node_support':True},
        field_qualification='Native convergence/voltage/mesh identity separately recorded. Vela strict carrier/global/KCL acceptance is independent; no failed point silently promoted.',
        new_native_sweeps=4,new_target_fields=84,production_changes=False,m82_released=False,m83_released=False))
    inputs.append(OUT/'native_contract.json');d.matrix.freeze(OUT/'native_freeze.json',inputs)
    with tarfile.open(LOCAL/'input.tgz','w:gz') as tar:
        tar.add(LOCAL/'bundle',arcname='bundle');tar.add(LOCAL/'run.sh',arcname='run.sh')
    print('Prepared and frozen 4 native sweeps / 84 target fields:',LOCAL/'input.tgz')


if __name__=='__main__':main()
