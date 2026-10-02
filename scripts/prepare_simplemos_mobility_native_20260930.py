"""Prepare eight Load/Plot-only decks for the existing M60 816 states."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile


def prepare(output):
    if output.exists():
        raise ValueError('Use a new output directory')
    output.mkdir(parents=True)
    source='/tmp/vela_simplemos_m60_fields_20260929/bundle'
    for n in range(17,25):
        device=f'n{n}'; folder=output/device;folder.mkdir()
        header=f'''File {{ Grid="{source}/{device}/input_fps.tdr" Plot="final.tdr" Current="field_only" Output="field_only.log" }}
Electrode {{ {{ Name="source" Voltage=0 }} {{ Name="drain" Voltage=0 }} {{ Name="gate" Voltage=0 }} {{ Name="substrate" Voltage=0 }} }}
Physics {{ EffectiveIntrinsicDensity(OldSlotboom) }}
Physics(Material="Silicon") {{ Mobility(PhuMob HighFieldSaturation Enormal) Recombination(SRH(DopingDependence)) }}
Plot {{ eDensity hDensity Potential eQuasiFermi hQuasiFermi eMobility/Element hMobility/Element Doping DonorConcentration AcceptorConcentration SRH }}
Math {{ Extrapolate RelErrControl Digits=8 ErrRef(Electron)=1e2 ErrRef(Hole)=1e2 Iterations=20 ExitOnFailure }}
Solve {{
'''
        for vd in ('0p05','1'):
            case=f'{device}_vd_{vd}'
            for i in range(51):
                header+=f' Load(FilePrefix="{source}/{device}/m60fields_{case}_state_{i:04d}")\n Plot(FilePrefix="{case}_vg_{i:03d}")\n'
        (folder/'native.cmd').write_text(header+'}\n',encoding='utf-8',newline='\n')
    run='''#!/bin/bash
set -u
cd "$(dirname "$0")" || exit 90
mkdir execution_lock || exit 91
sha256sum -c decks.sha256 > deck_verification.log || exit 92
sha256sum /tmp/vela_simplemos_m60_fields_20260929/bundle/n*/m60fields_*_state_*_des.tdr > native_input_before.sha256
run_device() {
  local device="$1"
  (cd "$device" && sdevice native.cmd > console.log 2>&1)
  local result=$?
  echo "$result" > "$device/job.exit"
  return "$result"
}
lane_a() { for device in n17 n19 n21 n23; do run_device "$device" || return $?; done; }
lane_b() { for device in n18 n20 n22 n24; do run_device "$device" || return $?; done; }
lane_a & pa=$!
lane_b & pb=$!
wait "$pa"; a=$?
wait "$pb"; b=$?
sha256sum -c native_input_before.sha256 > native_input_after.log; identity=$?
sha256sum n*/n*_vd_*_vg_*_des.tdr > output.sha256
if [ "$a" = 0 ] && [ "$b" = 0 ] && [ "$identity" = 0 ]; then echo 0 > job.exit; else echo 1 > job.exit; fi
date -Is > completed.txt
'''
    (output/'run.sh').write_text(run,encoding='utf-8',newline='\n')
    files=sorted(output.glob('n*/native.cmd'))+[output/'run.sh']
    hashes={p.relative_to(output).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (output/'decks.sha256').write_text(''.join(f'{v}  {k}\n' for k,v in hashes.items()),encoding='utf-8',newline='\n')
    (output/'manifest.json').write_text(json.dumps(dict(points=816,devices=8,new_dc_solves=0,hashes=hashes),indent=2)+'\n')
    with tarfile.open(output.with_suffix('.tgz'),'w:gz') as tar:
        for p in sorted(output.rglob('*')):
            if p.is_file():tar.add(p,arcname=p.relative_to(output).as_posix())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    prepare(parser.parse_args().output)
