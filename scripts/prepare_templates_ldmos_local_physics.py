"""Prepare isolated native Load/Plot diagnostics without Newton iterations.

Preserve the original D0 physics and finite contacts; add Auger to Plot only.
This prepares files locally and does not upload or execute them.
"""
import argparse,hashlib,json,shutil
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('native-full','control-bundle','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--single-process',action='store_true',help='Load/Plot every state in one native process')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    source=a.control_bundle/'vg8_closed.cmd';deck=source.read_text(encoding='utf-8')
    if deck.count('Solve {')!=1 or deck.count('SRH Band2Band')!=1:
        raise ValueError('Unexpected audited native deck')
    prefix=deck.split('Solve {')[0].replace('SRH Band2Band','SRH Auger Band2Band')
    manifest=dict(scope=__doc__,single_process=a.single_process,sources_sha256={},cases=[])
    def record(path):manifest['sources_sha256'][str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
    record(source)
    for name in ('n1_fps.tdr','sdevice.par','Siliconc100.par'):
        original=a.native_full/'bundle'/name;shutil.copy2(original,a.output/name);record(original)
    commands=['#!/bin/bash','set -eu']
    for gate in (4,8):
        for index in (0,1,10,30):
            name=f'local_vg{gate}_{index:02d}';seed=f'seed_vg{gate}_{index:02d}'
            original=a.native_full/'raw'/f'field_vg{gate}_{index:04d}_des.tdr'
            shutil.copy2(original,a.output/f'{seed}_des.tdr');record(original)
            text=prefix.replace('vg8_closed',name)+f'Solve {{\n Load(FilePrefix="{seed}")\n Plot(FilePrefix="{name}")\n}}\n'
            (a.output/f'{name}.cmd').write_text(text,encoding='utf-8',newline='\n')
            commands.append(f'/atctools/Synopsys/tcad/T-2022.03/bin/sdevice {name}.cmd > {name}.stdout 2>&1')
            manifest['cases'].append(dict(name=name,gate_V=gate,index=index,bias_V=index*40/30,seed=f'{seed}_des.tdr'))
    if a.single_process:
        solve='Solve {\n'+''.join(f' Load(FilePrefix="{c["seed"].removesuffix("_des.tdr")}")\n Plot(FilePrefix="{c["name"]}")\n' for c in manifest['cases'])+'}\n'
        (a.output/'local_batch.cmd').write_text(prefix.replace('vg8_closed','local_batch')+solve,encoding='utf-8',newline='\n')
        commands=commands[:2]+['/atctools/Synopsys/tcad/T-2022.03/bin/sdevice local_batch.cmd > local_batch.stdout 2>&1']
    commands+=['tar -czf results.tgz local_*.log local_*.stdout local_*_des.tdr','printf "completed\\n" > completed.txt']
    (a.output/'run.sh').write_text('\n'.join(commands)+'\n',encoding='utf-8',newline='\n')
    manifest['bundle_sha256']={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in a.output.iterdir() if path.is_file()}
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')

if __name__=='__main__':main()
