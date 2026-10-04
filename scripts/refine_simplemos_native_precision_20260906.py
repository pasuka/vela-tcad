"""Freeze and export 42 high-NWell native DC reclosures with tighter arithmetic."""
import argparse
from pathlib import Path
import shutil
import tarfile
from concurrent.futures import ThreadPoolExecutor
import prepare_simplemos_fullfield_native_20260906 as base
import export_simplemos_fullfield_native_20260906 as exporter

a=base.a;d=base.d
LOCAL=base.LOCAL/'native_precision';OUT=base.OUT/'native_precision'
REMOTE=base.REMOTE+'/precision'


def prepare():
    a.verify(base.OUT/'native_freeze.json');a.verify(base.OUT/'export_contract.json')
    original='Math { Extrapolate RelErrControl Digits=8 ErrRef(Electron)=1e2 ErrRef(Hole)=1e2 Iterations=20 ExitOnFailure CNormPrint }'
    improved='Math { ExtendedPrecision(128) Method=Super RelErrControl Digits=12 ErrRef(Electron)=1e-2 ErrRef(Hole)=1e-2 RhsMin=1e-20 Iterations=40 ExitOnFailure CNormPrint }'
    files=[Path(__file__).resolve(),base.OUT/'native_points.csv',base.OUT/'native_freeze.json',base.OUT/'export_contract.json']
    cases=[]
    for c in a.read(base.OUT/'native_contract.json')['cases']:
        if c['device']!='n23':continue
        old=base.LOCAL/'native_raw/bundle'/c['case'];dest=LOCAL/'bundle'/c['case'];dest.mkdir(parents=True,exist_ok=False)
        for source in [old/'input_fps.tdr']+sorted(old.glob('vg_*.sav')):
            shutil.copyfile(source,dest/source.name);files += [source,dest/source.name]
        text=(old/'native_des.cmd').read_text().split('Solve {',1)[0]
        assert text.count(original)==1
        text=text.replace(original,improved)+'Solve {\n'
        for index in range(21):
            text+=f''' Load(FilePrefix="vg_{index:03d}")
 NewCurrentPrefix="check_{index:03d}_"
 Coupled {{ Poisson Electron Hole }}
 Plot(FilePrefix="refined_{index:03d}")
 Save(FilePrefix="refined_{index:03d}")
'''
        text+='}\n';(dest/'native_des.cmd').write_text(text,newline='\n');files.append(dest/'native_des.cmd');cases.append(c)
    shell='''#!/bin/bash
set -u
cd "$(dirname "$0")"
for path in bundle/*; do
 (
  cd "$path"
  sdevice native_des.cmd > console.log 2>&1
  code=$?
  printf '%s\\n' "$code" > exit_code.txt
 )
done
tar -czf results.tgz bundle
printf 'complete\\n' > complete.txt
'''
    (LOCAL/'run.sh').write_text(shell,newline='\n');files.append(LOCAL/'run.sh')
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',cases=cases,remote=REMOTE,
        trigger='Original high-NWell native KCL/Id up to .005694 despite successful native Digits=8 convergence.',
        numerical_control=improved,physical_models_mesh_and_bias_unchanged=True,
        acceptance={'exit_zero':True,'bias_error_V':1e-10,'kcl_over_Id':1e-8},
        control='Independent reclosure from each saved native DC state; original 42 states and their KCL failures retained.',
        manual='Local T-2022.03 user guide pp.223-225: ExtendedPrecision(128), SUPER support, tighter Digits/RhsMin; runtime SP2 tested separately.',
        independent_target_states=42,production_changes=False))
    files.append(OUT/'contract.json');d.matrix.freeze(OUT/'freeze.json',files)
    with tarfile.open(LOCAL/'input.tgz','w:gz') as tar:
        tar.add(LOCAL/'bundle',arcname='bundle');tar.add(LOCAL/'run.sh',arcname='run.sh')
    print(LOCAL/'input.tgz',flush=True)


def unpack():
    a.verify(OUT/'freeze.json');root=(LOCAL/'native_raw').resolve()
    with tarfile.open(LOCAL/'results.tgz','r:gz') as tar:
        for member in tar.getmembers():
            target=(root/member.name).resolve()
            if not target.is_relative_to(root) or member.issym() or member.islnk():raise ValueError('Unsafe archive')
        tar.extractall(root,filter='data')
    for p in (LOCAL/'bundle').rglob('*'):
        if p.is_file():assert a.sha(p)==a.sha(root/'bundle'/p.relative_to(LOCAL/'bundle'))
    points=[];jobs=[]
    for c in a.read(OUT/'contract.json')['cases']:
        raw=root/'bundle'/c['case'];assert int((raw/'exit_code.txt').read_text())==0
        log=(raw/'console.log').read_text(errors='replace');assert 'Good Bye' in log and 'T-2022.03-SP2' in log
        for index,vg in enumerate(c['vg']):
            curve=exporter.pltrows(raw/f'check_{index:03d}_native_des.plt');assert len(curve)==1,(c['case'],index,len(curve))
            r=curve[0];assert abs(r['gate OuterVoltage']-vg)<1e-10 and abs(r['drain OuterVoltage']-c['vd'])<1e-10
            cc=[r[name+' TotalCurrent'] for name in ('drain','source','gate','substrate')]
            kcl=abs(sum(cc))/max(abs(cc[0]),1e-300)
            points.append(dict(case=c['case'],device=c['device'],vd=c['vd'],vg=vg,current_A_per_um=cc[0],kcl_A_per_um=sum(cc),kcl_over_Id=kcl,native_converged=True,native_kcl_qualified=kcl<=1e-8))
            jobs.append(dict(case=c['case'],index=index,tdr=str(raw/f'refined_{index:03d}_des.tdr'),export=str(LOCAL/'native_exports'/c['case']/f'vg_{index:03d}')))
    a.write_csv(OUT/'native_points.csv',points)
    a.write(OUT/'export_contract.json',dict(jobs=jobs,input_hashes={a.rel(LOCAL/'results.tgz'):a.sha(LOCAL/'results.tgz'),a.rel(Path(__file__).resolve()):a.sha(Path(__file__).resolve())}))
    print('Native reclosure outputs verified; KCL-qualified:',sum(r['native_kcl_qualified'] for r in points),flush=True)


def export():
    a.verify(OUT/'export_contract.json')
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(exporter.export_one,a.read(OUT/'export_contract.json')['jobs']))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','unpack','export'));globals()[p.parse_args().action]()
