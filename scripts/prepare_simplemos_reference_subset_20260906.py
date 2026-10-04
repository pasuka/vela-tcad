"""Freeze a cross-solver Masetti-only transport control using BJT M1 evidence."""
import argparse
from pathlib import Path
import shutil
import tarfile
from concurrent.futures import ThreadPoolExecutor
import math
import validate_simplemos_linear_refinement_20260906 as prior
import prepare_simplemos_fullfield_native_20260906 as n
import export_simplemos_fullfield_native_20260906 as exporter

a=prior.a;d=prior.d
LOCAL=d.REPO/'build-release/simplemos_reference_subset_20260906'
OUT=d.ROOT/'reference_subset_20260906'
REMOTE='/tmp/vela_simplemos_reference_subset_20260906'
POINTS=(16,20)
MATH='Math { ExtendedPrecision(128) Method=Super RelErrControl Digits=12 ErrRef(Electron)=1e-2 ErrRef(Hole)=1e-2 RhsMin=1e-20 Iterations=40 ExitOnFailure CNormPrint }'


def prepare():
    a.verify(prior.OUT/'validation_evidence.json');a.verify(n.OUT/'native_freeze.json')
    assert not (OUT/'native_freeze.json').exists()
    LOCAL.mkdir(parents=True,exist_ok=True)
    files=[Path(__file__).resolve(),prior.OUT/'validation_evidence.json',n.OUT/'native_freeze.json',
        d.REPO/'reference_tcad/genius_bjt_sentaurus2022/CASE_SUMMARY.md',
        d.REPO/'reference_tcad/genius_bjt_sentaurus2022/vela/configs/m1_collector_sweep.json',
        d.REPO/'include/vela/physics/MobilityModel.h',d.REPO/'src/physics/MobilityModel.cpp']
    jobs=[]
    for case in a.read(n.OUT/'native_contract.json')['cases']:
        original=n.LOCAL/'native_raw/bundle'/case['case'];header=(original/'native_des.cmd').read_text().split('Solve {',1)[0]
        oldmath='Math { Extrapolate RelErrControl Digits=8 ErrRef(Electron)=1e2 ErrRef(Hole)=1e2 Iterations=20 ExitOnFailure CNormPrint }'
        assert header.count(oldmath)==1 and header.count('Mobility(PhuMob HighFieldSaturation Enormal)')==1
        header=header.replace(oldmath,MATH)
        for model in ('baseline','masetti'):
            dest=LOCAL/'bundle'/model/case['case'];dest.mkdir(parents=True,exist_ok=True)
            saves=[p for i in POINTS for p in sorted(original.glob(f'vg_{i:03d}*.sav'))];assert len(saves)==4
            for p in [original/'input_fps.tdr']+saves:
                if (dest/p.name).exists():assert a.sha(p)==a.sha(dest/p.name)
                else:shutil.copyfile(p,dest/p.name)
                files += [p,dest/p.name]
            text=header
            if model=='masetti':text=text.replace('Mobility(PhuMob HighFieldSaturation Enormal)','Mobility(DopingDependence)')
            text+='Solve {\n'
            for i in POINTS:
                text+=f''' Load(FilePrefix="vg_{i:03d}")
 NewCurrentPrefix="check_{i:03d}_"
 Coupled {{ Poisson Electron Hole }}
 Plot(FilePrefix="result_{i:03d}")
 Save(FilePrefix="result_{i:03d}")
'''
            text+='}\n';(dest/'native_des.cmd').write_text(text,newline='\n');files.append(dest/'native_des.cmd')
            jobs.append(dict(case=case['case'],device=case['device'],vd=case['vd'],model=model,indices=list(POINTS)))
    shell='''#!/bin/bash
set -u
cd "$(dirname "$0")"
for path in bundle/*/*; do
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
    a.write(OUT/'native_contract.json',dict(status='frozen_before_execution',jobs=jobs,remote=REMOTE,
        question='Does the high-NWell current gap persist after removing field/surface mobility and using the Masetti bulk model covered by the accepted BJT M1 terminal/source chain?',
        axis='Original PhuMob+HFS+Enormal versus Masetti DopingDependence only, simultaneously in both solvers. This is a reduced-model composite control, not a single-mechanism attribution or a fully prequalified MOS model.',
        unchanged='Same original meshes, boron/dopants, contacts, bias, matched-ni/no-BGN Boltzmann and doping-dependent SRH; 300K. No material ni swap to the BJT Fermi/BGN convention.',
        native_math=MATH,native_runs=8,native_target_states=16,interpolation=False,
        gates=dict(native_exit_zero=True,bias_error_V=1e-10,native_kcl_over_Id=1e-8,
            vela_row_eps=1e-6,vela_all_active_rows=1814,vela_zero_scale_rows=0,vela_kcl_over_Id=1e-8,
            vela_global_tolerance=1e-6,vela_global_source_floor=1e-10),
        vela_arms=['baseline with existing corrected Jacobian and 4 linear refinements',
            'Masetti with experimental vector-chain correction disabled and no linear refinements',
            'Same Masetti with exactly 4 linear refinements, for convergence discrimination'],
        endpoint_selection='n19/n23 x Vd=.05/1 x Vg=.8/1, eight working points; independent saved-state reclosures, not full curves.',
        prior_acceptance_unchanged=True,production_changes=False,m82_released=False,m83_released=False))
    files.append(OUT/'native_contract.json');d.matrix.freeze(OUT/'native_freeze.json',files)
    with tarfile.open(LOCAL/'input.tgz','w:gz') as t:
        t.add(LOCAL/'bundle',arcname='bundle');t.add(LOCAL/'run.sh',arcname='run.sh')
    print('Frozen 8 native decks / 16 DC states:',LOCAL/'input.tgz',flush=True)


def unpack():
    a.verify(OUT/'native_freeze.json');root=(LOCAL/'native_raw').resolve();assert not root.exists()
    with tarfile.open(LOCAL/'results.tgz','r:gz') as t:
        for item in t.getmembers():
            p=(root/item.name).resolve()
            if not p.is_relative_to(root) or item.issym() or item.islnk():raise ValueError('Unsafe archive')
        t.extractall(root,filter='data')
    for p in (LOCAL/'bundle').rglob('*'):
        if p.is_file():assert a.sha(p)==a.sha(root/'bundle'/p.relative_to(LOCAL/'bundle'))
    points=[];exports=[]
    for job in a.read(OUT/'native_contract.json')['jobs']:
        path=root/'bundle'/job['model']/job['case'];code=int((path/'exit_code.txt').read_text());log=(path/'console.log').read_text(errors='replace')
        assert code==0 and 'T-2022.03-SP2' in log and 'Good Bye' in log,(job,code)
        for i in job['indices']:
            r=exporter.pltrows(path/f'check_{i:03d}_native_des.plt');assert len(r)==1
            r=r[0];vg=i*.05;assert abs(r['gate OuterVoltage']-vg)<=1e-10 and abs(r['drain OuterVoltage']-job['vd'])<=1e-10
            cc=[r[c+' TotalCurrent'] for c in ('drain','source','gate','substrate')];kcl=abs(math.fsum(cc))/max(abs(cc[0]),1e-300)
            points.append(dict(**job,index=i,vg=vg,Id_A_per_um=cc[0],kcl_over_Id=kcl,native_qualified=kcl<=1e-8))
            exports.append(dict(case=job['case'],model=job['model'],index=i,tdr=str(path/f'result_{i:03d}_des.tdr'),
                export=str(LOCAL/'native_exports'/job['model']/job['case']/f'vg_{i:03d}')))
    a.write_csv(OUT/'native_points.csv',points)
    a.write(OUT/'export_contract.json',dict(jobs=exports,input_hashes={a.rel(LOCAL/'results.tgz'):a.sha(LOCAL/'results.tgz'),a.rel(Path(__file__).resolve()):a.sha(Path(__file__).resolve())}))
    print('Verified native states:',len(points),'; KCL qualified:',sum(r['native_qualified'] for r in points),flush=True)


def export():
    a.verify(OUT/'export_contract.json')
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(exporter.export_one,a.read(OUT/'export_contract.json')['jobs']))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','unpack','export'));globals()[p.parse_args().action]()
