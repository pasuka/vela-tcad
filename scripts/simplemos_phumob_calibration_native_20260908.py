"""PhuMob cell/edge calibration on the qualified BGN geometry and biases.

Remote execution is explicit; no production model guard or acceptance changes.
"""
import argparse
import math
import shutil
import tarfile
from pathlib import Path
import simplemos_bgn_restore_native_20260908 as b
import prepare_simplemos_masetti_runtime_20260907 as runtime
import validate_simplemos_bgn_psi_scoped_20260908 as qualified

a,d,REPO=b.a,b.d,b.REPO
LOCAL=REPO/'build-release/simplemos_phumob_calibration_20260908'
OUT=REPO/'reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908'
REMOTE='/tmp/vela_simplemos_phumob_calibration_20260908'
MODELS=('masetti_control','phumob')
INDICES=(40,50)


def prepare():
    a.verify(qualified.OUT/'final_evidence.json')
    a.verify(qualified.OUT/'validation_freeze.json')
    a.verify(b.OUT/'native_evidence.json')
    files=[Path(__file__).resolve(),Path(runtime.__file__),runtime.prior.MANUAL,
           qualified.OUT/'final_evidence.json',b.OUT/'native_evidence.json']
    jobs=[]
    for c in a.read(b.OUT/'native_contract.json')['cases']:
        source=b.LOCAL/'native_raw/bundle/old_slotboom'/c['case']
        header=(source/'native_des.cmd').read_text().split('Plot {',1)[0]
        assert 'Mobility(DopingDependence)' in header
        for model in MODELS:
            for index in INDICES:
                key=f"{c['case']}_vg_{index:03d}"
                dest=LOCAL/'bundle'/model/key
                dest.mkdir(parents=True,exist_ok=False)
                for name in ('input_fps.tdr',f'result_{index:03d}_des.sav',f'result_{index:03d}_circuit_des.sav'):
                    shutil.copyfile(source/name,dest/name);files += [source/name,dest/name]
                text=header if model=='masetti_control' else header.replace('Mobility(DopingDependence)','Mobility(PhuMob)')
                text+='''Plot { eDensity hDensity Potential eQuasiFermi hQuasiFermi
 eMobility/Element hMobility/Element eCurrent/Vector/Element hCurrent/Vector/Element
 Doping DonorConcentration AcceptorConcentration SRHRecombination
 BandGapNarrowing EffectiveIntrinsicDensity }
CurrentPlot { Tcl(tcl="source runtime.tcl") }
'''+b.prior.prior.MATH+f'''
Solve {{ Load(FilePrefix="result_{index:03d}") Coupled {{ Poisson Electron Hole }}
 Plot(FilePrefix="final") Save(FilePrefix="final") }}
'''
                (dest/'native_des.cmd').write_text(text,newline='\n')
                (dest/'runtime.tcl').write_text(runtime.TCL,newline='\n')
                files += [dest/'native_des.cmd',dest/'runtime.tcl',source/'native_des.cmd']
                jobs.append(dict(**c,key=key,model=model,index=index,vg=index*.02))
    shell='''#!/bin/bash
set -u
cd "$(dirname "$0")"
count=0
for path in bundle/*/*; do
 (
  cd "$path"
  /atctools/Synopsys/tcad/T-2022.03/bin/sdevice native_des.cmd > console.log 2>&1
  printf '%s\\n' "$?" > exit_code.txt
 ) &
 count=$((count+1))
 if [ "$((count % 4))" -eq 0 ]; then wait; fi
done
wait
tar -czf results.tgz bundle
printf 'complete\\n' > complete.txt
'''
    (LOCAL/'run.sh').write_text(shell,newline='\n');files.append(LOCAL/'run.sh')
    a.write(OUT/'native_contract.json',dict(jobs=jobs,remote=REMOTE,target_states=16,
        scope='Eight PhuMob states and eight same-bias Masetti runtime-export controls, using qualified OldSlotboom saved seeds. Enormal/HFS remain off.',
        physics_axis='Mobility(DopingDependence) -> Mobility(PhuMob), default arsenic; total impurities and both live carrier populations retained.',
        gates=dict(native_exit_zero=True,runtime='T-2022.03-SP2',bias_error_V=1e-10,kcl_over_Id=1e-8,
                   masetti_current_reclosure_relative=1e-8,cell_mobility_relative=1e-7,edge_mobility_relative=1e-7,
                   native_SG_terminal_relative=1e-6,isolated_derivative_relative=1e-5),
        constitutive_candidates=['box-weighted vertex mobility at local n,p','box-weighted vertex mobility at cell arithmetic n,p','box-weighted vertex mobility at cell geometric n,p','mobility at box-averaged donors,acceptors,n,p'],
        source_and_state_gates='Unchanged inherited BGN contracts. Global SRH source-relative closure was inactive below its existing source floor; do not treat that as a passed small-source response calibration.',
        interpretation='Runtime edge/node mobility storage is not declared the solver operator before cell/edge/port checks. No fitting to drain current.',
        production_changed=False,acceptance_changed=False))
    d.matrix.freeze(OUT/'native_freeze.json',files+[OUT/'native_contract.json'])
    with tarfile.open(LOCAL/'input.tgz','w:gz') as tar:
        tar.add(LOCAL/'bundle',arcname='bundle');tar.add(LOCAL/'run.sh',arcname='run.sh')
    print('Prepared 16 independent native DC/runtime exports',flush=True)


def unpack():
    a.verify(OUT/'native_freeze.json')
    dest=LOCAL/'native_raw';assert not dest.exists()
    with tarfile.open(LOCAL/'results.tgz') as tar:
        for item in tar.getmembers():
            assert (dest/item.name).resolve().is_relative_to(dest.resolve()) and not item.issym() and not item.islnk()
        tar.extractall(dest,filter='data')
    for source in (LOCAL/'bundle').rglob('*'):
        if source.is_file():assert a.sha(source)==a.sha(dest/'bundle'/source.relative_to(LOCAL/'bundle'))
    rows=[]
    for job in a.read(OUT/'native_contract.json')['jobs']:
        path=dest/'bundle'/job['model']/job['key']
        code=int((path/'exit_code.txt').read_text())
        log=(path/'console.log').read_text(errors='replace')
        row=dict(**job,exit_code=code,native_qualified=False)
        plt=path/'native_des.plt'
        if plt.exists():
            points=b.exporter.pltrows(plt);assert len(points)==1,(path,len(points))
            p=points[0];currents=[p[n+' TotalCurrent'] for n in ('drain','source','gate','substrate')]
            kcl=abs(math.fsum(currents))/max(abs(currents[0]),1e-300)
            bias=max(abs(p['gate OuterVoltage']-job['vg']),abs(p['drain OuterVoltage']-job['vd']))
            snapshots=sorted(path.glob('runtime_*_vertices.csv'))
            row.update(Id_A_per_um=currents[0],kcl_over_Id=kcl,bias_error_V=bias,
                       runtime_prefix=snapshots[-1].name.removesuffix('_vertices.csv') if snapshots else '',
                       native_qualified=code==0 and 'T-2022.03-SP2' in log and 'Good Bye' in log and kcl<=1e-8 and bias<=1e-10 and bool(snapshots))
        rows.append(row)
    a.write_csv(OUT/'native_points.csv',rows)
    d.matrix.freeze(OUT/'native_evidence.json',[OUT/'native_freeze.json',OUT/'native_points.csv',LOCAL/'results.tgz']+[p for p in dest.rglob('*') if p.is_file()])
    print('Native qualified',sum(r['native_qualified'] for r in rows),'/',len(rows),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','unpack'))
    globals()[parser.parse_args().action]()
