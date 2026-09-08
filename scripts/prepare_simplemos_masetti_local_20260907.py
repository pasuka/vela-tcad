"""Freeze native element/box exports and eight qualified-state local operator audits."""
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import math
import shutil
import tarfile
import run_simplemos_reference_subset_20260906 as prev

a=prev.a;d=prev.d;REPO=d.REPO
LOCAL=REPO/'build-release/simplemos_masetti_local_20260907'
OUT=d.ROOT/'masetti_local_discretization_20260907'
REMOTE='/tmp/vela_simplemos_masetti_local_20260907'
MANUAL=REPO/'build-release/m79_research/sdevice_ug_local_2022.txt'


def native_prepare():
    a.verify(prev.OUT/'validation_evidence.json')
    files=[Path(__file__).resolve(),prev.OUT/'validation_evidence.json',MANUAL];jobs=[]
    for device in ('n19','n23'):
        case=f'm65_{device}_vd_0p050000_endpoint'
        source=prev.LOCAL/'native_raw/bundle/masetti'/case;dest=LOCAL/'bundle'/device;dest.mkdir(parents=True,exist_ok=False)
        for name in ('input_fps.tdr','result_020_des.sav','result_020_circuit_des.sav'):
            shutil.copyfile(source/name,dest/name);files += [source/name,dest/name]
        text=(source/'native_des.cmd').read_text().split('Plot {',1)[0]
        text+='''Plot {
 eDensity hDensity Potential eQuasiFermi hQuasiFermi
 eMobility/Element hMobility/Element
 DonorConcentration AcceptorConcentration Doping
 eCurrent/Vector/Element hCurrent/Vector/Element
 BM_ElementVolume EffectiveIntrinsicDensity
}
'''
        text+=prev.p.MATH.replace('Math {','Math { BoxMeasureFromFile(GrdNumbering)')
        text+='''
Solve {
 Load(FilePrefix="result_020")
 Coupled { Poisson Electron Hole }
 Plot(FilePrefix="local_export")
}
'''
        (dest/'native_des.cmd').write_text(text,newline='\n');files.append(dest/'native_des.cmd');jobs.append(dict(device=device,case=case,vd=.05,vg=1.))
    shell='''#!/bin/bash
set -u
cd "$(dirname "$0")"
for path in bundle/*; do
 (cd "$path"; sdevice native_des.cmd > console.log 2>&1; printf '%s\\n' "$?" > exit_code.txt)
done
tar -czf results.tgz bundle
printf 'complete\\n' > complete.txt
'''
    (LOCAL/'run.sh').write_text(shell,newline='\n');files.append(LOCAL/'run.sh')
    a.write(OUT/'native_contract.json',dict(status='frozen_before_execution',jobs=jobs,remote=REMOTE,
        question='Map native element mobility and box coefficients/measures on the unchanged Masetti model, then compare with local formula and Vela operators.',
        model_changes=False,math_change='BoxMeasureFromFile(GrdNumbering) with no external input file in the fresh directory, to write the native MeasureCoefficients.debug as in M34. Same 128-bit reclosure settings.',
        gates=dict(current_reclosure_relative=1e-8,bias_error_V=1e-10,kcl_over_Id=1e-8,coordinate_error_um=1e-12),
        scope='Two meshes at Vd=.05,Vg=1. Mobility is bias-independent in frozen Masetti; no full-current or arbitrary native edge-operator equivalence assumed.',
        manual_reference='T-2022.03 User Guide, Mobility Averaging p464 and box-method coefficient/measure export. Native node plots are not edge mobility.',
        raw_inputs_immutable=True,production_changes=False))
    files.append(OUT/'native_contract.json');d.matrix.freeze(OUT/'native_freeze.json',files)
    with tarfile.open(LOCAL/'input.tgz','w:gz') as t:t.add(LOCAL/'bundle',arcname='bundle');t.add(LOCAL/'run.sh',arcname='run.sh')
    print('Frozen 2 native element/box exports',flush=True)


def native_unpack():
    a.verify(OUT/'native_freeze.json');dest=(LOCAL/'native_raw').resolve();assert not dest.exists()
    with tarfile.open(LOCAL/'results.tgz') as t:
        for m in t.getmembers():assert (dest/m.name).resolve().is_relative_to(dest) and not m.issym() and not m.islnk()
        t.extractall(dest,filter='data')
    for p in (LOCAL/'bundle').rglob('*'):
        if p.is_file():assert a.sha(p)==a.sha(dest/'bundle'/p.relative_to(LOCAL/'bundle'))
    checks=[]
    for j in a.read(OUT/'native_contract.json')['jobs']:
        root=dest/'bundle'/j['device'];log=(root/'console.log').read_text(errors='replace')
        assert int((root/'exit_code.txt').read_text())==0 and 'Good Bye' in log and 'T-2022.03-SP2' in log
        points=prev.p.exporter.pltrows(root/'native_des.plt');assert len(points)==1
        row=points[0];assert abs(row['gate OuterVoltage']-j['vg'])<=1e-10 and abs(row['drain OuterVoltage']-j['vd'])<=1e-10
        cc=[row[c+' TotalCurrent'] for c in ('drain','source','gate','substrate')];kcl=abs(math.fsum(cc))/abs(cc[0])
        old=next(r for r in a.rows(prev.OUT/'native_points.csv') if r['model']=='masetti' and r['case']==j['case'] and int(r['index'])==20)
        drift=abs(cc[0]/float(old['Id_A_per_um'])-1)
        assert drift<=1e-8 and kcl<=1e-8 and (root/'MeasureCoefficients.debug').exists()
        checks.append(dict(**j,current_A_per_um=cc[0],relative_drift=drift,kcl_over_Id=kcl))
        prev.p.exporter.export_one(dict(case=j['case'],index=20,tdr=str(root/'local_export_des.tdr'),export=str(LOCAL/'native_exports'/j['device'])))
    a.write_csv(OUT/'native_checks.csv',checks);print('Verified/exported 2 native element/box states',flush=True)


def vela_prepare():
    a.verify(prev.OUT/'validation_evidence.json');cases=[];files=[Path(__file__).resolve(),prev.OUT/'validation_evidence.json',prev.RUNNER]
    for j in a.read(prev.OUT/'vela_contract.json')['jobs']:
        if j['arm']!='masetti_refined':continue
        assert a.read(Path(j['config']).parent/'result.json')['comparison_qualified']
        key=f"{j['case']}_vg{j['index']:03d}";base=Path(j['config']).parent;cfg=a.read(base/'config.json');cfg.pop('output_state_file');cfg['solver'].pop('local_update_diagnostics')
        case=dict(j,key=key,base=str(base),mapped=str(base/'initial.csv'),strict=str(base/'state.csv'));paths=[]
        for role,state in (('strict',base/'state.csv'),('mapped',base/'initial.csv')):
            for name,kind in (('functional','terminal_current_functional_probe'),('carrier','newton_carrier_term_probe'),('edges','sg_edge_flux_probe')):
                root=LOCAL/'vela'/key/role;deck=copy.deepcopy(cfg);deck.update(state_file=str(state),simulation_type=kind,contact='drain')
                if name=='functional':deck['residual_output_csv']=str(root/'residual.csv')
                else:deck['output_csv']=str(root/(name+'.csv'))
                if name=='carrier':deck['carrier_term_probe']={'solved_equation_terms':True}
                path=root/(name+'.json');a.write(path,deck);paths.append(str(path));files.append(path)
        root=LOCAL/'vela'/key;deck=copy.deepcopy(cfg);deck.update(state_file=str(base/'state.csv'),simulation_type='terminal_current_adjoint_probe',contact='drain',output_csv=str(root/'adjoint.csv'))
        path=root/'adjoint.json';a.write(path,deck);paths.append(str(path));files.append(path)
        case['probes']=paths;cases.append(case)
        files += [base/'config.json',base/'state.csv',base/'initial.csv',base/'result.json']+[Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file')]
    a.write(OUT/'vela_contract.json',dict(status='frozen_before_execution',cases=cases,read_only_probes=56,
        definitions=dict(mapped='Native potentials with coherent Vela densities as frozen in the previous run.',
            SG='Same-state Vela operator and common-geometry state secant; distinguish from actual native edge conductance.',
            poisson='Exact conditional decomposition of state difference through unchanged Vela K and physical charge volumes, with residual and boundary terms retained.',
            adjoint='Fresh Masetti full-DD weights; arbitrary finite gaps remain screening until direction/nonlinear response calibration.'),
        gates=dict(current_identity_relative=1e-8,edge_sum_relative=1e-8,adjoint_relative_residual=1e-10,poisson_reconstruction_relative=1e-5),
        production_changes=False,nonlinear_acceptance_changes=False))
    files.append(OUT/'vela_contract.json');d.matrix.freeze(OUT/'vela_freeze.json',files);print('Frozen 56 local operator probes',flush=True)


def vela_run():
    a.verify(OUT/'vela_freeze.json')
    def one(c):
        for path in c['probes']:
            s=prev.execute(Path(path),'masetti_plain');assert s['exit_code']==0,(path,s)
        print(c['key'],'7 probes complete',flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,a.read(OUT/'vela_contract.json')['cases']))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('native_prepare','native_unpack','vela_prepare','vela_run'));globals()[p.parse_args().action]()
