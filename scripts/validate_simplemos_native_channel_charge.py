"""Frozen native shallow-channel fixed-charge, two-amplitude signed FD pilot."""
import argparse
import json
import math
import os
import shutil
import subprocess
import tarfile
from pathlib import Path
import numpy as np
import compare_simplemos_qualified_charge_response as upstream
from run_simplemos_m81_output_amendment import rows as pltrows
from analyze_simplemos_m81_fixed_charge_fd import compare

old=upstream.old
REPO,ROOT=upstream.REPO,upstream.ROOT
LOCAL=REPO/'build-release/simplemos_native_channel_charge_fd'
OUT=ROOT/'native_channel_charge_fd'
CONTRACT=OUT/'contract.json'
EVIDENCE=OUT/'evidence.json'
REMOTE='/tmp/vela_simplemos_channel_fd_20260905'
SCRIPT=Path(__file__).resolve()
CASES=(('zero',0.),('plus_1e13',1e13),('minus_1e13',-1e13),('plus_5e12',5e12),('minus_5e12',-5e12))


def prepare():
    if CONTRACT.exists():raise FileExistsError(CONTRACT)
    upstream.verify()
    geo=upstream.m73.Geometry('n23');m,_,xy=old.m78.supports('n23',geo,.05)
    surface=float(xy[m['interface'],0].min());center=surface+.025
    target=m['channel'].copy();target[geo.contact_nodes]=False
    box=m['all_si'] & (abs(xy[:,0]-center)<.025000000001) & (abs(xy[:,1])<.125000000001)
    if not np.array_equal(box,target):raise ValueError('Window does not reproduce frozen node mask')
    mask=[{'node_id':i,'x_um':xy[i,0],'y_um':xy[i,1],'selected':bool(target[i])} for i in range(geo.count)]
    old.write_csv(OUT/'node_mask.csv',mask)
    src=upstream.native.LOCAL/'fixed_charge_bundle/n23/zero'
    for name,value in CASES:
        dest=LOCAL/'bundle'/name
        shutil.copytree(src,dest)
        path=dest/'m81_des.cmd';text=path.read_text()
        assert text.count('Traps(FixedCharge Conc=0)')==1
        text=text.replace('Traps(FixedCharge Conc=0)',
            f'Traps(FixedCharge Conc={value:.16g} SpatialShape=Uniform SpaceMid=({center:.17g} 0 0) SpaceSig=(0.025000000001 0.125000000001 1))')
        path.write_text(text,encoding='utf-8',newline='\n')
    shell='''#!/bin/bash
set -u
cd "$(dirname "$0")"
for label in zero plus_1e13 minus_1e13 plus_5e12 minus_5e12; do
  (
    cd "bundle/$label"
    sdevice m81_des.cmd > m81.console.log 2>&1
    result=$?
    printf '%s\\n' "$result" > exit_code.txt
    printf '%s %s\\n' "$label" "$result"
  )
done
tar -czf results.tgz bundle
'''
    (LOCAL/'run.sh').write_text(shell,encoding='utf-8',newline='\n')
    prediction=next(r for r in old.read_csv(upstream.OUT/'response_ledger.csv') if r['device']=='n23' and float(r['vd'])==.05 and r['support']=='channel_0p05um')
    pilot=next(r for r in old.read_json(upstream.native.OUT/'m81_native_pilot_result.json')['cases'] if r['device']=='n23')
    files=[SCRIPT,OUT/'node_mask.csv',upstream.CONTRACT,upstream.OUT/'response_ledger.csv',
           ROOT/'simplemos_strict_precision_closure_evidence.json',upstream.native.OUT/'m81_native_pilot_result.json',LOCAL/'run.sh',
           REPO/'build-release/sentaurus_import.exe',REPO/'build-release/m79_research/sdevice_ug_local_2022.txt']
    files += [p for p in (LOCAL/'bundle').rglob('*') if p.is_file()]
    old.write_json(CONTRACT,{'status':'frozen_before_execution','device':'n23','vd':.05,'vg':.9,
        'cases':[{'name':n,'charge_cm_3':v} for n,v in CASES], 'remote_root':REMOTE,
        'window':{'center_um':[center,0,0],'half_width_um':[.025000000001,.125000000001,1],'selected_nodes':int(target.sum()),
                  'scope':'Original shallow channel mask; x depth, y source-drain; no Dirichlet Poisson node included. 1e-12 um boundary margin changes no node membership.'},
        'manual_basis':'Local Sentaurus Device UG T-2022.03 p548 Eq523: Uniform trap spatial window, SpaceMid/SpaceSig in um; p543 signed FixedCharge Conc. No SingleTrap or SFactor, no mesh/dopant changes.',
        'native_green_prediction_A_per_um_at_1e13':float(prediction['sentaurus_green_delta_Id_A_per_um']),
        'vela_qualified_prediction_A_per_um_at_1e13':float(prediction['vela_delta_Id_A_per_um']),
        'pilot_current_A_per_um':pilot['dc_current_A_per_um'],
        'gates':{'all_exit_zero':True,'bias_error_V':1e-12,'zero_reclosure_dex':1e-5,
                 'fd_vs_green_relative':.05,'two_amplitude_derivative_relative':.05,
                 'even_nonlinear_fraction':.01,'minimum_signal_to_zero_drift':100,
                 'charge_mask_relative_tolerance':1e-6,'roundoff_multiplier':128},
        'mask_audit':'At every Si node compare Delta SpaceCharge - Delta p + Delta n - Delta ND + Delta NA to signed Conc * frozen mask, using zero control. SpaceCharge is exported number density. Tolerance max(|Conc|*1e-6,128*double_epsilon*max component magnitude).',
        'counts':{'dc_zero':1,'dc_perturbed':4,'bias_sweeps':0,'new_ac':0},
        'm82_released':False,'m83_released':False,
        'input_hashes':{old.portable(p):old.sha256(p) for p in sorted(set(files))}})
    with tarfile.open(LOCAL/'input.tgz','w:gz') as tar:
        tar.add(LOCAL/'bundle',arcname='bundle');tar.add(LOCAL/'run.sh',arcname='run.sh')
    print(json.dumps({'prepared':True,'selected_nodes':int(target.sum()),'archive':str(LOCAL/'input.tgz')}),flush=True)


def verify():
    for p,digest in old.read_json(CONTRACT)['input_hashes'].items():
        if old.sha256(REPO/p)!=digest:raise ValueError(f'Frozen input changed {p}')


def unpack():
    verify()
    with tarfile.open(LOCAL/'results.tgz','r:gz') as tar:
        root=(LOCAL/'raw').resolve()
        for member in tar.getmembers():
            path=(root/member.name).resolve()
            if not path.is_relative_to(root) or member.issym() or member.islnk():raise ValueError('Unsafe archive member')
        tar.extractall(root,filter='data')


def export():
    verify()
    env=os.environ.copy();env['PATH']='D:/msys64/ucrt64/bin'+os.pathsep+env['PATH']
    for name,_ in CASES:
        src=LOCAL/'raw/bundle'/name;dest=LOCAL/'exports'/name;dest.mkdir(parents=True,exist_ok=True)
        if (dest/'field_manifest.json').exists():continue
        p=subprocess.run([str(REPO/'build-release/sentaurus_import.exe'),'--tdr',str(src/'reclosed_des.tdr'),
            '--export-dir',str(dest),'--inventory-json',str(dest/'inventory.json')],env=env,capture_output=True,text=True)
        (dest/'export.log').write_text(p.stdout+p.stderr)
        if p.returncode:raise RuntimeError(p.stderr[-1000:])


def fields(name):
    path=LOCAL/'exports'/name
    return {k:upstream.m73.scalar(path/'fields'/f'{k}_region0.csv') for k in
            ('SpaceCharge','eDensity','hDensity','DonorConcentration','AcceptorConcentration')}


def analyze():
    verify();cfg=old.read_json(CONTRACT);g=cfg['gates'];currents={};case_rows=[];mask_rows=[]
    mask={int(r['node_id']):r['selected']=='True' for r in old.read_csv(OUT/'node_mask.csv')}
    base=fields('zero')
    for name,charge in CASES:
        src=LOCAL/'bundle'/name;raw=LOCAL/'raw/bundle'/name
        for p in src.iterdir():
            if old.sha256(p)!=old.sha256(raw/p.name):raise ValueError('Remote input changed')
        code=int((raw/'exit_code.txt').read_text())
        if code:raise ValueError(f'Sentaurus failed {name}: {code}')
        points=pltrows(raw/'m81_des.plt')
        if len(points)!=1:raise ValueError('Expected one DC point')
        r=points[0];current=r['drain TotalCurrent'];currents[name]=current
        log=(raw/'m81.console.log').read_text()
        if 'T-2022.03-SP2' not in log or 'Good Bye' not in log:raise ValueError('Runtime identity')
        assert abs(r['drain OuterVoltage']-.05)<=1e-12 and abs(r['gate OuterVoltage']-.9)<=1e-12
        case_rows.append({'name':name,'charge_cm_3':charge,'exit_code':code,'current_A_per_um':current})
        f=fields(name)
        for k in f:
            if set(f[k])!=set(base[k]):raise ValueError('Node support changed')
        for i in sorted(base['SpaceCharge']):
            delta={k:f[k][i]-base[k][i] for k in f}
            inferred=delta['SpaceCharge']-delta['hDensity']+delta['eDensity']-delta['DonorConcentration']+delta['AcceptorConcentration']
            expected=charge*mask[i]
            magnitude=max(abs(d[i]) for ff in (f,base) for d in ff.values())
            tolerance=max(abs(charge)*g['charge_mask_relative_tolerance'],g['roundoff_multiplier']*np.finfo(float).eps*magnitude,1.)
            mask_rows.append({'name':name,'node_id':i,'selected':mask[i],'inferred_charge_cm_3':inferred,
                'expected_charge_cm_3':expected,'error_cm_3':inferred-expected,'tolerance_cm_3':tolerance,'pass':abs(inferred-expected)<=tolerance})
    cal=[]
    for suffix,amp in (('1e13',1e13),('5e12',5e12)):
        r=compare(currents['plus_'+suffix],currents['minus_'+suffix],currents['zero'],
            cfg['native_green_prediction_A_per_um_at_1e13'],amp,cfg['pilot_current_A_per_um'])
        r['zero_reclosure_dex']=abs(math.log10(currents['zero']/cfg['pilot_current_A_per_um']))
        cal.append(r)
    change=abs((cal[1]['central_delta_A_per_um']/5e12)/(cal[0]['central_delta_A_per_um']/1e13)-1)
    support_pass=all(r['pass'] for r in mask_rows)
    for r in cal:
        r['two_amplitude_derivative_relative_change']=change
        r['pass']=bool(support_pass and r['relative_error']<=g['fd_vs_green_relative'] and r['same_sign'] and
            r['positive_perturbation_sign_pass'] and r['negative_perturbation_sign_pass'] and
            r['zero_reclosure_dex']<=g['zero_reclosure_dex'] and r['signal_to_zero_control_drift']>=100 and
            r['even_nonlinear_fraction']<=g['even_nonlinear_fraction'] and change<=g['two_amplitude_derivative_relative'])
    old.write_csv(OUT/'case_ledger.csv',case_rows);old.write_csv(OUT/'calibration.csv',cal)
    old.write_csv(OUT/'charge_support_audit.csv',mask_rows)
    result={'status':'passed_local_native_charge_fd' if all(r['pass'] for r in cal) else 'failed_local_native_charge_fd',
        'passed':all(r['pass'] for r in cal),'passed_amplitudes':sum(r['pass'] for r in cal),'native_dc_runs':5,
        'charge_support_pass':support_pass,'charge_support_failed_rows':sum(not r['pass'] for r in mask_rows),
        'maximum_fd_vs_green_relative':max(r['relative_error'] for r in cal),'two_amplitude_derivative_relative_change':change,
        'scope':'n23 Vd=.05 V Vg=.9 V shallow channel bulk fixed charge only; no inference to interface sheet charge, continuity sources, n19 or high Vd.',
        'm82_released':False,'m83_released':False}
    old.write_json(OUT/'result.json',result);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','verify','unpack','export','analyze'))
    {'prepare':prepare,'verify':verify,'unpack':unpack,'export':export,'analyze':analyze}[p.parse_args().action]()
