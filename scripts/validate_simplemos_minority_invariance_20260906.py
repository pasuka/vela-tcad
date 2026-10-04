"""Separately frozen all-row and initialization-invariance diagnostic."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from pathlib import Path
import copy
import json
import math
import subprocess
import numpy as np
import run_simplemos_fullfield_vela_20260906 as b
import audit_simplemos_surface_derivative_step_20260906 as surface
import refine_simplemos_native_precision_20260906 as native

a=b.a;d=b.d
LOCAL=b.LOCAL.parent/'simplemos_minority_invariance_20260906'
OUT=d.ROOT/'minority_invariance_20260906'


def env():
    e=surface.m.env();e['VELA_VALIDATE_TRANSPORT_STEP']='1e-7';return e


def prepare():
    a.verify(b.OUT/'followup_evidence.json');a.verify(surface.OUT/'freeze.json')
    files=[Path(__file__).resolve(),b.OUT/'followup_evidence.json',surface.OUT/'freeze.json',surface.RUNNER]
    points=a.rows(native.OUT/'comparison/points.csv');jobs=[]
    for point in points:
        vg=float(point['vg'])
        if vg not in (.35,1.):continue
        key=point['case'];idx=round(vg/.05);geo=d.matrix.spatial.m73.Geometry(point['device'])
        src=b.LOCAL/point['initialization']/key/f'vg_{idx:03d}'
        original=a.read(src/'config.json');state=d.ordered(src/'state.csv',geo.count)
        masks,_,_=d.matrix.spatial.old.m78.supports(point['device'],geo,.05)
        free=masks['all_si'].copy();free[geo.contact_nodes]=False
        export=(native.LOCAL if point['device']=='n23' else b.LOCAL)/'native_exports'/key/f'vg_{idx:03d}'
        psi,n,p,spread=d.fixed.upstream.coherent_state(export,geo);assert spread<1e-12
        en=d.matrix.spatial.m73.scalar(export/'fields/eQuasiFermiPotential_region0.csv');hp=d.matrix.spatial.m73.scalar(export/'fields/hQuasiFermiPotential_region0.csv')
        files += [src/'config.json',src/'state.csv']+list((export/'fields').glob('*.csv'))
        for label,offset in (('vela',0.),('hole_plus',.05),('hole_minus',-.05),('native',0.)):
            dest=LOCAL/key/f'vg_{idx:03d}'/label
            if label=='native':
                initial=[dict(node_id=i,psi=psi[i],phin=en.get(i,0.),phip=hp.get(i,0.),electrons_m3=n[i],holes_m3=p[i]) for i in range(geo.count)]
            else:
                initial=copy.deepcopy(state)
                for i,r in enumerate(initial):
                    if not free[i] or offset==0.:continue
                    r['hole_qf_increment_V']=str(Decimal(r['hole_qf_increment_V'])+Decimal(str(offset)))
                    r['phip']=str(Decimal(r['hole_qf_reference_V'])+Decimal(r['hole_qf_increment_V']))
                    r['holes_m3']=float(r['holes_m3'])*math.exp(offset/d.fixed.upstream.VT)
            a.write_csv(dest/'initial.csv',initial)
            cfg=copy.deepcopy(original);cfg.update(state_file=str(dest/'initial.csv'),output_state_file=str(dest/'state.csv'))
            # An additional diagnostic contract, not a change to the old gate.
            gate=cfg['solver']['carrier_row_convergence']
            gate.update(mode='enforce',eps_row=1e-6,scale_floor=0.,min_source_scale=0.,min_source_scale_fraction=0.,
                min_source_global_fraction=0.,min_carrier_density_m3=1e-300,min_flux_scale=0.,min_flux_scale_fraction=0.)
            cfg['solver']['diagnostics']=True
            a.write(dest/'config.json',cfg);files += [dest/'initial.csv',dest/'config.json']
            jobs.append(dict(case=key,device=point['device'],vd=float(point['vd']),vg=vg,initialization=label,config=str(dest/'config.json'),legacy_config=str(src/'config.json')))
    assert len(jobs)==32
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',jobs=jobs,
        new_gate='All positive-density active electron/hole continuity rows, no source/flux exclusion and zero row scale floor; same eps_row=1e-6. Zero-scale rows separately unresolved.',
        legacy_gate='Unmodified historical carrier-row/global/KCL gate also evaluated on each output; no historical qualification rewritten.',
        common_solver='Existing isolated vector-chain runner with independently tested local transport step coefficient 1e-7; SparseLU; original physical residual.',
        initializations='Existing Vela, coherent native, Vela hole phi +/-0.05 V on all free Si nodes with consistent p; original contacts unchanged.',
        gates={'all_row_relative':1e-6,'kcl_over_Id':1e-8,'initialization_phi_max_V':1e-6,'initialization_density_relative':1e-4,'initialization_Id_relative':1e-6},
        scope='8 target biases, 32 independent reclosures; do not assume a failed state has converged; no production or physical model change.',m82_released=False,m83_released=False))
    files.append(OUT/'contract.json');d.matrix.freeze(OUT/'freeze.json',files)
    print('Frozen 32 minority-field reclosures',flush=True)


def execute(path):
    output=path.with_suffix('.status.json')
    if output.exists():return a.read(output)
    r=subprocess.run([str(surface.RUNNER),'--config',str(path),'--log','off'],env=env(),capture_output=True,text=True)
    path.with_suffix('.stdout.txt').write_text(r.stdout);path.with_suffix('.stderr.txt').write_text(r.stderr)
    status=json.loads(r.stdout.strip().splitlines()[-1]);status['exit_code']=r.returncode;a.write(output,status);return status


def one(job):
    path=Path(job['config']);dest=path.parent;s=execute(path)
    result=dict(**job,exit_code=s['exit_code'],converged=s.get('converged',False),reason=s.get('convergence_reason',''),iterations=s.get('iterations',0),
        current_A_per_um=s.get('contact_currents_A_per_um',{}).get('drain',0.),all_row_qualified=False,legacy_qualified=False)
    if (dest/'state.csv').exists():
        audits={}
        for label,source in (('all_row',path),('legacy',Path(job['legacy_config']))):
            cfg=a.read(source);cfg.pop('output_state_file',None)
            cfg.update(simulation_type='newton_carrier_term_probe',state_file=str(dest/'state.csv'),output_csv=str(dest/(label+'.csv')),carrier_term_probe={'solved_equation_terms':True})
            cfg['solver']['carrier_row_convergence']['mode']='report'
            cfg['solver']['global_continuity_closure']={'mode':'enforce','tolerance':1e-6,'source_floor':1e-10}
            target=dest/(label+'.json')
            if not target.exists():a.write(target,cfg)
            audits[label]=execute(target)
        cc=s.get('contact_currents_A_per_um',{});kcl=abs(math.fsum(cc.values()))/max(abs(result['current_A_per_um']),1e-300)
        result['kcl_over_Id']=kcl
        for label,st in audits.items():
            result[label+'_violations']=st['carrier_row_convergence']['violation_count']
            result[label+'_max_ratio']=st['carrier_row_convergence']['max_ratio']
            result[label+'_checked_rows']=st['carrier_row_convergence']['qualified_row_count']
            result[label+'_qualified']=s['converged'] and s['exit_code']==0 and st['exit_code']==0 and st['carrier_row_convergence']['satisfied'] and st['global_continuity_closure']['satisfied'] and kcl<=1e-8
    a.write(dest/'result.json',result);print(job['case'],job['vg'],job['initialization'],result['reason'],result['all_row_qualified'],result.get('all_row_max_ratio'),flush=True)
    return result


def run():
    a.verify(OUT/'freeze.json')
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,a.read(OUT/'contract.json')['jobs']))
    a.write_csv(OUT/'runs.csv',rows)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','run'));globals()[p.parse_args().action]()
