"""Three frozen Vela arms on matching baseline/reduced-model native initial states."""
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import math
import json
import subprocess
import time
import prepare_simplemos_reference_subset_20260906 as p

a=p.a;d=p.d;LOCAL=p.LOCAL;OUT=p.OUT;RUNNER=p.prior.RUNNER
ARMS=('baseline_refined','masetti_plain','masetti_refined')


def env(arm,path):
    e=p.prior.env(arm.endswith('_refined'),path)
    if arm.startswith('masetti'):
        for k in ('VELA_VALIDATE_VECTOR_CHAIN_ENABLE','VELA_VALIDATE_VECTOR_CHAIN_STEP','VELA_VALIDATE_TRANSPORT_STEP'):e.pop(k,None)
    return e


def prepare():
    a.verify(OUT/'native_freeze.json');a.verify(OUT/'export_contract.json')
    native=a.rows(OUT/'native_points.csv');jobs=[];files=[Path(__file__).resolve(),RUNNER,OUT/'native_freeze.json',OUT/'export_contract.json',OUT/'native_points.csv']
    for point in native:
        for arm in ARMS:
            model='baseline' if arm.startswith('baseline') else 'masetti'
            if point['model']!=model:continue
            index=int(point['index']);vg=float(point['vg']);case=point['case'];device=point['device']
            src=LOCAL/'native_exports'/model/case/f'vg_{index:03d}';geo=d.matrix.spatial.m73.Geometry(device)
            psi,n,h,spread=d.fixed.upstream.coherent_state(src,geo);assert spread<1e-12
            en=d.matrix.spatial.m73.scalar(src/'fields/eQuasiFermiPotential_region0.csv');hp=d.matrix.spatial.m73.scalar(src/'fields/hQuasiFermiPotential_region0.csv')
            dest=LOCAL/'vela'/arm/case/f'vg_{index:03d}'
            initial=[dict(node_id=i,psi=psi[i],phin=en.get(i,0.),phip=hp.get(i,0.),electrons_m3=n[i],holes_m3=h[i]) for i in range(geo.count)]
            a.write_csv(dest/'initial.csv',initial)
            original=p.prior.v.m.LOCAL/case/'vg_020/native/config.json';cfg=a.read(original)
            for contact in cfg['contacts']:
                if contact['name']=='gate':contact['bias']=vg
            cfg.update(state_file=str(dest/'initial.csv'),output_state_file=str(dest/'state.csv'))
            if model=='masetti':cfg['solver']['mobility']=dict(model='masetti',doping_concentration_basis='total_impurity')
            assert cfg['solver']['carrier_row_convergence']['eps_row']==1e-6
            assert cfg['solver']['carrier_row_convergence']['scale_floor']==0
            assert cfg['solver']['global_continuity_closure']['mode']=='off'
            cfg['solver']['local_update_diagnostics']=dict(enabled=True,nodes=[41,355,453,949,973,1041],csv_file=str(dest/'updates.csv'),first_iterations=200,every_iterations=1)
            a.write(dest/'config.json',cfg)
            probe=copy.deepcopy(cfg);probe.pop('output_state_file');probe['solver'].pop('local_update_diagnostics')
            probe.update(simulation_type='newton_carrier_term_probe',state_file=str(dest/'state.csv'),output_csv=str(dest/'all_row.csv'),carrier_term_probe={'solved_equation_terms':True})
            probe['solver']['carrier_row_convergence']['mode']='report';probe['solver']['global_continuity_closure']=dict(mode='enforce',tolerance=1e-6,source_floor=1e-10)
            a.write(dest/'all_row.json',probe)
            jobs.append(dict(arm=arm,model=model,case=case,device=device,vd=float(point['vd']),vg=vg,index=index,config=str(dest/'config.json'),
                native_Id_A_per_um=float(point['Id_A_per_um']),native_qualified=point['native_qualified']=='True',native_kcl_over_Id=float(point['kcl_over_Id'])))
            files += [dest/n for n in ('initial.csv','config.json','all_row.json')]+[original]+list((src/'fields').glob('*.csv'))
            files += [Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file')]
    assert len(jobs)==24
    # Same reduced-model initial state/config in both numerical arms, apart from diagnostic destinations.
    for j in jobs:
        if j['arm']!='masetti_plain':continue
        other=next(k for k in jobs if k['arm']=='masetti_refined' and k['case']==j['case'] and k['index']==j['index'])
        assert a.sha(Path(j['config']).parent/'initial.csv')==a.sha(Path(other['config']).parent/'initial.csv')
        x=a.read(Path(j['config']));y=a.read(Path(other['config']))
        for c in (x,y):
            for key in ('state_file','output_state_file'):c.pop(key)
            c['solver'].pop('local_update_diagnostics')
        assert x==y
    a.write(OUT/'vela_contract.json',dict(status='frozen_before_execution',jobs=jobs,arms=ARMS,
        scope='24 independent native-coherent reclosures on 8 operating points; no full Id-Vg or initialization-invariance claim.',
        initialization='Preserve native psi and electron/hole quasi-Fermi potentials; recompute coherent initial densities with the frozen Vela matched ni. Use the corresponding model native solution; not a frozen-state current substitution.',
        numerical_axes='Baseline uses existing vector-chain patch and four linear refinements. Both Masetti arms disable experimental vector/transport FD controls; the only difference between Masetti arms is the fixed four-refinement switch.',
        preserved='All nonlinear tolerances, caps, scalar line-search, all-row and post global/KCL qualification unchanged from current strict diagnostic. Fixed material/geometry/contacts.',
        qualified_comparison='Require native qualification and Vela nonlinear success, all 1814 active rows, no zero-scale exclusion, global and KCL checks. Keep failed results as diagnostic only.',
        production_changes=False))
    files.append(OUT/'vela_contract.json');d.matrix.freeze(OUT/'vela_freeze.json',files)
    print('Frozen 24 Vela reclosures, three arms',flush=True)


def execute(path,arm):
    out=path.with_suffix('.status.json')
    if out.exists():return a.read(out)
    start=time.monotonic();r=subprocess.run([str(RUNNER),'--config',str(path),'--log','off'],env=env(arm,path.parent),capture_output=True,text=True)
    path.with_suffix('.stdout.txt').write_text(r.stdout);path.with_suffix('.stderr.txt').write_text(r.stderr)
    s=json.loads(r.stdout.strip().splitlines()[-1]);s.update(exit_code=r.returncode,elapsed_seconds=time.monotonic()-start);a.write(out,s);return s


def one(job):
    path=Path(job['config']);dest=path.parent;s=execute(path,job['arm']);assert (dest/'state.csv').exists(),s
    q=execute(dest/'all_row.json','masetti_plain')
    cc=s['contact_currents_A_per_um'];current=cc['drain'];kcl=abs(math.fsum(cc.values()))/max(abs(current),1e-300)
    gate=q['carrier_row_convergence'];geo=d.matrix.spatial.m73.Geometry(job['device'])
    mask=d.matrix.spatial.old.m78.supports(job['device'],geo,.05)[0]['all_si'].copy();mask[geo.contact_nodes]=False
    terms=d.ordered(dest/'all_row.csv',geo.count)
    zero=sum(max(float(t[c+'_flux_abs_sum']),abs(float(t[c+'_recombination'])),abs(float(t[c+'_impact'])))==0 for t,keep in zip(terms,mask) if keep for c in ('electron','hole'))
    assert gate['qualified_row_count']==2*int(sum(mask))==1814
    qualified=s['exit_code']==q['exit_code']==0 and s['converged'] and gate['satisfied'] and zero==0 and q['global_continuity_closure']['satisfied'] and kcl<=1e-8
    row=dict(**job,converged=s['converged'],exit_code=s['exit_code'],iterations=s['iterations'],failure=s['failure_reason'],
        Id_A_per_um=current,kcl_over_Id=kcl,all_row_violations=gate['violation_count'],all_row_max_ratio=gate['max_ratio'],zero_scale_rows=zero,
        electron_violations=sum(r['carrier']=='electron' for r in gate['violations']),hole_violations=sum(r['carrier']=='hole' for r in gate['violations']),
        global_satisfied=q['global_continuity_closure']['satisfied'],vela_qualified=qualified,comparison_qualified=qualified and job['native_qualified'],
        signed_Id_error_relative=current/job['native_Id_A_per_um']-1,Id_error_dex=math.log10(abs(current/job['native_Id_A_per_um'])),elapsed_seconds=s['elapsed_seconds'])
    a.write(dest/'result.json',row);print(job['arm'],job['case'],job['vg'],'qualified',row['comparison_qualified'],'rows',row['all_row_violations'],'error',row['signed_Id_error_relative'],flush=True)
    return row


def run():
    a.verify(OUT/'vela_freeze.json')
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,a.read(OUT/'vela_contract.json')['jobs']))
    a.write_csv(OUT/'vela_runs.csv',rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run'));globals()[parser.parse_args().action]()
