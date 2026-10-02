"""Explicit Enormal candidate: fixed-state replay, derivative and dual-DC gates."""
import argparse,copy,math,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import calibrate_simplemos_enormal_20260912 as p
a,d,c=p.a,p.d,p.c;v=c.v;V=v.V;R=p.R
L=R/'build-release/enormal_candidate_20260912/v1'
O=R/'reference_tcad/simplemos_sentaurus2022/enormal_candidate_20260912/v1'
RUNNER=L/'bin/vela_example_runner.exe'

def prepare_inputs():
    a.verify(p.O/'floor_position_evidence.json');assert not (O/'inputs.json').exists()
    accepted=a.rows(c.OUT/'continuation_attempts.csv');cases=[];files=[Path(__file__).resolve(),p.O/'floor_position_evidence.json']
    for j in a.read(p.O/'native_contract.json')['jobs']:
        if j['arm']!='baseline':continue
        key=j['name'].removesuffix('_baseline');dest=L/'inputs'/key
        previous=next(r for r in accepted if r['case']==j['case'] and int(r['index'])==j['index'] and r['qualified']=='True')
        phumob_seed=Path(previous['dest'])/'state.csv';cfg=a.read(Path(previous['dest'])/'config.json')
        cfg['solver']['mobility'].update(model='phumob_lombardi',surface=dict(discretization='element_distance_gradient',acoustic_factor=1.,roughness_factor=1.))
        cfg.pop('state_file');cfg.pop('output_state_file');a.write(dest/'template.json',cfg)
        cc=dict(case=j['case'],key=key,device=j['device'],vd=j['vd'],vg=j['vg'],index=j['index'],template=str(dest/'template.json'),phumob_seed=str(phumob_seed))
        src=p.L/'exports'/j['name'];geo,_=V.m.previous.prior.support(cc)
        psi,en,hp,spread=c.bgn.mapping.RAW_READER(src,geo);assert spread<=1e-12
        scalar=c.bgn.mapping.original.m73.scalar;fn=scalar(src/'fields/eQuasiFermiPotential_region0.csv');fp=scalar(src/'fields/hQuasiFermiPotential_region0.csv')
        raw=[dict(node_id=i,psi=psi[i],phin=fn.get(i,0.),phip=fp.get(i,0.),electrons_m3=en[i],holes_m3=hp[i]) for i in range(geo.count)]
        a.write_csv(dest/'native_raw.csv',raw)
        doping={int(r['node_id']):float(r['donors_cm3'])+float(r['acceptors_cm3']) for r in a.rows(Path(cfg['node_doping_file']))}
        ni=next(m['ni'] for m in a.read(Path(cfg['materials_file']))['materials'] if m['name']=='Si')*1e6
        for i in fn:
            effective=ni*math.exp(c.bgn.delta_eg(doping[i],'old_slotboom')/(2*c.bgn.VT))
            raw[i]['electrons_m3']=effective*math.exp((psi[i]-fn[i])/c.bgn.VT);raw[i]['holes_m3']=effective*math.exp((fp[i]-psi[i])/c.bgn.VT)
        a.write_csv(dest/'native_coherent.csv',raw);cc.update(native_seed=str(dest/'native_coherent.csv'),raw_native=str(dest/'native_raw.csv'))
        cases.append(cc);files += [phumob_seed,dest/'template.json',dest/'native_raw.csv',dest/'native_coherent.csv']+[Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file')]
    a.write(O/'inputs.json',dict(cases=cases,gates=v.GATES,scope='Restore Enormal only, alpha=0, original physical parameters and all acceptance gates. Raw native fields are for fixed constitutive replay; coherent native fields and prior qualified PhuMob states are independent DC starts.',weak_derivatives='Wide residual differences permit checking weak cross blocks; do not reuse old whole-double-residual exclusions as a passing result.'))
    d.matrix.freeze(O/'input_evidence.json',files+[O/'inputs.json'])
    print('Prepared 8 targets, raw physical replay and two independent starts.',flush=True)

def freeze():
    a.verify(O/'input_evidence.json');assert not (O/'candidate_freeze.json').exists()
    RUNNER.parent.mkdir(parents=True,exist_ok=False);shutil.copyfile(R/'build-release/vela_example_runner.exe',RUNNER)
    files=[RUNNER,O/'input_evidence.json',Path(__file__).resolve(),R/'build-release/libvela_core.a',R/'CMakeLists.txt',R/'docs/config_schema.md']
    files += [f for folder in ('src','include') for f in (R/folder).rglob('*') if f.suffix in ('.h','.cpp')]
    files += list((R/'scripts').glob('*simplemos*.py'))
    d.matrix.freeze(O/'candidate_freeze.json',files)

def config_case(cc):return copy.deepcopy(a.read(Path(cc['template'])))

def fixed():
    a.verify(O/'candidate_freeze.json');results=[]
    expected={(r['key'],r['carrier'],int(r['cell'])):r for r in a.rows(p.O/'cell_formation.csv') if r['mode']=='mathematical_floor'}
    for cc in a.read(O/'inputs.json')['cases']:
        dest=L/'fixed'/cc['key'];cfg=config_case(cc);cfg.update(simulation_type='element_box_probe',state_file=cc['raw_native'],output_csv=str(dest/'cells.csv'))
        a.write(dest/'cells.json',cfg);status=V.execute(dest/'cells.json',RUNNER,V.environment());assert status['exit_code']==0,status
        errors=[]
        for row in a.rows(dest/'cells.csv'):
            cid=int(row['cell_id'])
            for car,name in [('e','electron'),('h','hole')]:
                key=(cc['key'],car,cid)
                if key not in expected:continue
                errors.append(abs(float(row[name+'_mobility_m2_V_s'])*1e4/float(expected[key]['formula_mu'])-1))
        results.append(dict(key=cc['key'],samples=len(errors),max_relative=max(errors),qualified=max(errors)<=1e-11));print(results[-1],flush=True)
    a.write_csv(O/'fixed_cells.csv',results);d.matrix.freeze(O/'fixed_evidence.json',[O/'candidate_freeze.json',O/'fixed_cells.csv']+[f for f in (L/'fixed').rglob('*') if f.is_file()]);assert all(r['qualified'] for r in results)

def smoke():
    a.verify(O/'candidate_freeze.json');cc=a.read(O/'inputs.json')['cases'][0];dest=L/'smoke';cfg=config_case(cc)
    cfg.update(simulation_type='newton_jvp_probe',state_file=cc['native_seed'],output_csv=str(dest/'jvp.csv'),directions=[dict(name=f'{mode}_{h:.0e}',mode=mode,amplitude_V=h,node_ids=[1000],exclude_contacts=False) for mode in ('psi','phin','phip') for h in (1e-6,1e-20)])
    a.write(dest/'jvp.json',cfg);status=V.execute(dest/'jvp.json',RUNNER,V.environment());assert status['exit_code']==0,status
    rows=[]
    for row in a.rows(dest/'jvp.csv'):
        for block in ('psi','phin','phip'):
            an=float(row[f'analytic_{block}_norm']);fd=float(row[f'finite_difference_{block}_norm']);error=float(row[f'{block}_relative_error'])*max(1.,fd)/max(an,fd,1e-300)
            rows.append(dict(direction=row['direction'],block=block,relative=error,qualified=error<=1e-4))
    a.write_csv(O/'smoke_jvp.csv',rows);d.matrix.freeze(O/'smoke_evidence.json',[O/'candidate_freeze.json',O/'smoke_jvp.csv']+[f for f in dest.iterdir() if f.is_file()]);print(status,rows,flush=True);assert all(r['qualified'] for r in rows)

def configure():
    v.LOCAL=L;v.OUT=O;v.RUNNER=RUNNER;V.post_config=c.qualification.previous.post_config

def dc():
    a.verify(O/'fixed_evidence.json');a.verify(O/'smoke_evidence.json');a.verify(O/'candidate_freeze.json');configure()
    assert all(r['qualified']=='True' for r in a.rows(O/'fixed_cells.csv')) and all(r['qualified']=='True' for r in a.rows(O/'smoke_jvp.csv'))
    jobs=[(cc,arm,cc['phumob_seed' if arm=='phumob_seed' else 'native_seed']) for cc in a.read(O/'inputs.json')['cases'] for arm in ('phumob_seed','native_seed')]
    def run(job):cc,arm,seed=job;return v.solve_target(cc,cc['index'],arm,Path(seed))
    with ThreadPoolExecutor(max_workers=2) as pool:attempts=[row for result in pool.map(run,jobs) for row in result]
    v.csv_union(O/'attempts.csv',attempts)
    d.matrix.freeze(O/'dc_evidence.json',[O/'candidate_freeze.json',O/'attempts.csv']+[f for f in (L/'vela').rglob('*') if f.is_file()])

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare_inputs','freeze','fixed','smoke','dc'));globals()[parser.parse_args().action]()
