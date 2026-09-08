"""Frozen production migration, independent paths and unchanged DC gates."""
import argparse,copy,json,math,os,subprocess,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import simplemos_production_migration_20260907 as m
a=m.a;d=m.d;p=m.p;LOCAL=m.LOCAL;OUT=m.OUT;RUNNER=m.RUNNER

def environment():
    env=os.environ.copy();env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env['PATH']
    for key in list(env):
        if key.startswith(('VELA_DIAGNOSTIC_','VELA_VALIDATE_','VELA_MINORITY_','VELA_SIMPLEMOS_','VELA_TEST_','VELA_NATIVE_GEOMETRY_','VELA_CANDIDATE_')):env.pop(key)
    return env

def execute(path,runner=RUNNER,env=None):
    status=path.with_suffix('.status.json')
    if status.exists():return a.read(status)
    start=time.monotonic();r=subprocess.run([str(runner),'--config',str(path),'--log','off'],env=env or environment(),capture_output=True,text=True)
    path.with_suffix('.stdout.txt').write_text(r.stdout);path.with_suffix('.stderr.txt').write_text(r.stderr)
    try:s=json.loads(r.stdout.strip().splitlines()[-1])
    except (ValueError,IndexError):s=dict(failure_reason='No JSON',stderr=r.stderr[-2000:])
    s.update(exit_code=r.returncode,elapsed_seconds=time.monotonic()-start);a.write(status,s);return s

def config(c):
    cfg=a.read(Path(c['base'])/'config.json');cfg['solver'].pop('local_update_diagnostics',None)
    cfg['solver']['linear_refinement_iterations']=4
    ratios=[float(x.split()[1]) for x in (Path(c['ratios'])/'edges.txt').read_text().splitlines()]
    cfg['solver']['region_resolved_interface_assembly']=dict(poisson_charge_node_volume='signed_transport',transport_edge_coupling_ratios=ratios)
    return cfg

def probes(cfg,dest,state):
    result=[]
    for name,kind in [('functional','terminal_current_functional_probe'),('edges','sg_edge_flux_probe'),('terms','newton_carrier_term_probe'),('mobility','edge_mobility_probe'),('adjoint','terminal_current_adjoint_probe')]:
        deck=copy.deepcopy(cfg);deck.pop('output_state_file',None)
        deck.update(simulation_type=kind,state_file=str(state),contact='drain')
        deck['solver']['carrier_row_convergence']['mode']='report';deck['solver']['global_continuity_closure']['mode']='off'
        if name=='functional':deck.update(residual_output_csv=str(dest/'residual.csv'),contact_edge_output_csv=str(dest/'contact_edges.csv'))
        else:deck['output_csv']=str(dest/(name+'.csv'))
        if name=='terms':deck['carrier_term_probe']={'solved_equation_terms':True}
        path=dest/(name+'.json');a.write(path,deck);result.append(path)
    return result

def post_config(cfg,dest):
    q=copy.deepcopy(cfg);q.pop('output_state_file',None);q['solver'].pop('local_update_diagnostics',None)
    q.update(simulation_type='newton_carrier_term_probe',state_file=str(dest/'state.csv'),output_csv=str(dest/'all_row.csv'),carrier_term_probe={'solved_equation_terms':True})
    q['solver']['carrier_row_convergence']['mode']='report';q['solver']['global_continuity_closure']=dict(mode='enforce',tolerance=1e-6,source_floor=1e-10)
    a.write(dest/'all_row.json',q)

def prepare():
    m.verify_historical(m.previous.OUT/'validation_evidence.json')
    cases=m.previous.prior.cases();files=[Path(__file__).resolve(),Path(m.__file__).resolve(),RUNNER,OUT/'prechange_evidence.json']
    for c in cases:
        root=LOCAL/c['key'];cfg=config(c);c['probes']=[];c['jobs']=[]
        states=dict(joint=Path(c['root'])/'joint/replacement/state.csv',native_raw=Path(c['native_initial']),
            native_referenced=m.previous.LOCAL/c['key']/'native_referenced_joint/state.csv')
        for label,state in states.items():
            paths=probes(cfg,root/'fixed'/label,state)
            c['probes'] += [str(x) for x in paths if label=='joint' or x.stem in ('functional','edges')]
            files+=paths+[state]
        for label,state,refine in [('replay',states['joint'],4),('from_vela',Path(c['base'])/'state.csv',4),
                                    ('from_native',states['native_referenced'],4),('without_refinement',Path(c['base'])/'state.csv',0)]:
            dest=root/label;deck=copy.deepcopy(cfg);deck.update(state_file=str(state),output_state_file=str(dest/'state.csv'))
            deck['solver']['linear_refinement_iterations']=refine
            a.write(dest/'config.json',deck);post_config(deck,dest)
            c['jobs'].append(dict(label=label,dest=str(dest),refinement=refine));files += [state,dest/'config.json',dest/'all_row.json']
        files += [Path(cfg[k]) for k in ('mesh_file','node_doping_file','materials_file')]
    for name in a.read(OUT/'prechange.json')['files']:
        if name.startswith(('src/','include/')):files.append(p.REPO/name)
    files += list((p.REPO/'include/vela').rglob('*Refinement.h'))+[p.REPO/'include/vela/discretization/StableSGDerivative.h']
    gates=a.read(m.b.OUT/'contract.json')['gates']
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',cases=cases,DC=32,gates=gates,
        geometry='Explicit native Si transport-edge ratios and independent signed Si charge volume; original dielectric K and SRH volume.',
        numerical='Physical reference subtraction before scaling; stable unclipped equal-ni SG psi partials; optional MP100 linear defect, four corrections.',
        controls='Eight prior qualified joint-state replays; same eight replacements from former Vela and native states; refinement-off failures retained as controls.',
        acceptance_changed=False,defaults_changed=False))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json']);print('Frozen production 72 probes and 32 DC including 8 refinement-off controls',flush=True)

def qualify(c,dest):
    s=a.read(dest/'config.status.json');q=a.read(dest/'all_row.status.json');geo,mask=m.previous.prior.support(c)
    rr=[]
    for r,keep in zip(d.ordered(dest/'all_row.csv',geo.count),mask):
        if keep:
            for car in ('electron','hole'):
                scale=max(float(r[car+'_flux_abs_sum']),abs(float(r[car+'_recombination'])),abs(float(r[car+'_impact'])))
                rr.append(abs(float(r[car+'_residual']))/scale if scale else math.inf)
    assert len(rr)==q['carrier_row_convergence']['qualified_row_count']==1814
    bad=sum(x>1e-6 for x in rr);assert bad==q['carrier_row_convergence']['violation_count']
    cc=s.get('contact_currents_A_per_um',{});Id=cc.get('drain',math.nan);kcl=abs(math.fsum(cc.values()))/abs(Id)
    return dict(qualified=s['exit_code']==q['exit_code']==0 and s.get('converged',False) and bad==0 and q['global_continuity_closure']['satisfied'] and kcl<=1e-8,
        current_A_per_um=Id,max_row_ratio=max(rr),row_violations=bad,kcl_over_Id=kcl,iterations=s.get('iterations'),failure=s.get('failure_reason',''))

def run():
    a.verify(OUT/'freeze.json')
    def one(c):
        for path in c['probes']:
            s=execute(Path(path));assert s['exit_code']==0,(path,s)
        rows=[]
        for j in c['jobs']:
            dest=Path(j['dest']);s=execute(dest/'config.json');row=dict(key=c['key'],label=j['label'],device=c['device'],vg=c['vg'],vd=c['vd'],qualified=False,failure=s.get('failure_reason',''))
            if (dest/'state.csv').exists():
                execute(dest/'all_row.json');row.update(qualify(c,dest))
            rows.append(row);print(c['key'],j['label'],row['qualified'],row.get('max_row_ratio'),flush=True)
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/'dc.csv',[{k:r.get(k,'') for k in keys} for r in rows])

def analyze():
    a.verify(OUT/'freeze.json');rows=a.rows(OUT/'dc.csv');ports=[];paths=[];equivalence=[];dual=[]
    for c in a.read(OUT/'contract.json')['cases']:
        root=LOCAL/c['key'];geo,mask=m.previous.prior.support(c)
        for label in ('joint','native_raw','native_referenced'):
            dest=root/'fixed'/label;s=a.read(dest/'functional.status.json')
            err=abs(s['current_A_per_um']/s['contact_current_extractor_A_per_um']-1)
            ports.append(dict(key=c['key'],state=label,relative_error=err,qualified=err<=1e-8))
        edges=a.rows(root/'fixed/joint/edges.csv');mu=a.rows(root/'fixed/joint/mobility.csv');old=a.rows(m.previous.LOCAL/c['key']/'joint/edges.csv')
        for field in ('couple_m','electron_mobility_m2_V_s','hole_mobility_m2_V_s','electron_flux','hole_flux'):
            err=max(abs(float(x[field])-float(y[field]))/max(abs(float(x[field])),abs(float(y[field])),1e-300) for x,y in zip(edges,old))
            paths.append(dict(key=c['key'],check='isolated_'+field,relative_error=err,qualified=err<=1e-10))
        err=max(abs(float(x['couple_m'])-float(y['couple_m']))/max(abs(float(x['couple_m'])),1e-300) for x,y in zip(edges,mu))
        paths.append(dict(key=c['key'],check='ordinary_mobility_geometry',relative_error=err,qualified=err<=1e-14))
        states={}
        for label in ('replay','from_vela','from_native'):
            dest=root/label;info=qualify(c,dest);states[label]=d.ordered(dest/'state.csv',geo.count)
            oldstate=d.ordered(Path(c['root'])/'joint/replacement/state.csv',geo.count)
            diff=m.previous.prior.delta_states(states[label],oldstate,mask)
            oldId=a.read(Path(c['root'])/'joint/replacement/config.status.json')['contact_currents_A_per_um']['drain']
            err=abs(info['current_A_per_um']/oldId-1);ok=info['qualified'] and err<=1e-6 and max(diff[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and diff['density_max_relative']<=1e-4
            equivalence.append(dict(key=c['key'],label=label,Id_relative=err,**diff,qualified=ok))
        diff=m.previous.prior.delta_states(states['from_vela'],states['from_native'],mask)
        currents=[qualify(c,root/label)['current_A_per_um'] for label in ('from_vela','from_native')];err=abs(currents[0]/currents[1]-1)
        dual.append(dict(key=c['key'],Id_relative=err,**diff,qualified=err<=1e-6 and max(diff[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and diff['density_max_relative']<=1e-4))
    for name,data in [('ports',ports),('paths',paths),('equivalence',equivalence),('dual_init',dual)]:a.write_csv(OUT/(name+'.csv'),data)
    ok=all(r['qualified'] for r in ports+paths+equivalence+dual)
    a.write(OUT/'summary.json',dict(production_qualified=ok,DC=32,required_DC_qualified=sum(r['qualified']=='True' for r in rows if r['label']!='without_refinement'),
        refinement_off_qualified=sum(r['qualified']=='True' for r in rows if r['label']=='without_refinement')))
    files=[OUT/'freeze.json',OUT/'summary.json']+list(OUT.glob('*.csv'))+[x for c in a.read(OUT/'contract.json')['cases'] for x in (LOCAL/c['key']).rglob('*') if x.is_file()]
    d.matrix.freeze(OUT/'validation_evidence.json',files);print('Production migration qualified:',ok,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run','analyze'));globals()[parser.parse_args().action]()
