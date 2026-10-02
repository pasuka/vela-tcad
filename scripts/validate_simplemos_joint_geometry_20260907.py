"""Frozen 2-axis coefficient response and full self-consistent comparison."""
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import subprocess
import time
import numpy as np
import build_simplemos_joint_geometry_20260907 as b
import audit_simplemos_masetti_box_mobility_20260907 as m

p=b.p;a=p.a;d=p.d;LOCAL=b.LOCAL;OUT=b.OUT
AXES={'baseline':(0.,0.),'transport':(1.,0.),'poisson_volume':(0.,1.),'joint':(1.,1.)}


def prepare():
    a.verify(p.OUT/'validation_evidence.json');assert b.RUNNER.exists()
    files=[Path(__file__).resolve(),Path(b.__file__).resolve(),b.HEADER,b.RUNNER,p.OUT/'validation_evidence.json']
    files += [Path(x) for x in a.read(LOCAL/'build_command.json') if Path(x).is_file() and Path(x).is_relative_to(p.REPO)]
    old=a.read(p.OUT/'distributed_response/contract.json')['cases'];cases=[];geometry=[];predictions=[]
    factor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge']
    for device in ('n19','n23'):
        geo,xy,cells,vol,K,parts,info,modified=m.geom.geometry(device)
        c=next(c for c in old if c['device']==device);edges=a.rows(p.prior.LOCAL/'vela'/c['key']/'strict/edges.csv')
        root=LOCAL/'ratios'/device;root.mkdir(parents=True,exist_ok=False);lines=[]
        for e in edges:
            i,j=int(e['node0']),int(e['node1']);g0=float(e['couple_m'])/float(e['length_m'])
            gn=math.fsum(x['coefficient'] for x in parts[(i,j)] if x['material']=='Si')
            if g0<=0:assert gn==0;ratio=1.
            else:ratio=gn/g0
            assert ratio>=0
            lines.append(f"{e['edge_id']} {ratio:.17g}\n")
            geometry.append(dict(device=device,kind='edge',id=e['edge_id'],baseline=g0,native_Si=gn,ratio=ratio))
        (root/'edges.txt').write_text(''.join(lines),newline='\n');lines=[]
        for i,(v0,vn) in enumerate(zip(geo.volumes['all_cell'],vol['Si'])):
            assert v0>0 and vn>=0;ratio=float(vn/v0);lines.append(f'{i} {ratio:.17g}\n')
            geometry.append(dict(device=device,kind='node',id=i,baseline=float(v0),native_Si=float(vn),ratio=ratio))
        (root/'nodes.txt').write_text(''.join(lines),newline='\n');files += [root/'edges.txt',root/'nodes.txt']
    for oldc in old:
        c=copy.deepcopy(oldc);base=Path(c['base']);geo=d.matrix.spatial.m73.Geometry(c['device']);root=LOCAL/c['key']
        cfg=a.read(base/'config.json');c['ratios']=str(LOCAL/'ratios'/c['device']);c['native_initial']=str(base/'initial.csv')
        c['native_Id_A_per_um']=next(x['native_Id_A_per_um'] for x in a.read(p.OUT/'native_contract.json')['jobs'] if x['key']==c['key'])
        c['jobs']=[];c['probes']=[]
        eRat={int(x.split()[0]):float(x.split()[1]) for x in (Path(c['ratios'])/'edges.txt').read_text().splitlines()}
        nodeRat=np.array([float(x.split()[1]) for x in (Path(c['ratios'])/'nodes.txt').read_text().splitlines()])
        edges=a.rows(p.prior.LOCAL/'vela'/c['key']/'strict/edges.csv');oldR=d.ordered(p.prior.LOCAL/'vela'/c['key']/'strict/residual.csv',geo.count)
        state=d.ordered(base/'state.csv',geo.count);adj=d.ordered(p.prior.LOCAL/'vela'/c['key']/'adjoint.csv',geo.count)
        lam=np.array([[float(x[k]) for x in adj] for k in ('lambda_poisson','lambda_electron','lambda_hole')]);lam[:,geo.contact_nodes]=0
        rhsT=np.zeros((3,geo.count));absT=np.zeros_like(rhsT);rhsP=np.zeros_like(rhsT)
        mesh=a.read(Path(cfg['mesh_file']));drain=set(next(x['node_ids'] for x in mesh['contacts'] if x['name']=='drain'));direct=[]
        for e in edges:
            i,j=int(e['node0']),int(e['node1']);delta=eRat[int(e['edge_id'])]-1
            for block,car in [(1,'electron'),(2,'hole')]:
                f=float(e[car+'_flux'])*delta;rhsT[block,i]+=f;rhsT[block,j]-=f;absT[block,i]+=abs(f);absT[block,j]+=abs(f)
            direct.append(-d.fixed.Q*1e-6*(int(i in drain)-int(j in drain))*(float(e['electron_particle_line_flux_per_m_s'])-float(e['hole_particle_line_flux_per_m_s']))*delta)
        for i,(s,r) in enumerate(zip(state,oldR)):
            charge=float(s['electrons_m3'])-float(s['holes_m3'])-float(r['net_doping_m3'])
            rhsP[0,i]=factor*d.fixed.Q*charge*geo.volumes['all_cell'][i]*(nodeRat[i]-1)
        rhsT[:,geo.contact_nodes]=0;rhsP[:,geo.contact_nodes]=0
        for arm,(T,P) in AXES.items():
            rhs=T*rhsT+P*rhsP;feedback=-math.fsum(float(x*y) for x,y in zip(lam.flat,rhs.flat));di=T*math.fsum(direct)
            predictions.append(dict(key=c['key'],arm=arm,unit_feedback_A_per_um=feedback,unit_direct_A_per_um=di,unit_prediction_A_per_um=feedback+di))
            rows=[dict(node_id=i,psi_source=rhs[0,i],electron_source=rhs[1,i],hole_source=rhs[2,i],electron_source_abs=absT[1,i]*T,hole_source_abs=absT[2,i]*T) for i in range(geo.count)]
            source=root/'preflight'/arm/'expected.csv';a.write_csv(source,rows);files.append(source)
            for name,kind in [('functional','terminal_current_functional_probe'),('edges','sg_edge_flux_probe'),('jvp','joint_geometry_jvp')]:
                dest=root/'preflight'/arm;deck=copy.deepcopy(cfg);deck.pop('output_state_file');deck['solver'].pop('local_update_diagnostics')
                deck.update(simulation_type=kind,state_file=str(base/'state.csv'),contact='drain')
                deck['residual_output_csv' if name=='functional' else 'output_csv']=str(dest/('residual.csv' if name=='functional' else name+'.csv'))
                a.write(dest/(name+'.json'),deck);files.append(dest/(name+'.json'));c['probes'].append(dict(arm=arm,T=T,P=P,config=str(dest/(name+'.json'))))
            labels=[('zero',0.,'vela')] if arm=='baseline' else [(n,s,'vela') for n,s in [('plus_full',.001),('minus_full',-.001),('plus_half',.0005),('minus_half',-.0005),('replacement',1.)]]+[('replacement_native',1.,'native')]
            for label,scale,initial in labels:
                dest=root/arm/label;deck=copy.deepcopy(cfg);deck.update(state_file=str(base/('initial.csv' if initial=='native' else 'state.csv')),output_state_file=str(dest/'state.csv'))
                deck['solver']['local_update_diagnostics']['csv_file']=str(dest/'updates.csv');a.write(dest/'config.json',deck)
                q=copy.deepcopy(deck);q.pop('output_state_file');q['solver'].pop('local_update_diagnostics')
                q.update(simulation_type='newton_carrier_term_probe',state_file=str(dest/'state.csv'),output_csv=str(dest/'all_row.csv'),carrier_term_probe={'solved_equation_terms':True})
                q['solver']['carrier_row_convergence']['mode']='report';q['solver']['global_continuity_closure']=dict(mode='enforce',tolerance=1e-6,source_floor=1e-10);a.write(dest/'all_row.json',q)
                j=copy.deepcopy(q);j.pop('carrier_term_probe');j.update(simulation_type='joint_geometry_jvp',output_csv=str(dest/'jvp.csv'));a.write(dest/'jvp.json',j)
                fun=copy.deepcopy(j);fun.update(simulation_type='terminal_current_functional_probe',contact='drain',residual_output_csv=str(dest/'residual.csv'));fun.pop('output_csv');a.write(dest/'functional.json',fun)
                c['jobs'].append(dict(arm=arm,label=label,scale=scale,T=T*scale,P=P*scale,initialization=initial,config=str(dest/'config.json')))
                files += [dest/n for n in ('config.json','all_row.json','jvp.json','functional.json')]
        cases.append(c);files += [base/n for n in ('config.json','state.csv','initial.csv','result.json')]
        files += [p.prior.LOCAL/'vela'/c['key']/n for n in ('strict/edges.csv','strict/residual.csv','adjoint.csv')]
    a.write_csv(OUT/'geometry.csv',geometry);a.write_csv(OUT/'predictions.csv',predictions)
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',cases=cases,new_DC=sum(len(c['jobs']) for c in cases),
        axes=dict(transport='Both carriers use summed native Si element-edge geometry. Masetti edge mobility unchanged.',
            poisson_volume='Electron, hole and dopant Poisson charge volumes all use native Si node volume. Dielectric K, continuity/SRH source volumes unchanged.',
            joint='Both complete coefficient changes together. These coefficients remain in the residual and Jacobian throughout Newton.'),
        gates=dict(all_row=1e-6,kcl_over_Id=1e-8,base_current_drift=1e-8,prediction_relative=.001,two_amplitude_relative=.001,even_over_odd=.01,minimum_signal_over_zero_drift=100,zero_drift_dex=1e-5,
            initialization_phi_max_V=1e-6,initialization_density_relative=1e-4,initialization_Id_relative=1e-6,
            preflight_source_relative=1e-8,preflight_roundoff_baseline_flow_eps=64,preflight_absolute=1e-25,jvp_block_relative=1e-4),
        jvp='All three input blocks, sine direction on all non-contact nodes, steps 1e-5/5e-6/2.5e-6 V. Gate the seven potentially changed blocks without floor-1 normalization at the two smaller steps. The electron-hole and hole-electron SRH cross blocks must remain analytically bit-identical across fixed-state variants; retain their FD noise instead of certifying weak columns again. This is not an exhaustive column audit.',
        source_preflight='Compare coefficient replacement residual differences to independently formed edge/charge sources. Allow 64 machine eps times baseline absolute edge flow for residual evaluation roundoff; not a nonlinear acceptance gate.',
        response='Same four qualified Vg=1 V points; +/- .001 and .0005 for each parameter direction, common zero. Full replacement from Vela and native-coherent initial states.',
        promotion='All states qualified; high-NWell absolute relative errors and high/low log-current pairing improve at both Vd; low-NWell errors must not increase.',
        nonlinear_acceptance_changes=False,production_changes=False,full_native_mobility_replacement=False,full_native_Poisson_K_replacement=False))
    files += [OUT/n for n in ('geometry.csv','predictions.csv','contract.json')];d.matrix.freeze(OUT/'freeze.json',files)
    print('Frozen',sum(len(c['jobs']) for c in cases),'DC cases with original solver gates',flush=True)


def execute(path,c,T,P):
    dest=path.with_suffix('.status.json')
    if dest.exists():return a.read(dest)
    start=time.monotonic();r=subprocess.run([str(b.RUNNER),'--config',str(path),'--log','off'],env=b.env(path.parent,Path(c['ratios']),T,P),capture_output=True,text=True)
    path.with_suffix('.stdout.txt').write_text(r.stdout);path.with_suffix('.stderr.txt').write_text(r.stderr)
    try:s=json.loads(r.stdout.strip().splitlines()[-1])
    except (ValueError,IndexError):s={'failure_reason':'No JSON status','stderr_tail':r.stderr[-2000:]}
    s.update(exit_code=r.returncode,elapsed_seconds=time.monotonic()-start);a.write(dest,s);return s


def preflight():
    a.verify(OUT/'freeze.json')
    def one(c):
        for job in c['probes']:
            s=execute(Path(job['config']),c,job['T'],job['P']);assert s['exit_code']==0,(job,s)
        print('Preflight exported',c['key'],flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,a.read(OUT/'contract.json')['cases']))


def run():
    a.verify(OUT/'freeze.json');a.verify(OUT/'preflight_evidence.json')
    def one(c):
        result=[];geo=d.matrix.spatial.m73.Geometry(c['device']);mask=d.matrix.spatial.old.m78.supports(c['device'],geo,.05)[0]['all_si'].copy();mask[geo.contact_nodes]=False
        for j in c['jobs']:
            path=Path(j['config']);s=execute(path,c,j['T'],j['P']);dest=path.parent
            row=dict(key=c['key'],device=c['device'],vd=c['vd'],arm=j['arm'],label=j['label'],T=j['T'],P=j['P'],initialization=j['initialization'],qualified=False,exit_code=s['exit_code'],failure=s.get('failure_reason',''))
            if (dest/'state.csv').exists():
                q=execute(dest/'all_row.json',c,j['T'],j['P']);cc=s.get('contact_currents_A_per_um',{});ratios=[];zero=0
                for x,keep in zip(d.ordered(dest/'all_row.csv',geo.count),mask):
                    if not keep:continue
                    for car in ('electron','hole'):
                        scale=max(float(x[car+'_flux_abs_sum']),abs(float(x[car+'_recombination'])),abs(float(x[car+'_impact'])));zero+=scale==0;ratios.append(abs(float(x[car+'_residual']))/scale if scale else math.inf)
                assert len(ratios)==q['carrier_row_convergence']['qualified_row_count']==1814
                bad=sum(x>1e-6 for x in ratios);assert bad==q['carrier_row_convergence']['violation_count']
                current=cc.get('drain',float('nan'));kcl=abs(math.fsum(cc.values()))/abs(current)
                passed=s['exit_code']==q['exit_code']==0 and s.get('converged',False) and bad==zero==0 and q['global_continuity_closure']['satisfied'] and kcl<=1e-8
                row.update(current_A_per_um=current,signed_Id_relative_error=current/c['native_Id_A_per_um']-1,qualified=passed,iterations=s.get('iterations'),max_row_ratio=max(ratios),row_violations=bad,zero_scale_rows=zero,kcl_over_Id=kcl)
                if j['label'] in ('zero','replacement','replacement_native'):
                    execute(dest/'functional.json',c,j['T'],j['P']);execute(dest/'jvp.json',c,j['T'],j['P'])
            result.append(row);print(c['key'],j['arm'],j['label'],'qualified',row['qualified'],'Id',row.get('current_A_per_um'),'failure',row['failure'],flush=True)
        return result
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
    # Uniform keys preserve diagnostics even if a solver fails before writing state.
    keys=list(dict.fromkeys(k for row in rows for k in row));a.write_csv(OUT/'dc.csv',[{k:row.get(k,'') for k in keys} for row in rows])


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','preflight','run'));globals()[parser.parse_args().action]()
