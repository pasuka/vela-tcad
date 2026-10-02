"""Extend the frozen joint geometry experiment to four Vg=0.8 V controls."""
import argparse
import copy
import math
from pathlib import Path
import numpy as np
import build_simplemos_joint_stable_20260907 as b
import validate_simplemos_joint_geometry_20260907 as v
import check_simplemos_joint_geometry_20260907 as check

p=b.p; a=p.a; d=p.d
LOCAL=p.REPO/'build-release/simplemos_joint_eight_point_20260907'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/joint_eight_point_20260907'
v.b=b; v.LOCAL=LOCAL; v.OUT=OUT
v.AXES={'baseline':(0.,0.),'joint':(1.,1.)}
check.LOCAL=LOCAL; check.OUT=OUT


def prepare():
    for f in (b.OUT/'validation_evidence.json',p.prior.OUT/'validation_evidence.json'):a.verify(f)
    assert b.RUNNER.is_file()
    old_contract=a.read(b.OUT/'contract.json')
    files=[Path(__file__).resolve(),Path(v.__file__).resolve(),Path(check.__file__).resolve(),
           b.OUT/'validation_evidence.json',p.prior.OUT/'validation_evidence.json',b.RUNNER,b.HEADER,b.old.HEADER,
           p.prior.OUT/'vela_contract.json']
    cases=[]; predictions=[]
    factor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge']
    for old in a.read(p.prior.OUT/'vela_contract.json')['cases']:
        if old['vg']!=.8:continue
        c=copy.deepcopy(old); base=Path(c['base']); geo=d.matrix.spatial.m73.Geometry(c['device']); root=LOCAL/c['key']
        result=a.read(base/'result.json'); assert result['comparison_qualified'] and c['native_qualified']
        c['base_Id_A_per_um']=result['Id_A_per_um']; c['native_initial']=str(base/'initial.csv')
        c['ratios']=str(b.old.LOCAL/'ratios'/c['device'])
        # Stable build reuses the original, immutable geometry ratio tables.
        assert Path(c['ratios']).is_dir()
        cfg=a.read(base/'config.json'); c['jobs']=[]; c['probes']=[]
        eRat={int(x.split()[0]):float(x.split()[1]) for x in (Path(c['ratios'])/'edges.txt').read_text().splitlines()}
        nodeRat=np.array([float(x.split()[1]) for x in (Path(c['ratios'])/'nodes.txt').read_text().splitlines()])
        strict=p.prior.LOCAL/'vela'/c['key']/'strict'
        edges=a.rows(strict/'edges.csv'); oldR=d.ordered(strict/'residual.csv',geo.count)
        state=d.ordered(base/'state.csv',geo.count); adj=d.ordered(p.prior.LOCAL/'vela'/c['key']/'adjoint.csv',geo.count)
        lam=np.array([[float(x[k]) for x in adj] for k in ('lambda_poisson','lambda_electron','lambda_hole')]); lam[:,geo.contact_nodes]=0
        rhs=np.zeros((3,geo.count)); abssource=np.zeros_like(rhs)
        mesh=a.read(Path(cfg['mesh_file'])); drain=set(next(x['node_ids'] for x in mesh['contacts'] if x['name']=='drain')); direct=[]
        for e in edges:
            i,j=int(e['node0']),int(e['node1']); delta=eRat[int(e['edge_id'])]-1
            for block,car in [(1,'electron'),(2,'hole')]:
                f=float(e[car+'_flux'])*delta; rhs[block,i]+=f; rhs[block,j]-=f; abssource[block,i]+=abs(f); abssource[block,j]+=abs(f)
            direct.append(-d.fixed.Q*1e-6*(int(i in drain)-int(j in drain))*(float(e['electron_particle_line_flux_per_m_s'])-float(e['hole_particle_line_flux_per_m_s']))*delta)
        for i,(s,r) in enumerate(zip(state,oldR)):
            charge=float(s['electrons_m3'])-float(s['holes_m3'])-float(r['net_doping_m3'])
            rhs[0,i]=factor*d.fixed.Q*charge*geo.volumes['all_cell'][i]*(nodeRat[i]-1)
        rhs[:,geo.contact_nodes]=0
        for arm,(T,P) in v.AXES.items():
            source=rhs*T; feedback=-math.fsum(float(x*y) for x,y in zip(lam.flat,source.flat)); di=T*math.fsum(direct)
            predictions.append(dict(key=c['key'],arm=arm,unit_feedback_A_per_um=feedback,unit_direct_A_per_um=di,unit_prediction_A_per_um=feedback+di))
            path=root/'preflight'/arm/'expected.csv'
            a.write_csv(path,[dict(node_id=i,psi_source=source[0,i],electron_source=source[1,i],hole_source=source[2,i],electron_source_abs=abssource[1,i]*T,hole_source_abs=abssource[2,i]*T) for i in range(geo.count)]); files.append(path)
            for name,kind in [('functional','terminal_current_functional_probe'),('edges','sg_edge_flux_probe'),('jvp','joint_geometry_jvp')]:
                dest=root/'preflight'/arm; deck=copy.deepcopy(cfg); deck.pop('output_state_file'); deck['solver'].pop('local_update_diagnostics')
                deck.update(simulation_type=kind,state_file=str(base/'state.csv'),contact='drain')
                deck['residual_output_csv' if name=='functional' else 'output_csv']=str(dest/('residual.csv' if name=='functional' else name+'.csv'))
                a.write(dest/(name+'.json'),deck); files.append(dest/(name+'.json')); c['probes'].append(dict(arm=arm,T=T,P=P,config=str(dest/(name+'.json'))))
            labels=[('zero',0.,'vela')] if arm=='baseline' else [(n,s,'vela') for n,s in [('plus_full',.001),('minus_full',-.001),('plus_half',.0005),('minus_half',-.0005),('replacement',1.)]]+[('replacement_native',1.,'native')]
            for label,scale,initial in labels:
                dest=root/arm/label; deck=copy.deepcopy(cfg)
                deck.update(state_file=str(base/('initial.csv' if initial=='native' else 'state.csv')),output_state_file=str(dest/'state.csv'))
                deck['solver']['local_update_diagnostics']['csv_file']=str(dest/'updates.csv'); a.write(dest/'config.json',deck)
                q=copy.deepcopy(deck); q.pop('output_state_file'); q['solver'].pop('local_update_diagnostics')
                q.update(simulation_type='newton_carrier_term_probe',state_file=str(dest/'state.csv'),output_csv=str(dest/'all_row.csv'),carrier_term_probe={'solved_equation_terms':True})
                q['solver']['carrier_row_convergence']['mode']='report'; q['solver']['global_continuity_closure']=dict(mode='enforce',tolerance=1e-6,source_floor=1e-10); a.write(dest/'all_row.json',q)
                j=copy.deepcopy(q); j.pop('carrier_term_probe'); j.update(simulation_type='joint_geometry_jvp',output_csv=str(dest/'jvp.csv')); a.write(dest/'jvp.json',j)
                fun=copy.deepcopy(j); fun.update(simulation_type='terminal_current_functional_probe',contact='drain',residual_output_csv=str(dest/'residual.csv')); fun.pop('output_csv'); a.write(dest/'functional.json',fun)
                c['jobs'].append(dict(arm=arm,label=label,scale=scale,T=T*scale,P=P*scale,initialization=initial,config=str(dest/'config.json')))
                files += [dest/n for n in ('config.json','all_row.json','jvp.json','functional.json')]
        cases.append(c)
        files += [base/n for n in ('config.json','state.csv','initial.csv','result.json')]
        files += [strict/'edges.csv',strict/'residual.csv',p.prior.LOCAL/'vela'/c['key']/'adjoint.csv',Path(cfg['mesh_file'])]
        files += [Path(c['ratios'])/n for n in ('edges.txt','nodes.txt')]
    assert len(cases)==4 and sum(len(c['jobs']) for c in cases)==28
    a.write_csv(OUT/'predictions.csv',predictions)
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',cases=cases,new_DC=28,axes=old_contract['axes'],gates=old_contract['gates'],
        runner=str(b.RUNNER),prior_evidence=str(b.OUT/'validation_evidence.json'),
        response='Vg=0.8 V, n19/n23 x Vd 0.05/1 V: common zero, joint +/-0.001 and +/-0.0005, full joint replacement from Vela and native-coherent initialization.',
        promotion=old_contract['promotion'],jvp=old_contract['jvp'],source_preflight=old_contract['source_preflight'],
        near_threshold='Rank every active carrier row. At joint full replacement states with max row ratio >=5e-7, run one unchanged-gate restart from the accepted state. Retain original records; compare fields and Id with the frozen initialization-invariance gates. Include existing Vg=1 V states.',
        field_comparison='Eight baseline/joint endpoints, all 907 non-contact Si nodes, native Si volume weighted RMS and maximum. Physical quasi-Fermi potential is reference plus increment. Rank residual field maxima with coordinates, densities, current/source row terms and native densities; no density exclusion.',
        nonlinear_acceptance_changes=False,production_changes=False,full_native_mobility_replacement=False,full_native_Poisson_K_replacement=False))
    files += [OUT/'predictions.csv',OUT/'contract.json']; d.matrix.freeze(OUT/'freeze.json',files)
    print('Frozen four Vg=0.8 V cases, 28 DC jobs; reuse stable runner and identical joint coefficients.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('action',choices=('prepare','preflight','check','run'))
    action=parser.parse_args().action
    {'prepare':prepare,'preflight':v.preflight,'check':check.main,'run':v.run}[action]()
