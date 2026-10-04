"""Coefficient-path checks and local native/qualified-state flux ledgers."""
import argparse
import ast
from decimal import Decimal
import math
from pathlib import Path
import re
import subprocess
import numpy as np
import prepare_simplemos_production_consistency_20260907 as run
import audit_simplemos_masetti_box_mobility_20260907 as native

a=run.a; d=run.d; p=run.p; OUT=run.OUT; LOCAL=run.LOCAL; prior=run.prior
REPORT=p.REPO/'docs/validation/simplemos_production_consistency_and_local_fields_2026-09-07.md'


def analyze():
    for f in ('freeze.json','probe_evidence.json','representation_evidence.json'):a.verify(OUT/f)
    contract=a.read(OUT/'contract.json'); ports=[]; paths=[]; geometry=[]; sources=[]; localedges=[]; localnodes=[]; nativeports=[]; poisson=[]
    sourcefactor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge']
    cache={}; files=[Path(__file__).resolve()]
    for device in ('n19','n23'):
        geo,xy,cells,vol,K,parts,info,mods=native.geom.geometry(device)
        _,weights,checks=native.geometry(device); assert all(c['qualified'] for c in checks)
        cache[device]=(geo,vol,K,parts,weights)
        for name,other in [('legacy_all',geo.volumes['all_cell']),('builtin_transport_barycentric',geo.volumes['barycentric_si']),('signed_Si',geo.volumes['signed_si'])]:
            mask=vol['Si']>0; rel=np.abs(other[mask]/vol['Si'][mask]-1)
            geometry.append(dict(device=device,quantity='node_volume',variant=name,count=int(sum(mask)),different_at_relative_1e_10=int(sum(rel>1e-10)),max_relative=float(max(rel)),l1_relative=float(sum(abs(other[mask]-vol['Si'][mask]))/sum(vol['Si'][mask]))))
    for c in contract['cases']:
        geo,vol,K,parts,weights=cache[c['device']]; _,mask=prior.support(c); ids=np.flatnonzero(mask); root=LOCAL/c['key']
        meta=dict(key=c['key'],device=c['device'],vg=c['vg'],vd=c['vd'])
        modes=('legacy','builtin','joint','native_state_joint','native_referenced_joint')
        edges={mode:a.rows(root/mode/'edges.csv') for mode in modes}
        terms={mode:d.ordered(root/mode/'terms.csv',geo.count) for mode in modes}
        residual={mode:d.ordered(root/mode/'residual.csv',geo.count) for mode in modes}
        state=d.ordered(Path(c['root'])/'joint/replacement/state.csv',geo.count)
        mapped=d.ordered(root/'native_referenced_joint/state.csv',geo.count)
        mesh=a.read(Path(a.read(Path(c['base'])/'config.json')['mesh_file']))
        drain=set(next(x['node_ids'] for x in mesh['contacts'] if x['name']=='drain'))
        for mode in modes:
            s=a.read(root/mode/'functional.status.json'); error=abs(s['current_A_per_um']/s['contact_current_extractor_A_per_um']-1)
            edgeId=-d.fixed.Q*1e-6*math.fsum((int(int(e['node0']) in drain)-int(int(e['node1']) in drain))*(float(e['electron_particle_line_flux_per_m_s'])-float(e['hole_particle_line_flux_per_m_s'])) for e in edges[mode])
            ports.append(dict(**meta,mode=mode,functional_A_per_um=s['current_A_per_um'],extractor_A_per_um=s['contact_current_extractor_A_per_um'],
                SG_edge_sum_A_per_um=edgeId,functional_extractor_relative=error,edge_functional_relative=abs(edgeId/s['current_A_per_um']-1),qualified=error<=1e-8))
        for mode in ('builtin','joint'):
            mobility={int(e['edge_id']):e for e in a.rows(root/mode/'mobility.csv')}; different=0; maxdiff=0.; murerr=0.; active=0; drain_changed=0
            for e,old in zip(edges[mode],edges['legacy']):
                eid=int(e['edge_id']); assert eid==int(old['edge_id'])
                if float(e['electron_mobility_m2_V_s'])<=0:continue
                active+=1; g=float(e['couple_m']); rawg=float(mobility[eid]['couple_m']); rel=abs(g-rawg)/max(g,rawg,1e-300)
                different+=rel>1e-10; maxdiff=max(maxdiff,rel)
                if (int(e['node0']) in drain)!=(int(e['node1']) in drain):drain_changed+=abs(g-float(old['couple_m']))>1e-10*max(g,float(old['couple_m']))
                for car in ('electron','hole'):
                    assert e[car+'_mobility_m2_V_s']==old[car+'_mobility_m2_V_s']
                    murerr=max(murerr,abs(float(e[car+'_mobility_m2_V_s'])/float(mobility[eid][car+'_final_mobility_m2_V_s'])-1))
            paths.append(dict(**meta,mode=mode,active_edges=active,mobility_probe_coupling_different_edges=different,
                coupling_max_relative=maxdiff,mobility_max_relative=murerr,drain_support_changed_edges=drain_changed))
        for mode in ('legacy','builtin','joint'):
            diffs=[]; zero_native_nonzero=0
            for e in edges[mode]:
                i,j=int(e['node0']),int(e['node1']); gn=math.fsum(x['coefficient'] for x in parts[(i,j)] if x['material']=='Si'); gv=float(e['couple_m'])/float(e['length_m'])
                if gn>0:diffs.append(abs(gv/gn-1))
                elif float(e['electron_mobility_m2_V_s'])>0 and gv>0:zero_native_nonzero+=1
            geometry.append(dict(device=c['device'],quantity='edge_geometry',variant=mode,count=len(diffs),different_at_relative_1e_10=sum(x>1e-10 for x in diffs),
                max_relative=max(diffs),l1_relative='',native_zero_active_nonzero=zero_native_nonzero,key=c['key']))
        # Fixed-state Poisson source differences must use all three charge terms.
        oldR=np.array([float(r['psi_residual']) for r in residual['legacy']])
        oldsource=d.ordered(p.prior.LOCAL/'vela'/c['key']/'strict/residual.csv',geo.count)
        charge=np.array([float(s['electrons_m3'])-float(s['holes_m3'])-float(r['net_doping_m3']) for s,r in zip(state,oldsource)])
        for mode,target in [('builtin',geo.volumes['barycentric_si']),('joint',vol['Si'])]:
            expected=d.fixed.Q*sourcefactor*charge*(target-geo.volumes['all_cell'])
            actual=np.array([float(r['psi_residual']) for r in residual[mode]])-oldR
            err=float(np.linalg.norm((actual-expected)[mask])/max(np.linalg.norm(expected[mask]),1e-300))
            sources.append(dict(**meta,mode=mode,check='all_three_Poisson_volume_terms',relative_or_bound_error=err,qualified=err<=1e-8))
        # Independent oriented edge sums reproduce actual continuity row terms.
        for mode in modes:
            net=np.zeros((2,geo.count)); absolute=np.zeros_like(net)
            for e in edges[mode]:
                i,j=int(e['node0']),int(e['node1'])
                for k,car in enumerate(('electron','hole')):
                    f=float(e[car+'_flux']); net[k,i]+=f; net[k,j]-=f; absolute[k,i]+=abs(f); absolute[k,j]+=abs(f)
            worst=0.
            for k,car in enumerate(('electron','hole')):
                for i in ids:
                    err=abs(net[k,i]-float(terms[mode][i][car+'_flux']))
                    bound=max(128*np.finfo(float).eps*absolute[k,i],1e-300); worst=max(worst,err/bound)
            sources.append(dict(**meta,mode=mode,check='oriented_edge_sum_to_carrier_terms',relative_or_bound_error=worst,qualified=worst<=1))
        # Replay native qf SG with calibrated box-weighted cell mobility.
        case=Path(c['base']).parent.name; index=int(Path(c['base']).name.split('_')[1])
        assert c['key']==f'{case}_vg{index:03d}' and abs(c['vg']-index*.05)<1e-12
        src=p.prior.prev.LOCAL/'native_exports/masetti'/case/f"vg_{index:03d}"; fields=src/'fields'
        names={'psi':'ElectrostaticPotential','n':'eDensity','p':'hDensity','fn':'eQuasiFermiPotential','fp':'hQuasiFermiPotential','SRH':'srhRecombination'}
        raw={key:native.scalar(fields/(name+'_region0.csv')) for key,name in names.items()}
        files += [fields/(name+'_region0.csv') for name in names.values()]
        nativeF={}; nativeNet=np.zeros((2,geo.count)); nativeAbs=np.zeros_like(nativeNet); terminal=[]
        ratio=next(float(e['electron_particle_line_flux_per_m_s'])/float(e['electron_flux']) for e in edges['joint'] if float(e['electron_flux'])!=0)
        for en,em in zip(edges['joint'],edges['native_referenced_joint']):
            i,j=int(en['node0']),int(en['node1']); eid=int(en['edge_id']); w=weights[(i,j)]
            if i not in raw['n'] or j not in raw['n']:fn,fp=0.,0.
            else:
                fn,fp=native.sg_flux(raw['psi'][i],raw['psi'][j],raw['n'][i],raw['n'][j],raw['p'][i],raw['p'][j],raw['fn'][i],raw['fn'][j],raw['fp'][i],raw['fp'][j],'qf')
                fn=float(fn)*w['e']*100; fp=float(fp)*w['h']*100
            nativeF[eid]=(fn,fp); terminal.append(-native.Q*1e-6*(int(i in drain)-int(j in drain))*(fn-fp))
            for k,f in enumerate((fn,fp)):nativeNet[k,i]+=f; nativeNet[k,j]-=f; nativeAbs[k,i]+=abs(f); nativeAbs[k,j]+=abs(f)
            if not ({i,j}&set(contract['local_nodes'])) or i not in raw['n'] or j not in raw['n']:continue
            gn=math.fsum(x['coefficient'] for x in parts[(i,j)] if x['material']=='Si')
            if gn<=0:continue
            for k,car,field,rfield in [(0,'electron','phin','fn'),(1,'hole','phip','fp')]:
                Fv=float(en[car+'_particle_line_flux_per_m_s']); Fm=float(em[car+'_particle_line_flux_per_m_s']); Fn=(fn,fp)[k]
                dv=float(prior.physical(state[j],field)-prior.physical(state[i],field)); dn=float(Decimal(str(raw[rfield][j]))-Decimal(str(raw[rfield][i])))
                direction=1 if k==0 else -1; dv*=direction; dn*=direction
                Gv=Fv/dv if dv else math.nan; Gn=Fn/dn if dn else math.nan
                muV=float(en[car+'_mobility_m2_V_s'])*1e4; muN=w['e' if k==0 else 'h']/gn
                mobility_change=Fm*(muN/muV-1); convention=Fn-Fm-mobility_change
                dropPart=.5*(Gv+Gn)*(dn-dv) if dv and dn else math.nan
                conductancePart=.5*(dv+dn)*(Gn-Gv) if dv and dn else math.nan
                localedges.append(dict(**meta,edge_id=eid,node0=i,node1=j,carrier=car,native_geometry=gn,
                    Vela_edge_mu_cm2_Vs=muV,native_effective_mu_cm2_Vs=muN,mu_relative=muV/muN-1,
                    Vela_oriented_qf_drop_V=dv,native_oriented_qf_drop_V=dn,Vela_G_per_m_s_V=Gv,native_G_per_m_s_V=Gn,
                    Vela_flux_per_m_s=Fv,native_flux_per_m_s=Fn,mapped_native_flux_per_m_s=Fm,
                    fixed_native_state_mobility_change_per_m_s=mobility_change,fixed_native_state_remaining_convention_per_m_s=convention,
                    qf_drop_algebraic_part_per_m_s=dropPart,SG_conductance_algebraic_part_per_m_s=conductancePart,
                    algebraic_closure_relative=abs(dropPart+conductancePart-(Fn-Fv))/max(abs(Fn),abs(Fv),1e-300) if dv and dn else math.nan))
        nativeId=math.fsum(terminal); nativeerr=abs(nativeId/c['native_Id_A_per_um']-1)
        nativeports.append(dict(**meta,native_SG_Id_A_per_um=nativeId,target_Id_A_per_um=c['native_Id_A_per_um'],relative_error=nativeerr,qualified=nativeerr<=1e-6))
        # Local rows include physical native densities and the nodal SRH plot.
        for i in contract['local_nodes']:
            if not mask[i]:continue
            t=terms['joint'][i]; tm=terms['native_referenced_joint'][i]; s=state[i]
            record=dict(**meta,node_id=i,x_um=geo.coords[i][0],y_um=geo.coords[i][1],native_Si_volume_m2=vol['Si'][i],Vela_SRH_volume_m2=geo.volumes['all_cell'][i])
            for field,rfield in [('psi','psi'),('phin','fn'),('phip','fp')]:record[field+'_delta_V']=float(prior.physical(s,field)-Decimal(str(raw[rfield][i])))
            for k,car,field in [(0,'electron','electrons_m3'),(1,'hole','holes_m3')]:
                scale=max(float(t[car+'_flux_abs_sum']),abs(float(t[car+'_recombination'])),abs(float(t[car+'_impact'])))
                record.update({car+'_row_ratio':abs(float(t[car+'_residual']))/scale,car+'_native_net_edge_flux_per_m_s':nativeNet[k,i],
                    car+'_native_absolute_edge_flux_per_m_s':nativeAbs[k,i],car+'_Vela_net_edge_flux_per_m_s':float(t[car+'_flux'])*ratio,
                    car+'_Vela_SRH_integral_per_m_s':float(t[car+'_recombination'])*ratio,
                    car+'_mapped_residual_per_m_s':float(tm[car+'_residual'])*ratio})
                rawdensity=raw['n' if k==0 else 'p'][i]*1e6
                record.update({car+'_Vela_density_m3':float(s[field]),car+'_native_density_m3':rawdensity,car+'_density_relative':float(s[field])/rawdensity-1})
            record['native_plot_SRH_cm3_s']=raw['SRH'][i]
            record['Vela_SRH_cm3_s']=float(t['electron_recombination'])*ratio/geo.volumes['all_cell'][i]/1e6
            # Native nodal plot * box volume is a diagnostic quadrature only.
            record['native_plot_SRH_times_box_per_m_s']=raw['SRH'][i]*1e6*vol['Si'][i]
            localnodes.append(record)
        psi,n,h,spread=d.matrix.spatial.m73.sentaurus_state(src,geo); assert spread<=1e-12
        nd=native.scalar(fields/'DonorConcentration_region0.csv'); na=native.scalar(fields/'AcceptorConcentration_region0.csv')
        doping=np.array([(nd.get(i,0)-na.get(i,0))*1e6 for i in range(geo.count)])
        mq=np.array([float(s['electrons_m3'])-float(s['holes_m3']) for s in mapped])
        replay=K@psi+d.fixed.Q*(n-h-doping)*vol['Si']; dielectric=(geo.matrices['legacy']-K)@psi
        convention=d.fixed.Q*(mq-(n-h))*vol['Si']
        actual=np.array([float(r['psi_residual']) for r in residual['native_referenced_joint']])/sourcefactor
        error=float(np.linalg.norm((actual-replay-dielectric-convention)[mask])/max(np.linalg.norm(actual[mask]),1e-300))
        sources.append(dict(**meta,mode='native_referenced_joint',check='Poisson_native_geometry_density_replay',relative_or_bound_error=error,qualified=error<=1e-8))
        for i in contract['local_nodes']:
            if mask[i]:poisson.append(dict(**meta,node_id=i,actual_joint_Poisson_C_per_m=actual[i],native_geometry_charge_replay_C_per_m=replay[i],
                dielectric_geometry_difference_C_per_m=dielectric[i],density_convention_difference_C_per_m=convention[i],replay_remainder_C_per_m=actual[i]-replay[i]-dielectric[i]-convention[i]))
        print('Audited paths, native SG calibration and local ledger',c['key'],flush=True)
    outputs=[('port_paths',ports),('diagnostic_paths',paths),('geometry_comparison',geometry),('operator_replay_checks',sources),('local_edge_ledger',localedges),('local_node_ledger',localnodes),('native_SG_port_calibration',nativeports),('local_Poisson_ledger',poisson)]
    for name,rows in outputs:
        keys=list(dict.fromkeys(k for r in rows for k in r)); payload=[{k:r.get(k,'') for k in keys} for r in rows]
        dest=OUT/(name+'.csv')
        if dest.exists():assert a.rows(dest)==[{k:str(v) for k,v in r.items()} for r in payload],name
        else:a.write_csv(dest,payload)
    metrics=dict(read_only_probes=136,new_DC=0,production_changes=False,qualified_state_port_checks=sum(r['qualified'] for r in ports if r['mode'] in ('legacy','builtin','joint')),
        qualified_state_port_count=24,unreferenced_native_port_failures=sum(not r['qualified'] for r in ports if r['mode']=='native_state_joint'),
        referenced_native_port_qualified=sum(r['qualified'] for r in ports if r['mode']=='native_referenced_joint'),
        native_SG_port_qualified=sum(r['qualified'] for r in nativeports),native_SG_port_max_relative=max(r['relative_error'] for r in nativeports),
        operator_replay_qualified=sum(r['qualified'] for r in sources),operator_replay_count=len(sources),
        max_algebraic_flux_decomposition_closure=max(r['algebraic_closure_relative'] for r in localedges if math.isfinite(r['algebraic_closure_relative'])),
        production_joint_scheme_equivalent=False,M82_released=False,M83_released=False)
    metrics={k:v.item() if isinstance(v,np.generic) else v for k,v in metrics.items()}
    files += [p.REPO/x for x in ('src/post/StoredCharge.cpp','src/post/TerminalCharge.cpp')]
    files += [x for root in (LOCAL,OUT) for x in root.rglob('*') if x.is_file()]
    a.write(OUT/'analysis_evidence.json',dict(status='completed_production_path_and_local_physics_audit',metrics=metrics,input_hashes={a.rel(f):a.sha(f) for f in files}))
    print(metrics,flush=True)


def seal():
    for f in ('freeze.json','probe_evidence.json','representation_evidence.json','analysis_evidence.json'):a.verify(OUT/f)
    a.verify(prior.OUT/'validation_evidence.json')
    scripts=[Path(__file__).resolve(),Path(run.__file__).resolve(),p.REPO/'scripts/probe_simplemos_state_representation_20260907.py']
    for script in scripts:ast.parse(script.read_text(encoding='utf-8'))
    for link in re.findall(r'\]\(([^)]+)\)',REPORT.read_text(encoding='utf-8')):
        if '://' not in link:
            path=(REPORT.parent/link.split('#')[0]).resolve()
            if path.name!='validation_evidence.json':assert path.is_file(),path
    diff=subprocess.run(['git','-c','core.fsmonitor=false','diff','--numstat','--ignore-space-at-eol'],cwd=p.REPO,capture_output=True,text=True)
    assert diff.returncode==0 and diff.stdout.strip()=='127\t1\tsrc/tools/vela_example_runner.cpp'
    files=scripts+[REPORT,prior.OUT/'validation_evidence.json']+[x for root in (LOCAL,OUT) for x in root.rglob('*') if x.is_file()]
    a.write(OUT/'validation_evidence.json',dict(status='completed_bounded_audit',date='2026-09-07',metrics=a.read(OUT/'analysis_evidence.json')['metrics'],
        input_hashes={a.rel(f):a.sha(f) for f in sorted(set(files))},production_changes=False,acceptance_changes=False,
        limitations=['No production implementation or self-consistent candidate run.','Native local flux is reconstructed and terminal calibrated, not an exported native edge unknown.',
            'Local algebraic flux attribution is not a causal current-response calibration.','Native Poisson/SRH replay is not exported native row residual.','Current production Jacobian has not inherited isolated stable SG or linear-refinement changes.']))
    a.verify(OUT/'validation_evidence.json'); print('Audit sealed; original eight-point evidence unchanged.',flush=True)


if __name__=='__main__':
    q=argparse.ArgumentParser(); q.add_argument('action',choices=('analyze','seal')); globals()[q.parse_args().action]()
