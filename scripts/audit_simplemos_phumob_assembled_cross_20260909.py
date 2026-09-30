"""Frozen-state legacy PhuMob cross columns against a 100-digit chain rule.

SRH is disabled to isolate mobility. This is not a self-consistent PhuMob
qualification and does not enable the unsupported element_box mobility path.
"""
import argparse,copy,math
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal as D,localcontext
from pathlib import Path
import audit_simplemos_phumob_cross_derivatives_20260908 as h

p=h.p;a,d=p.a,p.d;q=p.qualified
LOCAL=p.LOCAL/'assembled_cross_20260909';OUT=p.OUT/'assembled_cross_20260909'
VT=1.380649e-23*300/1.602176634e-19
STEPS=(1e-5,3e-6)
NODES=h.NODES


def prepare():
    a.verify(q.OUT/'final_evidence.json')
    frozen=a.read(q.OUT/'validation_freeze.json')['input_hashes']
    runner=q.run.RUNNER
    assert a.sha(runner)==frozen[runner.resolve().relative_to(p.REPO).as_posix()]
    jobs=[];files=[Path(__file__).resolve(),Path(h.__file__),runner,q.OUT/'final_evidence.json']
    for row in a.rows(q.OUT/'comparison.csv'):
        if row['model']!='old_slotboom' or int(row['index'])!=40:continue
        base=q.LOCAL/'dc/old_slotboom'/row['case']/'vela/vg_040/attempt_0'
        cfg=a.read(base/'acceptance_edges.json');mesh=a.read(Path(cfg['mesh_file']))
        contacts={i for c in mesh['contacts'] for i in c['node_ids']}
        assert not contacts.intersection(NODES)
        edges=a.rows(base/'acceptance_edges.csv')
        neighbors=set(NODES)
        for e in edges:
            ij={int(e['node0']),int(e['node1'])}
            if ij.intersection(NODES):neighbors.update(ij)
        neighbors-=contacts
        dest=LOCAL/row['case'];cfg['solver']['mobility'].update(model='phumob',edge_averaging='legacy')
        cfg['solver']['recombination']=['none']
        cfg['solver']['srh_doping_dependence']['enabled']=False
        cfg['output_csv']=str(dest/'edges.csv');a.write(dest/'edges.json',cfg)
        jvp=copy.deepcopy(cfg);jvp.update(simulation_type='newton_jvp_probe',output_csv=str(dest/'jvp.csv'),row_output_csv=str(dest/'rows.csv'))
        jvp['directions']=[dict(name=f'{mode}_{node}_{step:.0e}',mode=mode,node_ids=[node],amplitude_V=step,exclude_contacts=True)
                           for node in NODES for mode in ('phin','phip') for step in STEPS]
        jvp['sample_rows']=[dict(block=block,node_id=node) for block in ('phin','phip') for node in sorted(neighbors)]
        a.write(dest/'jvp.json',jvp)
        jobs.append(dict(case=row['case'],device=row['device'],vd=row['vd'],dir=str(dest),contacts=sorted(contacts)))
        files += [dest/'edges.json',dest/'jvp.json',base/'acceptance_edges.json',base/'acceptance_edges.csv']+[Path(cfg[k]) for k in ('mesh_file','state_file','materials_file','node_doping_file')]
    assert len(jobs)==4
    a.write(OUT/'contract.json',dict(jobs=jobs,nodes=NODES,steps_V=STEPS,temperature_K=300,thermal_voltage_V=VT,
        models='PhuMob arsenic, OldSlotboom, legacy endpoint-state mobility; SRH/Enormal/HFS/avalanche off; previous qualified Masetti states held fixed.',
        geometry='Retain qualified element_box transport geometry, Delaunay transfer, cell-material permittivity, signed Si Poisson volume.',
        reference='100-digit analytic carrier chain rule at each edge endpoint-mean state, multiplied by independently exported base flux/mobility and scattered conservatively.',
        row_units='Raw assembler residual per physical volt; no row-scale weighting.',
        gates=dict(scalar_mu_relative=1e-12,base_flux_closure_over_abs_sum=1e-10,cross_derivative_relative=1e-5),
        interpretation='No absolute floor hides nonzero weak references. Report absolute entries and row sums as well as relative failures. This audit does not validate a new cell mobility stencil or self-consistent PhuMob states.'))
    files += [OUT/'contract.json',p.REPO/'src/equation/CoupledDDAssembler.cpp',p.REPO/'src/physics/MobilityModel.cpp',p.REPO/'include/vela/equation/AssemblerUtils.h',p.REPO/'include/vela/core/PhysicalConstants.h']
    d.matrix.freeze(OUT/'freeze.json',files)


def execute():
    a.verify(OUT/'freeze.json')
    def one(job):
        dest=Path(job['dir']);results=[]
        for name in ('edges','jvp'):
            status=q.run.V.execute(dest/(name+'.json'),q.run.RUNNER,q.run.V.environment())
            results.append(dict(case=job['case'],probe=name,**status))
            assert status['exit_code']==0,status
        return results
    with ThreadPoolExecutor(max_workers=2) as pool:
        results=[r for block in pool.map(one,a.read(OUT/'contract.json')['jobs']) for r in block]
    a.write(OUT/'execution.json',dict(results=results))
    d.matrix.freeze(OUT/'raw_evidence.json',[OUT/'freeze.json',OUT/'execution.json']+[x for x in LOCAL.rglob('*') if x.is_file()])


def analyze():
    a.verify(OUT/'raw_evidence.json')
    result=[];reference_edges=[];base_checks=[];scalar_checks=[]
    with localcontext() as ctx:
        ctx.prec=100
        for job in a.read(OUT/'contract.json')['jobs']:
            dest=Path(job['dir']);cfg=a.read(dest/'edges.json')
            doping={int(x['node_id']):x for x in a.rows(Path(cfg['node_doping_file']))}
            ref=defaultdict(lambda:D(0));fluxes=defaultdict(list)
            for e in a.rows(dest/'edges.csv'):
                i,j=int(e['node0']),int(e['node1'])
                for car,block in [('electron','phin'),('hole','phip')]:
                    flux=D(e[car+'_flux']);fluxes[block,i].append(flux);fluxes[block,j].append(-flux)
                if i not in NODES and j not in NODES:continue
                nd=sum(D(doping[n]['donors_cm3']) for n in (i,j))/2
                na=sum(D(doping[n]['acceptors_cm3']) for n in (i,j))/2
                densities={car:[D(e[f'{car}_density{k}_m3'])*D('1e-6') for k in range(2)] for car in ('electron','hole')}
                state=(nd,na,sum(densities['electron'])/2,sum(densities['hole'])/2)
                for ci,(car,block,mode,pop,sign) in enumerate([('electron','phin','phip','hole',1),('hole','phip','phin','electron',-1)]):
                    mu,deriv,clamped=h.hp(state,ci)
                    actual_mu=D(e[car+'_mobility_m2_V_s'])*D('1e4')
                    scalar_checks.append(dict(case=job['case'],edge=e['edge_id'],carrier=car,relative=float(abs(mu/actual_mu-1))))
                    derivative=deriv[1 if pop=='hole' else 0]
                    flux=D(e[car+'_flux'])
                    for k,node in enumerate((i,j)):
                        if node not in NODES:continue
                        value=flux/mu*derivative*densities[pop][k]/sum(densities[pop])*D(sign)/D(str(VT))
                        ref[mode,node,block,i]+=value;ref[mode,node,block,j]-=value
                        reference_edges.append(dict(case=job['case'],edge=e['edge_id'],input_mode=mode,input_node=node,output_block=block,flux=str(flux),mu=str(mu),dmu_dlog_population=str(derivative),dflux_dpotential=str(value),elasticity=str(derivative/mu),G_clamped=clamped))
            for row in a.rows(dest/'rows.csv'):
                mode,node_s,_=row['direction'].split('_');node=int(node_s);block=row['row_block'];rnode=int(row['row_node'])
                if row['direction_mode']==block:continue
                vals=fluxes[block,rnode];expected=sum(vals);den=sum(abs(x) for x in vals)
                base_checks.append(dict(case=job['case'],row=block,node=rnode,relative=float(abs(expected-D(row['base_residual']))/den) if den else 0))
                target=ref[mode,node,block,rnode];actual=D(row['analytic_derivative']);fd=D(row['finite_difference_derivative'])
                if not target:
                    assert actual==0 and fd==0,(row,target)
                    continue
                error=abs(actual/target-1);fderror=abs(fd/target-1)
                result.append(dict(case=job['case'],input_mode=mode,input_node=node,output_block=block,output_node=rnode,step_V=row['amplitude_V'],reference=str(target),production=str(actual),residual_fd=str(fd),absolute_error=str(abs(actual-target)),relative_error=float(error),fd_relative_error=float(fderror),production_zero=actual==0,fd_zero=fd==0,qualified=error<=D('1e-5')))
    a.write_csv(OUT/'cross_entries.csv',result);a.write_csv(OUT/'reference_edges.csv',reference_edges)
    a.write_csv(OUT/'base_closure.csv',base_checks);a.write_csv(OUT/'scalar_values.csv',scalar_checks)
    closure=max(x['relative'] for x in base_checks);muerror=max(x['relative'] for x in scalar_checks)
    preflight=closure<=1e-10 and muerror<=1e-12
    summary=dict(cases=4,nonzero_entry_checks=len(result),base_closure_max_relative=closure,scalar_mu_max_relative=muerror,reference_units_and_state_qualified=preflight,
        derivative_failures=sum(not x['qualified'] for x in result),production_zero_count=sum(x['production_zero'] for x in result),residual_fd_zero_count=sum(x['fd_zero'] for x in result),
        absolute_reference_min=min(abs(float(x['reference'])) for x in result),absolute_reference_max=max(abs(float(x['reference'])) for x in result),
        max_absolute_derivative_error=max(abs(float(x['absolute_error'])) for x in result),production_changed=False)
    a.write(OUT/'summary.json',summary)
    d.matrix.freeze(OUT/'evidence.json',[OUT/'raw_evidence.json']+[OUT/x for x in ('cross_entries.csv','reference_edges.csv','base_closure.csv','scalar_values.csv','summary.json')])
    print(summary,flush=True)
    assert preflight,'Reference state/units preflight failed; derivative gate is not interpretable.'


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','execute','analyze'))
    globals()[parser.parse_args().action]()
