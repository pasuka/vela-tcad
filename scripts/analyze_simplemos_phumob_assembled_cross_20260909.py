"""Analysis amendment: explicitly exclude transport-inactive zero-mobility edges.

The original inputs, raw results and failed preflight remain frozen unchanged.
No numerical gate changes; zero mobility is accepted here only with zero flux.
"""
from audit_simplemos_phumob_assembled_cross_20260909 import *
import audit_simplemos_phumob_assembled_cross_20260909 as original
OUT=original.OUT/'active_transport_analysis'

def analyze():
    a.verify(original.OUT/'raw_evidence.json')
    a.verify(OUT/'analysis_freeze.json')
    result=[];reference_edges=[];base_checks=[];scalar_checks=[]
    with localcontext() as ctx:
        ctx.prec=100
        for job in a.read(original.OUT/'contract.json')['jobs']:
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
                    if actual_mu==0:
                        assert D(e[car+'_flux'])==0, 'Zero mobility edge carries flux.'
                        continue
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
    d.matrix.freeze(OUT/'evidence.json',[original.OUT/'raw_evidence.json',OUT/'analysis_freeze.json']+[OUT/x for x in ('cross_entries.csv','reference_edges.csv','base_closure.csv','scalar_values.csv','summary.json')])
    print(summary,flush=True)
    assert preflight,'Reference state/units preflight failed; derivative gate is not interpretable.'


if __name__=='__main__':
    d.matrix.freeze(OUT/'analysis_freeze.json',[Path(__file__).resolve(),original.OUT/'raw_evidence.json'])
    analyze()
