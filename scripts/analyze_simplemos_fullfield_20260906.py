"""Compare all 84 target fields with explicit per-state strict qualification."""
from decimal import Decimal
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import run_simplemos_fullfield_vela_20260906 as b

a=b.a;d=b.d;LOCAL=b.LOCAL;OUT=b.OUT
REGIONS=('all_si','channel','gate_interface','source','drain','substrate')
FIELDS=('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential','eDensity','hDensity')


def main():
    for name in ('native_freeze','export_contract','vela_freeze','recovery_freeze'):
        a.verify(OUT/(name+'.json'))
    attempts=[];files=[Path(__file__).resolve(),OUT/'native_points.csv']
    for name,folder in (('vela_ascending','vela'),('vela_recovery','native_seed'),('vela_reclosure','reclosure')):
        path=OUT/(name+'.csv')
        if not path.exists():continue
        files.append(path)
        for row in a.rows(path):
            row['folder']=folder;attempts.append(row)
    native_points={(r['case'],round(float(r['vg']),10)):r for r in a.rows(OUT/'native_points.csv')}
    assert len(native_points)==84
    points=[];fields=[];sources=[];nodes=[]
    for c in a.read(OUT/'vela_contract.json')['cases']:
        key=c['case'];geo=d.matrix.spatial.m73.Geometry(c['device']);count=geo.count
        masks,_,xy=d.matrix.spatial.old.m78.supports(c['device'],geo,.05)
        edgepath=d.LOCAL/key/'strict/edges.csv';files.append(edgepath)
        edge=max(a.rows(edgepath),key=lambda r:abs(float(r['electron_flux'])))
        factor=float(edge['electron_particle_line_flux_per_m_s'])*d.fixed.Q*1e-6/float(edge['electron_flux'])
        for index in range(21):
            vg=round(.05*index,10);available=[r for r in attempts if r['case']==key and abs(float(r['vg'])-vg)<1e-10]
            assert available,(key,vg)
            qualified=[r for r in available if r['qualified']=='True']
            chosen=qualified[0] if qualified else available[-1]
            dest=LOCAL/chosen['folder']/key/f'vg_{index:03d}';native=native_points[key,vg]
            statepath=dest/'state.csv';termpath=dest/'carrier.csv';cfgpath=dest/'config.json'
            files += [statepath,termpath,cfgpath,dest/'config.status.json',dest/'acceptance.status.json']
            state=d.ordered(statepath,count);terms=d.ordered(termpath,count);cfg=a.read(cfgpath)
            biases={r['name']:r['bias'] for r in cfg['contacts']}
            assert abs(biases['gate']-vg)<1e-12 and abs(biases['drain']-c['vd'])<1e-12
            export=LOCAL/'native_exports'/key/f'vg_{index:03d}'
            assert d.matrix.spatial.m73.coordinate_error(export,geo)<1e-10
            nv={}
            for name in FIELDS+('srhRecombination','EffectiveIntrinsicDensity'):
                path=export/'fields'/(name+'_region0.csv');files.append(path)
                nv[name]=d.matrix.spatial.m73.scalar(path)
            ids=np.array(sorted(nv['eDensity']))
            assert set(ids)==set(np.flatnonzero(masks['all_si']))
            assert all(set(x)==set(ids) for x in nv.values())
            nv={name:np.array([x[i] for i in ids]) for name,x in nv.items()}
            vv={name:np.array([float(state[i][col])*scale for i in ids]) for name,col,scale in (
                ('ElectrostaticPotential','psi',1),('eDensity','electrons_m3',1e-6),('hDensity','holes_m3',1e-6))}
            for name,carrier in (('eQuasiFermiPotential','electron'),('hQuasiFermiPotential','hole')):
                vv[name]=np.array([float(Decimal(state[i][carrier+'_qf_reference_V'])+Decimal(state[i][carrier+'_qf_increment_V'])) for i in ids])
            doping=np.array([float(terms[i]['donors_m3'])+float(terms[i]['acceptors_m3']) for i in ids])
            ni=np.array([float(terms[i]['ni_eff_m3']) for i in ids]);tau={}
            # These legacy keys are internal cm^-3 for unit_scaling mode.
            for carrier in ('electron','hole'):
                q=cfg['solver']['srh_doping_dependence'][carrier]
                tau[carrier]=q['tau_min_s']+(q['tau_max_s']-q['tau_min_s'])/(1+(doping/q['reference_doping_m3'])**q['gamma'])
            split=np.array([float(Decimal(state[i]['hole_qf_reference_V'])+Decimal(state[i]['hole_qf_increment_V'])-Decimal(state[i]['electron_qf_reference_V'])-Decimal(state[i]['electron_qf_increment_V'])) for i in ids])
            rate=ni**2*np.expm1(split/d.fixed.upstream.VT)/(tau['hole']*(vv['eDensity']+ni)+tau['electron']*(vv['hDensity']+ni))
            vv['srhRecombination']=rate
            actual=np.array([float(terms[i]['electron_recombination']) for i in ids])*factor
            expected=rate*geo.volumes['all_cell'][ids]*d.fixed.Q
            rate_check=float(np.sum(abs(actual-expected))/max(np.sum(abs(actual)),1e-300))
            assert rate_check<1e-7,(key,vg,rate_check)
            common=dict(case=key,device=c['device'],vd=c['vd'],vg=vg,qualified=bool(qualified),initialization=chosen['folder'])
            si_area=geo.volumes['barycentric_si'][ids]
            for region in REGIONS:
                selected=masks[region][ids];weights=si_area[selected]
                for name in FIELDS:
                    left=vv[name][selected];right=nv[name][selected]
                    diff=np.log10(left/right) if 'Density' in name else left-right
                    peak=int(np.argmax(abs(diff)))
                    assert np.all(np.isfinite(diff))
                    fields.append(dict(**common,region=region,field=name,units='dex' if 'Density' in name else 'V',count=int(sum(selected)),
                        mean_signed=float(np.dot(weights,diff)/sum(weights)),weighted_rms=float(np.sqrt(np.dot(weights,diff**2)/sum(weights))),
                        max_abs=float(abs(diff[peak])),max_node=int(ids[selected][peak])))
                nr=nv['srhRecombination'][selected];vr=rate[selected]
                sources.append(dict(**common,region=region,native_integral_A_per_um=float(d.fixed.Q*np.dot(weights,nr)),
                    vela_integral_A_per_um=float(d.fixed.Q*np.dot(weights,vr)),delta_integral_A_per_um=float(d.fixed.Q*np.dot(weights,vr-nr)),
                    normalized_L1=float(np.dot(weights,abs(vr-nr))/max(np.dot(weights,abs(nr)),1e-300)),production_rate_reconstruction_L1=rate_check))
            vi=float(chosen['current_A_per_um']);sn=float(native['current_A_per_um'])
            points.append(dict(**common,native_A_per_um=sn,vela_A_per_um=vi,delta_A_per_um=vi-sn,relative_error_percent=100*(vi-sn)/sn,
                vela_kcl_over_Id=float(chosen['kcl_over_Id']),native_kcl_over_Id=float(native['kcl_over_Id']),local_violations=int(chosen['local_violations']),
                global_electron_source_qualified=chosen['global_electron_qualified'],global_hole_source_qualified=chosen['global_hole_qualified'],
                production_rate_reconstruction_L1=rate_check))
            for j,node in enumerate(ids):
                nodes.append(dict(**common,node_id=int(node),x_um=float(xy[node,0]),y_um=float(xy[node,1]),
                    **{name+'_native':float(values[j]) for name,values in nv.items()},
                    **{name+'_vela':float(values[j]) for name,values in vv.items()}))
    assert len(points)==84 and len(fields)==84*len(REGIONS)*len(FIELDS)
    for name,rows in (('points',points),('fields',fields),('srh',sources)):
        a.write_csv(OUT/(name+'.csv'),rows)
    a.write_csv(LOCAL/'comparison_nodes.csv',nodes)
    bins=[]
    for c in a.read(OUT/'vela_contract.json')['cases']:
        for low,high in ((0.,.5),(.55,1.),(0.,1.)):
            group=[r for r in points if r['case']==c['case'] and low<=r['vg']<=high and r['qualified']]
            worst=max(group,key=lambda r:abs(r['relative_error_percent']))
            bins.append(dict(case=c['case'],vg_low=low,vg_high=high,qualified_points=len(group),
                max_absolute_error_percent=abs(worst['relative_error_percent']),worst_vg=worst['vg'],
                median_absolute_error_percent=float(np.median([abs(r['relative_error_percent']) for r in group]))))
    a.write_csv(OUT/'curve_bins.csv',bins)
    fig,axes=plt.subplots(3,2,figsize=(12,11),sharex=True,layout='constrained')
    for col,vd in enumerate((.05,1.)):
        for device,color in (('n19','#3376b8'),('n23','#b33b35')):
            p=[r for r in points if r['device']==device and r['vd']==vd]
            x=[r['vg'] for r in p];ok=np.array([r['qualified'] for r in p])
            err=np.array([r['relative_error_percent'] for r in p]);err[~ok]=np.nan
            axes[0,col].plot(x,err,'o-',ms=3,label=device,color=color)
            if not all(ok):axes[0,col].plot(np.array(x)[~ok],[r['relative_error_percent'] for r in p if not r['qualified']],'x',color=color,label=device+' unqualified')
            for name,style,label in (('ElectrostaticPotential','-','psi'),('eQuasiFermiPotential','--','phi_n'),('hQuasiFermiPotential',':','phi_p')):
                f=[r for r in fields if r['device']==device and r['vd']==vd and r['region']=='channel' and r['field']==name]
                axes[1,col].plot([r['vg'] for r in f],[1e3*r['weighted_rms'] if r['qualified'] else np.nan for r in f],style,color=color,label=device+' '+label)
            sr=[r for r in sources if r['device']==device and r['vd']==vd and r['region']=='all_si']
            axes[2,col].plot([r['vg'] for r in sr],[100*r['normalized_L1'] if r['qualified'] else np.nan for r in sr],'o-',ms=3,color=color,label=device)
        axes[0,col].set_title(f'Vd = {vd:g} V');axes[2,col].set_xlabel('Vg (V)')
        for ax in axes[:,col]:ax.grid(alpha=.25);ax.legend(fontsize=8)
    axes[0,0].set_ylabel('Id error: (Vela / Sentaurus - 1) (%)')
    axes[1,0].set_ylabel('Channel weighted RMS difference (mV)')
    axes[2,0].set_ylabel('Si SRH normalized L1 difference (%)')
    fig.suptitle('Matched-ni SimpleMOS: strict state comparison, 0 to 1 V\nCrosses show excluded DC states; missing lines are qualification gaps',fontsize=13)
    fig.savefig(OUT/'fullfield_comparison.png',dpi=150);plt.close(fig)
    a.write(OUT/'analysis.json',dict(input_hashes={a.rel(p):a.sha(p) for p in sorted(set(files))},target_points=84,
        qualified_points=sum(r['qualified'] for r in points),unqualified_points=[r for r in points if not r['qualified']],
        attempted_vela_dc_states=len(attempts),native_dc_states=84,field_metrics='Vela minus Sentaurus; Si barycentric area weighted; densities log10(Vela/native).',
        native_acceptance='Native DC exit and bias checks; native KCL reported separately.',
        strict_acceptance='Unchanged carrier-row, global-floor and KCL gates; global source qualification flags retained separately.',
        srh='Configured generalized Boltzmann SRH checked at every state against production carrier probe; common Si support and volume for comparison.',
        unresolved_global_source_scale=True,production_changes=False,m82_released=False,m83_released=False))
    print('Compared 84 states, qualified:',sum(r['qualified'] for r in points),flush=True)


if __name__=='__main__':main()
