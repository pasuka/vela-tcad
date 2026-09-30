"""Full-curve native field comparisons, gated on all 204 dual-qualified points."""
import math
from decimal import Decimal,localcontext
from pathlib import Path
import numpy as np
import simplemos_enormal_curves_20260912 as curve
q=curve.q
a,d=q.a,q.d
O=curve.O/'fields'
D=lambda x:Decimal.from_float(float(x))

def main():
    a.verify(curve.O/'comparison_evidence.json');assert a.read(curve.O/'summary.json')['comparison_qualified']==204
    a.verify(curve.O/'supplemental_source_evidence.json')
    selected=a.rows(curve.O/'continuation_attempts.csv');fields=[];sources=[];checks=[];nodes=[];files=[Path(__file__).resolve(),curve.O/'comparison_evidence.json',curve.O/'export_evidence.json',curve.O/'supplemental_source_evidence.json']
    cases=[dict(**{k:v for k,v in c.items() if k!='vg'},key=c['case']+f'_vg_{i:03d}',vg=curve.n.GRID[i],index=i) for c in a.read(curve.O/'vela_contract.json')['cases'] for i in range(51)]
    for cc in cases:
        row=next(r for r in selected if r['case']==cc['case'] and int(r['index'])==cc['index'] and r['qualified']=='True')
        root=Path(row['dest']);cfg=a.read(root/'config.json')
        assert cfg['scaling']['mode']=='unit_scaling' and cfg['solver']['split_dd_state']
        assert cfg['solver']['srh_doping_dependence']['concentration_basis']=='total_impurity'
        geo,mask=q.V.m.previous.prior.support(cc);state=d.ordered(root/'state.csv',geo.count);terms=d.ordered(root/'all_row.csv',geo.count)
        native_root=curve.L/'native_exports'/cc['case']/f"vg_{cc['index']:03d}";native={}
        names=('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential','eDensity','hDensity','srhRecombination')
        for name in names:
            path=native_root/'fields'/(name+'_region0.csv');native[name]=q.p.scalar(native_root,name+'_region0.csv');files.append(path)
        # The qualification mask omits contacts; exported region fields include them.
        all_si=d.matrix.spatial.old.m78.supports(cc['device'],geo,.05)[0]['all_si']
        ids=np.array(sorted(native['eDensity']));assert set(ids)==set(np.flatnonzero(all_si))
        assert all(set(f)==set(ids) for f in native.values())
        vv={name:[] for name in names};density_errors=[]
        with localcontext() as ctx:
            ctx.prec=100
            for i in ids:
                s=state[i];t=terms[i];scale=D(s['packed_potential_scale_V'])
                psi=scale*(D(s['packed_psi'])+D(s['packed_psi_low']))
                fn=D(s['electron_qf_reference_V'])+scale*(D(s['packed_electron_qf_increment'])+D(s['packed_electron_qf_increment_low']))
                fp=D(s['hole_qf_reference_V'])+scale*(D(s['packed_hole_qf_increment'])+D(s['packed_hole_qf_increment_low']))
                vt=D(q.c.bgn.VT);assert scale==vt
                # Legacy probe *_m3 labels retain cm^-3 in this unit-scaling deck.
                ni=D(t['ni_eff_m3']);en=ni*((psi-fn)/vt).exp();hp=ni*((fp-psi)/vt).exp()
                density_errors += [abs(float(en)*1e6/float(s['electrons_m3'])-1),abs(float(hp)*1e6/float(s['holes_m3'])-1)]
                total=float(t['donors_m3'])+float(t['acceptors_m3']);tau={}
                for carrier in ('electron','hole'):
                    par=cfg['solver']['srh_doping_dependence'][carrier]
                    tau[carrier]=D(par['tau_min_s']+(par['tau_max_s']-par['tau_min_s'])/(1+(total/par['reference_doping_m3'])**par['gamma']))
                rate=ni*ni*(((fp-fn)/vt).exp()-1)/(tau['hole']*(en+ni)+tau['electron']*(hp+ni))
                for name,value in zip(names,(psi,fn,fp,en,hp,rate)):vv[name].append(float(value))
        weights=np.asarray(geo.volumes['barycentric_si'])[ids];assert np.all(weights>0)
        all_volume=np.asarray(geo.volumes['all_cell'])[ids]
        edges=a.rows(root/'acceptance_edges.csv');edge=max(edges,key=lambda r:abs(float(r['electron_flux'])))
        factor=float(edge['electron_particle_line_flux_per_m_s'])*1.602176634e-19*1e-6/float(edge['electron_flux'])
        actual=np.array([float(terms[i]['electron_recombination'])*factor for i in ids])
        expected=np.array(vv['srhRecombination'])*all_volume*1.602176634e-19
        denom=float(np.sum(abs(actual)));check=float(np.sum(abs(actual-expected)))/denom if denom else float(np.max(abs(expected)))
        checks.append(dict(key=cc['key'],density_reconstruction_max=max(density_errors),SRH_integral_reconstruction_L1=check,qualified=max(density_errors)<=1e-12 and check<=1e-7))
        for name in names[:-1]:
            nr=np.array([native[name][i] for i in ids]);vr=np.array(vv[name])
            diff=np.log10(vr/nr) if 'Density' in name else vr-nr
            peak=int(np.argmax(abs(diff)))
            fields.append(dict(key=cc['key'],device=cc['device'],vd=cc['vd'],vg=cc['vg'],scope='all_Si',field=name,units='dex' if 'Density' in name else 'V',mean_signed=float(np.dot(weights,diff)/sum(weights)),weighted_rms=float(np.sqrt(np.dot(weights,diff*diff)/sum(weights))),max_abs=float(abs(diff[peak])),max_node=int(ids[peak])))
        nr=np.array([native['srhRecombination'][i] for i in ids]);vr=np.array(vv['srhRecombination']);base=float(np.dot(weights,abs(nr)))
        sources.append(dict(key=cc['key'],device=cc['device'],vd=cc['vd'],vg=cc['vg'],common_volume='barycentric_Si_m2',native_charge_equivalent_A_per_um=float(1.602176634e-19*np.dot(weights,nr)),vela_charge_equivalent_A_per_um=float(1.602176634e-19*np.dot(weights,vr)),delta_charge_equivalent_A_per_um=float(1.602176634e-19*np.dot(weights,vr-nr)),normalized_L1=float(np.dot(weights,abs(vr-nr)))/base if base else None,native_max_abs_rate_cm3_s=float(max(abs(nr))),vela_max_abs_rate_cm3_s=float(max(abs(vr)))))
        for k,i in enumerate(ids):
            nodes.append(dict(key=cc['key'],node_id=int(i),**{name+'_native':native[name][i] for name in names},**{name+'_vela':vv[name][k] for name in names}))
        files += [root/'config.json',root/'state.csv',root/'all_row.csv',root/'acceptance_edges.csv']
        print(checks[-1],flush=True)
    for name,rows in [('fields',fields),('srh',sources),('reconstruction',checks),('nodes',nodes)]:a.write_csv(O/(name+'.csv'),rows);files.append(O/(name+'.csv'))
    summary=dict(points=len(checks),reconstruction_qualified=sum(r['qualified'] for r in checks),scope='Qualified Enormal Vg=0..1 V, .02 V lattice fields, comparing Vela with native, not comparing two Vela initializations.',units='Potentials V; densities and SRH rates cm^-3 and cm^-3 s^-1. Legacy all_row *_m3 headers are internal cm^-3 for this unit-scaling deck; state CSV density columns are SI m^-3.',SRH_meaning='Common positive barycentric Si weights isolate rate/state differences. These charge equivalents are not each solver own source integral, terminal current differences, or source-to-drain transfer predictions.',acceptance='Descriptive field metrics; no new acceptance threshold. The SRH formula is checked against the production integrated source before use.',all_reconstruction_qualified=all(r['qualified'] for r in checks))
    a.write(O/'summary.json',summary);files.append(O/'summary.json');d.matrix.freeze(O/'evidence.json',files)
    print(summary,flush=True);assert summary['all_reconstruction_qualified']

if __name__=='__main__':main()
