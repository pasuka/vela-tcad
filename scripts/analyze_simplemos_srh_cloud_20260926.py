"""Independent SRH rate/area reconstruction on qualified finite controls."""
import argparse,math
from pathlib import Path
from decimal import Decimal,localcontext
import simplemos_hfs_cloud_20260926 as h
from simplemos_srh_cloud_20260926 import signed_si

def analyze(root,baseline):
    selected=h.read(root/'finite/summary/selected.json');out=[];fields=[]
    contract=h.read(baseline/'inputs/contract.json')
    for r in selected:
        if r['arm']!='baseline':continue
        dest=Path(r['dest']);geo=h.read(baseline/'inputs'/r['case']/'geometry.json');cfg=h.read(dest/'config.json')
        native=h.read(baseline/'inputs'/r['case']/f"native/vg_{r['index']:03d}/fields.json")
        area=signed_si(h.read(baseline/'inputs'/r['case']/'mesh_file.json'));state=h.ordered(dest/'state.csv',geo['count']);terms=h.ordered(dest/'all_row.csv',geo['count']);ids=geo['all_si'];values={k:[] for k in ('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential','eDensity','hDensity','srhRecombination')};errs=[]
        D=lambda v:Decimal.from_float(float(v))
        with localcontext() as context:
            context.prec=100;vt=D(contract['VT'])
            for i in ids:
                s,t=state[i],terms[i];scale=D(s['packed_potential_scale_V']);assert scale==vt
                psi=scale*(D(s['packed_psi'])+D(s['packed_psi_low']))
                fn=D(s['electron_qf_reference_V'])+scale*(D(s['packed_electron_qf_increment'])+D(s['packed_electron_qf_increment_low']))
                fp=D(s['hole_qf_reference_V'])+scale*(D(s['packed_hole_qf_increment'])+D(s['packed_hole_qf_increment_low']))
                ni=D(t['ni_eff_m3']);n=ni*((psi-fn)/vt).exp();p=ni*((fp-psi)/vt).exp()
                errs.extend([abs(float(n)*1e6/float(s['electrons_m3'])-1),abs(float(p)*1e6/float(s['holes_m3'])-1)])
                dop=float(t['donors_m3'])+float(t['acceptors_m3']);tau={}
                for car in ('electron','hole'):
                    z=cfg['solver']['srh_doping_dependence'][car];tau[car]=D(z['tau_min_s']+(z['tau_max_s']-z['tau_min_s'])/(1+(dop/z['reference_doping_m3'])**z['gamma']))
                rate=ni*ni*(((fp-fn)/vt).exp()-1)/(tau['hole']*(n+ni)+tau['electron']*(p+ni))
                for k,v in zip(values,(psi,fn,fp,n,p,rate)):values[k].append(float(v))
        weights=[geo['barycentric_si'][i] for i in ids];label={k:r[k] for k in ('case','device','vd','index','vg')}
        for name in list(values)[:-1]:
            delta=[math.log10(v/native[name][str(i)]) if 'Density' in name else v-native[name][str(i)] for i,v in zip(ids,values[name])];k=max(range(len(ids)),key=lambda j:abs(delta[j]))
            fields.append(dict(**label,field=name,max_abs=abs(delta[k]),max_node=ids[k],weighted_rms=math.sqrt(math.fsum(w*d*d for w,d in zip(weights,delta))/math.fsum(weights))))
        edges=h.rows(dest/'acceptance_edges.csv');edge=max(edges,key=lambda e:abs(float(e['electron_flux'])));factor=float(edge['electron_particle_line_flux_per_m_s'])*1.602176634e-19*1e-6/float(edge['electron_flux'])
        actual=[float(terms[i]['electron_recombination'])*factor for i in ids];expected=[v*area[i]*1.602176634e-19 for i,v in zip(ids,values['srhRecombination'])]
        error=math.fsum(abs(a-b) for a,b in zip(actual,expected))/math.fsum(map(abs,actual))
        nr=[native['srhRecombination'][str(i)] for i in ids];vr=values['srhRecombination'];l1=math.fsum(w*abs(a-b) for w,a,b in zip(weights,vr,nr))/math.fsum(w*abs(a) for w,a in zip(weights,nr))
        out.append(dict(**label,density_reconstruction=max(errs),source_reconstruction=error,SRH_field_L1=l1,vela_integrated_SRH_A_per_um=math.fsum(expected),native_integrated_SRH_A_per_um=1.602192e-19*math.fsum(v*area[i] for i,v in zip(ids,nr)),qualified=max(errs)<=1e-12 and error<=1e-7))
    h.csvout(root/'finite_fields.csv',fields);h.csvout(root/'finite_sources.csv',out)
    ok=len(out)==16 and all(r['qualified'] for r in out)
    h.write(root/'field_review.json',dict(points=len(out),reconstruction_passed=ok,max_density_reconstruction=max(r['density_reconstruction'] for r in out),max_source_reconstruction=max(r['source_reconstruction'] for r in out),max_SRH_field_L1=max(r['SRH_field_L1'] for r in out),note='Field errors are measurements, not implied acceptance against absent native field thresholds.'))
    assert ok,'Independent field reconstruction failed'
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True);a=p.parse_args();analyze(a.root,a.baseline)
