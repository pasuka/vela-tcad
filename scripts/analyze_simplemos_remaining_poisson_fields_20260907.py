"""Physical field comparison on all 907 free silicon nodes, no new acceptance gate."""
from decimal import Decimal
from pathlib import Path
import math
import numpy as np
import validate_simplemos_remaining_dielectric_20260907 as t

a=t.a;d=t.d;l=t.l;OUT=t.OUT

def main():
    a.verify(OUT/'validation_evidence.json');rows=[];nodes=[];files=[Path(__file__).resolve(),OUT/'validation_evidence.json']
    for c in a.read(OUT/'contract.json')['cases']:
        geo,mask=t.v.m.previous.prior.support(c);ids=np.flatnonzero(mask);assert len(ids)==907
        native={}
        for field,name in [('psi','ElectrostaticPotential'),('phin','eQuasiFermiPotential'),('phip','hQuasiFermiPotential'),('electrons_m3','eDensity'),('holes_m3','hDensity'),('SRH','srhRecombination')]:
            path=l.raw_root(c)/'fields'/(name+'_region0.csv');files.append(path)
            native[field]={int(x['node_id']):Decimal(x['component0']) for x in a.rows(path)}
        for label,statepath,termpath,edgepath in [
            ('joint',t.v.LOCAL/c['key']/'from_vela/state.csv',t.v.LOCAL/c['key']/'from_vela/all_row.csv',l.LOCAL/c['key']/'joint/edges.csv'),
            ('native_mu',t.q.LOCAL/c['key']/'replacement/state.csv',t.q.LOCAL/c['key']/'replacement/all_row.csv',l.LOCAL/c['key']/'native_mu/edges.csv')]+[
            (axis,t.LOCAL/c['key']/axis/'replacement/state.csv',t.LOCAL/c['key']/axis/'replacement/all_row.csv',t.LOCAL/c['key']/axis/'replacement/post/edges.csv') for axis in t.AXES]:
            state=d.ordered(statepath,geo.count);terms=d.ordered(termpath,geo.count);edges=a.rows(edgepath)
            ratio=next(float(e['electron_particle_line_flux_per_m_s'])/float(e['electron_flux']) for e in edges if float(e['electron_flux'])!=0)
            for e in edges:
                if float(e['electron_flux'])!=0:
                    check=float(e['electron_particle_line_flux_per_m_s'])/float(e['electron_flux'])
                    assert abs(check/ratio-1)<=1e-12
            values={}
            for field in ('psi','phin','phip'):
                values[field]=np.array([float(t.v.m.previous.prior.physical(state[i],field)-native[field][i]) for i in ids])
                err=values[field];worst=int(ids[np.argmax(abs(err))])
                rows.append(dict(key=c['key'],state=label,field=field,units='V',nodes=907,signed_min=float(min(err)),signed_max=float(max(err)),rms=float(np.linalg.norm(err)/math.sqrt(len(err))),max_absolute=float(max(abs(err))),worst_node=worst,relative_l2='not_applicable'))
            for field in ('electrons_m3','holes_m3'):
                reference=np.array([float(native[field][i])*1e6 for i in ids]);actual=np.array([float(state[i][field]) for i in ids]);err=actual/reference-1;values[field]=err
                rows.append(dict(key=c['key'],state=label,field=field,units='relative',nodes=907,signed_min=float(min(err)),signed_max=float(max(err)),rms=float(np.linalg.norm(err)/math.sqrt(len(err))),max_absolute=float(max(abs(err))),worst_node=int(ids[np.argmax(abs(err))]),relative_l2=float(np.linalg.norm(actual-reference)/np.linalg.norm(reference))))
            rate=np.array([float(terms[i]['electron_recombination'])*ratio/geo.volumes['all_cell'][i]/1e6 for i in ids]);reference=np.array([float(native['SRH'][i]) for i in ids]);err=rate-reference
            rows.append(dict(key=c['key'],state=label,field='SRH',units='cm^-3 s^-1',nodes=907,signed_min=float(min(err)),signed_max=float(max(err)),rms=float(np.linalg.norm(err)/math.sqrt(len(err))),max_absolute=float(max(abs(err))),worst_node=int(ids[np.argmax(abs(err))]),relative_l2=float(np.linalg.norm(err)/np.linalg.norm(reference))))
            for pos,i in enumerate(ids):
                if i not in (320,324,338,343,1000,1009,1075,1077,1091,1092):continue
                nodes.append(dict(key=c['key'],state=label,node_id=int(i),x_um=geo.coords[i][0],y_um=geo.coords[i][1],
                    psi_delta_V=float(values['psi'][pos]),phin_delta_V=float(values['phin'][pos]),phip_delta_V=float(values['phip'][pos]),electron_density_relative=float(values['electrons_m3'][pos]),hole_density_relative=float(values['holes_m3'][pos]),
                    Vela_SRH_cm3_s=float(rate[pos]),native_plot_SRH_cm3_s=float(reference[pos]),SRH_relative=float(err[pos]/reference[pos]) if reference[pos]!=0 else 'undefined_zero_native'))
            files += [statepath,termpath,edgepath]
    a.write_csv(OUT/'full_Si_field_metrics.csv',rows);a.write_csv(OUT/'selected_node_fields.csv',nodes)
    # Reconcile this direct native-field calculation with the earlier local delta-chain calculation.
    checks=[]
    for old in a.rows(OUT/'local_fields.csv'):
        new=next(r for r in nodes if r['key']==old['key'] and r['state']==old['axis'] and r['node_id']==int(old['node_id']))
        error=max(abs(new[f+'_delta_V']-float(old[f+'_new_delta_V'])) for f in ('psi','phin','phip'))
        srh=abs(new['Vela_SRH_cm3_s']/float(old['new_SRH_cm3_s'])-1)
        checks.append(dict(key=old['key'],axis=old['axis'],node_id=old['node_id'],delta_chain_vs_direct_max_V=error,SRH_conversion_relative=srh,qualified=error<=1e-12 and srh<=1e-12))
    a.write_csv(OUT/'field_reconciliation.csv',checks)
    d.matrix.freeze(OUT/'fields_evidence.json',files+[OUT/'full_Si_field_metrics.csv',OUT/'selected_node_fields.csv',OUT/'field_reconciliation.csv'])
    print('Independent physical-field reconciliation',sum(x['qualified'] for x in checks),'/',len(checks),flush=True)

if __name__=='__main__':main()
