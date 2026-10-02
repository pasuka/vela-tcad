"""Independent native-field SG/SRH row replay at the selected source nodes.

This is a native exported-state formula audit, not native self-consistent FD.
Both qf and density SG reconstructions are retained, including weak signals.
"""
from pathlib import Path
from decimal import Decimal,localcontext
import math
import validate_simplemos_srh_finite_20260908 as t
import audit_simplemos_masetti_box_mobility_20260907 as native

a=t.a;d=t.d;OUT=t.OUT/'native_rows'

def main():
    a.verify(t.OUT/'freeze.json');scope=a.read(t.OUT/'contract.json');files=[Path(__file__).resolve(),Path(native.__file__),t.OUT/'freeze.json']
    fields=('ElectrostaticPotential','eDensity','hDensity','eQuasiFermiPotential','hQuasiFermiPotential','srhRecombination')
    for c in scope['cases']:
        root=t.c.s.prior.old.l.raw_root(c)/'fields';files += [root/(f+'_region0.csv') for f in fields]
        files.append(t.c.LOCAL/c['key']/'manual_constants/acceptance_edges.csv')
    files += [x for x in (native.prior.LOCAL/'native_exports').rglob('*') if x.is_file() and x.suffix=='.csv']
    files += [native.prior.LOCAL/'native_raw/bundle'/dev/'MeasureCoefficients.debug' for dev in ('n19','n23')]
    a.write(OUT/'contract.json',dict(scope='Two spatially matched nodes in four native Vg=.8 states. Native cell-weighted Masetti + native processed box g; no Vela state substituted.',
        residual='sum outward SG particle flux + native plotted SRH * selected volume, cm units. Inspect signed Si and original Vela SRH volume independently.',
        retained='Both density and qf SG representations. Tiny net source or rounded qf gradients may preclude a precise inferred volume; not a pass/fail native FD qualification.',native_new_runs=False))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json']);rows=[]
    for dev in ('n19','n23'):
        geo,weights,checks=native.geometry(dev);assert all(x['qualified'] for x in checks)
        full=t.c.s.prior.old.l.g.geometry(dev);vol=full[3]
        for c in (c for c in scope['cases'] if c['device']==dev):
            root=t.c.s.prior.old.l.raw_root(c)/'fields';raw={k:native.scalar(root/(f+'_region0.csv')) for k,f in zip(('psi','n','p','fn','fp','srh'),fields)}
            edges=a.rows(t.c.LOCAL/c['key']/'manual_constants/acceptance_edges.csv')
            for tag,record in c['mapped_nodes'].items():
                node=record['node'];source=raw['srh'][node];nativeV=float(vol['Si'][node])*1e4;oldV=float(geo.volumes['all_cell'][node])*1e4
                incident=[e for e in edges if int(e['node0'])==node or int(e['node1'])==node]
                for mode in ('qf','density'):
                    flux={'electron':[],'hole':[]}
                    for edge in incident:
                        i,j=int(edge['node0']),int(edge['node1']);sign=1 if i==node else -1
                        if i not in raw['n'] or j not in raw['n']:continue
                        en,hp=native.sg_flux(raw['psi'][i],raw['psi'][j],raw['n'][i],raw['n'][j],raw['p'][i],raw['p'][j],raw['fn'][i],raw['fn'][j],raw['fp'][i],raw['fp'][j],mode)
                        w=weights[tuple(sorted((i,j)))];flux['electron'].append(sign*float(en)*w['e']);flux['hole'].append(sign*float(hp)*w['h'])
                    for carrier,values in flux.items():
                        total=math.fsum(values);absolute=math.fsum(abs(x) for x in values);net=source*nativeV;scale=max(abs(total),abs(net),1e-300)
                        rows.append(dict(key=c['key'],device=dev,vd=c['vd'],tag=tag,node=node,carrier=carrier,mode=mode,incident_edges=len(incident),flux_per_cm_s=total,absolute_edge_flux_per_cm_s=absolute,native_SRH_cm3_s=source,signed_Si_source_per_cm_s=net,original_volume_source_per_cm_s=source*oldV,
                                         signed_Si_residual_over_source=abs(total+net)/scale,original_volume_residual_over_source=abs(total+source*oldV)/scale,
                                         signed_Si_residual_over_absolute_flux=abs(total+net)/max(absolute,abs(net),1e-300),inferred_volume_to_signed_Si=-total/net if net else 'zero_source',native_same_source_FD_qualified=False))
    a.write_csv(OUT/'row_replay.csv',rows);d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'row_replay.csv'])
    for r in rows:
        if r['device']=='n23' and r['vd']==1. and r['mode']=='qf':print({k:r[k] for k in ('node','carrier','signed_Si_residual_over_source','original_volume_residual_over_source','inferred_volume_to_signed_Si')},flush=True)

if __name__=='__main__':main()
