"""Expose which QF field differences lie outside the frozen carrier-row gate."""
from pathlib import Path
import numpy as np
import run_simplemos_fullfield_vela_20260906 as b
import refine_simplemos_native_precision_20260906 as precision

a=b.a;d=b.d;OUT=precision.OUT/'comparison/visibility'


def main():
    a.verify(precision.OUT/'comparison/analysis.json')
    files=[Path(__file__).resolve(),precision.OUT/'comparison/analysis.json',d.REPO/'src/solver/NewtonSolver.cpp']
    nodegroups={}
    for r in a.rows(precision.LOCAL/'comparison_nodes.csv'):nodegroups.setdefault((r['case'],round(float(r['vg']),10)),[]).append(r)
    summaries=[];peaks=[]
    for point in a.rows(precision.OUT/'comparison/points.csv'):
        key=point['case'];vg=round(float(point['vg']),10);geo=d.matrix.spatial.m73.Geometry(point['device'])
        masks,_,_=d.matrix.spatial.old.m78.supports(point['device'],geo,.05)
        dest=b.LOCAL/point['initialization']/key/f'vg_{round(vg/.05):03d}'
        terms=d.ordered(dest/'carrier.csv',geo.count);cfg=a.read(dest/'config.json')['solver']['carrier_row_convergence']
        status=a.read(dest/'acceptance.status.json')['carrier_row_convergence']
        assert cfg['min_carrier_density_m3']==0
        count=0;rows=nodegroups[key,vg];ids=np.array([int(r['node_id']) for r in rows]);area=geo.volumes['barycentric_si'][ids]
        for carrier,field,density in (('electron','eQuasiFermiPotential','eDensity'),('hole','hQuasiFermiPotential','hDensity')):
            flux=np.array([abs(float(r[carrier+'_flux_abs_sum'])) for r in terms])
            source=np.array([max(abs(float(r[carrier+'_recombination'])),abs(float(r[carrier+'_impact']))) for r in terms])
            scale=np.maximum(np.maximum(flux,source),cfg['scale_floor'])
            sq=(source>0)&(source>=cfg['min_source_scale'])&(source>=cfg['min_source_scale_fraction']*scale)&(source>=cfg['min_source_global_fraction']*max(source))
            fq=(flux>0)&(flux>=max(cfg['min_flux_scale'],cfg['min_flux_scale_fraction']*max(flux)))
            qualified=sq|fq;count+=int(sum(qualified))
            diff=np.array([float(r[field+'_vela'])-float(r[field+'_native']) for r in rows])
            for region in ('all_si','channel','gate_interface'):
                selected=masks[region][ids];selected_good=selected&qualified[ids];ignored=selected&~qualified[ids]
                energy=float(np.dot(area[selected],diff[selected]**2));peak=np.flatnonzero(selected)[np.argmax(abs(diff[selected]))]
                summaries.append(dict(case=key,device=point['device'],vd=point['vd'],vg=vg,carrier=carrier,region=region,
                    selected_nodes=int(sum(selected)),qualified_nodes=int(sum(selected_good)),ignored_nodes=int(sum(ignored)),
                    qualified_area_fraction=float(sum(area[selected_good])/sum(area[selected])),
                    full_rms_V=float(np.sqrt(energy/sum(area[selected]))),
                    qualified_rms_V=float(np.sqrt(np.dot(area[selected_good],diff[selected_good]**2)/sum(area[selected_good]))) if any(selected_good) else None,
                    ignored_squared_error_fraction=float(np.dot(area[ignored],diff[ignored]**2)/max(energy,1e-300))))
                peaks.append(dict(case=key,device=point['device'],vd=point['vd'],vg=vg,carrier=carrier,region=region,node_id=int(ids[peak]),
                    delta_V=float(diff[peak]),row_qualified=bool(qualified[ids[peak]]),
                    native_density_cm3=float(rows[peak][density+'_native']),vela_density_cm3=float(rows[peak][density+'_vela']),
                    normalized_flux_scale=float(flux[ids[peak]]),normalized_source_scale=float(source[ids[peak]])))
        assert count==status['qualified_row_count'],(key,vg,count,status['qualified_row_count'])
        files += [dest/'carrier.csv',dest/'config.json',dest/'acceptance.status.json']
    a.write_csv(OUT/'fields.csv',summaries);a.write_csv(OUT/'peaks.csv',peaks)
    a.write(OUT/'result.json',dict(input_hashes={a.rel(p):a.sha(p) for p in sorted(set(files))},
        checked_states=84,qualification_counts_match_production=True,
        definition='Exact frozen source/flux criteria; density criterion disabled in this deck. The reconstructed counts agree with the production acceptance probe at all states.',
        interpretation='Visibility of field differences to the existing Vela carrier gate; no claim that native states satisfy a Vela equation or that ignored fields are accurate.'))
    print('Carrier-row visibility checked at all 84 states',flush=True)


if __name__=='__main__':main()
