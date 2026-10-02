"""Record native precision sensitivity without replacing the original evidence."""
from pathlib import Path
import numpy as np
import refine_simplemos_native_precision_20260906 as p

a=p.a;d=p.d


def main():
    a.verify(p.OUT/'export_contract.json');a.verify(p.OUT/'comparison/analysis.json')
    original={(r['case'],r['vg']):r for r in a.rows(p.base.OUT/'native_points.csv')}
    points=a.rows(p.OUT/'native_points.csv');output=[]
    paths=[Path(__file__).resolve(),p.OUT/'export_contract.json',p.OUT/'comparison/analysis.json',p.base.OUT/'native_points.csv',p.OUT/'native_points.csv']
    for c in a.read(p.OUT/'contract.json')['cases']:
        geo=d.matrix.spatial.m73.Geometry('n23')
        for r in points:
            if r['case']!=c['case']:continue
            index=round(float(r['vg'])/.05)
            old=p.base.LOCAL/'native_exports'/c['case']/f'vg_{index:03d}'
            new=p.LOCAL/'native_exports'/c['case']/f'vg_{index:03d}'
            for name in ('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential','eDensity','hDensity','srhRecombination'):
                pair=[q/'fields'/(name+'_region0.csv') for q in (old,new)];paths+=pair
                values=[d.matrix.spatial.m73.scalar(q) for q in pair]
                ids=sorted(values[0]);assert ids==sorted(values[1])
                x,y=[np.array([q[i] for i in ids]) for q in values];delta=y-x;weights=geo.volumes['barycentric_si'][ids]
                output.append(dict(case=c['case'],vg=float(r['vg']),field=name,
                    max_abs_change=float(max(abs(delta))),weighted_rms_change=float(np.sqrt(np.average(delta**2,weights=weights))),
                    weighted_L1_relative=float(np.dot(weights,abs(delta))/max(np.dot(weights,abs(x)),1e-300))))
    changes=[dict(case=r['case'],vg=float(r['vg']),current_change_percent=100*(float(r['current_A_per_um'])/float(original[r['case'],r['vg']]['current_A_per_um'])-1),
                  original_kcl_over_Id=float(original[r['case'],r['vg']]['kcl_over_Id']),refined_kcl_over_Id=float(r['kcl_over_Id'])) for r in points]
    a.write_csv(p.OUT/'native_field_stability.csv',output);a.write_csv(p.OUT/'native_current_stability.csv',changes)
    a.write(p.OUT/'precision_result.json',dict(input_hashes={a.rel(q):a.sha(q) for q in sorted(set(paths))},
        points=42,kcl_qualified=sum(r['native_kcl_qualified']=='True' for r in points),
        maximum_kcl_over_Id=max(float(r['kcl_over_Id']) for r in points),
        maximum_absolute_current_change_percent=max(abs(r['current_change_percent']) for r in changes),
        maximum_potential_change_V=max(r['max_abs_change'] for r in output if 'Potential' in r['field']),
        limitation='Checks returned field precision and contact-current stability; not an independent proof of every native residual row.'))
    print('Native precision sensitivity summarized',flush=True)


if __name__=='__main__':main()
