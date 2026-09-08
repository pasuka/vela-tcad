"""Diagnostic nodal export check, not an assertion of native internal mobility semantics."""
from pathlib import Path
import math
import numpy as np
import prepare_simplemos_reference_subset_20260906 as p

a=p.a;d=p.d;sc=d.matrix.spatial.m73.scalar
# C++ Masetti defaults in cm^2/(V s), cm^-3 at 300 K, without fitting.
PARAMETERS={
    'e':(1417.,52.2,52.2,43.4,0.,9.68e16,3.43e20,.68,2.),
    'h':(470.5,44.9,0.,29.,9.23e16,2.23e17,6.10e20,.719,2.),
}


def formula(n,pars):
    mc,m1,m2,mx,pc,cr,cs,alpha,beta=pars
    return mc if n==0 else m1*math.exp(-pc/n)+(mc-m2)/(1+(n/cr)**alpha)-mx/(1+(cs/n)**beta)


def main():
    a.verify(p.OUT/'native_freeze.json');a.verify(p.OUT/'export_contract.json')
    summary=[];nodes=[];invariant={};files=[Path(__file__).resolve(),d.REPO/'include/vela/physics/MobilityModel.h',d.REPO/'src/physics/MobilityModel.cpp']
    for job in a.read(p.OUT/'export_contract.json')['jobs']:
        if job['model']!='masetti':continue
        src=Path(job['export']);device='n19' if 'n19' in job['case'] else 'n23';geo=d.matrix.spatial.m73.Geometry(device)
        manifest=a.read(src/'field_manifest.json')['fields']
        for name,unit in [('eMobility','cm^2*V^-1*s^-1'),('hMobility','cm^2*V^-1*s^-1'),('DonorConcentration','cm^-3'),('AcceptorConcentration','cm^-3')]:
            field=next(x for x in manifest if x['name']==name and x['region']==0)
            assert field['unit']==unit and field['support_kind']=='node' and field['mapping_status']=='complete'
            files.append(src/'fields'/field['csv_file'])
        nd=sc(src/'fields/DonorConcentration_region0.csv');na=sc(src/'fields/AcceptorConcentration_region0.csv')
        for c,pars in PARAMETERS.items():
            native=sc(src/f'fields/{c}Mobility_region0.csv');assert native.keys()==nd.keys()==na.keys()
            old=invariant.setdefault((device,c),native)
            assert old==native, 'Masetti native nodal export changed across bias'
            rr=[]
            for i,value in native.items():
                calc=formula(nd[i]+na[i],pars);error=calc/value-1
                row=dict(case=job['case'],index=job['index'],carrier=c,node=i,x_um=geo.coords[i][0],y_um=geo.coords[i][1],
                    total_impurity_cm3=nd[i]+na[i],local_formula_cm2_per_Vs=calc,native_nodal_export_cm2_per_Vs=value,
                    signed_relative_difference=error)
                rr.append(row);nodes.append(row)
            worst=max(rr,key=lambda x:abs(x['signed_relative_difference']))
            summary.append(dict(case=job['case'],index=job['index'],carrier=c,nodes=len(rr),
                median_absolute_relative_difference=float(np.median([abs(x['signed_relative_difference']) for x in rr])),
                p95_absolute_relative_difference=float(np.percentile([abs(x['signed_relative_difference']) for x in rr],95)),
                max_absolute_relative_difference=abs(worst['signed_relative_difference']),worst_node=worst['node'],
                worst_local_formula_cm2_per_Vs=worst['local_formula_cm2_per_Vs'],worst_native_export_cm2_per_Vs=worst['native_nodal_export_cm2_per_Vs']))
    a.write_csv(p.OUT/'mobility_nodal_diagnostic.csv',summary);a.write_csv(p.OUT/'mobility_nodal_diagnostic_rows.csv',nodes)
    a.write(p.OUT/'mobility_nodal_diagnostic.json',dict(status='diagnostic_only',parameters=PARAMETERS,bias_invariance='Exact equality of native nodal e/h mobility over four biases per device.',
        meaning='Compare the local C++-default formula at exported nodal total impurity with native node-plotted mobility. Native internal edge/cell averaging and node projection have not been mapped; these differences are NOT a calibrated formula-error or Id-attribution measure.',
        fitted_parameters=False,input_hashes={a.rel(f):a.sha(f) for f in files}))
    print('Recorded nodal diagnostic and confirmed bias invariance:',len(summary),'rows',flush=True)


if __name__=='__main__':main()
