"""Mesh-derived all-row scope for later original-matrix HFS validation.
Same physical acceptance equations and thresholds as the qualified four curves.
No global helper mutation and no 907-node assumption.
"""
import math
import numpy as np
import simplemos_hfs_curves_20260914 as e
a,d=e.a,e.d

def support(c):
    geo=d.matrix.spatial.m73.Geometry(c['device'])
    mask=d.matrix.spatial.old.m78.supports(c['device'],geo,.05)[0]['all_si'].copy()
    mask[geo.contact_nodes]=False
    expected=next(r for r in a.read(e.O/'original_mesh_scope.json')['rows'] if r['device']==c['device'])
    assert geo.count==expected['nodes'] and int(mask.sum())==expected['free_Si_nodes']
    return geo,mask

def qualify(c,dest,write=True):
    """Pure independent acceptance; runner/env chosen explicitly before this call."""
    s=a.read(dest/'config.status.json');ps=a.read(dest/'all_row.status.json');es=a.read(dest/'acceptance_edges.status.json')
    geo,mask=support(c);terms=d.ordered(dest/'all_row.csv',geo.count);edges=a.rows(dest/'acceptance_edges.csv');contacts=set(map(int,geo.contact_nodes))
    ratios=[];closure={}
    for r,keep in zip(terms,mask):
        if keep:
            for car in ('electron','hole'):
                scale=max(float(r[car+'_flux_abs_sum']),abs(float(r[car+'_recombination'])),abs(float(r[car+'_impact'])))
                ratios.append(abs(float(r[car+'_residual']))/scale if scale else math.inf)
    assert len(ratios)==2*int(np.count_nonzero(mask))
    for car in ('electron','hole'):
        flux=math.fsum((int(int(x['node0']) in contacts)-int(int(x['node1']) in contacts))*float(x[car+'_flux']) for x in edges)
        source=math.fsum(float(r[car+'_recombination'])+float(r[car+'_impact']) for i,r in enumerate(terms) if i not in contacts)
        ratio=abs(flux-source)/max(abs(flux),abs(source),1e-10)
        closure[car]=dict(contact_flux=flux,integrated_source=source,ratio=ratio,qualified=abs(source)>=1e-10,satisfied=abs(source)<1e-10 or ratio<=1e-6)
    cc=s.get('contact_currents_A_per_um',{});Id=cc.get('drain',math.nan);kcl=abs(math.fsum(cc.values()))/abs(Id);bad=sum(x>1e-6 for x in ratios)
    result=dict(qualified=s['exit_code']==ps['exit_code']==es['exit_code']==0 and s.get('converged',False) and bad==0 and all(x['satisfied'] for x in closure.values()) and kcl<=1e-8,
        current_A_per_um=Id,max_row_ratio=max(ratios),row_violations=bad,kcl_over_Id=kcl,iterations=s.get('iterations'),failure=s.get('failure_reason',''),
        global_electron_ratio=closure['electron']['ratio'],global_hole_ratio=closure['hole']['ratio'])
    record=dict(result=result,closure=closure);target=dest/'independent_acceptance.json'
    if target.exists():assert a.read(target)==record
    elif write:a.write(target,record)
    return result
