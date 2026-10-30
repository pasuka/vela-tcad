"""Mesh-independent cut and local-spacing measurements for the P6 corpus."""
import math
from collections import defaultdict


def source_point(ir,x,y):
    for op in ir['mesh_control'].get('operations',[]):
        if op['axis']=='x':x=abs(x)
        else:y=abs(y)
    return x,y


def inside(window,x,y,width=0.):
    shape=window['shape']
    if shape=='rectangle':
        a,b=window['lower_left'],window['upper_right']
        return a[0]-1e-12<=x<=b[0]+1e-12 and a[1]-1e-12<=y<=b[1]+1e-12
    if shape=='line':
        a,b=window['points'];dx,dy=b[0]-a[0],b[1]-a[1]
        t=max(0.,min(1.,((x-a[0])*dx+(y-a[1])*dy)/(dx*dx+dy*dy)))
        return math.hypot(x-a[0]-t*dx,y-a[1]-t*dy)<=width+1e-12
    points=window['points'];hit=False
    for a,b in zip(points,points[1:]+points[:1]):
        if (a[1]>y)!=(b[1]>y) and x<(b[0]-a[0])*(y-a[1])/(b[1]-a[1])+a[0]:hit=not hit
    return hit


def active(ir,p,size,regions,x,y):
    x,y=source_point(ir,x,y)
    names={r.removesuffix('_mirrored') for r in regions}
    if p['target_kind']=='region':return p['target'] in names
    if p['target_kind']=='material':
        return any(b['region'] in names and b['material']['sde']==p['target'] for b in ir['geometry']['bodies'])
    window=next(w for w in ir['mesh_control']['windows'] if w['name']==p['target'])
    return inside(window,x,y,max(size['max_x'],size['max_y']))


def multibox_spacing(size,window,x,y):
    values=[]
    for i,(axis,q) in enumerate(zip('xy',(x,y))):
        ratio=size['ratio_'+axis]
        if ratio==0:values.append(math.inf);continue
        if abs(ratio)==1:values.append(size['max_'+axis]);continue
        distance=q-window['lower_left'][i] if ratio>0 else window['upper_right'][i]-q
        values.append(min(size['max_'+axis],size['min_'+axis]+(abs(ratio)-1)*max(0.,distance)))
    return values


def local_minimum(ir,region,x,y,axis):
    control=ir['mesh_control'];sizes={s['name']:s for s in control['refinement_sizes']};limits=[]
    for p in control['refinement_placements']:
        size=sizes[p['refinement']]
        if active(ir,p,size,{region},x,y):limits.append(size.get('min_'+axis,size['max_'+axis]))
    x,y=source_point(ir,x,y)
    for p in control.get('multibox_placements',[]):
        window=next(w for w in control['windows'] if w['name']==p['window'])
        if inside(window,x,y):
            size=next(s for s in control['multibox_sizes'] if s['name']==p['multibox'])
            limits.append(multibox_spacing(size,window,x,y)['xy'.index(axis)])
    if not limits:raise ValueError('cut has no authored spacing minimum')
    return min(limits)


def cut_roots(model,fixed_axis,coordinate,region,species):
    fixed='xy'.index(fixed_axis);normal=1-fixed;roots=[]
    for name,ids in model['cells']:
        if name!=region:continue
        points=model['vertices']
        if min(points[n][fixed] for n in ids)>coordinate+1e-12 or max(points[n][fixed] for n in ids)<coordinate-1e-12:continue
        if len(ids)!=3:raise ValueError('cut qualification requires native Tri3 cells')
        hits=[]
        for a,b in zip(ids,ids[1:]+ids[:1]):
            pa,pb=points[a],points[b];va,vb=species(model,name,a)[2],species(model,name,b)[2]
            if abs(pa[fixed]-coordinate)<=1e-12:hits.append((pa[normal],va))
            if (pa[fixed]-coordinate)*(pb[fixed]-coordinate)<0:
                t=(coordinate-pa[fixed])/(pb[fixed]-pa[fixed]);hits.append((pa[normal]+t*(pb[normal]-pa[normal]),va+t*(vb-va)))
        hits=sorted(set(hits))
        if len(hits)<2:continue
        (lo,flo),(hi,fhi)=hits[0],hits[-1]
        if hi-lo<=1e-12 or flo*fhi>0:continue
        if flo==fhi==0:continue
        roots.append(lo if flo==0 else hi if fhi==0 else lo-flo*(hi-lo)/(fhi-flo))
    result=[]
    for root in sorted(roots):
        if not result or abs(root-result[-1])>1e-8:result.append(root)
    return result


def cuts(case,ir,native,own,species):
    kind=case['contract']['doping']['kind'];region='R.Substrate';fixed_axis='x'
    if kind=='bjt':coordinates=[1.625,2.6,3.5];region='R.Silicon'
    elif kind=='mos':coordinates=[-.125,.125];fixed_axis='y'
    elif kind=='bvmethods':coordinates=[-.3,.3]
    elif kind=='pn':
        bounds=case['contract']['bounds_um'];coordinates=[(bounds[1]+bounds[3])/2];region='R.Si';fixed_axis='y'
    else:return []
    rows=[]
    for coordinate in coordinates:
        a=cut_roots(native,fixed_axis,coordinate,region,species);b=cut_roots(own,fixed_axis,coordinate,region,species)
        comparisons=[]
        for ref,value in zip(a,b):
            x,y=(coordinate,(ref+value)/2) if fixed_axis=='x' else ((ref+value)/2,coordinate)
            limit=local_minimum(ir,region,x,y,'y' if fixed_axis=='x' else 'x')
            comparisons.append({'native_um':ref,'own_um':value,'error_um':abs(ref-value),'limit_um':limit,'passed':abs(ref-value)<=limit})
        rows.append({'fixed_axis':fixed_axis,'coordinate_um':coordinate,'native_count':len(a),'own_count':len(b),
                     'junctions':comparisons,'passed':bool(a) and len(a)==len(b) and all(r['passed'] for r in comparisons)})
    return rows


def sizing(ir,model,edge_map):
    control=ir['mesh_control'];sizes={s['name']:s for s in control['refinement_sizes']};points=model['vertices'];reports=[]
    samples=[((points[a][0]+points[b][0])/2,(points[a][1]+points[b][1])/2,
              abs(points[a][0]-points[b][0]),abs(points[a][1]-points[b][1]),owners)
             for (a,b),owners in edge_map.items()]
    for p in control['refinement_placements']:
        size=sizes[p['refinement']]
        selected=[(dx/size['max_x'],dy/size['max_y']) for x,y,dx,dy,owners in samples if active(ir,p,size,owners,x,y)]
        if not selected:
            removed=any(b['region']==p['target'] and b.get('removed') for b in ir['geometry']['bodies'])
            reports.append({'name':p['name'],'active':False,'passed':removed});continue
        p95=[sorted(v[i] for v in selected)[max(0,math.ceil(.95*len(selected))-1)] for i in (0,1)]
        reports.append({'name':p['name'],'active':True,'edges':len(selected),'p95_normalized_components':p95,'passed':max(p95)<=1.35+1e-9})
    for p in control.get('multibox_placements',[]):
        size=next(s for s in control['multibox_sizes'] if s['name']==p['multibox'])
        window=next(w for w in control['windows'] if w['name']==p['window'])
        live=False
        for _,ids in model['cells']:
            x,y=source_point(ir,sum(points[n][0] for n in ids)/len(ids),sum(points[n][1] for n in ids)/len(ids))
            lo,hi=window['lower_left'],window['upper_right']
            if lo[0]+1e-10<x<hi[0]-1e-10 and lo[1]+1e-10<y<hi[1]-1e-10:live=True;break
        if not live:reports.append({'name':p['name'],'active':False,'reason':'no surviving volume in multibox','passed':True});continue
        selected=[]
        for x,y,dx,dy,_ in samples:
            x,y=source_point(ir,x,y)
            if inside(window,x,y):
                hx,hy=multibox_spacing(size,window,x,y);selected.append(max(dx/hx,dy/hy))
        maximum=max(selected,default=math.inf)
        reports.append({'name':p['name'],'active':True,'edges':len(selected),'maximum_normalized_component':maximum,'passed':maximum<=1.35+1e-9})
    return reports
