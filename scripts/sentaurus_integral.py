"""Bounded C0 doping-integral convergence for Gaussian-containing meshes."""
import math

# Degree-five Dunavant triangle quadrature, with weights summing to one.
QUADRATURE=[((1/3,1/3,1/3),.225)]
for a,b,w in [(.059715871789770,.470142064105115,.132394152788506),
              (.797426985353087,.101286507323456,.125939180544827)]:
    QUADRATURE.extend([((a,b,b),w),((b,a,b),w),((b,b,a),w)])


def enabled(ir,target):
    return target>0 and any(p['type']=='gaussian' for p in ir['doping']['profiles'])


def analyse(ir,coordinates,by_surface,surface_owner,target,cache):
    from sentaurus_doping import _rectangles,_windows,_evaluate_m2_doping_for_owner
    rectangles,windows=_rectangles(ir),_windows(ir)
    def values(owner,point):
        key=(owner['id'],point[0],point[1])
        if key not in cache:
            cache[key]=_evaluate_m2_doping_for_owner(ir,owner,rectangles,windows,*point[:2])[:2]
        return cache[key]
    reference=[0.,0.];linear=[0.,0.];cells=[]
    for surface,triangles in by_surface.items():
        owner=surface_owner[surface]
        if owner['material']['vela'] not in {'Si','PolySilicon'}:continue
        for ids in triangles:
            p=[coordinates[n] for n in ids]
            area=abs((p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[1][1]-p[0][1])*(p[2][0]-p[0][0]))/2
            nodal=[values(owner,q) for q in p]
            fe=[area*sum(v[k] for v in nodal)/3 for k in (0,1)]
            quad=[0.,0.]
            for bary,w in QUADRATURE:
                q=tuple(sum(bary[i]*p[i][axis] for i in range(3)) for axis in (0,1))
                value=values(owner,q)
                for k in (0,1):quad[k]+=area*w*value[k]
            for k in (0,1):linear[k]+=fe[k];reference[k]+=quad[k]
            cells.append((ids,[fe[k]-quad[k] for k in (0,1)]))
    errors=[linear[k]-reference[k] for k in (0,1)]
    relative=[abs(errors[k])/max(abs(reference[k]),1.) for k in (0,1)]
    # A signed global error can hide cancellation or alternate sign after a
    # single large cell is split. Control the sum of absolute cell errors.
    absolute=[sum(abs(error[k]) for _,error in cells) for k in (0,1)]
    bounds=[absolute[k]/max(abs(reference[k]),1.) for k in (0,1)]
    selected=set()
    for k in (0,1):
        if bounds[k]<=target:continue
        aligned=sorted(((abs(error[k]),tuple(ids)) for ids,error in cells if error[k]!=0),reverse=True)
        # Bulk marking concentrates work on the largest local errors. Marking
        # the entire predicted deficit needlessly splits nearly every cell.
        covered=0.;wanted=.75*absolute[k]
        for amount,ids in aligned:
            selected.add(ids);covered+=amount
            if covered>=wanted:break
    marked=set()
    for ids in selected:
        edges=[tuple(sorted((ids[i],ids[(i+1)%3]))) for i in range(3)]
        marked.add(max(edges,key=lambda edge:math.dist(coordinates[edge[0]],coordinates[edge[1]])))
    return marked,{'relative_target':target,'relative_errors':relative,'relative_error_bounds':bounds,
                   'linear_integrals_cm3_um2':linear,'quadrature_integrals_cm3_um2':reference,
                   'mapping':'absolute cell-integral error bound versus degree-five quadrature; conforming longest-edge splits',
                   'applicability':'Gaussian-containing unstructured models'}
