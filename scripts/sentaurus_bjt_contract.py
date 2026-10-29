"""Independent BJT oracle from the frozen geometric doping contract."""
from __future__ import annotations
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ACCEPTANCE = REPO / 'configs/reference_bjt_acceptance.json'


def concentration(profiles, x, y):
    totals = {'donor': 0.0, 'acceptor': 0.0}
    for p in profiles:
        x0, y0, x1, y1 = p['window_um']
        value = float(p['peak_cm3'])
        if p['char_um']:
            dx = max(x0-x, 0.0, x-x1)/p['char_um'][0]
            dy = max(y0-y, 0.0, y-y1)/p['char_um'][1]
            value *= math.exp(-dx*dx-dy*dy)
        totals[p['species']] += value
    return totals['donor'], totals['acceptor']


def _axis_integral(lo, hi, wlo, whi, char):
    left = char*math.sqrt(math.pi)/2*(math.erf((min(hi,wlo)-wlo)/char)-math.erf((lo-wlo)/char)) if lo < wlo else 0.0
    center = max(0.0, min(hi,whi)-max(lo,wlo))
    right = char*math.sqrt(math.pi)/2*(math.erf((hi-whi)/char)-math.erf((max(lo,whi)-whi)/char)) if hi > whi else 0.0
    return left+center+right


def integrals(profiles, bounds):
    totals = {'donor': 0.0, 'acceptor': 0.0}
    for p in profiles:
        x0,y0,x1,y1 = p['window_um']
        a,b,c,d = bounds
        volume = (c-a)*(d-b) if p['char_um'] is None else (
            _axis_integral(a,c,x0,x1,p['char_um'][0])*_axis_integral(b,d,y0,y1,p['char_um'][1]))
        totals[p['species']] += p['peak_cm3']*volume
    return totals


def exact_roots(profiles, x):
    def net(y):
        d,a = concentration(profiles,x,y)
        return d-a
    roots=[]
    for i in range(2000):
        lo,hi = i/1000, (i+1)/1000
        flo,fhi=net(lo),net(hi)
        if flo*fhi>0:
            continue
        for _ in range(48):
            mid=(lo+hi)/2
            if flo*net(mid)<=0:
                hi=mid
            else:
                lo=mid;flo=net(lo)
        root=(lo+hi)/2
        if not roots or abs(root-roots[-1])>1e-8:
            roots.append(root)
    return roots


def mesh_roots(bundle, x):
    points={n['vela_id']:n for n in bundle['nodes']}
    values=bundle['fields']['DopingConcentration']['values']
    roots=[]
    for triangle in bundle['triangles']:
        ids=triangle['node_ids']
        intersections=[]
        for i in range(3):
            a,b=points[ids[i]],points[ids[(i+1)%3]]
            if abs(a['x']-b['x'])<1e-14:
                if abs(x-a['x'])<1e-12:
                    intersections += [(a['y'],values[ids[i]]),(b['y'],values[ids[(i+1)%3]])]
            elif min(a['x'],b['x'])-1e-12<=x<=max(a['x'],b['x'])+1e-12:
                t=(x-a['x'])/(b['x']-a['x'])
                intersections.append((a['y']+t*(b['y']-a['y']), values[ids[i]]+t*(values[ids[(i+1)%3]]-values[ids[i]])))
        if len(intersections)<2:
            continue
        a,b=min(intersections),max(intersections)
        if a[1]*b[1]<=0 and a[1]!=b[1]:
            roots.append(a[0]-a[1]*(b[0]-a[0])/(b[1]-a[1]))
    unique=[]
    for y in sorted(roots):
        if not unique or y-unique[-1]>1e-8:
            unique.append(y)
    return unique


def check(bundle, *, local=False):
    acceptance=json.loads(ACCEPTANCE.read_text())
    source=REPO/acceptance['source']
    if hashlib.sha256(source.read_bytes().replace(b'\r\n',b'\n')).hexdigest()!=acceptance['source_sha256_lf']:
        raise ValueError('independent BJT geometry contract changed')
    contract=json.loads(source.read_text())
    failures=[]
    points={n['vela_id']:n for n in bundle['nodes']}
    if len(points)<acceptance['minimum_total_nodes']:
        failures.append('minimum node count')
    max_node_error=0.0
    for n in points.values():
        nd,na=concentration(contract['profiles'],n['x'],n['y'])
        for name,exact in [('DonorConcentration',nd),('AcceptorConcentration',na),('DopingConcentration',nd-na)]:
            actual=bundle['fields'][name]['values'][n['vela_id']]
            error=abs(actual-exact)/max(abs(exact),1.0)
            max_node_error=max(max_node_error,error)
    if max_node_error>acceptance['pointwise_relative_tolerance']:
        failures.append('pointwise concentration')
    actual={'donor':0.0,'acceptor':0.0}
    for triangle in bundle['triangles']:
        ids=triangle['node_ids'];a,b,c=[points[i] for i in ids]
        area=((b['x']-a['x'])*(c['y']-a['y'])-(b['y']-a['y'])*(c['x']-a['x']))/2
        if area<=0:
            failures.append('triangle orientation');break
        for species,field in [('donor','DonorConcentration'),('acceptor','AcceptorConcentration')]:
            actual[species]+=area*sum(bundle['fields'][field]['values'][i] for i in ids)/3
    exact=integrals(contract['profiles'],contract['region']['bounds_um'])
    integral={k:{'analytic_cm3_um2':exact[k], 'mesh_linear_cm3_um2':actual[k],
                 'relative_error':abs(actual[k]-exact[k])/exact[k]} for k in exact}
    if any(v['relative_error']>acceptance['integral_relative_tolerance'] for v in integral.values()):
        failures.append('integrated concentration')
    # Frozen source sizing minima, independent of generated IR and mesh metric code.
    sizing=[([0,0,6,2],.008),([1,0,5,.85],.004),([0,1.45,6,2],.004)]
    if local:
        sizing += [([.75,.42,5.25,.88],.002),([2.52,0,2.98,.55],.0015),([4.02,0,4.48,.55],.0015)]
    cuts=[]
    for x in acceptance['cuts_x_um']:
        expected=exact_roots(contract['profiles'],x)
        observed=mesh_roots(bundle,x)
        rows=[]
        if len(expected)!=len(observed):
            failures.append(f'junction count at x={x}')
        for a,b in zip(expected,observed):
            limit=min(h for (x0,y0,x1,y1),h in sizing if x0<=x<=x1 and y0<=a<=y1)
            rows.append({'analytic_y_um':a,'mesh_linear_y_um':b,'absolute_error_um':abs(a-b),'limit_um':limit})
            if abs(a-b)>limit:
                failures.append(f'junction location at x={x}')
        cuts.append({'x_um':x,'analytic_count':len(expected),'mesh_count':len(observed),'junctions':rows})
    return {'schema':'vela.bjt_modeling_check.v1','status':'failed' if failures else 'passed',
            'acceptance_sha256':hashlib.sha256(ACCEPTANCE.read_bytes()).hexdigest(),
            'geometry_contract_sha256_lf':acceptance['source_sha256_lf'],
            'maximum_node_relative_error':max_node_error, 'integrals':integral, 'cuts':cuts,'failures':failures}
