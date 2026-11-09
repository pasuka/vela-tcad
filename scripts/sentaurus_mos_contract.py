"""Independent frozen MOS geometry/doping/reflection acceptance."""
import math
from collections import Counter,defaultdict


def check(bundle):
    nodes={n['vela_id']:n for n in bundle['nodes']}
    regions={r['id']:r['name'] for r in bundle['regions']}
    incidence=defaultdict(set);edges=Counter()
    for t in bundle['triangles']:
        for n in t['node_ids']:incidence[n].add(regions[t['region_id']])
        edges.update(tuple(sorted((t['node_ids'][i],t['node_ids'][(i+1)%3]))) for i in range(3))
    failures=[];max_error=0.
    for n in nodes.values():
        x,y=n['x'],abs(n['y']);expected_d=0.;expected_a=0.
        if 'R.Substrate' in incidence[n['vela_id']]:
            expected_a=1.5e18
            for peak,depth,factor,start in [(6e20,.05,.4,.085),(2e20,.0125,.25,.025)]:
                alpha=math.log(peak/1.5e18)
                expected_d+=peak*math.exp(-alpha*((x/depth)**2+(max(start-y,0)/(factor*depth))**2))
        if 'R.Polygate' in incidence[n['vela_id']]:expected_d+=1e20
        for field,expected in [('DonorConcentration',expected_d),('AcceptorConcentration',expected_a),('DopingConcentration',expected_d-expected_a)]:
            actual=bundle['fields'][field]['values'][n['vela_id']]
            max_error=max(max_error,abs(actual-expected)/max(abs(expected),1.))
    if max_error>1e-6:failures.append('analytic doping')
    index={(round(n['x'],12),round(n['y'],12)):n['vela_id'] for n in nodes.values()}
    reflected={}
    for i,n in nodes.items():
        other=index.get((round(n['x'],12),round(-n['y'],12)))
        if other is None:failures.append('missing mirror node');break
        reflected[i]=other
        for field in list(bundle['fields'].values())+list(bundle['fields']['species'].values()):
            if 'values' in field and field['values'][i]!=field['values'][other]:failures.append('mirror field parity')
    mirror=bundle['provenance'].get('mirror')
    if not mirror:failures.append('missing mirror provenance')
    if len(index)!=len(nodes):failures.append('duplicate merged axis node')
    triangles=Counter((regions[t['region_id']],tuple(sorted(t['node_ids']))) for t in bundle['triangles'])
    if len(reflected)==len(nodes):
        for t in bundle['triangles']:
            region=regions[t['region_id']]
            mirrored=region[:-9] if region.endswith('_mirrored') else region+'_mirrored' if region in {'R.Spacer','R.PolyReox'} else region
            if triangles[(mirrored,tuple(sorted(reflected[i] for i in t['node_ids'])))]!=1:
                failures.append('mirror triangle identity');break
    max_chord=0.;arc_edges=0
    for (a,b),count in edges.items():
        if count!=1:continue
        pa,pb=nodes[a],nodes[b]
        def in_arc(p):return -.1012-1e-10<=p['x']<=-.0512+1e-10 and .035-1e-10<=abs(p['y'])<=.085+1e-10
        if not(in_arc(pa) and in_arc(pb)):continue
        if pa['y']*pb['y']<0:continue
        cx,cy=-.0512,.035
        mx,my=(pa['x']+pb['x'])/2,(abs(pa['y'])+abs(pb['y']))/2
        defect=.05-math.hypot(mx-cx,my-cy)
        if defect<-1e-8:failures.append('fillet circle geometry')
        max_chord=max(max_chord,defect);arc_edges+=1
    if arc_edges<2 or max_chord>1e-5+1e-12:failures.append('chord error')
    return {'schema':'vela.mos_modeling_check.v1','status':'failed' if failures else 'passed',
            'maximum_node_relative_error':max_error,'maximum_chord_error_um':max_chord,'arc_edge_count':arc_edges,
            'node_count':len(nodes),'triangle_count':len(bundle['triangles']),'failures':sorted(set(failures))}
