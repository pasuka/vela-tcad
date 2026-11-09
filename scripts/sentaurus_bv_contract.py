"""Independent BVmethods source contract (um and cm^-3), frozen for P5."""
import math
from collections import Counter, defaultdict


def check(bundle):
    nodes = {n['vela_id']: n for n in bundle['nodes']}
    regions = {r['id']: r['name'] for r in bundle['regions']}
    incidence, edges = defaultdict(set), Counter()
    failures = []
    for triangle in bundle['triangles']:
        ids = triangle['node_ids']
        for i in ids:
            incidence[i].add(regions[triangle['region_id']])
        edges.update(tuple(sorted((ids[i],ids[(i+1)%3]))) for i in range(3))
        q = [nodes[i] for i in ids]
        if ((q[1]['x']-q[0]['x'])*(q[2]['y']-q[0]['y'])-
                (q[1]['y']-q[0]['y'])*(q[2]['x']-q[0]['x'])) <= 0:
            failures.append('triangle orientation')
        if any(-.065+1e-8 < n['x'] < .065-1e-8 and -.202+1e-8 < n['y'] < -.002-1e-8 for n in q):
            failures.append('removed gate volume')
    if any('polygate' in name.lower() for name in regions.values()):
        failures.append('removed gate region')
    maximum = 0.
    for i,n in nodes.items():
        donor = acceptor = 0.
        if 'R.Substrate' in incidence[i]:
            acceptor = 1e18
            for peak, depth, factor, start in [(5e20,.1,.4,.165),(1e19,.05,.8,.065)]:
                sigma = depth / math.sqrt(2*math.log(peak/1e18))
                # Analytic convolution of a finite unit top hat with erf tails.
                u = abs(n['x'])-start
                weight = .5*(math.erf(u/(factor*sigma))-math.erf((abs(n['x'])-1.13)/(factor*sigma)))
                donor += peak*math.exp(-.5*(n['y']/sigma)**2)*weight
        for field, expected in [('DonorConcentration',donor),('AcceptorConcentration',acceptor),('DopingConcentration',donor-acceptor)]:
            actual = bundle['fields'][field]['values'][i]
            maximum = max(maximum,abs(actual-expected)/max(abs(expected),1.))
    if maximum > 1e-6:
        failures.append('analytic doping')
    index = {(round(n['x'],12),round(n['y'],12)):i for i,n in nodes.items()}
    if len(index) != len(nodes):
        failures.append('duplicate axis node')
    reflected = {}
    for i,n in nodes.items():
        other = index.get((round(-n['x'],12),round(n['y'],12)))
        if other is None:
            failures.append('missing mirror node')
            break
        reflected[i] = other
        for field in list(bundle['fields'].values())+list(bundle['fields']['species'].values()):
            if 'values' in field and field['values'][i] != field['values'][other]:
                failures.append('mirror field parity')
    mirror = bundle['provenance'].get('mirror',{})
    if mirror.get('axis') != 'x':
        failures.append('mirror axis')
    triangles = Counter((regions[t['region_id']],tuple(sorted(t['node_ids']))) for t in bundle['triangles'])
    if len(reflected) == len(nodes):
        for t in bundle['triangles']:
            name = regions[t['region_id']]
            other = name[:-9] if name.endswith('_mirrored') else name+'_mirrored' if name in {'R.Spacer','R.PolyReox'} else name
            if triangles[other,tuple(sorted(reflected[i] for i in t['node_ids']))] != 1:
                failures.append('mirror triangle identity')
                break
    chord, arc_edges = 0., 0
    for (a,b),count in edges.items():
        if count != 1:
            continue
        pa,pb = nodes[a],nodes[b]
        def in_arc(p):
            return .085-1e-10 <= abs(p['x']) <= .165+1e-10 and -.202-1e-10 <= p['y'] <= -.122+1e-10
        if not (in_arc(pa) and in_arc(pb)) or pa['x']*pb['x'] < 0:
            continue
        for p in (pa,pb):
            if abs(math.hypot(abs(p['x'])-.085,p['y']+.122)-.08) > 1e-8:
                failures.append('fillet radius')
        mx,my = (abs(pa['x'])+abs(pb['x']))/2,(pa['y']+pb['y'])/2
        chord = max(chord,.08-math.hypot(mx-.085,my+.122))
        arc_edges += 1
    if arc_edges < 2 or chord > 1e-5+1e-12:
        failures.append('chord error')
    owners = bundle['provenance']['contact_owners'].get('gate',[])
    if set(owners) != {'R.Gateox','R.PolyReox','R.PolyReox_mirrored'}:
        failures.append('gate adjacency')
    inactive = bundle['provenance'].get('inactive_region_references',[])
    if not any(p['kind']=='doping' and p['name']=='PlaceCD.Gate' and p['target']=='R.Polygate' for p in inactive):
        failures.append('deleted gate profile audit')
    dag = bundle['provenance']['refinement_field_dag']
    gate = next((n for n in dag['nodes'] if n['id']=='multibox:MBPlace.Gate'),{})
    if gate.get('kind') != 'inactive-multibox':
        failures.append('removed gate multibox')
    return {'schema':'vela.bv_modeling_check.v1','status':'failed' if failures else 'passed',
            'maximum_node_relative_error':maximum,'maximum_chord_error_um':chord,
            'arc_edge_count':arc_edges,'node_count':len(nodes),'triangle_count':len(bundle['triangles']),
            'failures':sorted(set(failures))}
