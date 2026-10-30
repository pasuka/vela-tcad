"""Bounded, local MaxTransDiff remeshing; deliberately not native snmesh parity."""
from __future__ import annotations
import math

from sentaurus_doping import _rectangles, _windows, _evaluate_m2_doping_for_owner

DEFAULT_MAX_ITERATIONS = 12
DEFAULT_MAX_NODES = 200000


def validate(device_ir):
    from sentaurus_gmsh_backend import GmshBackendError
    control = device_ir['mesh_control']
    sizes = {s['name']: s for s in control['refinement_sizes']}
    functions = [f for f in control['refinement_functions'] if f['function'] == 'MaxTransDiff']
    for function in functions:
        placements = [p for p in control['refinement_placements'] if p['refinement'] == function['refinement']]
        size = sizes.get(function['refinement'], {})
        arguments = function.get('arguments', [])
        if (not placements or
                any(p['target_kind'] not in {'region','material','window'} for p in placements) or
                function.get('variable') != 'DopingConcentration' or len(arguments) != 1 or
                type(arguments[0]) not in (int,float) or not math.isfinite(arguments[0]) or arguments[0]<=0 or
                not all(k in size for k in ('min_x','min_y'))):
            raise GmshBackendError('MaxTransDiff requires a placed net doping function with positive Value and directional minima',
                                   {'stage':'capability','capability':'maxtransdiff-subset','function':function})


def _active(device_ir, placement, size, windows, x, y, owner):
    from sentaurus_gmsh_backend import _placement_active_at
    if owner is not None and placement['target_kind'] in {'region','material'}:
        value=owner['region'] if placement['target_kind']=='region' else owner['material']['sde']
        return value==placement['target']
    return _placement_active_at(device_ir,placement,size,windows,x,y)


def analyse(device_ir, coordinates, triangles, values, edge_marks=None, owner=None):
    """Pure edge qualification shared by remeshing and tests (coordinates in um)."""
    from sentaurus_gmsh_backend import _placement_active_at
    control = device_ir['mesh_control']
    sizes = {s['name']: s for s in control['refinement_sizes']}
    windows = {w['name']: w for w in control['windows']}
    functions = [f for f in control['refinement_functions'] if f['function'] == 'MaxTransDiff']
    placements=control['refinement_placements']
    edges = {tuple(sorted((t[i], t[(i+1) % 3]))) for t in triangles for i in range(3)}
    marked_nodes = set()
    stopped = 0
    maximum = 0.0
    for a, b in sorted(edges):
        pa, pb = coordinates[a], coordinates[b]
        x,y=(pa[0]+pb[0])/2,(pa[1]+pb[1])/2
        selected=[f for f in functions if any(p['refinement']==f['refinement'] and
                  _active(device_ir,p,sizes[p['refinement']],windows,x,y,owner) for p in placements)]
        if not selected:continue
        threshold=min(float(f['arguments'][0]) for f in selected)
        defect = abs(math.asinh(values[a]) - math.asinh(values[b]))
        maximum = max(maximum, defect)
        if defect <= threshold:
            continue
        x, y = (pa[0]+pb[0])/2, (pa[1]+pb[1])/2
        active = [sizes[p['refinement']] for p in control['refinement_placements']
                  if _active(device_ir,p,sizes[p['refinement']],windows,x,y,owner)]
        min_x = min(float(s['min_x']) for s in active if 'min_x' in s)
        min_y = min(float(s['min_y']) for s in active if 'min_y' in s)
        dx, dy = abs(pa[0]-pb[0]), abs(pa[1]-pb[1])
        if dx <= min_x and dy <= min_y:
            stopped += 1
            continue
        if edge_marks is not None:
            edge_marks.add((a, b))
        marked_nodes.update((a, b))
    return marked_nodes, {'edges': len(edges), 'max_transformed_difference': maximum,
                     'minimum_stopped_edges': stopped, 'refineable_nodes': len(marked_nodes)}


def split_triangles(coordinates, triangles, marked):
    """Conforming red/green midpoint split: each shared edge gets one vertex."""
    points = dict(coordinates)
    midpoint = {}
    next_tag = max(points)+1
    for edge in sorted(marked):
        a,b = edge
        points[next_tag] = tuple((x+y)/2 for x,y in zip(points[a], points[b]))
        midpoint[edge] = next_tag
        next_tag += 1
    def mid(a,b):
        return midpoint.get(tuple(sorted((a,b))))
    def length2(a,b):
        return sum((x-y)**2 for x,y in zip(points[a],points[b]))
    result = []
    for tri in triangles:
        a,b,c = tri
        ab,bc,ca = mid(a,b),mid(b,c),mid(c,a)
        count = sum(m is not None for m in (ab,bc,ca))
        if count == 0:
            result.append(list(tri))
        elif count == 3:
            result.extend([[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]])
        elif count == 1:
            if bc is not None:
                a,b,c,ab = b,c,a,bc
            elif ca is not None:
                a,b,c,ab = c,a,b,ca
            result.extend([[a,ab,c],[ab,b,c]])
        else:
            # Rotate the shared split corner to a; choose the shorter quad diagonal.
            if ca is None:
                a,b,c = b,c,a
            elif ab is None:
                a,b,c = c,a,b
            ab,ca = mid(a,b),mid(c,a)
            result.append([a,ab,ca])
            if length2(ab,c) <= length2(b,ca):
                result.extend([[ab,b,c],[ab,c,ca]])
            else:
                result.extend([[ab,b,ca],[b,c,ca]])
    return points,result,midpoint


def longest_edge_closure(coordinates, triangles, marked):
    """Split each affected triangle's longest edge to avoid repeated skinny greens."""
    from collections import defaultdict, deque
    incident = defaultdict(list)
    longest = []
    for index,tri in enumerate(triangles):
        edges = [tuple(sorted((tri[i],tri[(i+1)%3]))) for i in range(3)]
        longest.append(max(edges,key=lambda e: (sum((a-b)**2 for a,b in zip(coordinates[e[0]],coordinates[e[1]])),e)))
        for edge in edges:
            incident[edge].append(index)
    result = set(marked)
    queue = deque(sorted(marked))
    while queue:
        for triangle in incident[queue.popleft()]:
            edge = longest[triangle]
            if edge not in result:
                result.add(edge)
                queue.append(edge)
    return result


def _replace_mesh(gmsh, coordinates, triangles, midpoint, classification, lines, edge_surface):
    from collections import defaultdict
    # Preserve CAD entity classification and every curve, including non-contact edges.
    updated_lines = {}
    for curve, edges in lines.items():
        updated_lines[curve] = []
        for a,b in edges:
            m = midpoint.get(tuple(sorted((a,b))))
            if m is None:
                updated_lines[curve].append([a,b])
            else:
                classification[m] = (1,curve)
                updated_lines[curve].extend([[a,m],[m,b]])
    for edge,m in midpoint.items():
        classification.setdefault(m,(2,edge_surface[edge]))
    grouped = defaultdict(list)
    for tag,entity in classification.items():
        grouped[entity].append(tag)
    gmsh.model.mesh.clear()
    for (dimension,entity),tags in sorted(grouped.items()):
        gmsh.model.mesh.addNodes(dimension,entity,tags,[v for tag in tags for v in coordinates[tag]])
    for curve,edges in updated_lines.items():
        gmsh.model.mesh.addElementsByType(curve,1,[],[tag for e in edges for tag in e])
    for surface, cells in triangles.items():
        gmsh.model.mesh.addElementsByType(surface,2,[],[tag for t in cells for tag in t])
    return updated_lines


def remesh(gmsh, device_ir, field_dag, mesh_size, max_iterations=DEFAULT_MAX_ITERATIONS,
           max_nodes=DEFAULT_MAX_NODES, surface_owner=None, integral_relative_target=0.):
    from sentaurus_gmsh_backend import GmshBackendError, _triangles_for_surface, _lines_for_curve
    from sentaurus_integral import enabled,analyse as integral_analyse
    if not math.isfinite(integral_relative_target) or not 0<=integral_relative_target<1:
        raise GmshBackendError('invalid integral relative target')
    integral_active=enabled(device_ir,integral_relative_target)
    if not (any(f['function'] == 'MaxTransDiff' for f in device_ir['mesh_control']['refinement_functions'])
            or device_ir['mesh_control'].get('multibox_placements') or integral_active):
        return None
    validate(device_ir)
    if type(max_iterations) is not int or type(max_nodes) is not int or max_iterations < 0 or max_nodes < 3:
        raise GmshBackendError('invalid MaxTransDiff resource budget')
    report = {'schema': 'vela.gmsh_maxtransdiff.v1', 'status': 'running',
              'transform': 'asinh(N_net / (1 cm^-3))', 'edge_rule': 'absolute endpoint transformed difference > Value',
              'stop_rule': 'both absolute edge components <= smallest applicable directional minima',
              'mapping': 'conforming local red/green midpoint splits of authored BAMG mesh; shared-edge and longest-edge closure',
              'native_snmesh_equivalence': False, 'max_iterations': max_iterations,
              'max_nodes': max_nodes, 'rounds': [],
              'functions': [{'refinement': f['refinement'], 'variable': f['variable'],
                             'threshold': float(f['arguments'][0]),
                             'source_location': f['provenance']['source_location']}
                            for f in device_ir['mesh_control']['refinement_functions']
                            if f['function'] == 'MaxTransDiff']}
    body = device_ir['geometry']['bodies'][0]
    rectangles, windows = _rectangles(device_ir), _windows(device_ir)
    cache = {}
    integral_cache = {}
    classification = None
    lines = None
    surface = None
    for iteration in range(max_iterations+1):
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        coordinates = {int(t): (float(xyz[3*i]), float(xyz[3*i+1]), float(xyz[3*i+2])) for i, t in enumerate(tags)}
        if len(tags) > max_nodes:
            report.update(status='failed', reason='node_budget', observed_nodes=len(tags))
            raise GmshBackendError('MaxTransDiff node budget exhausted', report)
        by_surface={e:[t for _,t in _triangles_for_surface(gmsh,e)] for _,e in gmsh.model.getEntities(2)}
        triangles=[t for cells in by_surface.values() for t in cells]
        marked=set();reports=[]
        for entity,cells in by_surface.items():
            owner=surface_owner[entity] if surface_owner is not None else body
            values={}
            for tag in {n for t in cells for n in t}:
                point=coordinates[tag];key=(owner['id'],point)
                if key not in cache:
                    nd,na,_=_evaluate_m2_doping_for_owner(device_ir,owner,rectangles,windows,*point[:2])
                    cache[key]=nd-na
                values[tag]=cache[key]
            _,metrics=analyse(device_ir,coordinates,cells,values,marked,owner)
            reports.append(metrics)
        interface_marks=set()
        for node in field_dag.get('nodes',[]):
            if node['id'].startswith('maxlenint:') and node['kind']=='MathEvalAniso':
                for curve in node['curves']:
                    for _,edge in _lines_for_curve(gmsh,curve):
                        if math.dist(coordinates[edge[0]],coordinates[edge[1]])>float(node['value']):
                            interface_marks.add(tuple(sorted(edge)))
        marked.update(interface_marks)
        from sentaurus_multibox import marked_edges
        inactive={n['id'].removeprefix('multibox:') for n in field_dag.get('nodes',[])
                  if n['kind']=='inactive-multibox'}
        multibox_marks=marked_edges(device_ir,coordinates,triangles,mesh_size,inactive)
        marked.update(multibox_marks)
        integral_marks=set()
        if integral_active:
            integral_marks,integral_report=integral_analyse(device_ir,coordinates,by_surface,surface_owner,
                                                           integral_relative_target,integral_cache)
            marked.update(integral_marks)
            report['integral_control']=integral_report
        metrics={'edges':sum(m['edges'] for m in reports),
                 'max_transformed_difference':max(m['max_transformed_difference'] for m in reports),
                 'minimum_stopped_edges':sum(m['minimum_stopped_edges'] for m in reports),
                 'refineable_nodes':len({n for e in marked for n in e}), 'maxlenint_marked_edges':len(interface_marks), 'multibox_marked_edges':len(multibox_marks)}
        if integral_active:metrics['integral_marked_edges']=len(integral_marks)
        report['rounds'].append(dict(iteration=iteration,nodes=len(tags),triangles=len(triangles),marked_edges=len(marked),**metrics))
        gmsh.logger.write(f'MaxTransDiff round {iteration}: {len(tags)} nodes; {len(marked)} edges to split')
        if not marked:
            report.update(status='passed',remesh_iterations=iteration)
            return report
        if iteration == max_iterations:
            report.update(status='failed',reason='iteration_budget')
            raise GmshBackendError('MaxTransDiff iteration budget exhausted',report)
        marked = longest_edge_closure(coordinates,triangles,marked)
        report['rounds'][-1]['closure_edges'] = len(marked)
        predicted = len(coordinates)+len(marked)
        if predicted > max_nodes:
            report.update(status='failed',reason='node_budget',predicted_nodes=predicted)
            raise GmshBackendError('MaxTransDiff node budget exhausted',report)
        if classification is None:
            classification = {}
            for dim,entity in gmsh.model.getEntities():
                entity_nodes,_,_ = gmsh.model.mesh.getNodes(dim,entity)
                classification.update({int(tag):(dim,entity) for tag in entity_nodes})
            lines = {curve:[edge for _,edge in _lines_for_curve(gmsh,curve)] for _,curve in gmsh.model.getEntities(1)}
        edge_surface={tuple(sorted((t[i],t[(i+1)%3]))):entity for entity,cells in by_surface.items() for t in cells for i in range(3)}
        split={}
        for entity,cells in by_surface.items():
            points,split[entity],midpoint=split_triangles(coordinates,cells,marked)
        lines=_replace_mesh(gmsh,points,split,midpoint,classification,lines,edge_surface)
        field_dag['nodes'].append({'id':f'maxtransdiff:{iteration}','kind':'conforming-edge-split',
                                   'marked_edges':len(marked),'added_nodes':len(midpoint)})
