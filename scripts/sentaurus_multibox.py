"""Bounded Rectangle metrics (um), interpolating geometric layer sizes.

At d_n=h_min*(r**n-1)/(r-1), h(d)=h_min+(r-1)*d equals h_min*r**n.
The continuous metric does not reproduce snmesh topology.
"""
import math


def validate(device_ir):
    from sentaurus_gmsh_backend import GmshBackendError
    control = device_ir['mesh_control']
    sizes = control.get('multibox_sizes', [])
    windows = {w['name']: w for w in control['windows']}
    names = [s['name'] for s in sizes]
    if len(names) != len(set(names)):
        raise GmshBackendError('duplicate multibox definition')
    for size in sizes:
        for axis in 'xy':
            lo, hi, ratio = (size.get(k+'_'+axis) for k in ('min','max','ratio'))
            if (any(type(v) not in (int,float) or not math.isfinite(v) for v in (lo,hi,ratio))
                    or not 0 < lo <= hi or 0 < abs(ratio) < 1):
                raise GmshBackendError('invalid multibox sizes or ratios', {'stage':'refinement'})
    placements = control.get('multibox_placements', [])
    if len({p['name'] for p in placements}) != len(placements):
        raise GmshBackendError('duplicate multibox placement')
    for placement in placements:
        window = windows.get(placement['window'], {})
        if (placement['multibox'] not in names or window.get('shape') != 'rectangle'
                or any(a >= b for a,b in zip(window.get('lower_left',[]),window.get('upper_right',[])))):
            raise GmshBackendError('multibox requires a known definition and nonempty Rectangle window',
                                   {'stage':'refinement','placement':placement['name']})


def spacing(size, window, x, y, base):
    result = []
    for index, (axis,q) in enumerate(zip('xy',(x,y))):
        ratio = float(size['ratio_'+axis])
        low, high = float(size['min_'+axis]), float(size['max_'+axis])
        if ratio == 0:
            h = base
        elif abs(ratio) == 1:
            # Live U-2022.12 controls: ratio=1 equals max=max, not min=min.
            h = high
        else:
            lo, hi = window['lower_left'][index], window['upper_right'][index]
            distance = max(0., q-lo if ratio > 0 else hi-q)
            h = min(high, low + (abs(ratio)-1)*distance)
        result.append(h)
    return result


def _samples(device_ir, coordinates, triangles, base):
    control = device_ir['mesh_control']
    sizes = {s['name']:s for s in control.get('multibox_sizes',[])}
    windows = {w['name']:w for w in control['windows']}
    edges = {tuple(sorted((t[i],t[(i+1)%3]))) for t in triangles for i in range(3)}
    for placement in control.get('multibox_placements',[]):
        size, window = sizes[placement['multibox']], windows[placement['window']]
        x0,y0 = window['lower_left']; x1,y1 = window['upper_right']
        samples = []
        for a,b in edges:
            pa,pb = coordinates[a],coordinates[b]
            x,y = (pa[0]+pb[0])/2,(pa[1]+pb[1])/2
            if not (x0-1e-12 <= x <= x1+1e-12 and y0-1e-12 <= y <= y1+1e-12):
                continue
            h = spacing(size,window,x,y,base)
            ratios = [abs(pa[i]-pb[i])/h[i] if size['ratio_'+axis] else 0.
                      for i,axis in enumerate('xy')]
            samples.append(((a,b),ratios,(x,y),h))
        yield placement, size, window, samples


def marked_edges(device_ir, coordinates, triangles, base, inactive=()):
    return {edge for placement,_,_,samples in _samples(device_ir,coordinates,triangles,base)
            if placement['name'] not in inactive for edge,ratios,_,_ in samples if max(ratios) > 1.35 + 1e-10}


def qualify(device_ir, coordinates, triangles, field_dag):
    cells = [t['node_ids'] for t in triangles]
    dag = {n['id']:n for n in field_dag['nodes']}
    reports, failures = [], []
    for placement,size,window,samples in _samples(device_ir,coordinates,cells,1.):
        record = {'name':placement['name'],'ratio':[size['ratio_x'],size['ratio_y']],
                  'window_um':[window['lower_left'],window['upper_right']],
                  'edge_count':len(samples),'mapping':'continuous geometric-layer interpolation'}
        inactive = dag['multibox:'+placement['name']]['kind']=='inactive-multibox'
        if inactive:
            record.update(active=False, reason='window has no surviving volume',passed=True)
        else:
            maximum = max((max(r) for _,r,_,_ in samples), default=0.)
            record.update(active=True, maximum_normalized_component=maximum,
                          passed=bool(samples) and maximum <= 1.35+1e-9)
            if not record['passed']:
                failures.append('multibox '+placement['name']+' local spacing contract failed')
        reports.append(record)
    return {'placements':reports,'failures':failures}
