#!/usr/bin/env python3
"""Integrate reconstructed nodal current fields; never substitute for SG fluxes."""
from __future__ import annotations
import argparse
from pathlib import Path
from pn2d_variants import DEFAULT_ROOT, SPEC, read_json, write_json, tag
from compare_pn2d_variants import node_field
from compare_genius_bjt_transport_fields import read_vtk_point_data


def nodal_section(mesh, vectors, axis, cut_um, depth_um=1.0):
    """Exact line quadrature of the piecewise-linear normal current, in A/um.

    Coordinates are um and vectors A/cm2. An internal mesh edge is integrated
    once even when both adjacent triangles intersect along its entire length.
    This is a reconstruction diagnostic, not the simulator's discrete flux.
    """
    normal = {'x': 0, 'y': 1}[axis]
    tangent = 1-normal
    nodes = mesh['nodes']
    index = {n['id']: i for i, n in enumerate(nodes)}
    xy = [(n['x'], n['y']) for n in nodes]
    segments = {}
    for triangle in mesh['triangles']:
        ids = [index[i] for i in triangle['node_ids']]
        if cut_um < min(xy[i][normal] for i in ids)-1e-12 or cut_um > max(xy[i][normal] for i in ids)+1e-12:
            continue
        points = {}
        for a, b in zip(ids, ids[1:]+ids[:1]):
            da, db = xy[a][normal]-cut_um, xy[b][normal]-cut_um
            if abs(da) < 1e-12:
                points[round(xy[a][tangent], 12)] = (xy[a][tangent], vectors[a][normal])
            if da*db < 0 and abs(da) >= 1e-12 and abs(db) >= 1e-12:
                fraction = -da/(db-da)
                position = xy[a][tangent]+fraction*(xy[b][tangent]-xy[a][tangent])
                value = vectors[a][normal]+fraction*(vectors[b][normal]-vectors[a][normal])
                points[round(position, 12)] = (position, value)
        if len(points) < 2:
            continue
        low, high = min(points), max(points)
        p, q = points[low], points[high]
        segments[(low, high)] = (q[0]-p[0])*(p[1]+q[1])/2
    return sum(segments.values())*1e-8*depth_um


def analyze(root):
    data = read_json(root/'results.json')
    spec = read_json(SPEC)
    contract = read_json(root/'contract.json')
    rank = {m: i for i, m in enumerate(spec['meshes'])}
    preferred = {}
    for case in data['results']:
        if case['status'] != 'pass':
            continue
        if case['id'] not in preferred or rank[case['mesh']] > rank[preferred[case['id']]['mesh']]:
            preferred[case['id']] = case
    records = []
    for case in preferred.values():
        directory = root/case['id']/case['mesh']
        mesh = read_json(directory/'inputs/mesh.json')
        ids = [n['id'] for n in mesh['nodes']]
        for point in case['spatial']['points']:
            bias = point['bias_V']
            branch = 'reverse' if bias < 0 else 'forward'
            _, _, vectors = read_vtk_point_data(directory/point['vtk'])
            fields = directory/'fields'/f'{branch}_{tag(bias)}'/'fields'
            cathode = next(p for p in case['terminal']['metrics'] if p['branch'] == branch and p['bias_V'] == bias and p['contact'] == 'Cathode' and p['component'] == 'total')
            for component, native, vela in [('electron','eCurrentDensity','SentaurusElectronCurrentDensityVector'),
                                           ('hole','hCurrentDensity','SentaurusHoleCurrentDensityVector'),
                                           ('total','TotalCurrentDensity','SentaurusTotalCurrentDensityVector')]:
                reference = node_field(fields/f'{native}_region0.csv', ids, 2)
                for x in contract['conservation']['sections_x_um']:
                    item = {'id':case['id'], 'mesh':case['mesh'], 'bias_V':bias, 'x_um':x, 'component':component}
                    for tool, field in [('reference',reference), ('candidate',vectors[vela])]:
                        current = nodal_section(mesh, field, 'x', x, spec['depth_um'])
                        item[tool+'_nodal_A_per_um'] = current
                        if component == 'total':
                            target = -cathode[tool+'_A_per_um']
                            item[tool+'_minus_port_A_per_um'] = current-target
                            item[tool+'_relative_to_port'] = abs(current-target)/max(abs(target), contract['conservation']['absolute_floor_A_per_um'])
                    item['candidate_SG_A_per_um'] = next(s for s in point['sections'] if s['x_um'] == x)[component+'_A_per_um']
                    records.append(item)
    write_json(root/'nodal_section_diagnostics.json', {'status':'diagnostic_only',
        'method':'Exact piecewise-linear line integration of nodal reconstructed vectors, each shared edge counted once; comparison points are solved states, never bias interpolation.',
        'units':'A/um, positive toward +x; native discrete face flux is not exported',
        'records':records})
    return records


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=DEFAULT_ROOT)
    analyze(parser.parse_args().root.resolve())
