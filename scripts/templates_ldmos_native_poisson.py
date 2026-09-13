"""Explicit T-2022.03 LDMOS AverageBox Poisson geometry, in SI.

The audited local slot mapping is specific to this native export profile.
Both material-weighted edge coefficients and silicon charge measures move
together. Transport, recombination volumes and heat geometry are independent.
"""
import csv
import hashlib
import math
from pathlib import Path
from audit_templates_ldmos_averagebox_node4492 import (
    parse_debug_block, measure_in_tdr_vertex_order, coefficient_edge,
)


def native_poisson_geometry(mesh, coordinate_to_metres, debug, export):
    debug, export = Path(debug), Path(export)
    if coordinate_to_metres != 1e-6:
        raise ValueError('Audited native geometry requires mesh coordinates in micrometres')
    nodes, cells = mesh['nodes'], mesh['triangles']
    materials = {r['id']: r['material'] for r in mesh.get('regions', [])} or {0:'Si', 1:'SiO2'}
    with (export/'nodes.csv').open(encoding='utf-8-sig', newline='') as stream:
        native_nodes = list(csv.DictReader(stream))
    with (export/'elements.csv').open(encoding='utf-8-sig', newline='') as stream:
        native_cells = list(csv.DictReader(stream))
    if len(nodes) != len(native_nodes) or len(cells) != len(native_cells):
        raise ValueError('Native Poisson geometry topology size mismatch')
    for i, (node, native) in enumerate(zip(nodes, native_nodes)):
        if node['id'] != i or int(native['id']) != i:
            raise ValueError('Native Poisson node ordering mismatch')
        if any(not math.isfinite(node[k]) or not math.isfinite(float(native[k+'_um'])) or
               abs(node[k]-float(native[k+'_um'])) > 1e-12 for k in ('x', 'y')):
            raise ValueError('Native Poisson node coordinates mismatch')
    for cell, native in zip(cells, native_cells):
        if (cell['id'] != int(native['id']) or len(cell['node_ids']) != 3 or
            cell['node_ids'] != [int(native['node'+str(k)]) for k in range(3)]):
            raise ValueError('Native Poisson triangle connectivity mismatch')
        if native.get('material') != materials.get(cell['region_id']):
            raise ValueError('Native Poisson material identity mismatch')
    measures = parse_debug_block(debug, 'Measure')
    coefficients = parse_debug_block(debug, 'Coefficients')
    ids = {cell['id'] for cell in cells}
    if len(ids) != len(cells) or set(measures) != ids or set(coefficients) != ids:
        raise ValueError('Native Poisson coefficient cell set mismatch')
    area, edges = [0.]*len(nodes), {}
    eps0 = 8.8541878128e-12
    for cell in cells:
        cid, ns, region = cell['id'], cell['node_ids'], cell['region_id']
        if materials.get(region) != ('Si' if region == 0 else 'SiO2'):
            raise ValueError('Audited native Poisson requires silicon region 0 and silicon dioxide elsewhere')
        shares = measure_in_tdr_vertex_order(measures[cid])
        weights = coefficients[cid]
        if len(weights) != 3 or any(not math.isfinite(w) or w < 0. for w in weights):
            raise ValueError('Invalid native AverageBox edge coefficients')
        if region == 0:
            for node, measure in zip(ns, shares):
                area[node] += measure*coordinate_to_metres**2
        for slot, weight in enumerate(weights):
            a, b, _ = coefficient_edge(ns, slot)
            key = tuple(sorted((a, b)))
            edges[key] = edges.get(key, 0.) + eps0*(11.7 if region == 0 else 3.9)*weight
    silicon = {i for c in cells if c['region_id'] == 0 for i in c['node_ids']}
    if any((area[i] <= 0.) if i in silicon else (area[i] != 0.) for i in range(len(nodes))):
        raise ValueError('Native Poisson silicon charge support mismatch')
    sources = [Path(__file__), debug, export/'nodes.csv', export/'elements.csv']
    provenance = dict(policy='native_averagebox_poisson_v1',
        relative_permittivity={str(i):11.7 if m=='Si' else 3.9 for i,m in materials.items()}, measure_unit='um2',
        sources_sha256={str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    return area, edges, provenance


def apply_native_poisson(cfg, mesh, debug, export):
    """Return a new input; preserve every non-Poisson field and the input itself."""
    area, edges, provenance = native_poisson_geometry(mesh, cfg['coordinate_to_metres'], debug, export)
    if {tuple(e['nodes']) for e in cfg['edge_geometry']} != set(edges):
        raise ValueError('Native Poisson edge set mismatch')
    result = dict(cfg)
    result['silicon_area_m2'] = area
    # The old fallback source volume must not implicitly change with charge area.
    result['recombination_area_m2'] = list(cfg.get('recombination_area_m2', cfg['silicon_area_m2']))
    result['edge_geometry'] = [dict(e, poisson_F_per_m=edges[tuple(e['nodes'])]) for e in cfg['edge_geometry']]
    result['poisson_geometry_provenance'] = provenance
    return result
