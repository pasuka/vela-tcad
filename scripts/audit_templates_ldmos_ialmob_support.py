#!/usr/bin/env python3
"""Audit native element/node mobility support without changing the physics.

The supplied mesh must be the qualified LDMOS mesh in micrometers. Native
CSV fields come from sentaurus_import. Output directories must be new.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

try:
    from .audit_templates_ldmos_averagebox_node4492 import (
        measure_in_tdr_vertex_order, parse_debug_block,
    )
except ImportError:
    from audit_templates_ldmos_averagebox_node4492 import (
        measure_in_tdr_vertex_order, parse_debug_block,
    )


def area(points: list[tuple[float, float]]) -> float:
    a, b, c = points
    result = abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])) / 2
    if not math.isfinite(result) or result <= 0:
        raise ValueError("triangle must have finite positive area")
    return result


def project_to_nodes(points: dict, cells: list[dict], values: dict) -> dict:
    """Reconstruct this native Plot field using inverse element area weights."""
    numerator: dict = defaultdict(float)
    denominator: dict = defaultdict(float)
    for cell in cells:
        value = values[cell['id']]
        if not math.isfinite(value):
            raise ValueError("nonfinite element field")
        weight = 1 / area([points[n] for n in cell['node_ids']])
        for node in cell['node_ids']:
            numerator[node] += weight * value
            denominator[node] += weight
    return {n: value / denominator[n] for n, value in numerator.items()}


def stats(errors: list[float]) -> dict:
    if not errors or not all(math.isfinite(e) for e in errors):
        raise ValueError("expected nonempty finite error sample")
    absolute = sorted(abs(e) for e in errors)
    position = .95 * (len(absolute)-1)
    low, high = math.floor(position), math.ceil(position)
    return {'count':len(errors), 'p95_absolute_relative_error':
            absolute[low]+(position-low)*(absolute[high]-absolute[low]),
            'maximum_absolute_relative_error':absolute[-1]}


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--element-export', type=Path, required=True)
    parser.add_argument('--node-export', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--measure-coefficients', type=Path)
    parser.add_argument('--bulk-electron', type=Path)
    parser.add_argument('--bulk-hole', type=Path)
    parser.add_argument('--bulk-node-input', type=Path)
    args = parser.parse_args()
    bulk_args = (args.measure_coefficients,args.bulk_electron,args.bulk_hole,args.bulk_node_input)
    if any(bulk_args) and not all(bulk_args):
        parser.error('bulk audit requires all four bulk/measure arguments')
    sources = [args.mesh]
    mesh = json.loads(args.mesh.read_text())
    points = {n['id']:(n['x'],n['y']) for n in mesh['nodes']}
    si = {r['id'] for r in mesh['regions'] if r['material']=='Si'}
    cells = [c for c in mesh['triangles'] if c['region_id'] in si]
    if not cells:
        raise ValueError('mesh has no silicon triangles')
    if all(bulk_args):
        sources.extend(bulk_args)
        measures = parse_debug_block(args.measure_coefficients,'Measure')
        bulk = {s['id'] for s in json.loads(args.bulk_node_input.read_text())['states_SI']}
    args.output.mkdir(parents=True,exist_ok=False)
    report = {'scope':'Native output support and fixed-state bulk audit; not coupled qualification','results':{}}
    for prefix in ('e','h'):
        path = args.element_export/'fields'/f'{prefix}Mobility_region0_cells.csv'
        sources.append(path)
        element = {int(r['cell_id']):float(r['component0']) for r in read_csv(path)}
        path = args.node_export/'fields'/f'{prefix}Mobility_region0.csv'
        sources.append(path)
        native = {int(r['node_id']):float(r['component0']) for r in read_csv(path)}
        prediction = project_to_nodes(points,cells,element)
        if prediction.keys()!=native.keys() or any(v<=0 for v in native.values()):
            raise ValueError('native node support mismatch or nonpositive mobility')
        records = [dict(node_id=n,native=native[n],predicted=prediction[n],
                        relative_error=prediction[n]/native[n]-1) for n in sorted(native)]
        report['results'][prefix] = {'projection':stats([r['relative_error'] for r in records])}
        tables = {'nodes':records}
        if all(bulk_args):
            path = args.bulk_electron if prefix=='e' else args.bulk_hole
            local = {r['id']:r['mobility_m2_per_Vs']*1e4 for r in json.loads(path.read_text())['results']}
            records = []
            for cell in cells:
                i, ns = cell['id'], cell['node_ids']
                if not set(ns)<=bulk:
                    continue
                w = measure_in_tdr_vertex_order(measures[i])
                predicted = sum(local[n]*v for n,v in zip(ns,w))/sum(w)
                records.append(dict(cell_id=i,native=element[i],predicted=predicted,
                                    relative_error=predicted/element[i]-1))
            tables['bulk_cells'] = records
            report['results'][prefix]['bulk'] = stats([r['relative_error'] for r in records])
        for name, records in tables.items():
            with (args.output/f'{prefix}_{name}.csv').open('x',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    report['input_sha256'] = {str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['results'],indent=2))


if __name__=='__main__':
    main()
