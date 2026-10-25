#!/usr/bin/env python3
"""Compare independent BVMethods Poisson and coupled states by global node ID.

Potentials span all regions. Carrier metrics include every positive native
semiconductor node, including minority populations, without offset fitting.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def merged_scalar(folder, name):
    values = {}
    paths = sorted(Path(folder).glob(name + '_region*.csv'))
    if not paths:
        raise ValueError('missing native field: ' + name)
    for path in paths:
        rows = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
        for row in rows:
            node, value = int(row['node_id']), float(row['component0'])
            if not np.isfinite(value):
                raise ValueError('nonfinite native field')
            if node in values and not np.isclose(value, values[node], rtol=1e-12, atol=1e-15):
                raise ValueError('conflicting regional values at node ' + str(node))
            values[node] = value
    return values


def potential_metrics(candidate, reference):
    ids = np.array(sorted(reference), dtype=int)
    actual = np.array([candidate[i] for i in ids])
    expected = np.array([reference[i] for i in ids])
    if not len(ids) or not np.all(np.isfinite(actual)):
        raise ValueError('finite matching potentials required')
    delta = actual - expected
    return {'nodes': len(ids), 'max_absolute_V': float(max(abs(delta))),
            'rms_V': float(np.sqrt(np.mean(delta**2))),
            'worst_node': int(ids[np.argmax(abs(delta))]),
            'pass_1mV': bool(max(abs(delta)) <= 1e-3)}


def carrier_metrics(candidate, reference, semiconductor_ids):
    ids = np.array(sorted(semiconductor_ids), dtype=int)
    actual = np.array([candidate[i] for i in ids])
    expected = np.array([reference[i] for i in ids])
    if (not len(ids) or not np.all(np.isfinite(actual)) or np.any(actual <= 0)
            or not np.all(np.isfinite(expected)) or np.any(expected <= 0)):
        raise ValueError('all semiconductor populations must be finite and positive')
    error = np.log(actual) - np.log(expected)
    return {'nodes': len(ids), 'max_abs_ln_ratio': float(max(abs(error))),
            'rms_ln_ratio': float(np.sqrt(np.mean(error**2))),
            'max_relative': float(max(abs(actual / expected - 1))),
            'worst_node': int(ids[np.argmax(abs(error))]),
            'pass_1percent': bool(max(abs(error)) <= np.log(1.01))}


def compare(vela, poisson, coupled):
    result = {'independent_solver': json.loads((vela/'summary.json').read_text()),
              'fields': {}}
    for phase, folder in [('poisson', poisson), ('coupled', coupled)]:
        rows = np.atleast_1d(np.genfromtxt(vela/(phase+'.csv'), delimiter=',', names=True))
        if len(set(rows['node'])) != len(rows):
            raise ValueError('duplicate Vela node IDs')
        actual = {name: {int(row['node']): float(row[name]) for row in rows}
                  for name in ['psi_V', 'n_cm3', 'p_cm3']}
        psi = merged_scalar(folder, 'ElectrostaticPotential')
        if set(actual['psi_V']) != set(psi):
            raise ValueError('potential node sets differ')
        native = {name: merged_scalar(folder, field) for name, field in
                  [('n_cm3', 'eDensity'), ('p_cm3', 'hDensity')]}
        # Dielectric TDR fields may contain zero population. Define the mask
        # solely from native n AND p, then require full candidate coverage.
        ids = {i for i, value in native['n_cm3'].items() if value > 0}
        if ids != {i for i, value in native['p_cm3'].items() if value > 0}:
            raise ValueError('native electron/hole semiconductor masks disagree')
        metrics = {'potential': potential_metrics(actual['psi_V'], psi)}
        for name in native:
            metrics[name] = carrier_metrics(actual[name], native[name], ids)
        result['fields'][phase] = metrics
    p = merged_scalar(poisson, 'ElectrostaticPotential')
    c = merged_scalar(coupled, 'ElectrostaticPotential')
    result['native_P_to_DD'] = potential_metrics(c, p)
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--vela-dir', type=Path, required=True)
    ap.add_argument('--native-poisson-fields', type=Path, required=True)
    ap.add_argument('--native-coupled-fields', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    result = compare(a.vela_dir, a.native_poisson_fields, a.native_coupled_fields)
    a.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result['fields'], indent=2))
