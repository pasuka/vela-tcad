#!/usr/bin/env python3
"""Audit density updates against subsequent VAR states and saved TDR fields.

Uses the BVmethods 2026-10-09 short probe's documented iteration mapping.
Extrapolated states are deliberately excluded from direct update comparisons.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

from compare_sentaurus_debug_matrix import iterations, scalar_fields


def audit_update(old, delta, actual):
    if not old or set(old) != set(delta) or set(old) != set(actual):
        raise ValueError("density node mappings differ")
    if any(not math.isfinite(v) or v <= 0 for v in old.values()):
        raise ValueError("old densities must be positive and finite")
    bad = [k for k in old if old[k] + delta[k] <= 0]
    reconstructed = {k: old[k] + delta[k] if old[k] + delta[k] > 0 else .01*old[k] for k in old}
    if any(not math.isfinite(v) or v <= 0 for v in actual.values()):
        raise ValueError("actual densities must be positive and finite")
    return {"nodes": len(old), "nonpositive_raw_nodes": bad,
            "nonpositive_raw_count": len(bad),
            "max_reconstruction_relative_error": max(abs(actual[k]/reconstructed[k]-1) for k in old),
            "positive_below_one_percent_count": sum(0 < old[k]+delta[k] < .01*old[k] for k in old),
            "negative_update_ratios": [actual[k]/old[k] for k in bad],
            "max_relative_raw_update": max(abs(delta[k]/old[k]) for k in old),
            "min_actual_over_old": min(actual[k]/old[k] for k in old),
            "max_actual_over_old": max(actual[k]/old[k] for k in old)}


def vector_difference(a, b):
    if set(a) != set(b):
        raise ValueError("RHS mappings differ")
    numerator = math.fsum(abs(a[k]-b[k]) for k in a)
    denominator = math.fsum(abs(a[k])+abs(b[k]) for k in a)
    return numerator/denominator if denominator else 0.0


def independent_rhs(vector):
    """The eleven identical contact labels represent one circuit equation."""
    contacts = [v for (block,node),v in vector.items() if block == 1]
    if contacts and any(not math.isclose(v,contacts[0],abs_tol=0,rel_tol=1e-12) for v in contacts):
        raise ValueError('contact RHS aliases differ')
    return [v for (block,node),v in vector.items() if block != 1] + contacts[:1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--probe-root", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    with (args.probe_root/'native/debug').open(encoding='utf-8') as f:
        native = list(iterations(f))
    if [v['loop'] for v in native] != [0,1,2,3,4,5,6,0,0,0,1,2,0,0]:
        raise ValueError('unexpected diagnostic iteration chronology')
    comparisons = []
    # 8->9 uses the unextrapolated first QS starting state. 12->13 is excluded:
    # the new bias state is extrapolated and is not the previous Newton update.
    for i,j in [(7,8),(8,9),(9,10),(10,11)]:
        result = {'iteration': i, 'destination': f'VAR iteration {j}', 'carriers': {}}
        for block,name in [(2,'electron'),(3,'hole')]:
            old = {k[1]:v for k,v in native[i]['variables'].items() if k[0] == block}
            delta = {k:native[i]['solution'][block,k] for k in old}
            actual = {k:native[j]['variables'][block,k] for k in old}
            result['carriers'][name] = audit_update(old,delta,actual)
        result['rhs2_vs_next_rhs1_relative_l1'] = vector_difference(native[i]['rhs2'],native[j]['rhs1']) if (i,j) in [(9,10),(10,11)] else None
        comparisons.append(result)
    scale = 14824351000.0  # Independently established by the previous TDR audit.
    for i,label in [(7,'zero'),(11,'first'),(13,'second')]:
        result = {'iteration': i, 'destination': f'TDR {label}', 'carriers': {}}
        for block,name,field in [(2,'electron','eDensity'),(3,'hole','hDensity')]:
            old = {k[1]:v for k,v in native[i]['variables'].items() if k[0] == block}
            delta = {k:native[i]['solution'][block,k] for k in old}
            values = scalar_fields(args.probe_root/f'native_export_{label}'/'fields',field)
            actual = {k:values[k]/scale for k in old}
            result['carriers'][name] = audit_update(old,delta,actual)
        comparisons.append(result)
    mesh = json.loads((args.probe_root/'mesh.json').read_text(encoding='utf-8'))
    with (args.probe_root/'doping.csv').open(newline='',encoding='utf-8') as f:
        doping = {int(r['node_id']):r for r in csv.DictReader(f)}
    touched = sorted({k for c in comparisons for v in c['carriers'].values() for k in v['nonpositive_raw_nodes']})
    nodes = {n['id']:n for n in mesh['nodes']}
    result = {'comparisons':comparisons,'modified_density_node_union':touched,
              'modified_density_node_locations':[{'node':k,'x_um':nodes[k]['x'],'y_um':nodes[k]['y'],'donors_cm3':float(doping[k]['donors_cm3']),'acceptors_cm3':float(doping[k]['acceptors_cm3'])} for k in touched],
              'iteration_rhs':[{'index':i,'loop':v['loop'],'rhs1_l2':math.sqrt(math.fsum(x*x for x in independent_rhs(v['rhs1']))), 'rhs2_l2':math.sqrt(math.fsum(x*x for x in independent_rhs(v['rhs2']))), 'rhs1_inf':max(abs(x) for x in independent_rhs(v['rhs1'])), 'rhs2_inf':max(abs(x) for x in independent_rhs(v['rhs2']))} for i,v in enumerate(native) if i>=7],
              'excluded_transition':'12->13 includes bias extrapolation; not a direct update test',
              'native_density_scale_cm3':scale}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps([{'iteration':c['iteration'],'destination':c['destination'],'electron_negative':c['carriers']['electron']['nonpositive_raw_count'],'hole_negative':c['carriers']['hole']['nonpositive_raw_count'],'hole_error':c['carriers']['hole']['max_reconstruction_relative_error']} for c in comparisons]))


if __name__ == '__main__':
    main()
