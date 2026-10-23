#!/usr/bin/env python3
"""Audit a neutral continuity ledger without changing solver acceptance.

All residual, flux and source fields must share one row unit. Absolute floors
must be converted to that unit. The caller must bind the ledger to the actual
accepted state and supply flux_abs_sum AFTER merging contributions by edge.
This implements the source-qualified BVmethods profile, not every optional
density/flux qualification supported by NewtonSolver.
"""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def normalized_row_conversion(source_particles_per_m_s_per_row, target_continuity_scale,
                              target_current_density_A_m2_per_internal,
                              target_length_m_per_internal):
    """Convert particle-line residuals into the target normalized row unit.

    TCAD cm/um units have current-density factor 1e4 and length factor 1e-6;
    their product is 0.01. A target C0*D0 alone is not an SI particle scale.
    """
    values = (source_particles_per_m_s_per_row, target_continuity_scale,
              target_current_density_A_m2_per_internal, target_length_m_per_internal)
    if not all(math.isfinite(v) and v > 0 for v in values):
        raise ValueError('invalid physical row conversion')
    return values[0] / (values[1] * values[2] * values[3])


def audit_global_closure(contact_flux, integrated_source, *, source_floor=1e-14,
                         tolerance=0.01):
    """Evaluate the source-qualified Vela closure; flux is signed, not sum |I|."""
    if not all(math.isfinite(v) for v in (contact_flux, integrated_source, source_floor, tolerance)):
        raise ValueError('nonfinite closure input')
    if source_floor < 0 or tolerance <= 0:
        raise ValueError('invalid closure configuration')
    scale = max(abs(contact_flux), abs(integrated_source), source_floor)
    mismatch = contact_flux - integrated_source
    ratio = abs(mismatch)/scale if scale else 0.
    qualified = abs(integrated_source) >= source_floor
    return dict(contact_flux=contact_flux, integrated_source=integrated_source,
                mismatch=mismatch, ratio=ratio, qualified=qualified,
                satisfied=not qualified or ratio <= tolerance)


def audit_row(row, *, eps_row=1e-3, scale_floor=1e-30, min_source_scale=1e-18):
    values = {k: float(row[k]) for k in
              ('residual', 'flux_abs_sum', 'recombination', 'impact')}
    if not all(math.isfinite(v) for v in values.values()):
        raise ValueError('nonfinite row input')
    if values['flux_abs_sum'] < 0:
        raise ValueError('negative absolute flux sum')
    if not (math.isfinite(eps_row) and eps_row > 0 and
            math.isfinite(scale_floor) and scale_floor >= 0 and
            math.isfinite(min_source_scale) and min_source_scale >= 0):
        raise ValueError('invalid row-gate configuration')
    active = row.get('active', True)
    if not isinstance(active, bool):
        if str(active).lower() not in ('true', 'false', '1', '0'):
            raise ValueError('invalid active flag')
        active = str(active).lower() in ('true', '1')
    source = max(abs(values['recombination']), abs(values['impact']))
    scale = max(values['flux_abs_sum'], source, scale_floor)
    qualified = active and source > 0 and source >= min_source_scale
    # A zero-support nonzero residual is invalid evidence, not a passing zero.
    if scale == 0 and values['residual'] != 0:
        raise ValueError('nonzero residual has zero support')
    ratio = abs(values['residual']) / scale if scale else 0.
    return dict(row, source_scale=source, scale=scale, ratio=ratio,
                qualified=qualified, failed=qualified and ratio > eps_row)


def audit_rows(rows, **config):
    ledger = [audit_row(r, **config) for r in rows]
    groups = defaultdict(list)
    for r in ledger:
        groups[(str(r.get('step', '')), r['carrier'])].append(r)
    summary = []
    for (step, carrier), group in sorted(groups.items()):
        qualified = [r for r in group if r['qualified']]
        worst = max(qualified, key=lambda r: r['ratio']) if qualified else None
        summary.append(dict(step=step, carrier=carrier, rows=len(group),
                            qualified=len(qualified), failed=sum(r['failed'] for r in group),
                            maximum_ratio=None if worst is None else worst['ratio'],
                            worst_node=None if worst is None else worst['node']))
    return ledger, summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--eps-row', type=float, default=1e-3)
    p.add_argument('--scale-floor', type=float, default=1e-30)
    p.add_argument('--min-source-scale', type=float, default=1e-18)
    args = p.parse_args()
    with args.input.open(newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    ledger, summary = audit_rows(rows, eps_row=args.eps_row,
                                scale_floor=args.scale_floor,
                                min_source_scale=args.min_source_scale)
    args.output.write_text(json.dumps({'diagnostic_only': True, 'summary': summary,
                                      'ledger': ledger}, indent=2, allow_nan=False)+'\n',
                           encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
