#!/usr/bin/env python3
"""Summarize private matched-state source/transport precision diagnostics."""
import argparse
import csv
import json
import math
from pathlib import Path


def relative_error(actual, reference):
    if len(actual) != len(reference) or not actual:
        raise ValueError("nonempty matching vectors required")
    numerator = math.hypot(*(a-b for a,b in zip(actual,reference)))
    denominator = math.hypot(*reference)
    return numerator/denominator if denominator else (0.0 if not numerator else None)


def summarize(rows):
    required={'step','carrier','column_block','row_node','column_node','J_full',
              'J_difference','fd_total','fd_term_sum','fd_recomb'}
    if not rows or any(not required.issubset(r) for r in rows):
        raise ValueError("missing diagnostic fields")
    if any(not math.isfinite(float(r[key])) for r in rows for key in required):
        raise ValueError("nonfinite diagnostic field")
    result = []
    for step in sorted({float(r['step']) for r in rows}):
        for carrier in [1,2]:
            for block in [0,1,2]:
                selected=[r for r in rows if float(r['step'])==step and int(r['carrier'])==carrier and int(r['column_block'])==block]
                if not selected:
                    raise ValueError("missing diagnostic block")
                get=lambda key:[float(r[key]) for r in selected]
                diagonal=[r for r in selected if r['row_node']==r['column_node']]
                result.append({'step':step,'carrier':carrier,'column_block':block,
                    'full_vs_residual_fd_relative':relative_error(get('J_full'),get('fd_total')),
                    'full_vs_term_fd_relative':relative_error(get('J_full'),get('fd_term_sum')),
                    'source_difference_vs_direct_fd_relative':relative_error(get('J_difference'),get('fd_recomb')),
                    'source_fd_norm':math.hypot(*get('fd_recomb')),
                    'lost_nonzero_source_diagonals':[int(r['row_node']) for r in diagonal if float(r['J_difference'])==0 and float(r['fd_recomb'])!=0],
                    'focus790':[dict(r,source_over_full_ulp=abs(float(r['fd_recomb']))/math.ulp(float(r['J_full'])) if float(r['J_full']) else None) for r in diagonal if int(r['row_node'])==790]})
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args()
    with args.input.open(newline='',encoding='utf-8') as stream:
        result=summarize(list(csv.DictReader(stream)))
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__=='__main__':main()
