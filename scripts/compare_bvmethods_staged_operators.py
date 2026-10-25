#!/usr/bin/env python3
"""Compare same-state nonlinear Poisson operators and reduce DD constraints.

The Poisson block is differentiated at fixed quasi-Fermi levels. Full coupled
DD entries are NOT compared across density/QF coordinates here; only equivalent
constraint elimination and dimensions are checked for that system.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from compare_sentaurus_debug_matrix import iterations, super_iterations
from constraint_elimination import eliminate_constraints


def read_matrix(path, size):
    d = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
    return sparse.coo_matrix((d['value'], (d['row'].astype(int), d['col'].astype(int))),
                             shape=(size, size)).tocsr()


def row_l2(matrix):
    return np.sqrt(np.asarray(matrix.multiply(matrix).sum(axis=1)).ravel())


def poisson_metrics(a, b, ar, br, scale, reference_scale):
    an = sparse.diags(1/a.diagonal()) @ a
    arn = sparse.diags(1/ar.diagonal()) @ ar
    errors = row_l2(an-arn) / row_l2(arn)
    bv = b/a.diagonal()*scale
    br = np.asarray(br)/ar.diagonal()*reference_scale
    reference_norm = np.linalg.norm(br)
    return {'dimension': a.shape[0],
            'matrix_max_relative_row_L2': float(max(errors)),
            'matrix_worst_node': int(np.argmax(errors)),
            'rhs_diagonal_normalized_difference_L2_V': float(np.linalg.norm(bv-br)),
            'rhs_reference_L2_V': float(reference_norm),
            'rhs_relative_L2': float(np.linalg.norm(bv-br)/reference_norm) if reference_norm else None}


def compare(folder, debug, super_debug, candidate_potential_scale):
    with debug.open() as f:
        trace = list(iterations(f))
    native = list(super_iterations(super_debug.read_text()))
    first_dd = next(i for i, item in enumerate(trace) if len(item['names']) > 1)
    n = len(trace[0]['variables'])
    result = {'poisson': {}}
    for path in sorted(folder.glob('P*_matrix.csv')):
        index = int(path.name.split('_')[0][1:])
        if index >= first_dd or len(native[index]['rhs']) != n:
            raise ValueError('expected pure nonlinear Poisson reference')
        a = read_matrix(path, n)
        b = np.genfromtxt(folder/f'P{index}_rhs.csv', delimiter=',', names=True)['rhs']
        item = native[index]
        entries = item['matrix']
        ar = sparse.coo_matrix((list(entries.values()),
             ([r for r, _ in entries], [c for _, c in entries])), shape=(n, n)).tocsr()
        # Convert diagonal-normalized Newton RHS to volts. Row scaling and
        # packed potential scale must not be mistaken for physical mismatch.
        result['poisson'][str(index)] = poisson_metrics(
            a, b, ar, item['rhs'], candidate_potential_scale, .025851995266484913)
    full_size = 3*n+1
    a = read_matrix(folder/'D0_matrix.csv', full_size)
    b = np.genfromtxt(folder/'D0_rhs.csv', delimiter=',', names=True)['rhs']
    dd = trace[first_dd]
    blocks = {name: block for block, name in enumerate(dd['names'])}
    free_e = {node for block, node in dd['variables'] if block == blocks['eDensity']}
    free_h = {node for block, node in dd['variables'] if block == blocks['hDensity']}
    eliminated = [n+i for i in range(n) if i not in free_e]
    eliminated += [2*n+i for i in range(n) if i not in free_h]
    reduced, rb, kept, recover = eliminate_constraints(a, b, eliminated)
    v = np.sin(np.arange(len(kept))+.37)
    lifted = recover(v)
    full_defect = a@lifted-b
    reduced_defect = reduced@v-rb
    denominator = np.linalg.norm(abs(a)@abs(lifted)+abs(b))
    # Keep the numeric RHS correction even when prescribed rows are tiny.
    result['coupled_dimensions'] = {
        'vela_device': 3*n, 'vela_augmented': full_size,
        'native_labeled': len(dd['variables']),
        'native_labeled_blocks': {name: sum(k[0] == block for k in dd['variables'])
                                 for block, name in enumerate(dd['names'])},
        'native_internal': len(native[first_dd]['rhs']),
        'eliminated_carrier_constraints': len(eliminated),
        'vela_after_exact_elimination': len(kept),
        'matches_native_internal_dimension': len(kept) == len(native[first_dd]['rhs']),
        'recovery_constraint_relative_L2': float(np.linalg.norm(full_defect[eliminated])/denominator),
        'reduced_operator_identity_relative_L2': float(np.linalg.norm(full_defect[kept]-reduced_defect)/denominator),
        'numeric_matrix_equivalence_claimed': False}
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--vela-dir', type=Path, required=True)
    ap.add_argument('--native-debug', type=Path, required=True)
    ap.add_argument('--native-super-debug', type=Path, required=True)
    ap.add_argument('--candidate-potential-scale', type=float, default=.025851999786435535)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    result = compare(a.vela_dir, a.native_debug, a.native_super_debug, a.candidate_potential_scale)
    a.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
