"""Opt-in SLOT diagnostic linear guard; never a nonlinear acceptance test.

Internal coordinates must be potentials divided by potential_scale_V. Accept
|Ax-b|_i <= relative_tolerance*(|A||x|+|b|)_i +
              sum_j |A_ij|*absolute_direction_V/potential_scale_V,
in addition to block normwise and exact singleton-constraint checks. This is
a backward-residual allowance, NOT a forward-error or terminal-current bound.
"""
from pathlib import Path
import math
import numpy as np


def evaluate_coo(n, rows, columns, values, rhs, direction, *,
                 potential_scale_V, relative_tolerance=1e-6,
                 absolute_direction_V=1e-14, normwise_tolerance=1e-12,
                 block_size=None):
    policy = dict(relative_tolerance=relative_tolerance,
                  absolute_direction_V=absolute_direction_V,
                  normwise_tolerance=normwise_tolerance,
                  potential_scale_V=potential_scale_V)
    if any(not math.isfinite(v) or v <= 0 for v in policy.values()):
        raise ValueError('linear guard scales must be positive and finite')
    i, j = np.asarray(rows, int), np.asarray(columns, int)
    v, b, x = (np.asarray(a, dtype=np.longdouble) for a in (values, rhs, direction))
    if (len(i) != len(j) or len(j) != len(v) or b.shape != (n,) or
            x.shape != (n,) or np.any(i < 0) or np.any(i >= n) or
            np.any(j < 0) or np.any(j >= n)):
        raise ValueError('invalid matrix/vector dimensions or coordinates')
    if not all(np.isfinite(a).all() for a in (v, b, x)):
        return dict(accepted=False, reason='nonfinite_input', policy=policy)
    # Matrix Market input must have unique coordinates, as exported by Vela.
    if len(np.unique(i.astype(np.int64)*n+j)) != len(v):
        raise ValueError('duplicate matrix coordinates must be assembled first')
    product = np.zeros(n, dtype=np.longdouble)
    scale, row_sum = abs(b).copy(), np.zeros(n, dtype=np.longdouble)
    np.add.at(product, i, v*x[j])
    np.add.at(scale, i, abs(v)*abs(x[j]))
    np.add.at(row_sum, i, abs(v))
    residual = product-b
    allowance = relative_tolerance*scale + row_sum*absolute_direction_V/potential_scale_V

    def ratios(num, den):
        result = np.zeros_like(num)
        np.divide(num, den, out=result, where=den > 0)
        result[(den == 0) & (num > 0)] = np.inf
        return result

    raw = ratios(abs(residual), scale)
    mixed = ratios(abs(residual), allowance)
    counts = np.bincount(i[v != 0], minlength=n)
    diag = np.zeros(n, dtype=np.longdouble)
    on_diag = i == j
    diag[i[on_diag]] = v[on_diag]
    fixed = (counts == 1) & (diag != 0)
    constraint_V = (float(np.max(abs(x[fixed]-b[fixed]/diag[fixed]))*potential_scale_V)
                    if fixed.any() else 0.)
    width = n if block_size is None else int(block_size)
    if width <= 0 or n % width:
        raise ValueError('block_size must divide matrix dimension')
    block_norm = [float(ratios(np.array([max(abs(residual[k:k+width]))]),
                              np.array([max(scale[k:k+width])]))[0])
                  for k in range(0, n, width)]
    accepted = bool(np.all(mixed <= 1) and max(block_norm) <= normwise_tolerance
                    and constraint_V <= absolute_direction_V)
    # JSON-safe values: nonfinite diagnostics explicitly remain rejected/null.
    def finite(z):
        value = float(z)
        return value if math.isfinite(value) else None
    worst = int(np.argmax(mixed))
    return dict(accepted=accepted, reason='passed' if accepted else 'linear_quality_rejected',
                policy=policy, component_backward=finite(max(raw)),
                mixed_ratio=finite(max(mixed)), block_normwise=[finite(z) for z in block_norm],
                absolute_residual_max=finite(max(abs(residual))),
                raw_bad_count=int(sum(raw > relative_tolerance)),
                mixed_bad_count=int(sum(mixed > 1)),
                absolute_allowance_used_count=int(sum((raw > relative_tolerance) & (mixed <= 1))),
                eliminated_dofs=int(sum(fixed)), constraint_error_V=constraint_V,
                worst_mixed_row=worst, worst_mixed_residual=finite(residual[worst]),
                worst_mixed_allowance=finite(allowance[worst]))


def evaluate_matrix_market(path, rhs, direction, **kwargs):
    with Path(path).open() as f:
        if 'coordinate real general' not in f.readline():
            raise ValueError('expected general real coordinate Matrix Market')
        line = f.readline()
        while line.startswith('%'):
            line = f.readline()
        n, m, count = map(int, line.split())
        a = np.loadtxt(f, ndmin=2)
    if n != m or len(a) != count:
        raise ValueError('invalid square matrix or entry count')
    return evaluate_coo(n, a[:, 0].astype(int)-1, a[:, 1].astype(int)-1,
                        a[:, 2], rhs, direction, **kwargs)
