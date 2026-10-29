"""Exact diagnostic elimination of independent prescribed Newton constraints.

For A dx = b, retain K and eliminate E using the full constraint equation:
S = A_KK - A_KE A_EE^-1 A_EK; g = b_K - A_KE A_EE^-1 b_E.
This module never changes the production solver or its acceptance thresholds.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse


def eliminate_constraints(matrix, rhs, eliminated):
    matrix = sparse.csc_matrix(matrix)
    rhs = np.asarray(rhs, dtype=float)
    n = matrix.shape[0]
    ids = np.asarray(eliminated, dtype=int)
    if matrix.shape != (n, n) or rhs.shape != (n,):
        raise ValueError('square matrix and matching RHS required')
    if (len(set(ids)) != len(ids) or np.any(ids < 0) or np.any(ids >= n)
            or not np.all(np.isfinite(matrix.data)) or not np.all(np.isfinite(rhs))):
        raise ValueError('invalid constraint IDs or nonfinite system')
    kept = np.setdiff1d(np.arange(n), ids)
    diagonal_block = matrix[ids, :][:, ids]
    diagonal = diagonal_block.diagonal()
    remainder = diagonal_block - sparse.diags(diagonal)
    remainder.eliminate_zeros()
    if remainder.nnz or np.any(diagonal == 0):
        raise ValueError('eliminated constraints must have an invertible diagonal block')
    inverse = 1.0 / diagonal
    coupling = matrix[ids, :][:, kept]
    border = matrix[kept, :][:, ids]
    reduced = (matrix[kept, :][:, kept] - border @ sparse.diags(inverse) @ coupling).tocsc()
    reduced_rhs = rhs[kept] - border @ (inverse * rhs[ids])

    def recover(direction):
        direction = np.asarray(direction)
        if direction.shape != (len(kept),):
            raise ValueError('reduced direction has wrong dimension')
        full = np.empty(n)
        full[kept] = direction
        full[ids] = inverse * (rhs[ids] - coupling @ direction)
        return full

    return reduced, reduced_rhs, kept, recover
