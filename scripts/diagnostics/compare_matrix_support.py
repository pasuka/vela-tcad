"""Compare aligned Jacobian storage and numerical support without altering inputs."""
from __future__ import annotations

import numpy as np
from scipy import sparse


def compare_support(left, right, relative_thresholds=(0.0, 1e-14, 1e-12, 1e-10)):
    """Use the common row maximum for both masks; zero means strict nonzero.

    Matrices must already have identical variable ordering, units and constraint
    elimination. Stored explicit zeros are counted separately. No thresholded
    matrix is returned for use by a solver.
    """
    a, b = sparse.csr_matrix(left, copy=True), sparse.csr_matrix(right, copy=True)
    if a.shape != b.shape or a.shape[0] != a.shape[1]:
        raise ValueError("aligned square matrices required")
    for m in (a, b):
        m.sum_duplicates()
        if not np.all(np.isfinite(m.data)):
            raise ValueError("finite matrix coefficients required")
    result = {"dimension": a.shape[0], "left_stored": a.nnz, "right_stored": b.nnz,
              "left_strict_nonzero": int(np.count_nonzero(a.data)),
              "right_strict_nonzero": int(np.count_nonzero(b.data)), "thresholds": []}
    scale = np.maximum(abs(a).max(axis=1).toarray().ravel(),
                       abs(b).max(axis=1).toarray().ravel())
    positions = {}
    for threshold in relative_thresholds:
        if not np.isfinite(threshold) or threshold < 0:
            raise ValueError("finite nonnegative threshold required")
        masks = []
        for m in (a, b):
            c = m.tocoo()
            keep = abs(c.data) > threshold * scale[c.row]
            masks.append(set(zip(c.row[keep].tolist(), c.col[keep].tolist())))
        left_only, right_only = masks[0] - masks[1], masks[1] - masks[0]
        result["thresholds"].append({"relative": threshold, "left": len(masks[0]),
            "right": len(masks[1]), "shared": len(masks[0] & masks[1]),
            "left_only": len(left_only), "right_only": len(right_only)})
        positions[threshold] = {"left_only": sorted(left_only), "right_only": sorted(right_only)}
    return result, positions
