"""Reference/increment and element-vertex audit helpers; no solver gates."""
from __future__ import annotations
import numpy as np


def qf_drop_paths(reference, increment, edges):
    """Compare actual double endpoint anchoring with a wider direct drop.

    The wide result is a rounding diagnostic, not exact arithmetic or a new
    production transport rule. Callers must record longdouble's precision.
    """
    r, d = np.asarray(reference, float), np.asarray(increment, float)
    e = np.asarray(edges)
    if (r.ndim != 1 or r.shape != d.shape or e.ndim != 2 or e.shape[1] != 2
            or not np.issubdtype(e.dtype, np.integer) or np.any(e < 0)
            or np.any(e >= len(r)) or np.any(~np.isfinite(r))
            or np.any(~np.isfinite(d))):
        raise ValueError('finite node fields and valid integer endpoint pairs required')
    i, j = e[:, 0], e[:, 1]
    rw, dw = r.astype(np.longdouble), d.astype(np.longdouble)
    physical = r+d
    anchored_j = np.asarray((rw[j]-rw[i])+dw[j], float)
    return dict(absolute=physical[j]-physical[i],
                anchored=anchored_j-d[i],
                wide=np.asarray((rw[j]-rw[i])+(dw[j]-dw[i]), float))


def integrate_element_vertices(generation, measure, node_ids, element_ids,
                               local_ids, node_count):
    """Integrate every labeled element vertex once, without renormalization."""
    g, m = np.asarray(generation, float), np.asarray(measure, float)
    n, e, l = map(np.asarray, (node_ids, element_ids, local_ids))
    if (g.ndim != 1 or any(x.shape != g.shape for x in (m,n,e,l))
            or node_count <= 0 or np.any(~np.isfinite(g))
            or np.any(~np.isfinite(m)) or np.any(m < 0)
            or any(not np.issubdtype(x.dtype,np.integer) for x in (n,e,l))
            or np.any(n < 0) or np.any(n >= node_count)
            or np.any(e < 0) or np.any(l < 0)
            or len(set(zip(e.tolist(),l.tolist()))) != len(g)):
        raise ValueError('finite, uniquely labeled support and valid node IDs required')
    total = np.zeros(node_count);volume = np.zeros(node_count)
    np.add.at(total,n,g*m);np.add.at(volume,n,m)
    return total,volume
