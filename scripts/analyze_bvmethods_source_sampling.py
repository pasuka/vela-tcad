"""Fixed-state source probes; these helpers do not change solver physics."""
from __future__ import annotations

import numpy as np


def inverse_current_alpha(current_without_mobility, charge=1.602192e-19):
    """Diagnostic alpha that cancels a nonzero two-dimensional J/mu norm."""
    j = np.asarray(current_without_mobility, dtype=float)
    if j.ndim != 2 or j.shape[1] != 2 or not np.all(np.isfinite(j)) or not np.isfinite(charge) or charge <= 0:
        raise ValueError('finite two-dimensional currents and positive charge required')
    norm = np.hypot(j[:, 0], j[:, 1])
    alpha = np.zeros(len(j))
    np.divide(charge, norm, out=alpha, where=norm > 0)
    return alpha


def reconstruct_edge_mobility(directions, current_without_mobility, mobility):
    """Apply two edge mobilities before reconstructing a Cartesian current."""
    d = np.asarray(directions, dtype=float)
    j = np.asarray(current_without_mobility, dtype=float)
    mu = np.asarray(mobility, dtype=float)
    if d.shape != (2, 2) or j.shape != (2,) or mu.shape != (2,) or not all(np.all(np.isfinite(x)) for x in (d, j, mu)) or np.any(mu <= 0):
        raise ValueError('two finite directions, currents and positive mobilities required')
    return np.linalg.solve(d, j*mu)


def branch_mix_fraction(unit_source, step_source, swap_source):
    """Expose branch response and its linear complement without clipping.

    The fraction is undefined where the unit source is nonpositive. A value
    inside [0,1] alone does not prove a microscopic interpolation formula.
    """
    base, step, swap = (np.asarray(x, dtype=float) for x in (unit_source, step_source, swap_source))
    if base.ndim != 1 or base.shape != step.shape or base.shape != swap.shape or not all(np.all(np.isfinite(x)) for x in (base, step, swap)):
        raise ValueError('matching finite node source arrays required')
    fraction = np.full(len(base), np.nan)
    active = base > 0
    fraction[active] = step[active]/base[active]-1
    return dict(fraction=fraction, active=active, complement_error=step+swap-3*base)
