"""One-element fixed-state branch probes; no solver defaults or acceptance gates."""
from __future__ import annotations

import numpy as np


def observed_vo_transition(field, threshold=4e5):
    """Empirical T-2022.03-SP2 fixed-state VO weight, not a solver default.

    The 5% half-width and decimal pi are inferred from isolated unit/step
    controls. Apply to nonnegative Eparallel, not the signed electric field.
    Evidence and scope are in bvmethods_element_branch_sg_2026-10-12.md.
    """
    f = np.asarray(field, float)
    if (np.any(~np.isfinite(f)) or np.any(f < 0)
            or not np.isfinite(threshold) or threshold <= 0):
        raise ValueError('nonnegative finite force and positive threshold required')
    u = (f/threshold-0.95)/0.1
    interior = .5-.5*np.cos(3.141592654*np.clip(u, 0., 1.))
    return np.where(u <= 0, 0., np.where(u >= 1, 1., interior))


def observed_vo_alpha_300K(field, carrier, threshold=4e5):
    """Diagnostic alpha blend for the original BVmethods 300 K coefficients."""
    if carrier not in ['e', 'h']:
        raise ValueError('carrier must be e or h')
    f = np.asarray(field, float)
    w = observed_vo_transition(f, threshold)
    safe = np.where(f > 0, f, 1.)
    al, ah, bl, bh = ((7.03e5, 7.03e5, 1.231e6, 1.231e6) if carrier == 'e'
                      else (1.582e6, 6.71e5, 2.036e6, 1.693e6))
    return np.where(f > 0, (1-w)*al*np.exp(-bl/safe)+w*ah*np.exp(-bh/safe), 0.)


def single_element_alpha(unit_source, source, relative_tolerance=1e-12):
    """Read alpha from G/G(alpha=1), requiring uniform three-vertex sources.

    A zero unit source cannot identify alpha and returns None. This method
    applies to one isolated triangle, never to an integrated many-cell source.
    """
    unit, rate = (np.asarray(x, dtype=float) for x in (unit_source, source))
    if (unit.shape != (3,) or rate.shape != unit.shape
            or not np.isfinite(relative_tolerance) or relative_tolerance <= 0
            or not np.all(np.isfinite(unit)) or not np.all(np.isfinite(rate))
            or np.any(unit < 0) or np.any(rate < 0)):
        raise ValueError('three finite nonnegative vertex rates required')
    for x in (unit, rate):
        if np.max(x)-np.min(x) > relative_tolerance*np.max(x):
            raise ValueError('one-cell vertex source must be uniform')
    if np.max(unit) == 0:
        if np.max(rate) > 0:
            raise ValueError('positive generation with zero unit control')
        return None
    return float(np.mean(rate/unit))
