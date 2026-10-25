"""Conditional round-to-nearest radii from the actual printed VAR tokens.

The propagated |J| radius is a first-order state-serialization estimate. It
does not bound nonlinear remainder, assembly arithmetic, or material errors.
"""
from decimal import Decimal

import numpy as np

from compare_sentaurus_debug_matrix import VALUE, VARIABLE


def rounding_radius(token):
    value = Decimal(token)
    if not value.is_finite():
        raise ValueError("finite printed number required")
    return float(Decimal(5).scaleb(value.as_tuple().exponent - 1))


def variable_radii(lines, index):
    current = -1
    block = -1
    radii = {}
    for raw in lines:
        line = raw.strip()
        if line.startswith("Newton loop "):
            current += 1
            block = -1
            if current > index:
                break
        if current != index:
            continue
        if VARIABLE.fullmatch(line):
            block += 1
        match = VALUE.fullmatch(line)
        if match and match[1] == "VAR":
            key = (block, int(match[2]))
            if block < 0 or key in radii:
                raise ValueError("invalid VAR block or duplicate value")
            radii[key] = rounding_radius(match[3])
    if not radii:
        raise ValueError("selected native VAR iteration absent")
    return radii


def propagate_state_rounding(matrix, radius):
    radius = np.asarray(radius, dtype=float)
    if radius.shape != (matrix.shape[1],) or np.any(radius < 0) or not np.all(np.isfinite(radius)):
        raise ValueError("finite nonnegative state radii required")
    return np.asarray(abs(matrix) @ radius).ravel()
