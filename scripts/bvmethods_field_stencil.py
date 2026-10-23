"""Private fixed-state field/BTBT audit helpers; no production model defaults.

An identified stencil must be validated on independent potential increments.
These helpers alone do not qualify its locality, geometry, or native origin.
"""
import numpy as np


def apply_field_stencil(node_count, entries, potential_increment):
    """Apply dE/dpsi in V/cm/V to a physical potential increment in V."""
    delta = np.asarray(potential_increment, float)
    if node_count <= 0 or delta.shape != (node_count,) or not np.all(np.isfinite(delta)):
        raise ValueError('finite potential vector required')
    result = np.zeros((node_count, 2))
    seen = set()
    for row, col, vector in entries:
        a = np.asarray(vector, float)
        if (row, col) in seen or min(row, col) < 0 or max(row, col) >= node_count or a.shape != (2,) or not np.all(np.isfinite(a)):
            raise ValueError('invalid or duplicate stencil entry')
        seen.add((row, col))
        result[row] += a*delta[col]
    return result


def e2_potential_derivative(entries, field_partials):
    """Compose a fixed-box d(G V)/dE with the identified dE/dpsi.

    Returns sparse (row, column, derivative) entries. Field partial units are
    the caller's explicit source units divided by V/cm.
    """
    partials = np.asarray(field_partials, float)
    if partials.ndim != 2 or partials.shape[1] != 2 or not np.all(np.isfinite(partials)):
        raise ValueError('finite two-component field partials required')
    records = list(entries)
    apply_field_stencil(len(partials), records, np.zeros(len(partials)))
    return [(row, col, float(partials[row]@np.asarray(vector))) for row, col, vector in records]
