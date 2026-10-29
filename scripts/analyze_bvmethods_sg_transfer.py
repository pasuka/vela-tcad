"""Fixed-state audit primitives; none of these functions grants qualification.

The source conversion is for the labeled Sentaurus T-2022.03 BVmethods
normalization: G in cm^-3/s, 2-D box measure in um^2, ND in cm^-3.
It is not a universal conversion for another device or matrix convention.
"""
from __future__ import annotations

import math
import numpy as np


def source_rhs(generation_cm3_s, measure_um2, reference_density_cm3):
    """Convert a positive pair-generation rate to the native RHS magnitude."""
    g, m = np.broadcast_arrays(np.asarray(generation_cm3_s, float),
                                np.asarray(measure_um2, float))
    if (not math.isfinite(reference_density_cm3) or reference_density_cm3 <= 0
            or np.any(~np.isfinite(g)) or np.any(~np.isfinite(m))
            or np.any(g < 0) or np.any(m < 0)):
        raise ValueError("finite nonnegative generation/measure and positive ND required")
    return g * m / (reference_density_cm3 * 1e8)


def subtraction_resolution(full, disabled, ulps=64):
    """Return every difference and its resolution marker; do not filter rows."""
    a, b = np.broadcast_arrays(np.asarray(full, float), np.asarray(disabled, float))
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)) or ulps <= 0:
        raise ValueError("finite RHS values and a positive ULP multiplier required")
    limit = ulps * np.maximum(np.abs(np.spacing(a)), np.abs(np.spacing(b)))
    delta = a-b
    return delta, limit, np.abs(delta) > limit


def protected_qf_inverse(reference, required, active, contacts):
    """Apply an absolute-inverse diagnostic only on noncontact active nodes.

    This deliberately exposes inverse-map sensitivity; it is not the proposed
    relative finite-step Newton recovery and must not be used as one.
    """
    ref = np.asarray(reference, float)
    req = np.asarray(required, float)
    active, contacts = np.asarray(active, bool), np.asarray(contacts, bool)
    if any(x.shape != ref.shape for x in (req, active, contacts)):
        raise ValueError("incompatible node masks or QF arrays")
    if np.any(~np.isfinite(ref)) or np.any(~np.isfinite(req[active & ~contacts])):
        raise ValueError("nonfinite participating QF")
    out = ref.copy()
    out[active & ~contacts] = req[active & ~contacts]
    return out


def block_comparison(full, transport, avalanche, other, native_full,
                     native_transport, native_avalanche, native_other):
    """Keep full and component errors separate, with no arbitrary unit floor."""
    values = [np.asarray(v, float) for v in (full, transport, avalanche, other,
              native_full, native_transport, native_avalanche, native_other)]
    if any(v.shape != values[0].shape or np.any(~np.isfinite(v)) for v in values):
        raise ValueError("incompatible/nonfinite residual components")
    f,t,a,o,n,nt,na,no = values
    norm = float(np.linalg.norm(f))
    closure = float(np.linalg.norm(f-t-a-o))
    return dict(rows=f.size, full_difference_L2=float(np.linalg.norm(f-n)),
                transport_difference_L2=float(np.linalg.norm(t-nt)),
                avalanche_difference_L2=float(np.linalg.norm(a-na)),
                other_source_difference_L2=float(np.linalg.norm(o-no)),
                internal_closure_relative=closure/norm if norm else (0. if closure == 0 else None))
