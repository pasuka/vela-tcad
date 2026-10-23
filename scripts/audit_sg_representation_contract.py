"""Audit current SG values and increments without introducing a fitted offset.

Inputs are already mapped, equally normalized row values. ``compatibility``
and ``canonical`` must be evaluated independently at each current state.
Anchoring is a rejected-or-tested diagnostic, never an automatic correction.
The interval result tests a specified secant/Jv contract, not differentiability
of a rounded function or correctness of the supplied Jacobian.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from pathlib import Path


def number(value):
    # Preserve binary64 inputs exactly; decimal strings retain their low bits.
    if isinstance(value, bool):
        raise ValueError("A numeric operand is required")
    try:
        return Fraction(value)
    except (ValueError, TypeError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError("Finite numeric operands are required") from exc


def text(value):
    with localcontext() as ctx:
        ctx.prec = 90
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def relative(a, b):
    return abs(a - b) / (abs(a) + abs(b)) if a or b else Fraction(0)


def directional_interval(jv, tolerance):
    """Exact interval satisfying |fd-jv|/(|fd|+|jv|) <= tolerance."""
    jv, tolerance = number(jv), number(tolerance)
    if not 0 <= tolerance < 1:
        raise ValueError("Directional tolerance must be in [0, 1)")
    return sorted((jv * (1-tolerance)/(1+tolerance),
                   jv * (1+tolerance)/(1-tolerance)))


def audit_row(samples, jv, absolute_tolerance="1e-8", directional_tolerance="1e-5"):
    a, t, jv = map(number, (absolute_tolerance, directional_tolerance, jv))
    if a < 0 or not 0 <= t < 1:
        raise ValueError("Invalid tolerances")
    states = {}
    for sample in samples:
        h, c, m, b = (number(sample[k]) for k in
                      ("h", "compatibility", "canonical", "budget"))
        if h in states or b < 0:
            raise ValueError("Unique steps and nonnegative budgets required")
        states[h] = (c, m, b)
    if 0 not in states or len(states) < 3:
        raise ValueError("Reference and signed nonzero controls required")
    if any(-h not in states for h in states):
        raise ValueError("Every step requires its opposite control")
    c0, m0, _ = states[0]
    gap0 = m0-c0
    absolute = []
    for h, (c, m, b) in sorted(states.items()):
        gap = m-c
        # A(x;r)=C(r)+M(x)-M(r). Replacing r by x changes A by gap(r)-gap(x).
        drift = gap-gap0
        absolute.append({"h": text(h), "representation_gap": text(gap),
            "anchor_drift": text(drift),
            "canonical_absolute_pass": abs(gap) <= a*b,
            "anchored_absolute_pass": abs(drift) <= a*b,
            "rebase_invariant_exact": drift == 0,
            "budget": text(b)})
    directions = []
    for h in sorted(h for h in states if h > 0):
        cp, mp, bp = states[h]; cm, mm, bm = states[-h]
        fc, fm = (cp-cm)/(2*h), (mp-mm)/(2*h)
        change = ((mp-cp)-(mm-cm))/(2*h)
        assert fm == fc+change
        # Every pair passing both absolute endpoint gates must have a secant
        # inside this interval. Intersect with the exact symmetric Jv gate.
        radius = a*(bp+bm)/(2*h)
        allowed = (fc-radius, fc+radius)
        required = directional_interval(jv, t)
        feasible = max(allowed[0], required[0]) <= min(allowed[1], required[1])
        directions.append({"h": text(h), "compatibility_fd": text(fc),
            "canonical_fd": text(fm), "jv": text(jv),
            "representation_change_fd": text(change),
            "compatibility_direction_error": text(relative(fc,jv)),
            "canonical_direction_error": text(relative(fm,jv)),
            "compatibility_direction_pass": relative(fc,jv) <= t,
            "canonical_direction_pass": relative(fm,jv) <= t,
            "nonzero_direction": bool(fm or fc or jv),
            "endpoint_absolute_and_Jv_intervals_overlap": feasible,
            "absolute_secant_interval": list(map(text,allowed)),
            "Jv_secant_interval": list(map(text,required))})
    return {"states": absolute, "directions": directions,
            "reference_offset_is_diagnostic_only": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    result = [{"id": row["id"], **audit_row(row["samples"], row["jv"])}
              for row in data["rows"]]
    args.output.write_text(json.dumps({"rows": result}, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
