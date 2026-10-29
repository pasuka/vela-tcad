"""Bind a directional response to actual coordinates, with explicit uncertainty.

All arithmetic is exact for binary64 inputs and decimal strings. A printed
matrix action is not an autonomous Jacobian, and an unresolved response is not
a passing derivative check.
"""
from fractions import Fraction


def number(value):
    if isinstance(value, bool):
        raise ValueError("Finite numeric coordinates required")
    try:
        return Fraction(value)
    except (ValueError, TypeError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError("Finite numeric coordinates required") from exc


def bind_increment(reference, current, reference_id, current_reference_id):
    if not reference_id or reference_id != current_reference_id:
        raise ValueError("The same reference state is required")
    if not reference or reference.keys() != current.keys():
        raise ValueError("Complete matching coordinate maps required")
    return {k: number(current[k]) - number(reference[k]) for k in reference}


def linear_action(rows, increment):
    delta = {k: number(v) for k, v in increment.items()}
    result = {}
    for row, columns in rows.items():
        if not columns.keys() <= delta.keys():
            raise ValueError("Jacobian references an absent coordinate")
        result[row] = sum((number(v) * delta[k] for k, v in columns.items()), Fraction())
    return result


def compare_response(residual_change, action, residual_bound=0, action_bound=0,
                     tolerance="1e-5"):
    r, j, rb, jb, tol = map(number, (residual_change, action, residual_bound,
                                     action_bound, tolerance))
    if rb < 0 or jb < 0 or not 0 <= tol < 1:
        raise ValueError("Nonnegative uncertainty and tolerance in [0,1) required")
    relative = lambda a, b: abs(a-b)/(abs(a)+abs(b)) if a or b else Fraction()
    lo, hi, jl, jh = r-rb, r+rb, j-jb, j+jb
    values = [relative(a, b) for a in (lo, hi) for b in (jl, jh)]
    minimum = Fraction() if max(lo, jl) <= min(hi, jh) else min(values)
    maximum = max(values)
    if r == j == rb == jb == 0:
        status = "zero"
    elif (rb > 0 and abs(r) <= rb) or (jb > 0 and abs(j) <= jb):
        status = "unresolved"
    elif maximum <= tol:
        status = "pass"
    elif minimum > tol:
        status = "fail"
    else:
        status = "unresolved"
    return {"status": status, "relative_error": relative(r, j),
            "minimum_error": minimum, "maximum_error": maximum}


def endpoint_intersection(left, right, left_budget, right_budget, action,
                          absolute_tolerance="1e-8", directional_tolerance="1e-5"):
    """Intersect endpoint absolute gates with J(x0)*(x_right-x_left).

    No symmetric steps are assumed. The caller must independently validate
    that this Jacobian action describes the chosen finite chord.
    """
    l, r, lb, rb, j, a, t = map(number, (left, right, left_budget, right_budget,
                                       action, absolute_tolerance, directional_tolerance))
    if lb < 0 or rb < 0 or a < 0 or not 0 <= t < 1:
        raise ValueError("Nonnegative budgets and valid tolerances required")
    radius = a*(lb+rb)
    absolute_interval = (r-l-radius, r-l+radius)
    required = tuple(sorted((j*(1-t)/(1+t), j*(1+t)/(1-t))))
    return {"absolute_change_interval": absolute_interval,
            "direction_change_interval": required,
            "overlap": max(absolute_interval[0], required[0]) <=
                       min(absolute_interval[1], required[1])}
