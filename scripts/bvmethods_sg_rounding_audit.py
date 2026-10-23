"""Fixed exported-operand SG rounding controls; not native-state error bounds."""
from itertools import product
import math
from analyze_bvmethods_live_sg_source import sg_edge_current


def adjacent_values(value, *, positive=False, fixed=False):
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError('finite operand and positive population required')
    if fixed:
        return (value,)
    values = (math.nextafter(value, -math.inf), value,
              math.nextafter(value, math.inf))
    if any(not math.isfinite(v) or (positive and v <= 0) for v in values):
        raise ValueError('adjacent values outside finite positive domain')
    return values


def edge_operand_envelope(n0, n1, q0, q1, sign, *, fixed_q0=False,
                          fixed_q1=False, vt=.02585199526648486):
    """Enumerate +/- one representable neighbour per independent operand.

    This deliberately exceeds half-ULP conversion rounding. It audits this
    explicit-QF kernel only, with fixed geometry/mobility; it is not a bound on
    an unknown internal statistical transform or an entire nonlinear solver.
    """
    if sign not in (-1, 1) or not math.isfinite(vt) or vt <= 0:
        raise ValueError('valid carrier sign and thermal voltage required')
    operands = [adjacent_values(n0, positive=True),
                adjacent_values(n1, positive=True),
                adjacent_values(q0, fixed=fixed_q0),
                adjacent_values(q1, fixed=fixed_q1)]
    nominal = sg_edge_current(n0, n1, q1-q0, 1., sign, vt=vt)
    values = [sg_edge_current(a, b, d-c, 1., sign, vt=vt)
              for a, b, c, d in product(*operands)]
    return dict(nominal=nominal, minimum=min(values), maximum=max(values),
                max_deviation=max(abs(v-nominal) for v in values),
                combinations=len(values))
