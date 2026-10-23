"""Private-assembly controls for a parabolic Fermi density integral.

This module is not connected to the solver. Supplied exponentials and logarithms
keep transcendental evaluation distinct from rounded arithmetic. The rounding
model covers normal operands with an unbounded exponent, not hardware underflow.
"""
from fractions import Fraction as F

from bvmethods_sg_argument_audit import round_significand as rounding


def integral_from_exponentials(weights, exp_minus_nodes, exp_minus_eta, normalizer):
    """Return density, its eta derivative, and unnormalized partial sums.

    Weights include the generalized Laguerre x**(1/2) measure. Normalize only
    after summing, and retain the derivative product order f*f, then exp(-eta),
    then weight. Caller-provided transcendental results may be autonomous or
    captured controls; this function makes neither claim about their origin.
    """
    weights, exponentials = tuple(map(F, weights)), tuple(map(F, exp_minus_nodes))
    exponent, normalizer = F(exp_minus_eta), F(normalizer)
    if (not weights or len(weights) != len(exponentials) or exponent <= 0
            or normalizer <= 0 or min(weights) <= 0 or min(exponentials) < 0):
        raise ValueError('positive weights, exponent and normalizer required; lengths must match')
    value = derivative = F(0)
    partials = []
    for weight, node_exp in zip(weights, exponentials):
        occupation = rounding(1 / rounding(node_exp + exponent))
        value = rounding(value + rounding(weight * occupation))
        term = rounding(rounding(rounding(occupation * occupation) * exponent) * weight)
        derivative = rounding(derivative + term)
        partials.append((value, derivative))
    return rounding(normalizer * value), rounding(normalizer * derivative), tuple(partials)


def valley_sum(values, derivatives, weights):
    """Accumulate valleys in the supplied order without renormalizing weights."""
    values, derivatives, weights = (tuple(map(F, x)) for x in (values, derivatives, weights))
    if (not weights or len(values) != len(weights) or len(derivatives) != len(weights)
            or min(weights) < 0):
        raise ValueError('nonnegative weights and matching nonempty valley arrays required')
    value = derivative = F(0)
    for v, d, w in zip(values, derivatives, weights):
        value = rounding(value + rounding(w * v))
        derivative = rounding(derivative + rounding(w * d))
    return value, derivative


def eta_step(eta, value, derivative, log_value_over_target, step_limit, tolerance):
    """One log-density Newton step with explicit acceptance operands.

    The logarithm input must be R80(value/target); it is evaluated externally.
    Acceptance uses the clipped step and the updated eta, with tolerance added
    to its absolute scale. A successful step is not a full nonlinear solve.
    """
    eta, value, derivative, logarithm, limit, tolerance = map(
        F, (eta, value, derivative, log_value_over_target, step_limit, tolerance))
    if value <= 0 or derivative <= 0 or limit <= 0 or tolerance <= 0:
        raise ValueError('positive value, derivative, step limit and tolerance required')
    raw_step = rounding(rounding(logarithm / derivative) * value)
    step = max(-limit, min(limit, raw_step))
    updated = rounding(eta - step)
    criterion = rounding(abs(step) / rounding(abs(updated) + tolerance))
    return {'raw_step': raw_step, 'step': step, 'updated_eta': updated,
            'criterion': criterion, 'accepted': tolerance > criterion}
