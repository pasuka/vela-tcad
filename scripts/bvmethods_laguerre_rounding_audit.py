"""Finite-significand Laguerre polynomial and weight dependency controls.

No native tables, root seeds or transcendental implementation are included.
The rounding model covers normal operands with an unbounded exponent.
"""
from fractions import Fraction as F
from bvmethods_sg_argument_audit import round_significand as rounding


def polynomial_pair(order, alpha, x):
    """Return L_n, L_(n-1), and L_n' with separate operation rounding."""
    alpha, x = F(alpha), F(x)
    if not isinstance(order, int) or not 1 <= order <= 512 or alpha <= -1 or x <= 0:
        raise ValueError('positive x, alpha > -1 and order in [1, 512] required')
    previous, current = F(0), F(1)
    for j in range(1, order + 1):
        coefficient = rounding(rounding(2*(j-1) + alpha) + 1)
        coefficient = rounding(coefficient - x)
        numerator = rounding(rounding(coefficient*current)
                             - rounding(rounding(j-1 + alpha)*previous))
        previous, current = current, rounding(numerator/j)
    derivative = rounding(rounding(rounding(order*current)
                                   - rounding(rounding(order+alpha)*previous))/x)
    return current, previous, derivative


def newton_step(order, alpha, x):
    """Keep the polynomial evaluation point distinct from the updated root."""
    x = F(x)
    value, previous, derivative = polynomial_pair(order, alpha, x)
    if not derivative:
        raise ArithmeticError('zero Laguerre derivative')
    delta = rounding(value/derivative)
    return {'evaluation_point': x, 'updated_root': rounding(x-delta),
            'step': delta, 'previous': previous, 'derivative': derivative}


def weight_from_evaluation(order, factor, previous, derivative):
    """Use Gamma(n+alpha)/Gamma(n) and the supplied evaluation operands.

    Supplying operands from before the final root update deliberately retains
    that dependency; this function does not silently refresh or normalize them.
    """
    factor, previous, derivative = map(F, (factor, previous, derivative))
    if not isinstance(order, int) or order < 1 or factor <= 0:
        raise ValueError('positive order and normalization factor required')
    denominator = rounding(rounding(order*previous)*derivative)
    if not denominator:
        raise ArithmeticError('zero quadrature weight denominator')
    weight = rounding(-factor/denominator)
    if weight <= 0:
        raise ArithmeticError('nonpositive quadrature weight')
    return weight
