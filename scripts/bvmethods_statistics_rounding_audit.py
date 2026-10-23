"""Independent finite-significand statistics controls, outside the solver.

These controls cover normal operands with an unbounded exponent. Mathematical
log and exp use Decimal; no vendor transcendental implementation is assumed.
"""
from fractions import Fraction

from bvmethods_sg_argument_audit import logarithm_control, round_significand


def log_gamma_control(population, dos, eta, *, logarithm=logarithm_control):
    """Separate population division, logarithm evaluation and eta subtraction.

    An explicit logarithm callable permits a captured-result control without
    confusing it with independently generated statistics. Eta is supplied,
    not inverted from a defining Fermi integral by this function.
    """
    population, dos, eta = map(Fraction, (population, dos, eta))
    if population <= 0 or dos <= 0:
        raise ValueError('positive population and DOS required')
    ratio = round_significand(population / dos)
    log_value = Fraction(logarithm(ratio))
    return round_significand(log_value - eta)


def bernoulli_pair_series(value):
    """Return B(x), B'(x), B(-x), B'(-x) for |x| <= 1/2.

    Sum S=(exp(t)-1)/t and -S' at t=abs(x), stopping when both sums
    cease changing at the controlled significand. Form the opposite pair
    through B(-t)=B(t)+t and B'(-t)=-1-B'(t), with separate rounding.
    This is a diagnostic arithmetic candidate, not a universal error bound.
    """
    value = Fraction(value)
    t = abs(value)
    if t > Fraction(1, 2):
        raise ValueError('series control requires abs(x) <= 1/2')
    rounding = round_significand
    term, total, negative_derivative = Fraction(1), Fraction(1), Fraction(-1, 2)
    for n in range(1, 256):
        term = rounding(term * rounding(t / n))
        next_total = rounding(total + rounding(term / (n + 1)))
        next_derivative = rounding(negative_derivative - rounding(term / (n + 2)))
        unchanged = (next_total == total and next_derivative == negative_derivative)
        total, negative_derivative = next_total, next_derivative
        if unchanged:
            break
    else:
        raise ArithmeticError('series control did not stabilize')
    positive = rounding(1 / total)
    derivative = rounding(negative_derivative * rounding(positive * positive))
    negative = rounding(positive + t)
    reverse_derivative = rounding(-derivative - 1)
    result = (positive, derivative, negative, reverse_derivative)
    return result if value >= 0 else (negative, reverse_derivative, positive, derivative)
