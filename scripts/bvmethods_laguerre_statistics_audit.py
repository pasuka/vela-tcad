"""Independent generalized Gauss-Laguerre Fermi-half diagnostic.

NumPy eigenvalues seed roots only. Decimal refines each polynomial root and
generates normalized weights without importing native tables or fitting data.
This finite-order quadrature is distinct from the defining-integral audit and
is not connected to the production carrier statistics.
"""
from decimal import Decimal as D, localcontext
from functools import lru_cache
import numpy as np


def half_rule(order):
    """Generate nodes and weights for x**(1/2)*exp(-x)/Gamma(3/2)."""
    if not isinstance(order, int) or not 2 <= order <= 512:
        raise ValueError('audit order must be an integer in [2, 512]')
    alpha = D('0.5')
    index = np.arange(order, dtype=float)
    diagonal = 2 * index + 1.5
    off = np.sqrt((index[1:]) * (index[1:] + .5))
    seeds = np.linalg.eigvalsh(np.diag(diagonal) + np.diag(off, 1) + np.diag(off, -1))
    with localcontext() as context:
        context.prec = 100
        def polynomial(x):
            previous, current = D(0), D(1)
            for k in range(1, order + 1):
                previous, current = current, ((2*k - 1 + alpha - x)*current
                                               - (k - 1 + alpha)*previous) / k
            derivative = (order*current - (order + alpha)*previous) / x
            return current, derivative
        normalization = D(1)
        for k in range(1, order + 1):
            normalization *= 1 + alpha / k
        rule = []
        for seed in seeds:
            x = D.from_float(float(seed))
            for _ in range(12):
                value, derivative = polynomial(x)
                step = value / derivative
                x -= step
                if abs(step) < D('1e-80') * max(D(1), abs(x)):
                    break
            else:
                raise ArithmeticError('Laguerre root did not converge')
            value, derivative = polynomial(x)
            weight = normalization / (x * derivative**2)
            if x <= 0 or weight <= 0:
                raise ArithmeticError('nonpositive quadrature operand')
            rule.append((+x, +weight))
        if any(a[0] >= b[0] for a, b in zip(rule, rule[1:])):
            raise ArithmeticError('quadrature roots not distinct and ordered')
        if abs(sum(w for _, w in rule) - 1) > D('1e-70'):
            raise ArithmeticError('normalized quadrature mass failed')
        return tuple(rule)


@lru_cache(maxsize=4)
def _prepared(rule):
    with localcontext() as context:
        context.prec = 100
        return tuple(((-x).exp(), weight) for x, weight in rule)


def evaluate(rule, eta):
    """Normalized F_1/2 and its eta derivative using a finite-order rule."""
    eta = D(eta)
    if not eta.is_finite() or not -150 <= eta <= 150:
        raise ValueError('audit eta outside [-150, 150]')
    with localcontext() as context:
        context.prec = 100
        e = (-eta).exp()
        value, derivative = D(0), D(0)
        for exponential, weight in _prepared(tuple(rule)):
            denominator = exponential + e
            value += weight / denominator
            derivative += weight * e / denominator**2
        return +value, +derivative


def inverse(rule, population_ratio):
    """Invert the same finite quadrature with a tight log-density Newton step."""
    ratio = D(population_ratio)
    if not ratio.is_finite() or not D('1e-40') <= ratio <= 100:
        raise ValueError('audit population ratio outside [1e-40, 100]')
    with localcontext() as context:
        context.prec = 100
        eta = ratio.ln() if ratio < 1 else ratio**(D(2)/3)
        for _ in range(50):
            value, derivative = evaluate(rule, eta)
            step = (value / ratio).ln() * value / derivative
            step = max(D(-10), min(D(10), step))
            eta -= step
            if abs(step) < D('1e-70'):
                value, derivative = evaluate(rule, eta)
                if abs(value/ratio - 1) > D('1e-65'):
                    raise ArithmeticError('quadrature inverse residual failed')
                return +eta, +value, +derivative
        raise ArithmeticError('quadrature inverse did not converge')
