"""Independent finite-significand seed and parabolic valley controls.

These helpers are diagnostic only. Material inputs, approximation coefficients,
and transcendental providers are explicit; no native seeds or weights are read.
The arithmetic model has an unbounded exponent and covers normal operands.
"""
from fractions import Fraction as F
from bvmethods_sg_argument_audit import round_significand as rounding


def parabolic_valley_weights(principal_masses, degeneracies, *, square_root):
    """Normalize g*sqrt(m1*m2*m3) in valley order at the controlled precision.

    Masses are relative to the free-electron mass. This is an unstressed,
    parabolic-band control; it does not model nonparabolic or strain effects.
    """
    masses = tuple(tuple(map(F, m)) for m in principal_masses)
    degeneracies = tuple(map(F, degeneracies))
    if (not masses or len(masses) != len(degeneracies)
            or any(len(m) != 3 or min(m) <= 0 for m in masses)
            or min(degeneracies) <= 0):
        raise ValueError('three positive principal masses and positive degeneracy per valley required')
    factors = []
    total = F(0)
    for m, g in zip(masses, degeneracies):
        determinant = rounding(rounding(m[0] * m[1]) * m[2])
        factor = rounding(g * F(square_root(determinant)))
        if factor <= 0:
            raise ArithmeticError('nonpositive density-of-states factor')
        factors.append(factor)
        total = rounding(total + factor)
    return tuple(rounding(f / total) for f in factors)


def polynomial_density_seed(ratio, coefficients, *, logarithm):
    """Evaluate log(ratio) + ratio*(c0+ratio*(c1+ratio*(c2+ratio*c3)))."""
    ratio = F(ratio)
    coefficients = tuple(map(F, coefficients))
    if ratio <= 0 or len(coefficients) != 4:
        raise ValueError('positive population/DOS ratio and four coefficients required')
    term = rounding(coefficients[3] * ratio)
    for coefficient in reversed(coefficients[:3]):
        term = rounding(rounding(term + coefficient) * ratio)
    return rounding(F(logarithm(ratio)) + term)


def asymptotic_density_seed(ratio, factors, exponent, correction, shift, *, power, square_root):
    """Evaluate sqrt((ratio*f0*f1)**exponent-correction)-shift in order."""
    ratio = F(ratio)
    factors = tuple(map(F, factors))
    exponent, correction, shift = map(F, (exponent, correction, shift))
    if ratio <= 0 or len(factors) != 2 or min(factors) <= 0 or exponent <= 0:
        raise ValueError('positive ratio, two positive factors and positive exponent required')
    base = rounding(rounding(ratio * factors[0]) * factors[1])
    radicand = rounding(F(power(base, exponent)) - correction)
    if radicand < 0:
        raise ValueError('asymptotic seed outside its valid density range')
    return rounding(F(square_root(radicand)) - shift)


def valley_seed_shift(seed, weights, offsets, *, exponential, logarithm):
    """Retain the rounded valley partition sum in the initial seed correction."""
    weights, offsets = tuple(map(F, weights)), tuple(map(F, offsets))
    if not weights or len(weights) != len(offsets) or min(weights) <= 0:
        raise ValueError('matching nonempty valley arrays with positive weights required')
    total = F(0)
    for weight, offset in zip(weights, offsets):
        factor = F(exponential(offset))
        if factor <= 0:
            raise ArithmeticError('nonpositive exponential')
        total = rounding(total + rounding(weight * factor))
    return rounding(F(seed) - F(logarithm(total)))
