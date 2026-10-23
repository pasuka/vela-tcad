"""Scalar isothermal SG argument audit; no production solver integration.

Inputs are normalized band drive, DOS and log(F_half(eta)/exp(eta)).
The significand control covers finite, normal diagnostic operands; it does
not emulate IEEE exponent limits, underflow, exceptions or a vendor logl.
"""
from decimal import Decimal, localcontext
from fractions import Fraction


def round_significand(value, bits=64):
    """Round a rational to nearest, ties to even, with an unbounded exponent."""
    if not isinstance(bits, int) or bits < 2:
        raise ValueError('at least two significand bits required')
    value = Fraction(value)
    if not value:
        return value
    magnitude = abs(value)
    exponent = magnitude.numerator.bit_length() - magnitude.denominator.bit_length()
    power = lambda n: Fraction(2**n) if n >= 0 else Fraction(1, 2**-n)
    if magnitude < power(exponent):
        exponent -= 1
    step = power(exponent - bits + 1)
    scaled = magnitude / step
    integer, remainder = divmod(scaled.numerator, scaled.denominator)
    if (2 * remainder > scaled.denominator or
            (2 * remainder == scaled.denominator and integer % 2)):
        integer += 1
    return (1 if value > 0 else -1) * integer * step


def logarithm_control(value):
    """100-decimal-digit defining log rounded to a 64-bit significand."""
    value = Fraction(value)
    if value <= 0:
        raise ValueError('positive logarithm operand required')
    with localcontext() as context:
        context.prec = 100
        result = (Decimal(value.numerator) / Decimal(value.denominator)).ln()
        return round_significand(result)


def scalar_fermi_argument(band_drive, dos0, dos1, log_gamma0, log_gamma1,
                         carrier):
    """Compose the scalar 300 K audit chain with separate operation rounding.

    This reconstructs an argument from supplied statistics; it does not
    generate log-gamma, qualify nonuniform temperature or imply that an
    imported native intermediate is an autonomous material model.
    """
    if carrier not in ('electron', 'hole'):
        raise ValueError('electron or hole carrier required')
    d, n0, n1, g0, g1 = map(Fraction,
                             (band_drive, dos0, dos1, log_gamma0, log_gamma1))
    if n0 <= 0 or n1 <= 0:
        raise ValueError('positive DOS required')
    sign = -1 if carrier == 'electron' else 1
    dos_term = logarithm_control(round_significand(n1 / n0))
    gamma_term = round_significand(g1 - g0)
    result = round_significand(d + sign * dos_term)
    return sign * round_significand(result + sign * gamma_term)
