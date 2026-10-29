"""Independent Decimal reference for the explicit BVMethods JD/SG contracts.

This module does not run a solver or change any physical model. Input floats
are converted to their exact binary64 values before reference evaluation.
"""
from __future__ import annotations

from decimal import Decimal, localcontext

from derive_joyce_dixon_coefficients import derive

PI = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128")


def exact(value):
    return Decimal.from_float(value) if isinstance(value, float) else Decimal(value)


class JdReference:
    """Real JD law at high precision, with separately selectable constants.

    ideal: independently reverted coefficients and analytic pi constants.
    binary_constants: binary64 coefficients, norm, c and 2/3 exponent, with
    high-precision evaluation. The latter separates stored constants from
    runtime evaluation, without fitting a native value.
    """

    def __init__(self, digits=80, constants="ideal"):
        if not 40 <= digits <= 120:
            raise ValueError("Reference precision must be between 40 and 120 digits")
        self.digits = digits
        with localcontext() as ctx:
            ctx.prec = digits
            if constants == "ideal":
                self.coefficients = [Decimal(x) for x in derive(digits + 10)["inverse_log_coefficients"]]
                self.norm = Decimal(4) / (3 * PI.sqrt())
                self.c = PI * PI / 6
                self.exponent = Decimal(2) / 3
            elif constants == "binary_constants":
                self.coefficients = [exact(x) for x in [.3535533905932738, -.004950089729875255,
                    .00014838577128872333, -4.425630118996707e-6]]
                self.norm = exact(.75225277806367504926)
                self.c = exact(1.64493406684822643647)
                self.exponent = exact(2. / 3.)
            else:
                raise ValueError("Unknown constants policy")
            lo, hi = Decimal(8), Decimal(9)
            for _ in range(4 * digits):
                mid = (lo + hi) / 2
                if mid == lo or mid == hi:
                    break
                if self._low(mid)[1] > self._high(mid)[1]:
                    lo = mid
                else:
                    hi = mid
            self.join = (lo + hi) / 2
            self.shift = self._high(self.join)[0] - self._low(self.join)[0]

    def _low(self, r):
        a, b, c, d = self.coefficients
        polynomial = r * (a + r * (b + r * (c + r * d)))
        derivative = a + r * (2 * b + r * (3 * c + 4 * d * r))
        return r.ln() + polynomial, 1 / r + derivative, -polynomial, -derivative

    def _high(self, r):
        u = ((r / self.norm).ln() * self.exponent).exp()
        h = (u * u - self.c).sqrt()
        return h, self.exponent * u * u / (r * h)

    def evaluate(self, ratio):
        with localcontext() as ctx:
            ctx.prec = self.digits
            r = exact(ratio)
            if not r.is_finite() or r <= 0:
                raise ValueError("JD reference requires a positive finite ratio")
            if r < self.join:
                eta, slope, gamma, gamma_slope = self._low(r)
            else:
                h, slope = self._high(r)
                eta = h - self.shift
                gamma = r.ln() - eta
                gamma_slope = 1 / r - slope
            return {"eta": +eta, "eta_prime": +slope, "g": +gamma,
                    "g_prime": +gamma_slope, "r_g_prime": +(r * gamma_slope)}


def bernoulli(argument):
    """Caller selects Decimal context. Stable positive pair, B(-a)=B(a)+a."""
    x = exact(argument)
    a = abs(x)
    if a == 0:
        positive = Decimal(1)
    elif a < Decimal("0.01"):
        term = total = Decimal(1)
        for k in range(1000):
            term *= a / (k + 1)
            update = total + term / (k + 2)
            if update == total:
                break
            total = update
        positive = 1 / total
    else:
        positive = a / (a.exp() - 1)
    return (positive, positive + a) if x >= 0 else (positive + a, positive)


def sg_difference(population_from, population_to, argument):
    """Unscaled c_i B(-a)-c_j B(a), without a population floor."""
    b, opposite = bernoulli(argument)
    return exact(population_from) * opposite - exact(population_to) * b
