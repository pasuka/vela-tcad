"""Density-coordinate derivative controls for fixed-state SG audits.

This module does not invert a Fermi integral or modify the solver. Callers must
provide F'(eta) evaluated at the accepted eta, and the population/DOS state at
that same point. No QF-coordinate Jacobian qualification is implied.
"""
from decimal import Decimal


def statistics_density_derivatives(population, dos, forward_derivative):
    """Return d eta/d n and d logGamma/d n for rho=n/DOS, fixed DOS."""
    n, N, slope = map(Decimal, (population, dos, forward_derivative))
    if not all(x.is_finite() and x > 0 for x in (n, N, slope)):
        raise ValueError('finite positive population, DOS and Fprime required')
    eta_slope = 1 / (N * slope)
    return eta_slope, 1 / n - eta_slope


def sg_density_gradient(n0, n1, dos0, dos1, fprime0, fprime1, pair):
    """Gradient of n0*B(-a)-n1*B(a) at fixed bands, DOS and mobility.

    a=c*band_drive+log(DOS1/DOS0)+Gamma1-Gamma0. Pair order is
    B(a), B'(a), B(-a), B'(-a). Products here use the caller's Decimal
    context to audit the analytic constitutive derivative separately from
    binary80 rounding and finite-step numerical differentiation.
    """
    n0, n1 = Decimal(n0), Decimal(n1)
    _, gamma0 = statistics_density_derivatives(n0, dos0, fprime0)
    _, gamma1 = statistics_density_derivatives(n1, dos1, fprime1)
    B, Bp, Bm, Bmp = map(Decimal, pair)
    if not all(x.is_finite() for x in (B, Bp, Bm, Bmp)):
        raise ValueError('finite Bernoulli pair required')
    da0, da1 = -gamma0, gamma1
    dflux_da = -n0 * Bmp - n1 * Bp
    return Bm + dflux_da * da0, -B + dflux_da * da1
