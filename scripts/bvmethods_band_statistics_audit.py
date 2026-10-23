"""Decimal band/eta edge audit, independent of exported QF and Vela Fermi code."""
from decimal import Decimal as D


def qf_from_band(band, eta, vt, carrier):
    if carrier not in ('electron','hole') or not all(x.is_finite() for x in (band,eta,vt)) or vt<=0:
        raise ValueError('carrier and positive thermal voltage required')
    return -band+(-1 if carrier=='electron' else 1)*vt*eta


def edge_driver(band0, band1, eta0, eta1, vt, carrier):
    if carrier not in ('electron','hole') or not all(x.is_finite() for x in (band0,band1,eta0,eta1,vt)) or vt<=0:
        raise ValueError('carrier and positive thermal voltage required')
    band_term=-(band1-band0)
    statistics_term=(-1 if carrier=='electron' else 1)*vt*(eta1-eta0)
    return {'band_term_V':band_term,'statistics_term_V':statistics_term,
            'qf_drop_V':band_term+statistics_term}


def explicit_driver_flux(n0,n1,drop,vt,carrier):
    """Current without mobility, per cm; caller supplies independent eta/bands.

    Decimal operations retain cancellation in the driver supplied by the
    caller. They cannot restore bits already lost in a TDR band or population.
    """
    if carrier not in ('electron','hole') or not all(x.is_finite() for x in (n0,n1,drop,vt)) or min(n0,n1,vt)<=0:
        raise ValueError('carrier, populations and thermal voltage required')
    sign=1 if carrier=='electron' else -1
    delta=sign*drop/vt
    if delta==0:return D(0)
    arg=(n1/n0).ln()+delta
    def bernoulli(x):return D(1) if x==0 else x/(x.exp()-1)
    if delta>0:value=n0*bernoulli(-arg)*(1-(-delta).exp())
    else:value=n1*bernoulli(arg)*(delta.exp()-1)
    return -sign*D('1.602192e-19')*vt*value
