"""Independent constitutive properties, with resolvable finite steps."""
import sys
import unittest
from decimal import Decimal as D, localcontext
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from bvmethods_statistics_derivative_audit import statistics_density_derivatives, sg_density_gradient


def pair(a):
    if a == 0:
        return D(1), -D('.5'), D(1), -D('.5')
    e = a.exp()
    b = a / (e-1)
    bp = ((e-1)-a*e) / (e-1)**2
    return b, bp, b+a, -1-bp


class DerivativeProperties(unittest.TestCase):
    def test_boltzmann_cancels_gamma_slope(self):
        with localcontext() as ctx:
            ctx.prec = 70
            n, N = D('7'), D('11')
            eta, gamma = statistics_density_derivatives(n, N, n/N)
            self.assertLess(abs(eta-1/n), D('1e-68'))
            self.assertLess(abs(gamma), D('1e-68'))

    def test_accepted_point_matters(self):
        fresh = statistics_density_derivatives('8', '1', '3')[0]
        stale = statistics_density_derivatives('8', '1', '2.99')[0]
        self.assertNotEqual(fresh, stale)

    def test_independent_power_density_finite_direction(self):
        # Defining constitutive law rho=eta**(3/2), independent of the SG helper.
        with localcontext() as ctx:
            ctx.prec = 70
            n0, n1, N0, N1 = map(D, ('9', '15', '2', '3'))
            def eta(n, N):
                return ((n/N).ln()*D(2)/3).exp()
            def flux(x, y):
                g0, g1 = (x/N0).ln()-eta(x, N0), (y/N1).ln()-eta(y, N1)
                a = D('.37')+(N1/N0).ln()+g1-g0
                b, _, bm, _ = pair(a)
                return x*bm-y*b
            e0, e1 = eta(n0, N0), eta(n1, N1)
            fp0, fp1 = D('1.5')*e0.sqrt(), D('1.5')*e1.sqrt()
            a = D('.37')+(n1/N1).ln()-e1-(n0/N0).ln()+e0+(N1/N0).ln()
            d0, d1 = sg_density_gradient(n0, n1, N0, N1, fp0, fp1, pair(a))
            v0, v1, h = D('.4'), D('-.7'), D('1e-7')
            fd = (flux(n0+h*v0, n1+h*v1)-flux(n0-h*v0, n1-h*v1))/(2*h)
            self.assertLess(abs(fd-(d0*v0+d1*v1)), D('1e-16'))

    def test_endpoint_reversal(self):
        with localcontext() as ctx:
            ctx.prec = 60
            a = D('.21')
            f = sg_density_gradient('2', '3', '5', '7', '.3', '.4', pair(a))
            r = sg_density_gradient('3', '2', '7', '5', '.4', '.3', pair(-a))
            self.assertLess(abs(f[0]+r[1])+abs(f[1]+r[0]), D('1e-57'))

    def test_common_density_unit_scale(self):
        with localcontext() as ctx:
            ctx.prec = 60
            a = pair(D('.2'))
            f = sg_density_gradient('2', '3', '5', '7', '.3', '.4', a)
            r = sg_density_gradient('200', '300', '500', '700', '.3', '.4', a)
            self.assertLess(sum(abs(x-y) for x, y in zip(f, r)), D('1e-57'))

    def test_equilibrium_nondegenerate_gradient(self):
        self.assertEqual(sg_density_gradient('2', '2', '4', '4', '.5', '.5', pair(D(0))), (D(1), D(-1)))

    def test_invalid_statistical_domain(self):
        for n, N, fp in [('0','1','1'),('1','-1','1'),('1','1','0'),('NaN','1','1')]:
            with self.assertRaises(ValueError):
                statistics_density_derivatives(n, N, fp)


if __name__ == '__main__':
    unittest.main()
