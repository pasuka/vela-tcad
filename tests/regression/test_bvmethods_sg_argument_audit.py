from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from bvmethods_sg_argument_audit import (
    logarithm_control, round_significand, scalar_fermi_argument,
)


class SgArgumentAuditTests(unittest.TestCase):
    def test_uniform_equilibrium_and_reference_shift(self):
        for carrier in ['electron', 'hole']:
            self.assertEqual(scalar_fermi_argument(0, 7, 7, -5, -5, carrier), 0)
            # A common statistical reference cancels before composition.
            a = scalar_fermi_argument(F(1, 16), 7, 7, -5, F(-39, 8), carrier)
            b = scalar_fermi_argument(F(1, 16), 7, 7, 95, F(761, 8), carrier)
            self.assertEqual(a, b)

    def test_dilute_nonuniform_dos_matches_density_ratio(self):
        # DOS doubles and the Boltzmann exponent increases by log(2):
        # the population ratio is four, with gamma=1 at both endpoints.
        log2 = logarithm_control(2)
        for carrier, sign in [('electron', -1), ('hole', 1)]:
            a = scalar_fermi_argument(sign * log2, 3, 6, 0, 0, carrier)
            with localcontext() as context:
                context.prec = 80
                error = abs(Decimal(a.numerator) / Decimal(a.denominator) - Decimal(4).ln())
                self.assertLess(error, Decimal(2)**-62)

    def test_endpoint_reversal_in_equal_dos_control(self):
        for carrier in ['electron', 'hole']:
            a = scalar_fermi_argument(F(1, 32), 9, 9, F(-17, 8), -2, carrier)
            b = scalar_fermi_argument(F(-1, 32), 9, 9, -2, F(-17, 8), carrier)
            self.assertEqual(a, -b)

    def test_resolvable_signed_band_direction(self):
        step = F(1, 2**40)
        for carrier, sign in [('electron', -1), ('hole', 1)]:
            plus = scalar_fermi_argument(F(1, 8) + step, 2, 2, -3, -3, carrier)
            minus = scalar_fermi_argument(F(1, 8) - step, 2, 2, -3, -3, carrier)
            self.assertEqual((plus - minus) / (2 * step), sign)

    def test_rounding_ties_even_and_odd_symmetry(self):
        half_step = F(1, 2**64)
        self.assertEqual(round_significand(1 + half_step), 1)
        self.assertEqual(round_significand(1 + 3 * half_step), 1 + F(1, 2**62))
        self.assertEqual(round_significand(-1 - half_step), -1)
        for x in [F(7, 3), F(1, 10**30), F(10**20, 3)]:
            self.assertEqual(round_significand(-x), -round_significand(x))

    def test_invalid_domains_fail(self):
        for dos in [0, -1]:
            with self.assertRaises(ValueError):
                scalar_fermi_argument(0, dos, 1, 0, 0, 'electron')
        with self.assertRaises(ValueError):
            scalar_fermi_argument(0, 1, 1, 0, 0, 'other')
        with self.assertRaises(ValueError):
            logarithm_control(0)


if __name__ == '__main__':
    unittest.main()
