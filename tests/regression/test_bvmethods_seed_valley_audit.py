from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from bvmethods_seed_valley_audit import (
    parabolic_valley_weights, polynomial_density_seed, asymptotic_density_seed, valley_seed_shift)
from bvmethods_sg_argument_audit import round_significand as rounding


def math_control(operation, *args):
    with localcontext() as context:
        context.prec = 100
        values = [D(F(x).numerator)/D(F(x).denominator) for x in args]
        return rounding(F(operation(*values)))


sqrt = lambda x: math_control(lambda v: v.sqrt(), x)
log = lambda x: math_control(lambda v: v.ln(), x)
exp = lambda x: math_control(lambda v: v.exp(), x)
power = lambda x, y: math_control(lambda v, p: v**p, x, y)


class SeedValleyAuditTests(unittest.TestCase):
    def test_isotropic_density_of_states_and_degeneracy(self):
        weights = parabolic_valley_weights([(1, 1, 1), (4, 4, 4)], [2, 1], square_root=sqrt)
        self.assertLess(abs(weights[0]-F(1, 5)), F(1, 2**63))
        self.assertLess(abs(weights[1]-F(4, 5)), F(1, 2**63))

    def test_equivalent_ellipsoids_and_common_mass_scale(self):
        original = [(1, 2, 8), (4, 4, 1)]
        equal = parabolic_valley_weights(original, [1, 1], square_root=sqrt)
        self.assertEqual(equal, (F(1, 2), F(1, 2)))
        scaled = [tuple(4*x for x in m) for m in original]
        self.assertEqual(parabolic_valley_weights(scaled, [1, 1], square_root=sqrt), equal)

    def test_dilute_seed_has_log_limit_and_first_density_correction(self):
        r = F(1, 2**20)
        seed = polynomial_density_seed(r, [F(1, 3), F(-1, 100), 0, 0], logarithm=log)
        self.assertLess(abs((seed-log(r))/r-F(1, 3)), F(1, 10**7))

    def test_asymptotic_seed_is_monotone_and_has_fermi_scaling(self):
        evaluate = lambda r: asymptotic_density_seed(r, [F(3, 4), 1], F(4, 3), 1, 0,
                                                    power=power, square_root=sqrt)
        self.assertGreater(evaluate(1000), evaluate(100))
        # In the high-density limit eta scales as density**(2/3).
        self.assertLess(abs(evaluate(8000000)/evaluate(1000000)-4), F(1, 10**7))

    def test_partition_shift_retains_nonunit_mass_and_energy_offset(self):
        self.assertEqual(valley_seed_shift(3, [1, 1], [0, 0], exponential=exp, logarithm=log),
                         rounding(F(3)-log(2)))
        # A common valley energy offset translates the seed oppositely.
        shifted = valley_seed_shift(3, [F(1, 2), F(1, 2)], [1, 1], exponential=exp, logarithm=log)
        self.assertLess(abs(shifted-2), F(1, 2**61))

    def test_resolvable_log_seed_direction(self):
        r, h = F(1, 16), F(1, 2**20)
        evaluate = lambda x: polynomial_density_seed(x, [0, 0, 0, 0], logarithm=log)
        direction = (evaluate(r+h)-evaluate(r-h))/(2*h)
        self.assertLess(abs(direction*r-1), F(1, 10**9))

    def test_invalid_domains(self):
        for masses, g in [([], []), ([(1, 2)], [1]), ([(1, 0, 1)], [1]), ([(1, 1, 1)], [0])]:
            with self.assertRaises(ValueError):
                parabolic_valley_weights(masses, g, square_root=sqrt)
        with self.assertRaises(ValueError):
            polynomial_density_seed(0, [0]*4, logarithm=log)
        with self.assertRaises(ValueError):
            asymptotic_density_seed(1, [1, 1], 1, 2, 0, power=power, square_root=sqrt)
        with self.assertRaises(ValueError):
            valley_seed_shift(1, [1], [], exponential=exp, logarithm=log)


if __name__ == '__main__':
    unittest.main()
