from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from bvmethods_integral_update_audit import integral_from_exponentials, valley_sum, eta_step
from bvmethods_sg_argument_audit import round_significand as rounding


class IntegralUpdateAuditTests(unittest.TestCase):
    def test_single_level_half_occupation(self):
        value, derivative, _ = integral_from_exponentials([2], [1], 1, 3)
        self.assertEqual(value, 3)
        self.assertEqual(derivative, F(3, 2))

    def test_density_derivative_resolvable_direction(self):
        with localcontext() as context:
            context.prec = 90
            nodes, weights = [D('.2'), D('1.5'), D('4')], [F(1, 8), F(1, 2), F(3, 8)]
            h, eta = D('1e-5'), D('.7')
            exps = [rounding(F((-x).exp())) for x in nodes]
            evaluate = lambda z: integral_from_exponentials(weights, exps, rounding(F((-z).exp())), 1)
            plus, minus, center = evaluate(eta+h), evaluate(eta-h), evaluate(eta)
            finite_difference = (plus[0]-minus[0])/F(2*h)
            self.assertLess(abs(finite_difference/center[1]-1), F(1, 10**10))
            self.assertGreater(center[0], 0)
            self.assertGreater(center[1], 0)

    def test_boltzmann_limit(self):
        value, derivative, _ = integral_from_exponentials([F(1, 4), F(3, 4)], [1, F(1, 8)], 2**40, 1)
        self.assertLess(abs(derivative/value-1), F(1, 2**39))

    def test_valley_mass_is_not_implicitly_normalized(self):
        self.assertEqual(valley_sum([3, 3, 3], [2, 2, 2], [1, 2, 3]), (18, 12))
        self.assertEqual(valley_sum([3, 3], [2, 2], [F(1, 4), F(3, 4)]), (3, 2))

    def test_log_density_step_reduces_exponential_residual(self):
        # For F(eta)=exp(eta), F/F'=1 and log(F/target)=eta-log(target).
        result = eta_step(3, 8, 8, F(1, 2), 10, F(1, 10**9))
        self.assertEqual(result['updated_eta'], F(5, 2))
        self.assertFalse(result['accepted'])
        final = eta_step(F(5, 2), 8, 8, 0, 10, F(1, 10**9))
        self.assertTrue(final['accepted'])

    def test_clip_and_acceptance_use_updated_state(self):
        result = eta_step(4, 2, 1, 10, 1, F(1, 8))
        self.assertEqual(result['step'], 1)
        self.assertEqual(result['updated_eta'], 3)
        self.assertEqual(result['criterion'], rounding(F(8, 25)))
        self.assertFalse(result['accepted'])

    def test_invalid_domains(self):
        for args in [([], [], 1, 1), ([1], [], 1, 1), ([1], [1], 0, 1), ([-1], [1], 1, 1)]:
            with self.assertRaises(ValueError):
                integral_from_exponentials(*args)
        with self.assertRaises(ValueError):
            valley_sum([1], [1, 2], [1])
        with self.assertRaises(ValueError):
            eta_step(1, 1, 0, 0, 1, 1)


if __name__ == '__main__':
    unittest.main()
