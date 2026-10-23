from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from bvmethods_laguerre_rounding_audit import polynomial_pair, newton_step, weight_from_evaluation


class LaguerreRoundingAuditTests(unittest.TestCase):
    def test_second_order_closed_form(self):
        # L_2^(1/2)(x) = x^2/2 - 5x/2 + 15/8.
        for x in [F(1, 2), F(3, 2), F(4)]:
            value, _, derivative = polynomial_pair(2, F(1, 2), x)
            self.assertEqual(value, x*x/2 - 5*x/2 + F(15, 8))
            self.assertEqual(derivative, x-F(5, 2))

    def test_derivative_matches_resolvable_direction(self):
        x, step = F(3, 2), F(1, 2**20)
        minus = polynomial_pair(16, F(1, 2), x-step)[0]
        plus = polynomial_pair(16, F(1, 2), x+step)[0]
        derivative = polynomial_pair(16, F(1, 2), x)[2]
        self.assertLess(abs((plus-minus)/(2*step)/derivative-1), F(1, 10**10))

    def test_updated_root_and_weight_evaluation_remain_distinct(self):
        with localcontext() as context:
            context.prec = 80
            root = F(D('2.5')-D('2.5').sqrt())
        initial = root+F(1, 10**6)
        result = newton_step(2, F(1, 2), initial)
        fresh = polynomial_pair(2, F(1, 2), result['updated_root'])
        self.assertLess(abs(fresh[0]), abs(polynomial_pair(2, F(1, 2), initial)[0])/10**5)
        factor = F(1)  # Common positive factor cancels in this dependency test.
        stale_weight = weight_from_evaluation(2, factor, result['previous'], result['derivative'])
        fresh_weight = weight_from_evaluation(2, factor, fresh[1], fresh[2])
        self.assertGreater(stale_weight, 0)
        self.assertGreater(fresh_weight, 0)
        self.assertGreater(abs(stale_weight/fresh_weight-1), F(1, 10**8))
        self.assertEqual(result['evaluation_point'], initial)

    def test_invalid_domains(self):
        for n, alpha, x in [(0, 0, 1), (513, 0, 1), (2, -1, 1), (2, 0, 0)]:
            with self.assertRaises(ValueError):
                polynomial_pair(n, alpha, x)
        with self.assertRaises(ArithmeticError):
            newton_step(2, F(1, 2), F(5, 2))
        with self.assertRaises(ArithmeticError):
            weight_from_evaluation(2, 1, 0, 1)


if __name__ == '__main__':
    unittest.main()
