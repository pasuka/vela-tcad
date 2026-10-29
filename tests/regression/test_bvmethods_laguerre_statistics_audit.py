from decimal import Decimal as D, localcontext
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from bvmethods_laguerre_statistics_audit import half_rule, evaluate, inverse


class LaguerreStatisticsAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rule = half_rule(24)

    def test_positive_rule_and_exact_polynomial_moments(self):
        self.assertTrue(all(x > 0 and w > 0 for x, w in self.rule))
        with localcontext() as context:
            context.prec = 90
            moment = D(1)
            for k in range(9):
                if k:
                    moment *= D(k) + D('.5')
                observed = sum(w*x**k for x, w in self.rule)
                self.assertLess(abs(observed/moment - 1), D('1e-70'))

    def test_dilute_normalization(self):
        with localcontext() as context:
            context.prec = 90
            value, derivative = evaluate(self.rule, -40)
            z = D(-40).exp()
            self.assertLess(abs(value/z - (1 - z/(2*D(2).sqrt()))), D('1e-28'))
            self.assertGreater(value, derivative)
            self.assertGreater(derivative, 0)

    def test_analytic_derivative_and_monotonicity(self):
        with localcontext() as context:
            context.prec = 90
            step = D('1e-7')
            for eta in [D(-2), D(0), D(8)]:
                minus = evaluate(self.rule, eta-step)[0]
                value, derivative = evaluate(self.rule, eta)
                plus = evaluate(self.rule, eta+step)[0]
                self.assertTrue(minus < value < plus)
                self.assertLess(abs((plus-minus)/(2*step)/derivative - 1), D('1e-13'))

    def test_inverse_uses_same_operator(self):
        with localcontext() as context:
            context.prec = 90
            results = [inverse(self.rule, D(s)) for s in ['.03', '1', '18']]
            self.assertTrue(all(a[0] < b[0] for a, b in zip(results, results[1:])))
            for row, target in zip(results, [D('.03'), D(1), D(18)]):
                self.assertLess(abs(row[1]/target-1), D('1e-65'))

    def test_invalid_domains(self):
        for order in [0, 1, 513, 3.5]:
            with self.assertRaises(ValueError):
                half_rule(order)
        for eta in [D('Infinity'), D(151)]:
            with self.assertRaises(ValueError):
                evaluate(self.rule, eta)
        for ratio in [D(0), D(-1), D(101)]:
            with self.assertRaises(ValueError):
                inverse(self.rule, ratio)


if __name__ == '__main__':
    unittest.main()
