from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from bvmethods_sg_argument_audit import logarithm_control, round_significand
from bvmethods_statistics_rounding_audit import bernoulli_pair_series, log_gamma_control


class StatisticsRoundingAuditTests(unittest.TestCase):
    def test_zero_limits_and_endpoint_reversal(self):
        self.assertEqual(bernoulli_pair_series(0), (1, F(-1, 2), 1, F(-1, 2)))
        for x in [F(1, 10**12), F(7, 1000), F(1, 2)]:
            p = bernoulli_pair_series(x)
            self.assertEqual(bernoulli_pair_series(-x), (p[2], p[3], p[0], p[1]))
            self.assertGreater(p[0], 0)
            self.assertLess(p[1], 0)
            self.assertLess(p[3], 0)

    def test_values_and_derivatives_against_defining_math(self):
        with localcontext() as context:
            context.prec = 100
            for x in [F(1, 10**12), F(1, 10**4), F(7, 1000), F(1, 2)]:
                actual = bernoulli_pair_series(x)
                d = Decimal(x.numerator) / Decimal(x.denominator)
                for sign, offset in [(1, 0), (-1, 2)]:
                    z = sign * d
                    e = z.exp()
                    expected = [z / (e - 1), ((e - 1) - z * e) / (e - 1)**2]
                    for a, b in zip(actual[offset:offset + 2], expected):
                        error = abs(Decimal(a.numerator) / Decimal(a.denominator) - b)
                        self.assertLess(error, Decimal(2)**-61)

    def test_paired_identities_with_explicit_rounding(self):
        for x in [F(1, 10**20), F(1, 300), F(1, 2)]:
            b, db, bm, dbm = bernoulli_pair_series(x)
            self.assertEqual(bm, round_significand(b + x))
            self.assertEqual(dbm, round_significand(-db - 1))

    def test_near_zero_resolvable_direction(self):
        step = F(1, 2**35)
        plus, minus = bernoulli_pair_series(step)[0], bernoulli_pair_series(-step)[0]
        self.assertLess(abs((plus - minus) / (2 * step) + F(1, 2)), F(1, 2**27))

    def test_log_input_and_subtraction_are_separate(self):
        seen = []
        def supplied_log(ratio):
            seen.append(ratio)
            return F(3, 2) + F(1, 2**63)
        result = log_gamma_control(1, 3, F(3, 2), logarithm=supplied_log)
        self.assertEqual(seen, [round_significand(F(1, 3))])
        self.assertEqual(result, F(1, 2**63))
        self.assertEqual(log_gamma_control(8, 8, -3), 3)
        self.assertEqual(log_gamma_control(2, 1, 0), logarithm_control(2))

    def test_invalid_domains(self):
        for n, d in [(0, 1), (-1, 1), (1, 0), (1, -1)]:
            with self.assertRaises(ValueError):
                log_gamma_control(n, d, 0)
        for x in [F(1), F(-1)]:
            with self.assertRaises(ValueError):
                bernoulli_pair_series(x)


if __name__ == '__main__':
    unittest.main()
