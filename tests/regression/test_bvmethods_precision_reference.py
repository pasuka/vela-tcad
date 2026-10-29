"""Numerical invariants for the independent BVMethods precision reference."""
from decimal import Decimal, localcontext
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "diagnostics"))
from bvmethods_precision_reference import JdReference, bernoulli, sg_difference


class PrecisionReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = JdReference(80)

    def test_statistical_tangents_against_independent_finite_differences(self):
        with localcontext() as ctx:
            ctx.prec = 80
            for text in ["1e-17", ".5", "4", "8", "9", "17.8"]:
                r = Decimal(text)
                h = r * Decimal("1e-22")
                center = self.reference.evaluate(r)
                lo, hi = self.reference.evaluate(r-h), self.reference.evaluate(r+h)
                for value, tangent in [("eta", "eta_prime"), ("g", "g_prime")]:
                    finite = (hi[value]-lo[value])/(2*h)
                    self.assertLess(abs((finite-center[tangent])/center[tangent]), Decimal("1e-35"))

    def test_nondegenerate_limit_retains_nonzero_gamma_tangent(self):
        with localcontext() as ctx:
            ctx.prec = 80
            r = Decimal("1e-30")
            value = self.reference.evaluate(r)
            a1 = self.reference.coefficients[0]
            self.assertLess(abs(value["g"]/r+a1), Decimal("1e-30"))
            self.assertLess(abs(value["g_prime"]+a1), Decimal("1e-30"))
            self.assertGreater(value["eta_prime"], 0)

    def test_join_is_continuous_in_value_and_first_derivative(self):
        with localcontext() as ctx:
            ctx.prec = 80
            t = self.reference.join
            a = self.reference.evaluate(t*(1-Decimal("1e-25")))
            b = self.reference.evaluate(t*(1+Decimal("1e-25")))
            self.assertLess(abs(a["eta"]-b["eta"]), Decimal("1e-23"))
            self.assertLess(abs(a["eta_prime"]-b["eta_prime"]), Decimal("1e-23"))

    def test_sg_balance_and_reversal_without_current_cutoffs(self):
        with localcontext() as ctx:
            ctx.prec = 80
            n = Decimal("1e19")
            for text in ["-100", "-.1", "-1e-20", "0", "1e-20", ".1", "100"]:
                a = Decimal(text)
                p = n*a.exp()
                balanced = sg_difference(n, p, a)
                self.assertLess(abs(balanced)/max(n,p), Decimal("1e-60"))
                forward = sg_difference(n, p*Decimal("1.01"), a)
                reverse = sg_difference(p*Decimal("1.01"), n, -a)
                self.assertEqual(forward, -reverse)
                positive, negative = bernoulli(a)
                if a != 0:
                    self.assertLess(abs(positive*(a.exp()-1)/a-1), Decimal("1e-50"))
                else:
                    self.assertEqual(positive, 1)
                    self.assertEqual(negative, 1)


if __name__ == "__main__":
    unittest.main()
