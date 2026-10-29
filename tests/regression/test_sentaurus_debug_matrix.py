"""Parsing and sign/storage properties of actual native debug formats."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location("debug_matrix", Path(__file__).resolve().parents[2] / "scripts/compare_sentaurus_debug_matrix.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class NativeMatrixTests(unittest.TestCase):
    def test_global_backward_error_does_not_certify_weak_rows(self):
        item = {"matrix": {(0, 0): 1.0, (1, 1): 1e-30},
                "rhs": [1.0, 2e-30], "solution": [1.0, 1.0]}
        result = module.super_closure(item)
        self.assertLess(result["linear_global_l1_backward_error"], 1e-29)
        self.assertAlmostEqual(result["linear_componentwise_max_backward_error"], 1/3)
        self.assertEqual(result["linear_componentwise_worst_row"], 1)
        # Positive row scaling preserves each componentwise relative defect.
        scaled = {"matrix": {(0, 0): 1e-20, (1, 1): 1e10},
                  "rhs": [1e-20, 2e10], "solution": [1.0, 1.0]}
        self.assertAlmostEqual(module.super_closure(scaled)[
            "linear_componentwise_max_backward_error"], 1/3)

    def test_exact_zero_row_has_zero_defect(self):
        result = module.super_closure({"matrix": {}, "rhs": [0.0], "solution": [0.0]})
        self.assertEqual(result["linear_componentwise_max_backward_error"], 0.0)

    def test_native_negative_density_rule_does_not_floor_small_positive(self):
        old = {1: 2.0, 2: 2.0, 3: 2.0}
        delta = {1: -3.0, 2: -1.999, 3: 0.5}
        accepted = {1: 0.02*1e10, 2: (2-1.999)*1e10, 3: 2.5*1e10}
        result = module.density_update_audit(old, delta, accepted, 1e10)
        self.assertEqual(result["raw_update_nonpositive_count"], 1)
        self.assertEqual(result["observed_rule_max_relative_error"], 0)
        self.assertEqual(result["negative_update_accepted_over_old"], [0.01])

    def test_small_rows_have_no_unit_error_floor(self):
        item = {"matrix": {(0,0): 1e-30}, "rhs": [2e-30], "solution": [1.0]}
        self.assertAlmostEqual(module.super_closure(item)["linear_componentwise_l1_backward_error"], 1/3)

    def test_labeled_asymmetric_matrix_and_native_sign(self):
        fixture = '''Newton loop 0
VAR:Variable ElectrostaticPotential:
VAR:#7 = 0
VAR:#9 = 0
RHS1:Block 0
RHS1:#7 : 12
RHS1:#9 : 13
JAC(0,0)d:#7 : 2
JAC(0,0)d:#9 : 3
JAC(0,0)u:#(7,9) : 5
JAC(0,0)l:#(9,7) : 7
Solution:
Block 0
SOL:#7 : 1
SOL:#9 : 2
RHS2:Block 0
RHS2:#7 : 0
RHS2:#9 : 0
'''
        parsed = list(module.iterations(fixture.splitlines()))
        check = module.closure(parsed[0])["ElectrostaticPotential"]
        self.assertEqual(check["J_SOL_minus_RHS1_relative"], 0)
        self.assertGreater(check["J_SOL_plus_RHS1_relative"], 0.5)
        with self.assertRaises(ValueError):
            list(module.iterations((fixture + 'JAC(0,0)d:#7 : 2\n').splitlines()))
        parsed[0]["solution"].pop((0,9))
        with self.assertRaises(ValueError):
            module.closure(parsed[0])

    def test_super_msr_preserves_both_triangles(self):
        fixture = '''new matrix:
ja:
size = 4
4
5
5
2
matrix:
size = 5
2
3
1
5
7
rhs:
size = 2
12
13
result:
size = 2
1
2
'''
        item = next(module.super_iterations(fixture))
        self.assertEqual(item["matrix"][0,1], 5)
        self.assertEqual(item["matrix"][1,0], 7)
        self.assertEqual(module.super_closure(item)["linear_componentwise_l1_backward_error"], 0)
        item["solution"][0] += 0.1
        self.assertGreater(module.super_closure(item)["linear_componentwise_l1_backward_error"], 0.01)
        with self.assertRaises(ValueError):
            list(module.super_iterations(fixture.replace('size = 5', 'size = 6')))
        with self.assertRaises(ValueError):
            list(module.super_iterations(fixture.replace('5\n7\nrhs:', 'nan\n7\nrhs:')))


if __name__ == "__main__":
    unittest.main()
